
import numpy as np

# Patch deprecated NumPy types for compatibility with older libraries
if not hasattr(np, 'bool'):
    np.bool = bool
if not hasattr(np, 'object'):
    np.object = object
if not hasattr(np, 'int'):
    np.int = int
if not hasattr(np, 'float'):
    np.float = float
if not hasattr(np, 'str'):
    np.str = str

import pandas as pd
import re
import requests
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
import time

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_FILES = {
    "s1": "study.txt",
    "s1s2": "study.txt"
}
OUTPUT_DIR = "/path/to/input_dir"
CACHE_FILE = "ec_to_metacyc_cache.json"
MAX_WORKERS = 8
RETRY_LIMIT = 3
RETRY_DELAY = 2  # seconds (for exponential backoff)
# -----------------------------

def extract_ko_ec(gene_str):
    match = re.match(r'/path/to/input_dir]', gene_str)
    if match:
        return match.group(1), match.group(2)
    return None, None

def load_data(files):
    dataframes = {}
    all_ec_set = set()
    for tag, path in files.items():
        df = pd.read_csv(path, sep='\t')
        first_col = df.columns[0]
        df = df.rename(columns={first_col: "Gene_Family"})
        df[['KO_ID', 'EC_ID']] = df['Gene_Family'].apply(lambda x: pd.Series(extract_ko_ec(x)))
        dataframes[tag] = df
        all_ec_set.update(df['EC_ID'].dropna().unique())
    return dataframes, sorted(all_ec_set)

def load_cache(cache_path):
    if os.path.exists(cache_path):
        with open(cache_path, 'r') as f:
            return json.load(f)
    return {}

def save_cache(cache, cache_path):
    with open(cache_path, 'w') as f:
        json.dump(cache, f, indent=2)

def get_metacyc_pathways(ec, max_retries=RETRY_LIMIT, delay=RETRY_DELAY):
    url = f"https://biocyc.org/META/NEW-IMAGE?object=EC-{ec}"
    for attempt in range(max_retries):
        try:
            res = requests.get(url, timeout=10)
            res.raise_for_status()
            soup = BeautifulSoup(res.content, "html.parser")
            pathways = sorted({a.text for a in soup.find_all("a") if "PWY" in a.text})
            return ec, ", ".join(pathways)
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(delay * (2 ** attempt))  # exponential backoff
            else:
                print(f"[ERROR] Failed EC {ec} after {max_retries} attempts: {e}")
                return ec, ""

def fetch_ecs_parallel(ec_list, existing_cache, max_workers=8):
    new_ecs = [ec for ec in ec_list if ec not in existing_cache]
    results = {}

    print(f"Querying {len(new_ecs)} EC numbers in parallel with {max_workers} threads...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(get_metacyc_pathways, ec): ec for ec in new_ecs}
        for future in as_completed(futures):
            ec, pathways = future.result()
            results[ec] = pathways

    return results

def apply_and_save(dataframes, ec_to_pathways, outdir):
    for tag, df in dataframes.items():
        df['MetaCyc_Pathways'] = df['EC_ID'].map(ec_to_pathways)
        ordered_cols = ['Gene_Family', 'KO_ID', 'EC_ID', 'MetaCyc_Pathways'] + [c for c in df.columns if c not in ['Gene_Family', 'KO_ID', 'EC_ID', 'MetaCyc_Pathways']]
        df = df[ordered_cols]
        output_path = os.path.join(outdir, f"study.tsv")
        df.to_csv(output_path, sep='\t', index=False)
        print(f"[DONE] Saved: {output_path}")

# -------- MAIN WORKFLOW --------
if __name__ == "__main__":
    dataframes, ec_list = load_data(INPUT_FILES)
    cache = load_cache(CACHE_FILE)
    new_results = fetch_ecs_parallel(ec_list, cache, max_workers=MAX_WORKERS)
    cache.update(new_results)
    save_cache(cache, CACHE_FILE)
    apply_and_save(dataframes, cache, OUTPUT_DIR)
