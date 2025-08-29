# ============================================================================
# Snakefile — Robust sample→file mapping + first real step (fastp on .gz)
# Strategies supported:
#   A) Manifest: samples.tsv (explicit paths; allows comma-separated lanes)
#   B) Convention: raw_dir + suffixes (_R1/_R2)
#   C) Auto-discovery: lane_glob (e.g., *_L001_R1.fastq.gz)
# The pipeline merges multiple lanes by gzip concatenation, then runs fastp.
# ============================================================================

import os, csv, glob
from pathlib import Path
configfile: "config/config.yaml"

# ---------- Load sample IDs ---------------------------------------------------
if config.get("use_manifest", False):
    # When using a manifest, we still define SAMPLES from that file to drive DAG
    with open(config["manifest"], newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        SAMPLES = [row["sample"].strip() for row in rdr if row.get("sample", "").strip()]
else:
    with open(config["sample_list"]) as fh:
        SAMPLES = [s.strip() for s in fh if s.strip()]

RESULTS = config["results_dir"]
WORK = config["work_dir"]
RAW_DIR = config.get("raw_dir", ".")
LAYOUT = config.get("layout", "paired").lower()
FASTP_THREADS = int(config.get("fastp_threads", 8))
FASTP_OPTS = config.get("fastp_opts", "")
FASTP_MOD = config.get("fastp_module", "")

# ---------- Manifest loader: supports comma-separated lanes -------------------
def load_manifest(path):
    d = {}
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            sid = row["sample"].strip()
            r1s = [p for p in row["r1_path"].split(",") if p.strip()] if row.get("r1_path") else []
            r2s = [p for p in row.get("r2_path", "").split(",") if p.strip()]
            d[sid] = {"r1": r1s, "r2": r2s}
    return d

MANIFEST = load_manifest(config["manifest"]) if config.get("use_manifest", False) else None

# ---------- Map sample ID → list of input files per read ----------------------
def fastq_lists(sample):
    """Return dict with keys 'r1' and 'r2' mapping to lists of .fastq.gz paths.
    Priority:
      1) Manifest entries (if enabled)
      2) Auto-discovered lanes (lane_glob)
      3) Convention (raw_dir + suffix)
    """
    if MANIFEST:
        return MANIFEST[sample]

    # Auto-discovery
    patt_r1 = os.path.join(RAW_DIR, config["lane_glob"].format(sample=sample, read=1))
    r1_lanes = sorted(glob.glob(patt_r1))
    r2_lanes = []
    if LAYOUT == "paired":
        patt_r2 = os.path.join(RAW_DIR, config["lane_glob"].format(sample=sample, read=2))
        r2_lanes = sorted(glob.glob(patt_r2))
    if r1_lanes:
        return {"r1": r1_lanes, "r2": r2_lanes}

    # Convention
    r1 = os.path.join(RAW_DIR, f"{sample}{config['r1_suffix']}")
    r2 = os.path.join(RAW_DIR, f"{sample}{config['r2_suffix']}") if LAYOUT == "paired" else None
    return {"r1": [r1], "r2": [r2] if r2 else []}

# ---------- Make merged (or linked) inputs for fastp --------------------------
# We always create a *single* merged .gz per read:
#   - If multiple lanes: cat them into merged.gz
#   - If a single file: symlink to merged.gz to avoid copying
# fastp consumes these merged files directly.

rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/fastp/report.json", sample=SAMPLES)

rule prepare_inputs:
    input:
        lambda wc: fastq_lists(wc.sample)
    output:
        r1=f"{WORK}/{{sample}}/raw/{{sample}}_R1.merged.fastq.gz",
        r2=f"{WORK}/{{sample}}/raw/{{sample}}_R2.merged.fastq.gz" if LAYOUT == "paired" else temp("/dev/null")
    threads: 1
    resources:
        mem_mb=2000, runtime=20, cpus=1
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.r1})
        # Merge/Link R1
        set -- {input[r1]}
        if [ "$#" -gt 1 ]; then
          cat "$@" > {output.r1}
        else
          ln -sf "$1" {output.r1}
        fi
        # Merge/Link R2 if paired
        if [ "''' + LAYOUT + '''" = "paired" ]; then
          set -- {input[r2]}
          if [ "$#" -gt 1 ]; then
            cat "$@" > {output.r2}
          elif [ "$#" -eq 1 ] && [ -n "$1" ]; then
            ln -sf "$1" {output.r2}
          fi
        fi
        '''

rule fastp_trim:
    input:
        r1=rules.prepare_inputs.output.r1,
        r2=rules.prepare_inputs.output.r2 if LAYOUT == "paired" else None
    output:
        r1=f"{RESULTS}/{{sample}}/fastp/{{sample}}_R1.trim.fastq.gz",
        r2=f"{RESULTS}/{{sample}}/fastp/{{sample}}_R2.trim.fastq.gz" if LAYOUT == "paired" else temp("/dev/null"),
        report_json=f"{RESULTS}/{{sample}}/fastp/report.json",
        report_html=f"{RESULTS}/{{sample}}/fastp/report.html"
    threads: FASTP_THREADS
    resources:
        mem_mb=8000, runtime=90, cpus=FASTP_THREADS
    shell:
        r'''
        set -euo pipefail
        # Load fastp via module if provided (else assume it is in PATH)
        if [ -n "''' + FASTP_MOD + '''" ]; then module purge || true; module load ''' + FASTP_MOD + ''' || true; fi
        mkdir -p $(dirname {output.r1})
        if [ "''' + LAYOUT + '''" = "paired" ] && [ -s "{input.r2}" ]; then
          fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
                -w {threads} {FASTP_OPTS} -h {output.report_html} -j {output.report_json}
        else
          fastp -i {input.r1} -o {output.r1} -w {threads} {FASTP_OPTS} \
                -h {output.report_html} -j {output.report_json}
        fi
        '''