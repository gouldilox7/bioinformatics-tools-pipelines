#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, wilcoxon, kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from skbio import DistanceMatrix
from skbio.diversity import beta_diversity
from skbio.stats.distance import anosim, permanova
from scipy.spatial.distance import pdist, squareform
import os
import logging


def _normalize_index(idx):
    return pd.Index([str(x).strip() for x in idx])


logging.basicConfig(
    filename='diversity_analysis.log',
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    force=True
)


def parse_args():
    parser = argparse.ArgumentParser(description="Perform alpha & beta diversity statistical analysis.")
    parser.add_argument('--alpha', required=True, help='Alpha diversity table (samples as rows)')
    parser.add_argument('--beta', required=True, help='Distance matrix, mothur .dist, or feature table (shared file allowed)')
    parser.add_argument('--metadata', required=True, help='Metadata file study.tsv or .txt)')
    parser.add_argument('--out', required=True, help='Output Excel file name')
    parser.add_argument('--out_matrix', default='study.tsv', help='Optional: Output Bray-Curtis distance matrix if generated')
    parser.add_argument('--nperm', type=int, default=999, help='Number of permutations for ANOSIM/PERMANOVA [default: 999]')
    return parser.parse_args()


def read_tables(alpha_fp, beta_fp, metadata_fp):
    alpha = pd.read_csv(alpha_fp, sep='\t', index_col=0)
    metadata = pd.read_csv(metadata_fp, sep='\t', index_col=0)

    # Normalize indices (strip whitespace, ensure strings)
    alpha.index = _normalize_index(alpha.index)
    metadata.index = _normalize_index(metadata.index)

    # If beta is a mothur lower-triangle file, pass through as filename
    if isinstance(beta_fp, str) and beta_fp.endswith('.lt.dist'):
        logging.info(f"Beta input detected as mothur lower-triangle .lt.dist: {beta_fp}")
        return alpha, beta_fp, metadata

    try:
        beta = pd.read_csv(beta_fp, sep='\t', header=0, index_col=0)
        # Normalize DF beta labels too
        beta.index = _normalize_index(beta.index)
        beta.columns = _normalize_index(beta.columns)
        logging.info("Loaded beta table as DataFrame")
        return alpha, beta, metadata
    except Exception:
        logging.info(f"Loaded beta input as filename (likely .dist): {beta_fp}")
        return alpha, beta_fp, metadata


def merge_alpha_metadata(alpha, metadata):
    shared_samples = alpha.index.intersection(metadata.index)
    alpha = alpha.loc[shared_samples]
    metadata = metadata.loc[shared_samples]
    merged = metadata.merge(alpha, left_index=True, right_index=True)
    merged.insert(0, 'Sample', merged.index)
    logging.info(f"Merged alpha and metadata on {len(shared_samples)} shared samples")
    return merged.reset_index(drop=True)


def test_alpha_diversity(merged_df, metadata_cols, alpha_cols):
    results = []
    for var in metadata_cols:
        groups = merged_df[var].dropna().unique()
        n_groups = len(groups)
        for metric in alpha_cols:
            data = merged_df[[var, metric]].dropna()
            if n_groups == 2:
                g1 = data[data[var] == groups[0]][metric]
                g2 = data[data[var] == groups[1]][metric]
                logging.info(f"Alpha test for {var}-{metric}: group sizes = {len(g1)}, {len(g2)}")
                if len(g1) == 0 or len(g2) == 0:
                    logging.warning(f"Skipping {var}-{metric}: One of the groups is empty.")
                    continue
                # t-test
                try:
                    t_stat, t_p = ttest_ind(g1, g2)
                except Exception as e:
                    logging.warning(f"t-test failed for {var}-{metric}: {e}")
                    t_stat, t_p = np.nan, np.nan
                results.append([var, metric, 't-test', t_stat, t_p])
                # Non-parametric: Wilcoxon if paired (equal sizes), else Mann-Whitney U
                try:
                    if len(g1) == len(g2):
                        w_stat, w_p = wilcoxon(g1, g2)
                        results.append([var, metric, 'Wilcoxon', w_stat, w_p])
                    else:
                        mw_stat, mw_p = mannwhitneyu(g1, g2, alternative='two-sided')
                        results.append([var, metric, 'Mann-Whitney U', mw_stat, mw_p])
                except Exception as e:
                    logging.warning(f"Non-parametric test failed for {var}-{metric}: {e}")
                    results.append([var, metric, 'Wilcoxon/MWU', np.nan, np.nan])
            elif n_groups > 2:
                try:
                    grouped_data = [data[data[var] == g][metric] for g in groups]
                    k_stat, k_p = kruskal(*grouped_data)
                    results.append([var, metric, 'Kruskal-Wallis', k_stat, k_p])
                except Exception as e:
                    logging.warning(f"Alpha Kruskal-Wallis failed for {var}-{metric}: {e}")
                    results.append([var, metric, 'Kruskal-Wallis', np.nan, np.nan])
    df = pd.DataFrame(results, columns=['Variable', 'Metric', 'Test', 'Statistic', 'P-value'])
    # FDR correction across all alpha tests (ignore NaNs)
    if not df.empty:
        pvals = df['P-value'].astype(float)
        valid = pvals.notna()
        if valid.sum() > 0:
            _, adj_p, _, _ = multipletests(pvals[valid].values, method='fdr_bh')
            df['Adj_P-value (FDR)'] = np.nan
            df.loc[valid, 'Adj_P-value (FDR)'] = adj_p
        else:
            df['Adj_P-value (FDR)'] = np.nan
    return df


def preprocess_beta_table(beta):
    # Handle mothur distance files by extension
    if isinstance(beta, str):
        if beta.endswith('.lt.dist'):
            return parse_mothur_lt_dist(beta)
        if beta.endswith('.dist'):
            return parse_mothur_dist(beta)

    # Handle DataFrame feature tables or square matrices
    if isinstance(beta, pd.DataFrame):
        if {'label', 'Group', 'numOtus'}.issubset(beta.columns):
            beta = beta.drop(columns=['label', 'numOtus'])
            beta = beta.set_index('Group')
        non_numeric_cols = beta.columns[~beta.dtypes.apply(np.issubdtype, args=(np.number,))].tolist()
        if non_numeric_cols:
            logging.warning(f"Dropping non-numeric columns from beta table: {non_numeric_cols}")
            beta = beta.drop(columns=non_numeric_cols)
    return beta


def parse_mothur_dist(dist_fp):
    """Parse mothur pairwise (3-column) .dist files: s1, s2, value per line."""
    dist_dict = {}
    all_ids = set()
    with open(dist_fp, 'r') as f:
        for line in f:
            parts = [p.strip() for p in line.strip().split('\t')]
            if len(parts) != 3:
                continue
            s1, s2, dist = parts
            try:
                dist = float(dist)
            except ValueError:
                continue
            s1 = s1.strip(); s2 = s2.strip()
            all_ids.update([s1, s2])
            dist_dict[(s1, s2)] = dist
            dist_dict[(s2, s1)] = dist
            dist_dict[(s1, s1)] = 0.0
            dist_dict[(s2, s2)] = 0.0
    all_ids = sorted(list(all_ids))
    dist_matrix = np.zeros((len(all_ids), len(all_ids)))
    for i, id1 in enumerate(all_ids):
        for j, id2 in enumerate(all_ids):
            dist_matrix[i, j] = dist_dict.get((id1, id2), 0.0)
    logging.info(f"Parsed Mothur .dist (3-column) file with {len(all_ids)} samples")
    return DistanceMatrix(dist_matrix, ids=all_ids)


def parse_mothur_lt_dist(dist_fp):
    """Parse mothur lower-triangle .lt.dist files.
    Format:
      First non-empty line may be an integer (number of samples). Following lines:
      <sample_id> [d( current, id1 )] [d( current, id2 )] ... in lower-triangle order.
    """
    ids = []
    rows = []
    with open(dist_fp, 'r') as f:
        lines = [ln.rstrip('\n') for ln in f if ln.strip() != '']
    # If first token is an integer, drop it
    first_tokens = [t.strip() for t in lines[0].split('\t')] if lines else []
    if first_tokens and len(first_tokens) == 1:
        try:
            _ = int(first_tokens[0])
            lines = lines[1:]
        except ValueError:
            pass
    for ln in lines:
        parts = [p.strip() for p in ln.split('\t')]
        if not parts:
            continue
        sid = parts[0].strip()
        vals = [p for p in parts[1:] if p != '']
        try:
            dists = [float(x) for x in vals]
        except ValueError:
            continue
        ids.append(sid)
        rows.append(dists)
    n = len(ids)
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j, val in enumerate(rows[i][:i]):
            M[i, j] = val
            M[j, i] = val
    logging.info(f"Parsed Mothur .lt.dist (lower triangle) with {n} samples")
    return DistanceMatrix(M, ids=ids)


def ensure_distance_matrix(beta):
    beta = preprocess_beta_table(beta)
    if isinstance(beta, DistanceMatrix):
        return beta, None
    if isinstance(beta, pd.DataFrame):
        if beta.shape[0] == beta.shape[1] and (beta.columns == beta.index).all():
            return DistanceMatrix(beta.values, ids=beta.index.tolist()), None
        else:
            bc = beta_diversity('braycurtis', beta.values, ids=beta.index.tolist())
            return bc, pd.DataFrame(bc.data, index=bc.ids, columns=bc.ids)
    raise ValueError("Unable to process beta input into a distance matrix.")


# ---- New helper to compute PERMANOVA ANOVA-like fields (Df, SumOfSqs, R2, F) ----
def _permanova_effects(dm: DistanceMatrix, groups: pd.Series):
    """Compute Anderson (2001) one-factor PERMANOVA components to mirror vegan::adonis2.

    Returns dict with keys: 'Df', 'Df_within', 'SumOfSqs_among', 'SumOfSqs_within',
    'SumOfSqs_total', 'R2', 'F'.
    """
    ids = [i for i in dm.ids if i in groups.index]
    groups = groups.loc[ids]

    D = pd.DataFrame(dm.filter(ids, strict=False).data, index=ids, columns=ids)
    D2 = D.values ** 2

    N = len(ids)
    k = int(groups.nunique())

    # total SS_T = (1/N) * sum_{i<j} d_ij^2
    tri = np.triu_indices(N, 1)
    ss_total = D2[tri].sum() / N if N > 0 else np.nan

    # within-group SS_W = sum_g (1/n_g) * sum_{i<j in g} d_ij^2
    ss_within = 0.0
    id_pos = {sid: pos for pos, sid in enumerate(ids)}
    for _, idx_vals in groups.groupby(groups).groups.items():
        idx_list = list(idx_vals)
        n_g = len(idx_list)
        if n_g < 2:
            continue
        pos = [id_pos[sid] for sid in idx_list]
        sub = D2[np.ix_(pos, pos)]
        tri_g = np.triu_indices(n_g, 1)
        ss_within += sub[tri_g].sum() / n_g

    ss_among = ss_total - ss_within

    df_among = k - 1
    df_within = N - k

    ms_among = ss_among / df_among if df_among > 0 else np.nan
    ms_within = ss_within / df_within if df_within > 0 else np.nan

    F = (ms_among / ms_within) if (ms_among == ms_among and ms_within > 0) else np.nan
    R2 = (ss_among / ss_total) if (ss_total and ss_total > 0) else np.nan

    return {
        'Df': int(df_among) if pd.notna(df_among) else np.nan,
        'Df_within': int(df_within) if pd.notna(df_within) else np.nan,
        'SumOfSqs_among': float(ss_among) if pd.notna(ss_among) else np.nan,
        'SumOfSqs_within': float(ss_within) if pd.notna(ss_within) else np.nan,
        'SumOfSqs_total': float(ss_total) if pd.notna(ss_total) else np.nan,
        'R2': float(R2) if pd.notna(R2) else np.nan,
        'F': float(F) if pd.notna(F) else np.nan,
    }


def test_beta_diversity(dm, metadata, nperm=999):
    """Run ANOSIM and PERMANOVA and return an R-like table with:
    Variable, Test, R_statistic, P_value, Df, SumOfSqs, R2, F, Permutations
    """
    rows = []
    metadata = metadata.copy()
    metadata.index = _normalize_index(metadata.index)

    ids = [i for i in dm.ids if i in metadata.index]
    missing_in_meta = [i for i in dm.ids if i not in metadata.index]
    logging.info(f"Beta diversity: {len(ids)} shared samples between distance matrix and metadata")
    if missing_in_meta:
        logging.warning(f"Samples in distance matrix but missing in metadata (showing up to 10): {missing_in_meta[:10]}")
    if len(ids) < 2:
        return pd.DataFrame(columns=['Variable', 'Test', 'R_statistic', 'P_value', 'Df', 'SumOfSqs', 'R2', 'F', 'Permutations'])

    sub_dm = dm.filter(ids, strict=False)

    for var in metadata.columns:
        vec = metadata.loc[ids, var]
        vec = vec.dropna()
        if vec.nunique() < 2 or len(vec) < 2:
            logging.warning(f"Skipping {var}: <2 groups or <2 samples after NA drop (n={len(vec)})")
            continue

        use_ids = vec.index.tolist()
        vdm = sub_dm.filter(use_ids, strict=False)

        # ANOSIM
        try:
            a = anosim(vdm, vec, permutations=nperm)
            # skbio may return an object with attrs or a Series-like
            a_stat = getattr(a, 'statistic', np.nan)
            a_p = getattr(a, 'p_value', np.nan)
            a_perm = getattr(a, 'permutations', None)
            # Series-like fallback
            if isinstance(a, pd.Series):
                a_stat = a.get('R', a.get('test statistic', a_stat))
                a_p = a.get('p-value', a.get('p value', a_p))
                a_perm = a.get('number of permutations', a_perm)

            rows.append({
                'Variable': var,
                'Test': 'ANOSIM',
                'R_statistic': float(a_stat) if pd.notna(a_stat) else np.nan,
                'P_value': float(a_p) if pd.notna(a_p) else np.nan,
                'Df': np.nan,
                'SumOfSqs': np.nan,
                'R2': np.nan,
                'F': np.nan,
                'Permutations': int(a_perm) if a_perm is not None else int(nperm)
            })
        except Exception as e:
            logging.warning(f"ANOSIM failed for {var}: {e}")
            rows.append({'Variable': var, 'Test': 'ANOSIM', 'R_statistic': np.nan, 'P_value': np.nan,
                         'Df': np.nan, 'SumOfSqs': np.nan, 'R2': np.nan, 'F': np.nan, 'Permutations': int(nperm)})

        # PERMANOVA
        try:
            p = permanova(vdm, vec, permutations=nperm)
            p_stat = getattr(p, 'statistic', np.nan)
            p_p = getattr(p, 'p_value', np.nan)
            p_perm = getattr(p, 'permutations', None)
            if isinstance(p, pd.Series):
                p_stat = p.get('F', p.get('test statistic', p_stat))
                p_p = p.get('p-value', p.get('p value', p_p))
                p_perm = p.get('number of permutations', p_perm)

            eff = _permanova_effects(vdm, vec)

            rows.append({
                'Variable': var,
                'Test': 'PERMANOVA',
                'R_statistic': np.nan,
                'P_value': float(p_p) if pd.notna(p_p) else np.nan,
                'Df': eff['Df'],
                'SumOfSqs': eff['SumOfSqs_among'],
                'R2': eff['R2'],
                'F': float(p_stat) if pd.notna(p_stat) else eff['F'],
                'Permutations': int(p_perm) if p_perm is not None else int(nperm)
            })
        except Exception as e:
            logging.warning(f"PERMANOVA failed for {var}: {e}")
            rows.append({'Variable': var, 'Test': 'PERMANOVA', 'R_statistic': np.nan, 'P_value': np.nan,
                         'Df': np.nan, 'SumOfSqs': np.nan, 'R2': np.nan, 'F': np.nan, 'Permutations': int(nperm)})

    out = pd.DataFrame(rows, columns=[
        'Variable', 'Test', 'R_statistic', 'P_value', 'Df', 'SumOfSqs', 'R2', 'F', 'Permutations'
    ])
    return out


def write_excel(alpha_merged, alpha_results, beta_results, output_file):
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        alpha_merged.to_excel(writer, sheet_name='Alpha Diversity', index=False)
        start_col = alpha_merged.shape[1] + 5
        alpha_results.to_excel(writer, sheet_name='Alpha Diversity', index=False, startcol=start_col)
        beta_results.to_excel(writer, sheet_name='Beta Diversity', index=False)
    logging.info(f"Excel workbook written to {output_file}")


def main():
    args = parse_args()
    alpha, beta, metadata = read_tables(args.alpha, args.beta, args.metadata)

    alpha_merged = merge_alpha_metadata(alpha, metadata)
    metadata_cols = metadata.columns.tolist()
    alpha_cols = alpha.columns.tolist()
    alpha_results = test_alpha_diversity(alpha_merged, metadata_cols, alpha_cols)

    dist_matrix, bc_df = ensure_distance_matrix(beta)
    beta_results = test_beta_diversity(dist_matrix, metadata, nperm=args.nperm)

    write_excel(alpha_merged, alpha_results, beta_results, args.out)

    if bc_df is not None:
        bc_df.to_csv(args.out_matrix, sep='\t')
        logging.info(f"Bray-Curtis distance matrix saved to {args.out_matrix}")


if __name__ == '__main__':
    main()
