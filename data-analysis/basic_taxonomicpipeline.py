#!/usr/bin/env python3
# =========================================================
# taxonomic_pipeline.py — Full Taxonomic Analysis Pipeline
# =========================================================
# Purpose
# -------
# Build Taxonomy/Genus tables from a Mothur-style .shared and .taxonomy,
# join metadata, compute relative-abundance tables (no rounding), and
# write consolidated mean-relative-abundance tables plus 100% stacked
# column charts to an Excel workbook. Optionally style charts via JSON
# (exported from another system) without changing the underlying values.
#
# Key behaviors
# -------------
# • No rounding of percentages. Values are written as raw floats; any
#   decimals you see are Excel *display* formatting, not data changes.
# • "Less abundant genera" is a row-wise sum of all genera beyond Top10
#   (Top10 determined by global totals across groups).
# • Keep all genus columns in the consolidated tables; charts only use
#   Top10 + "Less abundant genera".
# • Hierarchical consolidation supports ANY number of metadata variables.
# • Optional: read a chart-only JSON to control chart aesthetics.
#
# Outputs
# -------
# <output>.xlsx with sheets:
#   - Metadata             (as provided)
#   - Taxonomy             (OTU→ranks with sample counts)
#   - Genus                (Counts, Relative %, Transposed+Metadata)
#   - Charts               (Flat groupings + full hierarchical combo)
#   - All Consolidated     (optional; all size-2..m-1 hier combos)
# ======================================================================

from __future__ import annotations

# --- std deps ------------------------------------------------------------
import argparse
import json
from itertools import combinations

# --- data stack ----------------------------------------------------------
import numpy as np
np.bool = bool   # compatibility for older code that references np.bool
np.object = object
import pandas as pd

# =========================
# I/O and core transforms
# =========================

def read_shared(shared_file: str) -> pd.DataFrame:
    """Read a Mothur .shared file and return wide table with OTU rows
    and Sample columns (fill missing with 0). Drops 'label'/'numOtus'
    if present; pivots Group→columns.
    """
    df = pd.read_csv(shared_file, sep='	')
    df = df.drop(columns=['label', 'numOtus'], errors='ignore')
    # Long → wide (OTU rows, Sample columns)
    df_long = df.melt(id_vars='Group', var_name='OTU', value_name='Count')
    df_pivot = df_long.pivot_table(index='OTU', columns='Group', values='Count', fill_value=0)
    df_pivot = df_pivot.reset_index()  # col 0 = OTU
    return df_pivot


def read_taxonomy(taxonomy_file: str) -> pd.DataFrame:
    """Read taxonomy map and split the semicolon string into ranks.
    Expects columns: [OTU, Taxonomy]; trailing ';' is tolerated.
    """
    df = pd.read_csv(taxonomy_file, sep='	')
    taxonomy_levels = df['Taxonomy'].str.rstrip(';').str.split(';', expand=True)
    taxonomy_levels.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    df_final = pd.concat([df[['OTU']], taxonomy_levels], axis=1)
    return df_final


def build_taxonomy_sheet(shared_df: pd.DataFrame, taxonomy_df: pd.DataFrame) -> pd.DataFrame:
    """Merge counts with taxonomy labels; move Genus near front; rename
    'OTU'→'Group' per preference.
    """
    merged = pd.merge(shared_df, taxonomy_df, on='OTU', how='left')
    sample_cols = [c for c in merged.columns if c not in ['OTU', 'Genus', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']]
    final_cols = ['OTU', 'Genus'] + sample_cols + ['Kingdom', 'Phylum', 'Class', 'Order', 'Family']
    merged = merged[final_cols]
    merged = merged.rename(columns={'OTU': 'Group'})
    return merged


def build_genus_sheet(taxonomy_sheet: pd.DataFrame, metadata_df: pd.DataFrame | None = None):
    """Compute genus-level counts, relative %, and a transposed table.

    Returns
    -------
    abundance : DataFrame
        Genus × Sample (counts), sorted by total abundance (desc).
    rel : DataFrame
        Genus × Sample (%) — no rounding and no row-sum forcing.
    rel_T : DataFrame
        Sample × Genus (%) merged with metadata (if provided).
    genus_cols : list[str]
        Genus columns present in rel_T (to drive consolidation).
    meta_cols : list[str]
        Metadata columns found in rel_T (excludes 'Sample').
    """
    sample_cols = [c for c in taxonomy_sheet.columns if c not in ['Group', 'Genus', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']]

    # Counts by Genus
    genus_df = taxonomy_sheet.groupby('Genus')[sample_cols].sum().reset_index()

    # Sort genera by total counts (desc) BEFORE computing relative
    genus_df['__total__'] = genus_df[sample_cols].sum(axis=1)
    genus_df = genus_df.sort_values('__total__', ascending=False).drop(columns='__total__')

    # Relative abundance per sample (no rounding)
    col_sums = genus_df[sample_cols].sum(axis=0)
    rel = genus_df.copy()
    for s in sample_cols:
        rel[s] = (rel[s] / col_sums[s] * 100.0) if col_sums[s] > 0 else 0.0

    # Transpose for Table 3 and merge metadata
    rel_T = rel.set_index('Genus').transpose().reset_index().rename(columns={'index': 'Sample'})
    if metadata_df is not None:
        rel_T = rel_T.merge(metadata_df, on='Sample', how='left')
        meta_cols = [c for c in metadata_df.columns if c != 'Sample']
        genus_cols = [c for c in rel_T.columns if c not in ['Sample'] + meta_cols]
        rel_T = rel_T[['Sample'] + meta_cols + genus_cols]
    else:
        meta_cols = []
        genus_cols = [c for c in rel_T.columns if c != 'Sample']

    abundance = genus_df
    return abundance, rel, rel_T, genus_cols, meta_cols


def build_grouped_table_no_rounding(relative_df_T: pd.DataFrame, group_col: str, genus_cols: list[str]):
    """Flat consolidation by a single metadata column.

    - Group by `group_col` and compute mean of genus percentages.
    - Determine Top10 genera by totals across groups (for series ordering).
    - Insert "Less abundant genera" immediately after the 10th genus; value is
      the row-wise sum over all non-Top10 genera (by that global ordering).
    - **Do not** round or force sums to 100.
    - **Retain every genus column** (chart uses Top10+Less).
    """
    grouped = relative_df_T.groupby(group_col)[genus_cols].mean()
    totals = grouped.sum(axis=0).sort_values(ascending=False)
    top10 = list(totals.head(10).index)
    rest = [g for g in totals.index if g not in top10]

    full = pd.concat([grouped[top10], grouped[rest]], axis=1).reset_index()

    less_col = grouped[rest].sum(axis=1)
    insert_at = 1 + len(top10)  # +1 for group column at index 0
    full.insert(insert_at, 'Less abundant genera', less_col.values)

    return full, top10


def build_hierarchical_table_no_rounding(relative_df_T: pd.DataFrame, hier_cols: list[str], genus_cols: list[str]):
    """Hierarchical consolidation by composite of metadata columns.

    Accepts ANY number of metadata variables (len(hier_cols) >= 2) and ANY
    number of levels per variable; the GroupCombo key reflects all combos
    present in the data. Missing keys are filled with 'NA'.
    """
    df = relative_df_T.copy()
    missing = [c for c in hier_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing hierarchical metadata columns: {missing}")

    df[hier_cols] = df[hier_cols].fillna('NA').astype(str)
    df['GroupCombo'] = df[hier_cols].agg('_'.join, axis=1)

    grouped = df.groupby('GroupCombo')[genus_cols].mean()

    totals = grouped.sum(axis=0).sort_values(ascending=False)
    top10 = list(totals.head(10).index)
    rest = [g for g in totals.index if g not in top10]

    full = pd.concat([grouped[top10], grouped[rest]], axis=1).reset_index()

    less_col = grouped[rest].sum(axis=1)
    insert_at = 1 + len(top10)
    full.insert(insert_at, 'Less abundant genera', less_col.values)

    return full, top10


# =========================
# Excel writing (xlsxwriter)
# =========================

import numpy as np
import pandas as pd

def write_excel_xlsxwriter(
    metadata_df: pd.DataFrame | None,
    taxonomy_sheet: pd.DataFrame,
    genus_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str], list[str]],
    charts_specs: list[dict],
    output_file: str,
    all_combos_specs: list[dict] | None = None,
    chart_style_raw: dict | None = None,
) -> None:
    """
    Parameters
    ----------
    charts_specs : list of dict
        Each dict has keys:
          - title (str): grouping column name (e.g., 'Treatment' or 'GroupCombo')
          - df (DataFrame): consolidated table (all genera + Less abundant)
          - top10 (list[str]): Top10 genera names (series ordering)
          - display_title (optional str): pretty label in chart title (e.g., 'Treatment + Timepoint')

    all_combos_specs : list of dict, optional
        When provided, these blocks are written to a separate sheet
        named **'All Consolidated'**.

    chart_style_raw : dict, optional
        Chart-only JSON (as dict). We'll normalize it to what xlsxwriter
        expects and apply when building each chart.
    """

    # ---------- Helper: ensure Excel-safe values (Inf/NaN -> blank cell) ----------
    def make_excel_safe(df: pd.DataFrame) -> pd.DataFrame:
        # Convert ±inf to NaN, then NaN to None (xlsxwriter writes None as a blank cell)
        df = df.replace([np.inf, -np.inf], np.nan)
        return df.where(pd.notna(df), None)

    abundance, rel, rel_T, genus_cols, meta_cols = genus_tables

    # --- Normalize chart style (JSON -> internal dict) ------------------
    def rgb_to_hex(rgb):
        if isinstance(rgb, dict):
            r, g, b = rgb.get('r', 0), rgb.get('g', 0), rgb.get('b', 0)
        else:
            r, g, b = rgb
        return f"#{int(r):02X}{int(g):02X}{int(b):02X}"

    CHART_TYPE_MAP = {
        51: ('column', None),              # clustered column
        52: ('column', 'stacked'),
        54: ('column', 'percent_stacked'),
        53: ('column', 'percent_stacked'), # fallback
    }
    LEGEND_POS_MAP = {
        -4152: 'right',
        -4131: 'left',
        -4160: 'top',
        -4107: 'bottom',
        -4161: 'right',
    }

    def normalize_chart_style(raw):
        if not raw:
            return {}
        ctype = raw.get('ChartType')
        chart_type, chart_subtype = CHART_TYPE_MAP.get(ctype, ('column', 'percent_stacked'))
        has_title  = bool(raw.get('HasTitle', True))
        has_legend = bool(raw.get('HasLegend', True))
        legend_pos = LEGEND_POS_MAP.get(raw.get('LegendPosition'), 'right')
        axes = raw.get('Axes', {})
        val_ax = axes.get('Value', {}) if isinstance(axes, dict) else {}
        y_title = val_ax.get('TitleText', 'Mean Relative Abundance (%)')
        y_grid  = bool(val_ax.get('MajorGridlines', False))
        series_colors = {}
        for s in (raw.get('Series', []) or []):
            name = s.get('Name')
            frgb = s.get('FillRGB')
            if name and frgb:
                series_colors[name] = rgb_to_hex(frgb)
        return {
            'chart_type': chart_type,
            'chart_subtype': chart_subtype,
            'has_title': has_title,
            'has_legend': has_legend,
            'legend_position': legend_pos,
            'y_axis': {'name': y_title, 'major_gridlines': y_grid},
            'series_colors': series_colors,
        }

    style = normalize_chart_style(chart_style_raw or {})

    # ---------- SANITIZE dataframes up-front ----------
    abundance = make_excel_safe(abundance.copy())
    rel       = make_excel_safe(rel.copy())
    rel_T     = make_excel_safe(rel_T.copy())
    if metadata_df is not None:
        metadata_df = make_excel_safe(metadata_df.copy())
    taxonomy_sheet = make_excel_safe(taxonomy_sheet.copy())

    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        wb = writer.book
        header_fmt = wb.add_format({'bold': True})

        # Metadata
        if metadata_df is not None:
            metadata_df.to_excel(writer, sheet_name='Metadata', index=False)
            ws_meta = writer.sheets['Metadata']
            for col, name in enumerate(metadata_df.columns):
                ws_meta.write(0, col, name, header_fmt)

        # Taxonomy
        taxonomy_sheet.to_excel(writer, sheet_name='Taxonomy', index=False)
        ws_tax = writer.sheets['Taxonomy']
        for col, name in enumerate(taxonomy_sheet.columns):
            ws_tax.write(0, col, name, header_fmt)

        # Genus sheet (3 tables stacked)
        ws_gen = writer.book.add_worksheet('Genus')
        writer.sheets['Genus'] = ws_gen
        r = 0
        # Table 1: Abundance (counts)
        ws_gen.write_row(r, 0, abundance.columns, header_fmt); r += 1
        for _, row in abundance.iterrows():
            ws_gen.write_row(r, 0, row.tolist()); r += 1
        r += 5
        # Table 2: Relative (%) — values are NOT rounded
        ws_gen.write_row(r, 0, rel.columns, header_fmt); r += 1
        for _, row in rel.iterrows():
            ws_gen.write_row(r, 0, row.tolist()); r += 1
        r += 5
        # Table 3: Transposed + metadata
        ws_gen.write_row(r, 0, rel_T.columns, header_fmt); r += 1
        for _, row in rel_T.iterrows():
            ws_gen.write_row(r, 0, row.tolist()); r += 1

        # Charts sheet: write tables and add charts for DEFAULT specs
        ws_ch = writer.book.add_worksheet('Charts')
        writer.sheets['Charts'] = ws_ch
        row_cursor = 0

        def write_block(ws, start_row, spec, sheet_name_for_series):
            title = spec['title']
            # sanitize the block’s DF before writing
            df_safe = make_excel_safe(spec['df'].copy())
            display_title = spec.get('display_title', title)

            ws.write_row(start_row, 0, df_safe.columns, header_fmt)
            for ridx in range(len(df_safe)):
                ws.write_row(start_row + 1 + ridx, 0, df_safe.iloc[ridx].tolist())

            la_col_name = 'Less abundant genera'
            la_idx = df_safe.columns.get_loc(la_col_name)
            first_series_col = 1
            last_series_col = la_idx

            n_rows = len(df_safe)
            ctype = style.get('chart_type', 'column')
            csub = style.get('chart_subtype', 'percent_stacked')
            chart = wb.add_chart({'type': ctype, 'subtype': csub})

            if style.get('has_title', True):
                chart.set_title({'name': f"Mean Relative Abundance by {display_title}"})
            else:
                chart.set_title({'name': ''})

            if style.get('has_legend', True):
                chart.set_legend({'position': style.get('legend_position', 'right')})
            else:
                chart.set_legend({'none': True})

            y_cfg = style.get('y_axis', {})
            y_axis = {'name': y_cfg.get('name', 'Mean Relative Abundance (%)')}
            if 'major_gridlines' in y_cfg:
                y_axis['major_gridlines'] = {'visible': bool(y_cfg['major_gridlines'])}
            chart.set_y_axis(y_axis)
            chart.set_x_axis({'name': title})

            cats_first_row = start_row + 1
            cats_last_row = start_row + n_rows
            cats = [sheet_name_for_series, cats_first_row, 0, cats_last_row, 0]

            name_colors = style.get('series_colors', {})

            for c in range(first_series_col, last_series_col + 1):
                series_kwargs = {
                    'name':       [sheet_name_for_series, start_row, c],
                    'categories': cats,
                    'values':     [sheet_name_for_series, start_row + 1, c, start_row + n_rows, c],
                }
                header_name = df_safe.columns[c]
                hex_color = name_colors.get(header_name)
                if hex_color:
                    series_kwargs['fill'] = {'color': hex_color}
                    series_kwargs['border'] = {'color': hex_color}
                chart.add_series(series_kwargs)

            chart.set_size({'width': 600, 'height': 380})
            chart_anchor_row = start_row + n_rows + 2
            ws.insert_chart(chart_anchor_row, 0, chart)
            return chart_anchor_row + 20

        # Default blocks → 'Charts'
        for spec in charts_specs:
            row_cursor = write_block(ws_ch, row_cursor, spec, 'Charts')

        # Optional ALL-COMBOS → separate sheet
        if all_combos_specs:
            ws_all = writer.book.add_worksheet('All Consolidated')
            writer.sheets['All Consolidated'] = ws_all
            row_cursor_all = 0
            for spec in all_combos_specs:
                row_cursor_all = write_block(ws_all, row_cursor_all, spec, 'All Consolidated')



# =========================
# CLI
# =========================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Taxonomic Analysis: build Genus tables and Excel charts "
            "(xlsxwriter) — No rounding, auto-grouping, optional JSON styling"
        )
    )
    parser.add_argument('--shared', required=True, help='Path to .shared file (tab-delimited)')
    parser.add_argument('--taxonomy', required=True, help='Path to .taxonomy file (tab-delimited)')
    parser.add_argument('--metadata', required=True, help="Path to metadata (.txt, tab-delimited) with 'Sample' column")
    parser.add_argument('--output', default='Taxonomic_Analysis.xlsx', help='Output Excel file name')
    parser.add_argument('--all_hier_combos', action='store_true',
                        help=('Also build hierarchical consolidated tables for EVERY combination '
                              'of metadata columns (size >= 2) and write them to a separate sheet '
                              '(All Consolidated).'))
    parser.add_argument('--chart_style_json', help='Path to a JSON file with chart styling')

    args = parser.parse_args()

    # ---- Load optional chart style JSON ----
    chart_style_raw: dict = {}
    if args.chart_style_json:
        with open(args.chart_style_json, 'r') as fh:
            chart_style_raw = json.load(fh) or {}

    # ---- Load inputs ----
    shared_df = read_shared(args.shared)
    taxonomy_df = read_taxonomy(args.taxonomy)

    # Metadata is tab‑delimited; ensure 'Sample' column exists
    metadata_df = pd.read_csv(args.metadata, sep='	')
    if 'Sample' not in metadata_df.columns:
        raise ValueError("Metadata must contain a 'Sample' column (case-sensitive).")

    # ---- Build base tables ----
    taxonomy_sheet = build_taxonomy_sheet(shared_df, taxonomy_df)
    abundance, rel, rel_T, genus_cols, meta_cols = build_genus_sheet(taxonomy_sheet, metadata_df)

    # Discover metadata columns from rel_T (everything except Sample and genus cols)
    meta_cols = [c for c in rel_T.columns if c not in ['Sample'] + genus_cols]

    charts_specs: list[dict] = []            # DEFAULT sheet: Charts
    all_combos_specs: list[dict] = []        # OPTIONAL sheet: All Consolidated

    # 1) Flat grouping for EACH metadata variable automatically
    for col in meta_cols:
        tbl, top10 = build_grouped_table_no_rounding(rel_T, col, genus_cols)
        charts_specs.append({'title': col, 'df': tbl, 'top10': top10})

    # 2) Multi‑level consolidation (default: ALL metadata variables together)
    if len(meta_cols) >= 2:
        tbl_all, top10_all = build_hierarchical_table_no_rounding(rel_T, meta_cols, genus_cols)
        charts_specs.append({
            'title': 'GroupCombo',
            'df': tbl_all,
            'top10': top10_all,
            'display_title': ' + '.join(meta_cols),
        })

        # Optional: write ALL SIZE‑2..(m‑1) combinations to the separate sheet
        if args.all_hier_combos:
            m = len(meta_cols)
            for k in range(2, m):  # exclude size m (already included in default Charts sheet)
                for cols in combinations(meta_cols, k):
                    tbl_k, top10_k = build_hierarchical_table_no_rounding(rel_T, list(cols), genus_cols)
                    all_combos_specs.append({
                        'title': 'GroupCombo',
                        'df': tbl_k,
                        'top10': top10_k,
                        'display_title': ' + '.join(cols),
                    })

    # ---- Write Excel ----
    write_excel_xlsxwriter(
        metadata_df=metadata_df,
        taxonomy_sheet=taxonomy_sheet,
        genus_tables=(abundance, rel, rel_T, genus_cols, meta_cols),
        charts_specs=charts_specs,
        output_file=args.output,
        all_combos_specs=all_combos_specs if args.all_hier_combos else None,
        chart_style_raw=chart_style_raw,
    )

    print(f"✅ Pipeline complete. Output saved to: {args.output}")


if __name__ == '__main__':
    main()
