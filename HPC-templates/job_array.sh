#!/usr/bin/env bash
# ============================================================================
# job_array.sh — Per-sample SLURM array with robust mapping & fastp
# Submit with: sbatch --array=1-$(wc -l < samples.txt) job_array.sh samples.txt
# Inputs : samples.txt (one sample ID per line)
# Mapping: choose via MODE env var — convention | manifest | autodiscover
#          - convention: RAW_DIR/<ID>_R1.fastq.gz (+ _R2.fastq.gz)
#          - manifest  : samples.tsv with columns sample r1_path r2_path (comma‑sep lanes OK)
#          - autodiscover: RAW_DIR/<ID>_L*_R{1,2}.fastq.gz (merged by cat)
# First step: fastp reads .gz directly, outputs trimmed .gz + HTML/JSON reports
# ============================================================================

#SBATCH --job-name=reads_qc
#SBATCH --output=logs/reads_qc_%A_%a.out
#SBATCH --error=logs/reads_qc_%A_%a.err
#SBATCH --partition=your_partition
#SBATCH --account=your_account
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --array=1-1
#SBATCH --requeue

set -euo pipefail

# ---------- Inputs & environment ---------------------------------------------
SAMPLE_LIST=${1:?"Provide samples.txt"}
SAMPLE_LIST=$(readlink -f "$SAMPLE_LIST")
MODE=${MODE:-"convention"}         # convention | manifest | autodiscover
LAYOUT=${LAYOUT:-"paired"}         # paired | single
RAW_DIR=${RAW_DIR:-"/path/to/raw/data"}
MANIFEST=${MANIFEST:-"samples.tsv"}
FASTP_MOD=${FASTP_MOD:-"fastp/0.23.4"}
THREADS=${THREADS:-$SLURM_CPUS_PER_TASK}
RESULTS_BASE=${RESULTS_BASE:-"/home/gould209/..."}
WORK_BASE=${WORK_BASE:-"/scratch.global/gould209/..."}

mkdir -p logs "$RESULTS_BASE" "$WORK_BASE"

# ---------- Resolve the sample for this array index ---------------------------
NUM_SAMPLES=$(wc -l < "$SAMPLE_LIST")
: "${SLURM_ARRAY_TASK_ID:?Submit with --array=1-${NUM_SAMPLES}}"
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
[ -n "$SAMPLE_ID" ] || { echo "No sample for index $SLURM_ARRAY_TASK_ID" >&2; exit 3; }

echo "[array] Task $SLURM_ARRAY_TASK_ID/$NUM_SAMPLES  sample=$SAMPLE_ID  MODE=$MODE"

# ---------- Helpers to produce merged .gz inputs per read ---------------------
WORK_DIR="$WORK_BASE/$SAMPLE_ID/raw"; mkdir -p "$WORK_DIR"
MERGED_R1="$WORK_DIR/${SAMPLE_ID}_R1.merged.fastq.gz"
MERGED_R2="$WORK_DIR/${SAMPLE_ID}_R2.merged.fastq.gz"

link_or_cat() {
  # Usage: link_or_cat <merged_out.gz> <file1.gz> [file2.gz ...]
  local out=$1; shift
  if (( $# == 0 )); then return 0; fi
  if (( $# == 1 )); then ln -sf "$1" "$out"; else cat "$@" > "$out"; fi
}

get_paths_convention() {
  R1=("$RAW_DIR/${SAMPLE_ID}_R1.fastq.gz")
  if [[ "$LAYOUT" == "paired" ]]; then R2=("$RAW_DIR/${SAMPLE_ID}_R2.fastq.gz"); else R2=(); fi
}

get_paths_manifest() {
  # Expect TSV with fields: sample \t r1_path \t r2_path(optional, comma-separated)
  local line
  line=$(awk -v s="$SAMPLE_ID" 'BEGIN{FS="\t"} $1==s{print $0}' "$MANIFEST") || true
  [[ -n "$line" ]] || { echo "[manifest] Sample $SAMPLE_ID not found" >&2; exit 4; }
  local r1_col r2_col
  r1_col=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $2}')
  r2_col=$(echo "$line" | awk 'BEGIN{FS="\t"}{print $3}')
  IFS=',' read -r -a R1 <<< "$r1_col"
  if [[ "$LAYOUT" == "paired" && -n "$r2_col" ]]; then IFS=',' read -r -a R2 <<< "$r2_col"; else R2=(); fi
}

get_paths_autodiscover() {
  shopt -s nullglob
  R1=( "$RAW_DIR/${SAMPLE_ID}"_L*_R1.fastq.gz )
  if [[ "$LAYOUT" == "paired" ]]; then R2=( "$RAW_DIR/${SAMPLE_ID}"_L*_R2.fastq.gz ); else R2=(); fi
  (( ${#R1[@]} > 0 )) || { echo "[auto] No lane files for $SAMPLE_ID in $RAW_DIR" >&2; exit 5; }
}

case "$MODE" in
  convention)    get_paths_convention;;
  manifest)      get_paths_manifest;;
  autodiscover)  get_paths_autodiscover;;
  *) echo "Unknown MODE=$MODE" >&2; exit 6;;
esac

# Merge or link inputs to produce single .gz per read
link_or_cat "$MERGED_R1" "${R1[@]}"
if [[ "$LAYOUT" == "paired" ]]; then link_or_cat "$MERGED_R2" "${R2[@]}"; fi

echo "[inputs] R1: ${R1[*]}"
[[ "$LAYOUT" == "paired" ]] && echo "[inputs] R2: ${R2[*]}"

echo "[merged] R1→ $MERGED_R1"
[[ "$LAYOUT" == "paired" ]] && echo "[merged] R2→ $MERGED_R2"

# ---------- Run fastp (reads .gz directly) -----------------------------------
module purge || true; module load "$FASTP_MOD" || true
OUT_DIR="$RESULTS_BASE/$SAMPLE_ID/fastp"; mkdir -p "$OUT_DIR"

if [[ "$LAYOUT" == "paired" && -s "$MERGED_R2" ]]; then
  fastp -i "$MERGED_R1" -I "$MERGED_R2" \
        -o "$OUT_DIR/${SAMPLE_ID}_R1.trim.fastq.gz" \
        -O "$OUT_DIR/${SAMPLE_ID}_R2.trim.fastq.gz" \
        -w "$THREADS" -h "$OUT_DIR/report.html" -j "$OUT_DIR/report.json"
else
  fastp -i "$MERGED_R1" -o "$OUT_DIR/${SAMPLE_ID}_R1.trim.fastq.gz" \
        -w "$THREADS" -h "$OUT_DIR/report.html" -j "$OUT_DIR/report.json"
fi

echo "status=0" > "$OUT_DIR/.last_status"
echo "[done] $SAMPLE_ID"