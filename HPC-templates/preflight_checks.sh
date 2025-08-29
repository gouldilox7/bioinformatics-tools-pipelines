#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Purpose  : Validate directories, quotas, and node visibility before running.
# Usage    : bash preflight_checks.sh
# Behavior : Exits non-zero on failure; prints actionable messages.
# -----------------------------------------------------------------------------

set -euo pipefail

# ===(1) User-configurable paths ===============================================
SCRATCH="/scratch.global/gould209"
PERM="/home/gould209"
RUN_DIR="$SCRATCH/..."
RESULTS_DIR="$PERM/..."

# ===(2) Basic path checks =====================================================
for d in "$SCRATCH" "$PERM"; do
  if [[ ! -d $d ]]; then
    echo "[preflight] ERROR: Missing directory: $d" >&2
    exit 2
  fi
  if [[ ! -w $d ]]; then
    echo "[preflight] ERROR: Directory not writable: $d" >&2
    exit 3
  fi
  echo "[preflight] OK path: $d"
done

mkdir -p "$RUN_DIR" "$RESULTS_DIR"

echo "[preflight] Run dir    : $RUN_DIR"
echo "[preflight] Results dir: $RESULTS_DIR"

# ===(3) Quota sanity (site-specific) ==========================================
# If your site exposes 'quota' or 'lfs quota', add checks so jobs don’t die
# mid-run due to ENOSPC. These commands vary by system—examples shown only.
if command -v lfs >/dev/null 2>&1; then
  # Lustre example; harmless if not available.
  lfs df -h "$SCRATCH" || true
fi

if command -v quota >/dev/null 2>&1; then
  quota -s || true
fi

# ===(4) Node visibility / scheduler reachability ==============================
# Quick sanity: can we see SLURM and submit? Avoid immediate failures later.
if ! command -v sbatch >/dev/null 2>&1; then
  echo "[preflight] ERROR: sbatch not found. Load your slurm module." >&2
  exit 4
fi

sinfo || true

echo "[preflight] All checks passed."