#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Purpose  : Glue everything together. Run checks, then launch Snakemake with
#            the SLURM profile, using config from config.yaml.
# Usage    : bash 07_submit_workflow.sh
# -----------------------------------------------------------------------------

set -euo pipefail

# (1) Prepare environment
source env_setup.sh

# (2) Preflight checks
bash preflight_checks.sh

# (3) Launch Snakemake DAG with SLURM executor and profile
snakemake \
  -j 100 \
  --profile profiles/snakemake_slurm_profile.yaml \
  --configfile config/config.yaml \
  --rerun-triggers mtime \
  --restart-times 2 \
  --keep-going \
  --scheduler greedy

# Notes:
# - --rerun-triggers mtime: if inputs changed, recompute targets.
# - --restart-times: auto-retry transient failures per-rule.
# - --scheduler greedy: starts runnable jobs ASAP to improve throughput.