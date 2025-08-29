#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Purpose  : Create or activate an isolated software environment reproducibly.
#            This script is idempotent: run it multiple times safely.
# Usage    : source env_setup.sh       # so the environment persists in shell
# Notes    : Adjust module names/versions for your HPC. MSI examples included.
# -----------------------------------------------------------------------------

set -euo pipefail

# ===(1) Load site-provided compilers/interpreters via modules =================
# Modules guarantee consistent software stacks across login and compute nodes.
# At MSI, 'module purge' can sometimes remove needed baseline modules; prefer 'module reset' 
module reset

# Example module loads (customize to your site). Keep versions explicit for
# reproducibility; avoid floating versions where possible.
module load python3
module load conda
module load R/4.4.0-openblas-rocky8

# ===(2) Choose environment location ==========================================
# Use a persistent location in your home or group directory, not scratch,
# so the env survives scratch cleanups. Keep name pinned to pipeline version.
ENV_DIR="$HOME/.conda/envs/microbiome-pipeline-1.0"
ENV_NAME="microbiome-pipeline-1.0"

# ===(3) Create the environment if it doesn't exist ============================
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[env] Creating environment: $ENV_NAME"
  # Provide a lockfile or explicit versions to maximize reproducibility.
  # If you have 'environment.yml', prefer that; otherwise specify packages.
  conda create -y -p "$ENV_DIR" python=3.8.3 \
    snakemake=8.14 \
    pandas=2.2 \
    numpy=2.0 \
    pip

  # Example: pip installs pinned to a requirements file.
  # If present, uncomment the next two lines:
  # pip install -r requirements.txt
fi

# ===(4) Activate the environment =============================================
# 'source' ensures the activation modifies the current shell environment.
source activate "$ENV_DIR"

# ===(5) Record provenance =====================================================
# Helpful for debugging and paper methods sections.
python -V
snakemake --version || true
conda list | sed -n '1,40p' || true

echo "[env] Environment ready: $ENV_NAME"