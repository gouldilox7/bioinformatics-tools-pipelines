#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Purpose  : Mirror key outputs from scratch to permanent storage efficiently.
# Usage    : bash rsync.sh
# -----------------------------------------------------------------------------

set -euo pipefail

SRC="/scratch.global/gould209/..."
DST="/home/gould209/.../archive_$(date +%Y%m%d_%H%M)"

# Create a time-stamped destination so we can roll back if needed.
mkdir -p "$DST"

# -a: archive (preserves perms/times)
# -v: verbose
# -h: human readable
# --partial: keep partially transferred files (useful on interruption)
# --info=progress2: nice progress output
rsync -avh --partial --info=progress2 "$SRC/" "$DST/"

echo "[sync] Copied $(du -sh "$DST" | cut -f1) to $DST"