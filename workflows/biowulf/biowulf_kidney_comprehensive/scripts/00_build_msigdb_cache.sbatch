#!/bin/bash
#SBATCH --job-name=passage_kid_msigdb
#SBATCH --cpus-per-task=1
#SBATCH --mem=12g
#SBATCH --time=01:30:00
#SBATCH --gres=lscratch:20
#SBATCH --output=logs/00_msigdb.%j.out
#SBATCH --error=logs/00_msigdb.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
COLLECTIONS="${2:-H,C2,C3,C5,C6,C7,C8}"
MIN_GENES="${3:-15}"
MAX_GENES="${4:-250}"
mkdir -p "$ROOT/logs" "$ROOT/refs/msigdb" "$ROOT/rlib"
mkdir -p "$ROOT/cache" "$ROOT/tmp"
export R_USER_CACHE_DIR="$ROOT/cache"
export XDG_CACHE_HOME="$ROOT/cache"
if [[ -n "${LSCRATCH:-}" && -d "$LSCRATCH" ]]; then
  export TMPDIR="$LSCRATCH"
else
  export TMPDIR="$ROOT/tmp"
fi
Rscript "$ROOT/scripts/00_build_msigdb_cache.R" "$ROOT" "$COLLECTIONS" "$MIN_GENES" "$MAX_GENES"
