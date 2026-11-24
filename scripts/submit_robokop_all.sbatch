#!/usr/bin/env bash
# Submit with:
#   ROBO_RUN_DIR=/work/your_group/your_user/robokop_all \
#     sbatch scripts/submit_robokop_all.sbatch
#
# hov is listed first to prefer it when your account is allowed to use that
# partition.  Remove it if sbatch reports that your account lacks access.
# This is intentionally CPU-only: do not add --gres or a GPU partition.

#SBATCH --job-name=robokop-all
#SBATCH --partition=hov,general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --requeue
#SBATCH --output=slurm-robokop-%j.out
#SBATCH --error=slurm-robokop-%j.err

set -Eeuo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_PREFIX="${ROBO_CONDA_ENV:-${PROJECT_DIR}/.conda/robokop-proteome}"

# The all-protein export is very large.  Place output and SQLite state on a
# high-capacity project/work filesystem, not a home directory.  Explicitly set
# ROBO_RUN_DIR at submission time as shown above.
if [[ -z "${ROBO_RUN_DIR:-}" ]]; then
  printf '%s\n' "ERROR: Set ROBO_RUN_DIR to a high-capacity Longleaf work/project path." >&2
  exit 2
fi

mkdir -p "${ROBO_RUN_DIR}" "${ROBO_RUN_DIR}/tmp"
cd "${PROJECT_DIR}"

if ! command -v conda >/dev/null 2>&1; then
  module purge
  module load miniforge 2>/dev/null || module load anaconda 2>/dev/null
fi
CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_PREFIX}"

export PYTHONUNBUFFERED=1
export PYTHONHASHSEED=0
export TMPDIR="${ROBO_RUN_DIR}/tmp"
export SQLITE_TMPDIR="${TMPDIR}"

printf 'Started: %s\n' "$(date --iso-8601=seconds)"
printf 'Host: %s\n' "$(hostname)"
printf 'Allocated memory: %s\n' "${SLURM_MEM_PER_NODE:-all node memory}"
printf 'Allocated CPUs: %s\n' "${SLURM_CPUS_ON_NODE:-unknown}"

# The state DB is committed after each successful query batch.  If Longleaf
# requeues or the job reaches its time limit, submit the same command again:
# completed proteins will be skipped safely.
python "${PROJECT_DIR}/robokop_extract_and_split.py" \
  --all-proteins \
  --out-dir "${ROBO_RUN_DIR}" \
  --state-db "${ROBO_RUN_DIR}/.robokop_extract_state.sqlite3" \
  --timeout 300 \
  --batch-size 25 \
  --protein-page-size 500

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
