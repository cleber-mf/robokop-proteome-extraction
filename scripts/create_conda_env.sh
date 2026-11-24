#!/usr/bin/env bash
# Create the self-contained Conda environment on Longleaf.
# Run once from the repository directory:
#   bash scripts/create_conda_env.sh

set -Eeuo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_PREFIX="${ROBO_CONDA_ENV:-${PROJECT_DIR}/.conda/robokop-proteome}"

# Start with an already available Conda installation.  If none is active,
# try the common Longleaf module names; module avail will show the current
# site-specific alternatives if both names change.
if ! command -v conda >/dev/null 2>&1; then
  module purge
  module load miniforge 2>/dev/null || module load anaconda 2>/dev/null
fi

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if [[ -d "${ENV_PREFIX}" ]]; then
  conda env update --prefix "${ENV_PREFIX}" --file "${PROJECT_DIR}/environment.yml" --prune
else
  conda env create --prefix "${ENV_PREFIX}" --file "${PROJECT_DIR}/environment.yml"
fi

conda run --prefix "${ENV_PREFIX}" python --version
conda run --prefix "${ENV_PREFIX}" python -m py_compile "${PROJECT_DIR}/robokop_extract_and_split.py"
printf '\nEnvironment ready: %s\n' "${ENV_PREFIX}"
