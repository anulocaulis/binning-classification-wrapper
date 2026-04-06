#!/bin/bash
#SBATCH --job-name=gtdbtk_db_setup
#SBATCH --output=logs/gtdbtk_db_setup_%j.out
#SBATCH --error=logs/gtdbtk_db_setup_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --partition=math-alderaan

set -euo pipefail

REPO_DIR="/storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper"
DB_PARENT_DIR="/storage/biology/projects/miller-lowry/beitner/databases"
TARBALL="${DB_PARENT_DIR}/gtdbtk_data.tar.gz"
DB_DIR="${DB_PARENT_DIR}/gtdbtk_db"
TMP_DIR="${DB_PARENT_DIR}/.gtdbtk_extract_${SLURM_JOB_ID:-manual}"
CONTAINER="${REPO_DIR}/containers/classification.sif"
URL="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz"

mkdir -p "${REPO_DIR}/logs"
mkdir -p "${DB_PARENT_DIR}"

if ! command -v singularity >/dev/null 2>&1; then
  echo "[$(date)] ERROR: singularity not found on PATH" >&2
  exit 1
fi

if [ ! -f "${CONTAINER}" ]; then
  echo "[$(date)] ERROR: container not found: ${CONTAINER}" >&2
  exit 1
fi

if [ ! -f "${TARBALL}" ]; then
  echo "[$(date)] GTDB-Tk tarball not found; downloading..."
  wget -c "${URL}" -O "${TARBALL}"
else
  echo "[$(date)] Found existing tarball: ${TARBALL}"
fi

echo "[$(date)] Preparing extraction directory: ${TMP_DIR}"
rm -rf "${TMP_DIR}" "${DB_DIR}.new"
mkdir -p "${TMP_DIR}" "${DB_DIR}.new"

echo "[$(date)] Extracting GTDB-Tk data with singularity..."
singularity exec -B "${DB_PARENT_DIR}:${DB_PARENT_DIR}" "${CONTAINER}" bash -lc "set -euo pipefail; tar -xzf '${TARBALL}' -C '${TMP_DIR}'"

echo "[$(date)] Normalizing extracted directory layout..."
shopt -s nullglob dotglob
entries=("${TMP_DIR}"/*)

if [ ${#entries[@]} -eq 1 ] && [ -d "${entries[0]}" ]; then
  mv "${entries[0]}"/* "${DB_DIR}.new/"
else
  mv "${TMP_DIR}"/* "${DB_DIR}.new/"
fi

if [ ! -d "${DB_DIR}.new/metadata" ] || [ ! -d "${DB_DIR}.new/taxonomy" ] || [ ! -d "${DB_DIR}.new/masks" ]; then
  echo "[$(date)] ERROR: Extracted DB missing expected directories (metadata/taxonomy/masks)." >&2
  exit 1
fi

if [ -d "${DB_DIR}" ]; then
  echo "[$(date)] Removing old DB directory: ${DB_DIR}"
  rm -rf "${DB_DIR}"
fi
mv "${DB_DIR}.new" "${DB_DIR}"
rm -rf "${TMP_DIR}"

echo "[$(date)] Verifying GTDB-Tk installation against new DB..."
singularity exec -B "${DB_PARENT_DIR}:${DB_PARENT_DIR}" "${CONTAINER}" bash -lc "set -euo pipefail; export GTDBTK_DATA_PATH='${DB_DIR}'; gtdbtk check_install"

echo "[$(date)] GTDB-Tk database ready: ${DB_DIR}"
echo "Set config.yaml gtdbtk_db to: ${DB_DIR}"
