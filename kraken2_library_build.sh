#!/bin/bash
#SBATCH --job-name=kraken2_microbe_db
#SBATCH --output=logs/kraken2_microbe_db_%j.out
#SBATCH --error=logs/kraken2_microbe_db_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --partition=math-alderaan

set -euo pipefail

DB_DIR="/storage/biology/projects/miller-lowry/beitner/databases"
TARBALL="${DB_DIR}/k2_standard_20240112.tar.gz"
UNPACK_DIR="${DB_DIR}/k2_standard"

mkdir -p "${DB_DIR}"
mkdir -p "${UNPACK_DIR}"

if [ ! -f "${TARBALL}" ]; then
  echo "[$(date)] ERROR: tarball not found: ${TARBALL}" >&2
  echo "Run: wget -c 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz' -O ${TARBALL}"
  exit 1
fi

echo "[$(date)] Unpacking Kraken2 standard DB tarball..."
tar -xzf "${TARBALL}" -C "${UNPACK_DIR}"

# Validate
for f in hash.k2d taxo.k2d opts.k2d; do
  if [ ! -f "${UNPACK_DIR}/${f}" ]; then
    echo "[$(date)] ERROR: missing expected file: ${UNPACK_DIR}/${f}" >&2
    exit 1
  fi
done

echo "[$(date)] Kraken2 DB ready at: ${UNPACK_DIR}"
echo "Set config.yaml kraken2_db to: ${UNPACK_DIR}"

set -euo pipefail

REPO_DIR="/storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper"
DB_PARENT_DIR="${REPO_DIR}/databases"
DB_NAME="microbe_db"
DB_DIR="${DB_PARENT_DIR}/${DB_NAME}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
CONTAINER="${REPO_DIR}/containers/qc_binning_annotation.sif"

# Optional: point this to a prebuilt Kraken2 DB directory if available on shared storage.
# If set and valid, the script will copy/sync it and skip network downloads.
PREBUILT_DB_SOURCE=""

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

cd "${DB_PARENT_DIR}"

echo "[$(date)] Building Kraken2 DB at: ${DB_DIR}"
echo "[$(date)] Threads: ${THREADS}"

if [ -n "${PREBUILT_DB_SOURCE}" ]; then
  echo "[$(date)] PREBUILT_DB_SOURCE set: ${PREBUILT_DB_SOURCE}"
  if [ ! -d "${PREBUILT_DB_SOURCE}" ]; then
    echo "[$(date)] ERROR: PREBUILT_DB_SOURCE does not exist: ${PREBUILT_DB_SOURCE}" >&2
    exit 1
  fi
  mkdir -p "${DB_DIR}"
  rsync -a --delete "${PREBUILT_DB_SOURCE}/" "${DB_DIR}/"
  if [ ! -f "${DB_DIR}/hash.k2d" ] || [ ! -f "${DB_DIR}/taxo.k2d" ] || [ ! -f "${DB_DIR}/opts.k2d" ]; then
    echo "[$(date)] ERROR: Copied prebuilt DB missing one of hash.k2d/taxo.k2d/opts.k2d" >&2
    exit 1
  fi
  echo "[$(date)] Prebuilt Kraken2 DB staged successfully: ${DB_DIR}"
  exit 0
fi

# 1) Taxonomy
# NOTE: Some Kraken2 versions break against current NCBI FTP layout.
	kraken2-build --download-library bacteria --use-ftp --db "${DB_NAME}"
singularity exec -B "${REPO_DIR}:${REPO_DIR}" "${CONTAINER}" \
	kraken2-build --download-library archaea --use-ftp --db "${DB_NAME}"
singularity exec -B "${REPO_DIR}:${REPO_DIR}" "${CONTAINER}" \
	kraken2-build --download-library viral --use-ftp --db "${DB_NAME}"

# 3) Final build (hash.k2d)
singularity exec -B "${REPO_DIR}:${REPO_DIR}" "${CONTAINER}" \
	kraken2-build --build --threads "${THREADS}" --db "${DB_NAME}"

echo "[$(date)] Kraken2 DB build complete: ${DB_DIR}"
echo "Check files: ${DB_DIR}/hash.k2d ${DB_DIR}/taxo.k2d ${DB_DIR}/opts.k2d"
echo "Set config.yaml kraken2_db to: ${DB_DIR}"
