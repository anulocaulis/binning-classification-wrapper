#!/usr/bin/env bash
set -euo pipefail

# Download and install the CheckM2 DIAMOND database via the workflow container.
# Usage:
#   ./checkm2_data_setup.sh [DEST_DIR]

DEST_DIR="${1:-databases/checkm2_db}"
CONTAINER="containers/qc_binning_annotation.sif"
CHECKM2_BIN="/opt/conda/envs/checkm2/bin/checkm2"

if [[ ! -f "$CONTAINER" ]]; then
  echo "ERROR: container not found: $CONTAINER" >&2
  exit 1
fi

mkdir -p "$DEST_DIR"

download_if_needed() {
  local existing_dmnd
  existing_dmnd="$(find "$DEST_DIR" -type f -name '*.dmnd' | head -n1 || true)"
  if [[ -n "$existing_dmnd" ]]; then
    echo "CheckM2 database already present at: $existing_dmnd"
    return
  fi

  singularity exec "$CONTAINER" "$CHECKM2_BIN" database --download --path "$DEST_DIR" --no_write_json_db
}

download_if_needed

DB_FILE="$(find "$DEST_DIR" -type f -name '*.dmnd' | head -n1 || true)"
if [[ -z "$DB_FILE" ]]; then
  echo "ERROR: CheckM2 database download did not produce a .dmnd file under $DEST_DIR" >&2
  exit 1
fi

echo "CheckM2 database installed at: $DB_FILE"
echo "Set this in config.yaml:"
echo "  checkm2_db: \"$DB_FILE\""