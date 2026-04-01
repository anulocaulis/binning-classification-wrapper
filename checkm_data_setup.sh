#!/usr/bin/env bash
set -euo pipefail

# Download and install the CheckM data bundle.
# Usage:
#   ./checkm_data_setup.sh [DEST_DIR]
# Optional env var:
#   CHECKM_DB_URL=<url-to-checkm-data-tar.gz>

DEST_DIR="${1:-databases/checkm_data}"
CHECKM_DB_URL="${CHECKM_DB_URL:-https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz}"
DEST_PARENT="$(dirname "$DEST_DIR")"
DEST_BASENAME="$(basename "$DEST_DIR")"
TMP_ARCHIVE="${TMPDIR:-/tmp}/checkm_data_$$.tar.gz"

mkdir -p "$DEST_PARENT"

if [[ -f "$DEST_DIR/hmms/phylo.hmm" ]]; then
  echo "CheckM data already present at: $DEST_DIR"
  exit 0
fi

download_with() {
  local url="$1"
  local out="$2"

  if command -v curl >/dev/null 2>&1; then
    curl -fL "$url" -o "$out"
    return
  fi

  if command -v wget >/dev/null 2>&1; then
    wget -O "$out" "$url"
    return
  fi

  echo "ERROR: neither curl nor wget is available for download." >&2
  exit 1
}

echo "Downloading CheckM DB from: $CHECKM_DB_URL"
download_with "$CHECKM_DB_URL" "$TMP_ARCHIVE"

TOP_DIR="$(tar -tzf "$TMP_ARCHIVE" | head -n1 | cut -d/ -f1)"
tar -xzf "$TMP_ARCHIVE" -C "$DEST_PARENT"
rm -f "$TMP_ARCHIVE"

# Normalize extracted directory name to DEST_DIR.
if [[ "$TOP_DIR" != "$DEST_BASENAME" ]]; then
  if [[ -d "$DEST_DIR" ]]; then
    rm -rf "$DEST_DIR"
  fi
  mv "$DEST_PARENT/$TOP_DIR" "$DEST_DIR"
fi

if [[ ! -f "$DEST_DIR/hmms/phylo.hmm" ]]; then
  echo "ERROR: installation incomplete; expected file missing: $DEST_DIR/hmms/phylo.hmm" >&2
  exit 1
fi

echo "CheckM data installed at: $DEST_DIR"
echo "Set this when running CheckM-enabled tools:"
echo "  export CHECKM_DATA_PATH=$DEST_DIR"
