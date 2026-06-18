#!/usr/bin/env bash
# Thin wrapper around bin/preflight_gtdb.py for the case
# "I just want to know which GTDB release is installed and what's the latest".
# Doesn't download anything by itself — print-only.
#
# Usage:
#   bash bin/check_gtdb_release.sh
#   HIFITAX_DB=/path/to/db bash bin/check_gtdb_release.sh

set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$HERE/.." && pwd)"
DB="${HIFITAX_DB:-$PROJ/db}"

if ! command -v python3 >/dev/null; then
    echo "ERROR: python3 not on PATH" >&2
    exit 1
fi

# Run the existing preflight in dry-run mode (--check-only). If the user
# wants to actually trigger the download, they re-run with --assume-yes.
python3 "$PROJ/bin/preflight_gtdb.py" \
    --db-dir "$DB" \
    --blca-db "$DB/gtdb_ssu_BLCAparsed.fasta" \
    --blca-tax "$DB/gtdb_ssu_BLCAparsed.taxonomy" \
    --assume-no \
    --release latest
