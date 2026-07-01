#!/usr/bin/env bash
# Build a BLCA-formatted UNITE reference for fungal ITS, the fungal analogue of
# HiFiTaxa's bin/build_gtdb_blca_db.sh. Use this when params.marker == 'ITS'.
#
#   usage: build_unite_blca_db.sh <unite_fasta> <db-dir> [min-len]
#   e.g.   bash bin/build_unite_blca_db.sh \
#              ~/Downloads/sh_general_release_dynamic_19.02.2025.fasta \
#              db_unite 0
#
#   <db-dir> is relative to the HiFiTaxa projectDir; the pipeline expects the
#   outputs under ${projectDir}/db_unite (params.unite_db_dir).
#
#   min-len: minimum sequence length to keep. DEFAULT 0 (keep all). ITS is short
#            (UNITE median ~526 bp); do NOT reuse the 1000 bp 16S threshold here.
#
# Produces in <db-dir>:
#   unite_BLCAparsed.fasta        bare-id headers (+ BLAST index .n*)
#   unite_BLCAparsed.taxonomy     id <TAB> superkingdom:Fungi;...;species:Genus epithet;
#   placeholder_species_ids.txt   records whose species is a UNITE placeholder
#   UNITE_PARSE_STATS.json        parse diagnostics
#   UNITE_VERSION.txt             release tag
#
# When marker==ITS the pipeline resolves blca_db -> <db-dir>/unite_BLCAparsed.fasta
# and blca_tax -> <db-dir>/unite_BLCAparsed.taxonomy (see nextflow.config).
#
# Needs: python3, makeblastdb (BLAST+) on PATH.
set -euo pipefail

UNITE_FASTA="${1:?UNITE FASTA required}"
DBDIR="${2:?db dir required}"
MINLEN="${3:-0}"
VERSION_TAG="${UNITE_VERSION:-sh_general_release_dynamic_19.02.2025}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTA="$DBDIR/unite_BLCAparsed.fasta"

command -v python3     >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
command -v makeblastdb >/dev/null || { echo "ERROR: makeblastdb (BLAST+) not on PATH" >&2; exit 1; }
[ -f "$UNITE_FASTA" ] || { echo "ERROR: $UNITE_FASTA not found" >&2; exit 1; }

mkdir -p "$DBDIR"

echo "[build] 1/2 parsing UNITE -> BLCA FASTA + taxonomy (min-length=${MINLEN})"
python3 "$HERE/parse_unite_to_blca.py" \
    --in-fasta   "$UNITE_FASTA" \
    --out-dir    "$DBDIR" \
    --min-length "$MINLEN"

echo "[build] 2/2 makeblastdb (-parse_seqids) on the parsed reference"
# -parse_seqids is REQUIRED for BLCA: it lets blastdbcmd retrieve hit sequences
# for the per-query alignment step (same requirement as the GTDB build).
rm -f "$FASTA".n*
makeblastdb -in "$FASTA" -dbtype nucl -parse_seqids -out "$FASTA" 2>&1 | tail -5

echo "$VERSION_TAG" > "$DBDIR/UNITE_VERSION.txt"
echo "[build] done: $(grep -c '^>' "$FASTA") sequences, UNITE ${VERSION_TAG}"
echo "[build] BLAST index:"; ls -lh "$FASTA".n* 2>/dev/null | head
