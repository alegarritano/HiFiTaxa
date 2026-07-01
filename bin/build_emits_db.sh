#!/usr/bin/env bash
# Stage the raw UNITE general-release FASTA as the EMITS minimap2 target. EMITS
# is the ITS read-level EM classifier (the fungal analogue of Emu); it aligns
# reads to the UNITE sequences with minimap2 and reads the species label directly
# from the UNITE headers, so its reference is just the UNITE FASTA itself — no
# build step, only a copy into place. Use this when params.marker == 'ITS'.
#
#   usage: build_emits_db.sh <unite_fasta> <db-dir>
#   e.g.   bash bin/build_emits_db.sh \
#              ~/Downloads/sh_general_release_dynamic_19.02.2025.fasta \
#              db_unite
#
#   <db-dir> is relative to the HiFiTaxa projectDir.
#
# Produces in <db-dir>:
#   unite.fasta            the raw UNITE FASTA, staged as the minimap2 target
#   EMITS_DB_VERSION.txt   release tag + provenance
#
# The pipeline points params.emits_db at <db-dir>/unite.fasta and runs:
#   minimap2 -cx map-hifi --secondary=yes -N 10 -p 0.95 ${params.emits_db} reads > aln.paf
#   emits run --input aln.paf --preset pacbio-hifi --rank species|genus
# (see modules/taxonomy_emits.nf). EMITS reads the species directly from the
# UNITE headers, so the headers must NOT be stripped — stage the file verbatim.
#
# Needs: nothing beyond coreutils (cp/gzip).
set -euo pipefail

UNITE_FASTA="${1:?UNITE FASTA required}"
DBDIR="${2:?db dir required}"
VERSION_TAG="${UNITE_VERSION:-sh_general_release_dynamic_19.02.2025}"

OUT="$DBDIR/unite.fasta"

[ -f "$UNITE_FASTA" ] || { echo "ERROR: $UNITE_FASTA not found" >&2; exit 1; }

mkdir -p "$DBDIR"

echo "[emits-db] staging raw UNITE FASTA -> $OUT (verbatim headers, minimap2 target)"
case "$UNITE_FASTA" in
    *.gz) gzip -dc "$UNITE_FASTA" > "$OUT" ;;
    *)    cp -f    "$UNITE_FASTA"   "$OUT" ;;
esac

{
    echo "Staged from: $UNITE_FASTA"
    echo "UNITE:       $VERSION_TAG"
    echo "Staged at:   $(date -u +%FT%TZ)"
} > "$DBDIR/EMITS_DB_VERSION.txt"

echo "[emits-db] done: $(grep -c '^>' "$OUT") sequences, UNITE ${VERSION_TAG}"
echo "[emits-db] EMITS target: $OUT  (set params.emits_db to this path)"
