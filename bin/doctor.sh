#!/usr/bin/env bash
# `hifitax doctor` — quick install health check. Run after a fresh install to
# confirm every component is in place: tools on PATH, conda/container envs,
# GTDB database, optional Emu DB and NB classifier.
#
# Exits non-zero if any required component is missing; warnings about optional
# components (Emu DB, NB classifier) don't cause a failing exit.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$HERE/.." && pwd)"
DB="${HIFITAX_DB:-$PROJ/db}"
EMU_DB="${HIFITAX_EMU_DB:-$PROJ/db_emu}"
NB_DIR="${HIFITAX_NB_DB_DIR:-$PROJ/db_nb}"
NB_GENUS="$NB_DIR/gtdb_ssu_dada2_genus.fa.gz"
NB_SPECIES="$NB_DIR/gtdb_ssu_dada2_species.fa.gz"

OK="\033[1;32m✓\033[0m"
WARN="\033[1;33m!\033[0m"
BAD="\033[1;31m✗\033[0m"
NL=$'\n'
fail=0

check_tool() {
    local label="$1" cmd="$2" required="$3"
    if command -v "$cmd" >/dev/null 2>&1; then
        printf "  ${OK} %-25s %s${NL}" "$label" "$(command -v "$cmd")"
    else
        if [ "$required" = "required" ]; then
            printf "  ${BAD} %-25s missing — REQUIRED${NL}" "$label"
            fail=1
        else
            printf "  ${WARN} %-25s missing (optional)${NL}" "$label"
        fi
    fi
}

check_file() {
    local label="$1" path="$2" required="$3"
    if [ -e "$path" ]; then
        local size
        size="$(du -sh "$path" 2>/dev/null | cut -f1)"
        printf "  ${OK} %-25s %s (%s)${NL}" "$label" "$path" "$size"
    else
        if [ "$required" = "required" ]; then
            printf "  ${BAD} %-25s missing — REQUIRED${NL}" "$label"
            fail=1
        else
            printf "  ${WARN} %-25s missing (optional)${NL}" "$label"
        fi
    fi
}

echo
echo "💻  HiFiTaxa doctor  🧬"
echo
echo "[tools on PATH]"
check_tool "nextflow"   nextflow   required
check_tool "python3"    python3    required
check_tool "blastn"     blastn     optional
check_tool "qiime"      qiime      optional
check_tool "emu"        emu        optional
check_tool "clustalo"   clustalo   optional
check_tool "muscle"     muscle     optional

echo
echo "[reference databases]"
check_file "GTDB BLCA fasta"       "$DB/gtdb_ssu_BLCAparsed.fasta"    required
check_file "GTDB BLCA taxonomy"    "$DB/gtdb_ssu_BLCAparsed.taxonomy" required
check_file "Emu DB dir"            "$EMU_DB"                          optional
check_file "DADA2 NB genus ref"    "$NB_GENUS"                        optional
check_file "DADA2 NB species ref"  "$NB_SPECIES"                      optional

echo
echo "[conda envs]"
for env in blca_nf qiime2-amplicon-2024.10 hifi_emu; do
    if [ -d "$HOME/miniconda3/envs/$env" ]; then
        printf "  ${OK} %-25s %s${NL}" "$env" "$HOME/miniconda3/envs/$env"
    else
        printf "  ${WARN} %-25s missing (will be set up on first run with -profile conda)${NL}" "$env"
    fi
done

echo
if [ "$fail" -eq 0 ]; then
    echo "✓ HiFiTaxa install looks healthy."
    echo "  Run the example end-to-end with:"
    echo "      python3 bin/run_pipeline.py --test"
    exit 0
else
    echo "✗ One or more REQUIRED components is missing — see the report above."
    echo "  Re-run the welcome wizard with:"
    echo "      python3 bin/run_pipeline.py"
    exit 1
fi
