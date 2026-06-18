#!/usr/bin/env bash
# Build a GTDB-formatted Emu database from the existing BLCA-parsed GTDB SSU
# files. Keeps the two classifiers (BLCA, Emu) anchored on the SAME reference
# release so their assignments are directly comparable.
#
# Usage:
#   bash bin/build_gtdb_emu_db.sh <blca_fasta> <blca_tax> <out_dir>
#
# Example:
#   bash bin/build_gtdb_emu_db.sh \
#       db/gtdb_ssu_BLCAparsed.fasta \
#       db/gtdb_ssu_BLCAparsed.taxonomy \
#       db_emu
#
# Requires:
#   - the BLCA-parsed GTDB DB already built (bin/build_gtdb_blca_db.sh)
#   - `emu` on PATH (conda env hifitax_emu, or run inside the Emu container)
#   - python3, awk

set -eo pipefail

# --------------------------------------------------------------------------- #
# Self-activate a conda env with `emu` if one isn't already on PATH. Lets the
# script Just Work whether called by the launcher or directly. set -u is
# deferred until AFTER conda activate (activate.d hooks trip nounset).
# --------------------------------------------------------------------------- #
if ! command -v emu >/dev/null 2>&1; then
    for conda_sh in \
        "$HOME/miniconda3/etc/profile.d/conda.sh" \
        "$HOME/anaconda3/etc/profile.d/conda.sh" \
        "$HOME/.conda/etc/profile.d/conda.sh" \
        "/opt/conda/etc/profile.d/conda.sh" \
        "/srv/scratch/$USER/miniconda3/etc/profile.d/conda.sh" \
        "/scratch/$USER/miniconda3/etc/profile.d/conda.sh" ; do
        [ -f "$conda_sh" ] && source "$conda_sh" && break
    done
    if type conda >/dev/null 2>&1; then
        for env in hifitax_emu hifi_emu emu_nf emu; do
            conda activate "$env" 2>/dev/null && command -v emu >/dev/null && break
            conda deactivate 2>/dev/null || true
        done
        if ! command -v emu >/dev/null; then
            for env in $(conda env list 2>/dev/null | awk '/^[^#]/{print $1}'); do
                [ -z "$env" ] && continue
                conda activate "$env" 2>/dev/null && command -v emu >/dev/null && break
                conda deactivate 2>/dev/null || true
            done
        fi
    fi
fi
set -u

BLCA_FASTA="${1:?missing arg: BLCA-parsed FASTA}"
BLCA_TAX="${2:?missing arg: BLCA taxonomy TSV}"
OUT_DIR="${3:?missing arg: output Emu DB dir}"

[ -f "$BLCA_FASTA" ] || { echo "ERROR: $BLCA_FASTA not found" >&2; exit 1; }
[ -f "$BLCA_TAX"   ] || { echo "ERROR: $BLCA_TAX not found"   >&2; exit 1; }
command -v emu     >/dev/null || { echo "ERROR: 'emu' not on PATH. The launcher (bin/run_pipeline.py) auto-installs one; or run: mamba create -n hifitax_emu -c bioconda -c conda-forge 'emu>=3.6.2' 'minimap2>=2.24' 'samtools>=1.17' -y" >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUT_DIR"
STAGE="${OUT_DIR}/_emu_stage"
mkdir -p "$STAGE"
# Absolutise STAGE: the emu build below cd's into OUT_DIR, after which a relative
# STAGE would resolve wrong (realpath returns '' -> emu gets an empty path and
# crashes). Absolute STAGE is cwd-proof, so a relative <out_dir> works too.
STAGE="$(cd "$STAGE" && pwd)"

echo "[emu-db] 1/3 preparing Emu-format inputs from GTDB BLCA-parsed files"
python3 "$HERE/_prepare_emu_inputs.py" \
    --blca-fasta "$BLCA_FASTA" \
    --blca-tax   "$BLCA_TAX" \
    --out-dir    "$STAGE"

echo "[emu-db] 2/3 running 'emu build-database' (silent step; growth-watcher below shows disk-write progress)"
DB_NAME="$(basename "$OUT_DIR")"

# Background growth watcher: `emu build-database` emits no progress, so we sample
# the OUT_DIR size every 10s and print growth + elapsed time. Killed once the
# build returns.
_watch_emu_build() {
    local target="$1"
    local started=$SECONDS
    local prev_kb=0
    while :; do
        sleep 10
        local kb=0
        if [ -d "$target" ]; then
            kb=$(du -sk "$target" 2>/dev/null | awk '{print $1}')
            kb=${kb:-0}
        fi
        local elapsed=$(( SECONDS - started ))
        local mb=$(( kb / 1024 ))
        local delta_mb=$(( (kb - prev_kb) / 1024 ))
        printf "[emu-db]   building... %5d MB written  (+%d MB last 10s, %ds elapsed)\n" \
               "$mb" "$delta_mb" "$elapsed"
        prev_kb=$kb
    done
}
_watch_emu_build "$OUT_DIR" &
WATCHER_PID=$!
# Make sure the watcher dies even on errors below
trap 'kill "$WATCHER_PID" 2>/dev/null || true' EXIT INT TERM

# Emu writes into a subdir named DB_NAME inside the cwd. Work in OUT_DIR so it
# lands in place, then flatten if Emu created a nested copy.
(
    cd "$OUT_DIR"
    emu build-database "$DB_NAME" \
        --sequences     "$(realpath "$STAGE/sequences.fasta")" \
        --seq2tax       "$(realpath "$STAGE/seq2tax.map")" \
        --taxonomy-list "$(realpath "$STAGE/taxonomy.tsv")"
)

# Build done -> stop the watcher
kill "$WATCHER_PID" 2>/dev/null || true
wait "$WATCHER_PID" 2>/dev/null || true
trap - EXIT INT TERM

# If Emu produced "<OUT_DIR>/<DB_NAME>/..." flatten it up one level so the
# pipeline can point --emu_db_dir at <OUT_DIR> directly.
if [ -d "$OUT_DIR/$DB_NAME" ]; then
    shopt -s dotglob nullglob
    mv "$OUT_DIR/$DB_NAME"/* "$OUT_DIR"/ 2>/dev/null || true
    rmdir "$OUT_DIR/$DB_NAME" 2>/dev/null || true
fi

# Optional cleanup of the staging dir; comment out to keep for debugging
rm -rf "$STAGE"

echo "[emu-db] 3/3 writing version stamp"
gtdb_release="$(cat "$(dirname "$BLCA_FASTA")/GTDB_VERSION.txt" 2>/dev/null || echo unknown)"
{
    echo "Built from: $BLCA_FASTA"
    echo "Taxonomy:   $BLCA_TAX"
    echo "GTDB:       $gtdb_release"
    echo "Built at:   $(date -u +%FT%TZ)"
    echo "Emu:        $(emu --version 2>&1 | head -1 || echo unknown)"
} > "$OUT_DIR/EMU_DB_VERSION.txt"

echo "[emu-db] done. Database at: $OUT_DIR"
echo "[emu-db] If --emu_db_dir is not already set to this path in nextflow.config,"
echo "[emu-db] pass --emu_db_dir $OUT_DIR on the launcher command line."
