#!/usr/bin/env bash
# Interactive one-stop GTDB database builder for HiFiTaxa.
#
# Downloads and parses the GTDB SSU reference ONCE, then asks which reference
# formats to build (BLCA, NB, Emu) and builds only the ones you choose. All three
# formats derive from the same parsed FASTA + taxonomy, so the download/parse
# happens a single time.
#
# Use this to prepare the database by hand -- e.g. on an HPC login (head) node
# before an offline/interactive run. The launcher (bin/run_pipeline.py) builds
# these automatically on first use; this is the manual, interactive equivalent.
#
#   usage: bash bin/build_gtdb_db.sh [release] [db-dir] [nb-dir] [emu-dir]
#          release  GTDB release number          (default 232)
#          db-dir   BLCA-parsed DB + BLAST index  (default db)
#          nb-dir   DADA2 NB references           (default db_nb)
#          emu-dir  Emu database                  (default db_emu)
#
# Emu note: the Emu build needs `emu`. If it is on PATH it is used directly;
# otherwise the build runs inside the Emu container (singularity); otherwise the
# Emu step is skipped with instructions. Run on a node WITH internet.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # bin/
PROJ="$(cd "$HERE/.." && pwd)"                          # repo root

REL="${1:-232}"
DBDIR="${2:-db}"
NBDIR="${3:-db_nb}"
EMUDIR="${4:-db_emu}"

ask() {  # ask "question"  -> 0 = yes (default), 1 = no
  local ans
  read -r -p "  $1 [Y/n] " ans </dev/tty || ans=""
  case "${ans,,}" in n|no) return 1 ;; *) return 0 ;; esac
}

echo "=== HiFiTaxa GTDB database builder — release r${REL} ==="
echo

if ask "Remove reference sequences shorter than 1000 bp? (recommended for full-length 16S)"; then
  MINLEN=1000
else
  MINLEN=0
fi
ask "Format for BLCA?              (BLAST index)"                        && DO_BLCA=1 || DO_BLCA=0
ask "Format for the NB classifier? (DADA2 genus + species references)"  && DO_NB=1   || DO_NB=0
ask "Format for Emu?"                                                    && DO_EMU=1  || DO_EMU=0
echo

if [ "$DO_BLCA$DO_NB$DO_EMU" = "000" ]; then
  echo "Nothing selected — no formats to build. Exiting."
  exit 0
fi

FASTA="$DBDIR/gtdb_ssu_BLCAparsed.fasta"
TAX="$DBDIR/gtdb_ssu_BLCAparsed.taxonomy"

# 1) Shared download + parse. Build the BLAST index only if BLCA was requested;
#    NB and Emu only need the parsed FASTA + taxonomy.
echo ">>> downloading + parsing GTDB SSU r${REL} into ${DBDIR}/"
if [ "$DO_BLCA" = 1 ]; then
  bash "$HERE/build_gtdb_blca_db.sh" "$REL" "$DBDIR" "$MINLEN"
else
  SKIP_BLAST=1 bash "$HERE/build_gtdb_blca_db.sh" "$REL" "$DBDIR" "$MINLEN"
fi

# 2) NB references (python only — no container, no internet).
if [ "$DO_NB" = 1 ]; then
  echo ">>> building NB (DADA2) references into ${NBDIR}/"
  mkdir -p "$NBDIR"
  bash "$HERE/build_gtdb_dada2_db.sh" "$FASTA" "$TAX" \
    "$NBDIR/gtdb_ssu_dada2_genus.fa.gz" \
    "$NBDIR/gtdb_ssu_dada2_species.fa.gz"
fi

# 3) Emu database (needs `emu`: PATH -> container -> skip).
if [ "$DO_EMU" = 1 ]; then
  echo ">>> building Emu database into ${EMUDIR}/"
  mkdir -p "$EMUDIR"
  if command -v emu >/dev/null 2>&1; then
    bash "$HERE/build_gtdb_emu_db.sh" "$FASTA" "$TAX" "$EMUDIR"
  elif command -v singularity >/dev/null 2>&1; then
    EMU_IMG="$(grep -E "emu_container[[:space:]]*=" "$PROJ/nextflow.config" | sed -E "s/.*'([^']+)'.*/\1/")"
    echo "    emu not on PATH; building inside the Emu container"
    echo "    ($EMU_IMG)"
    singularity exec --bind "$PWD" --bind "$PROJ" "docker://$EMU_IMG" \
      bash "$HERE/build_gtdb_emu_db.sh" "$FASTA" "$TAX" "$EMUDIR"
  else
    echo "    SKIPPED: neither 'emu' nor 'singularity' is available."
    echo "    Activate an env with emu, or run by hand:"
    echo "      bash bin/build_gtdb_emu_db.sh $FASTA $TAX $EMUDIR"
  fi
fi

echo
echo "=== done: GTDB r${REL} references built in ${DBDIR}/$([ "$DO_NB" = 1 ] && echo " + ${NBDIR}/")$([ "$DO_EMU" = 1 ] && echo " + ${EMUDIR}/") ==="
