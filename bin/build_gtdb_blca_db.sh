#!/usr/bin/env bash
# Build a BLCA-formatted GTDB SSU database for a given release.
#   usage: build_gtdb_blca_db.sh <release-number> <db-dir> [min-len]   e.g. 232 ./db 1000
#   min-len: minimum SSU sequence length to keep (default 1000). Use 0 to keep ALL
#            sequences (no <1000bp removal). Can also be set via env MINLEN.
#
# Produces in <db-dir>:
#   ssu_all_r<REL>.fna.gz            raw GTDB SSU download (kept)
#   gtdb_ssu_BLCAparsed.fasta        seqs with bare accession headers (+ blast index .n*)
#   gtdb_ssu_BLCAparsed.taxonomy     accession <TAB> superkingdom:..;...;species:..;
#   GTDB_VERSION.txt                 the release number
#
# env SKIP_BLAST=1: parse only -> write the FASTA + taxonomy but skip makeblastdb.
#            NB and Emu only need the parsed FASTA + taxonomy, not the BLAST index;
#            bin/build_gtdb_db.sh sets this when BLCA formatting was not requested.
#
# Needs: curl, gzip, awk, makeblastdb (BLAST+) on PATH.
set -euo pipefail

REL="${1:?release number required, e.g. 232}"
DBDIR="${2:?db dir required}"
MINLEN="${3:-${MINLEN:-1000}}"   # min SSU length to keep; 0 = keep all (no <1000bp removal)

# colours + live-progress only when attached to a terminal (suppressed in logs/CI)
if [ -t 1 ] && [ -t 2 ]; then
  C_INFO=$'\033[1;36m'; C_OK=$'\033[1;32m'; C_WARN=$'\033[1;33m'; C_RST=$'\033[0m'; TTY=1
else
  C_INFO=''; C_OK=''; C_WARN=''; C_RST=''; TTY=0
fi

mkdir -p "$DBDIR"
RAW="$DBDIR/ssu_all_r${REL}.fna.gz"
FASTA="$DBDIR/gtdb_ssu_BLCAparsed.fasta"
TAX="$DBDIR/gtdb_ssu_BLCAparsed.taxonomy"
# GTDB serves the file version-pinned under the release dir, and un-versioned
# (ssu_all.fna.gz) under releases/latest/. Try the pinned URL first, then latest.
URL_REL="https://data.gtdb.ecogenomic.org/releases/release${REL}/${REL}.0/genomic_files_all/ssu_all_r${REL}.fna.gz"
URL_LATEST="https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_all/ssu_all.fna.gz"

echo "${C_INFO}[build] obtaining GTDB r${REL} SSU${C_RST}"
if [ -s "$RAW" ]; then
    echo "${C_INFO}[build] reusing existing $RAW (skip download)${C_RST}"
elif curl -fSL --retry 3 -o "$RAW" "$URL_REL"; then
    echo "${C_OK}[build] got $URL_REL${C_RST}"
elif curl -fSL --retry 3 -o "$RAW" "$URL_LATEST"; then
    echo "${C_OK}[build] got $URL_LATEST${C_RST}"
else
    echo "${C_WARN}[build] ERROR: could not download GTDB SSU for r${REL}${C_RST}" >&2
    exit 1
fi

if [ "$MINLEN" -gt 0 ]; then
    echo "${C_INFO}[build] parsing -> BLCA FASTA (removing sequences < ${MINLEN} bp) + taxonomy${C_RST}"
else
    echo "${C_INFO}[build] parsing -> BLCA FASTA (keeping ALL sequences, no length filter) + taxonomy${C_RST}"
fi
# Single pass: emit cleaned FASTA and the accession->taxonomy map. Live counts
# (read / kept) go to stderr so the user sees the parse advancing.
gzip -dc "$RAW" | awk -v FASTA="$FASTA" -v TAX="$TAX" -v MINLEN="$MINLEN" \
                      -v CINFO="$C_INFO" -v CRST="$C_RST" -v TTY="$TTY" '
  function flush(   nm,i,rank,val,name,out) {
    if (id=="" || length(seq) < MINLEN) return
    print ">" id "\n" seq > FASTA
    kept++
    # taxstr like: d__Bacteria;p__..;c__..;o__..;f__..;g__..;s__Escherichia coli
    nlev = split(taxstr, L, ";")
    out = ""
    for (i=1;i<=nlev;i++) {
      rank = substr(L[i],1,1); val = substr(L[i],4)   # strip "x__"
      name = (rank=="d")?"superkingdom":(rank=="p")?"phylum":(rank=="c")?"class": \
             (rank=="o")?"order":(rank=="f")?"family":(rank=="g")?"genus":(rank=="s")?"species":""
      if (name!="" && val!="") out = out name ":" val ";"
    }
    print id "\t" out > TAX
  }
  /^>/ {
    flush()
    n++
    if (n % 100000 == 0) {
      if (TTY) printf "\r%s[build] parsing GTDB SSU… %d read, %d kept%s", CINFO, n, kept, CRST > "/dev/stderr"
      else     printf "%s[build] parsing GTDB SSU… %d read, %d kept%s\n", CINFO, n, kept, CRST > "/dev/stderr"
    }
    line = substr($0,2)
    sp = index(line," ")
    id = (sp>0)? substr(line,1,sp-1) : line
    rest = (sp>0)? substr(line,sp+1) : ""
    b = index(rest," [")
    taxstr = (b>0)? substr(rest,1,b-1) : rest
    seq = ""
    next
  }
  { seq = seq $0 }
  END {
    flush()
    if (TTY) printf "\r%s[build] parsing GTDB SSU… %d read, %d kept%s\n", CINFO, n, kept, CRST > "/dev/stderr"
    else     printf "%s[build] parsing GTDB SSU… done: %d read, %d kept%s\n", CINFO, n, kept, CRST > "/dev/stderr"
  }
'

if [ "${SKIP_BLAST:-0}" = 1 ]; then
  echo "${C_INFO}[build] SKIP_BLAST=1 — parse only: wrote FASTA + taxonomy, skipping makeblastdb (no BLAST index)${C_RST}"
else
echo "${C_INFO}[build] makeblastdb (-parse_seqids) — building BLAST index…${C_RST}"
mblog=$(mktemp)
makeblastdb -in "$FASTA" -dbtype nucl -parse_seqids -out "$FASTA" >"$mblog" 2>&1 &
mbpid=$!
start=$SECONDS; spin='|/-\'; i=0
while kill -0 "$mbpid" 2>/dev/null; do
  if [ "$TTY" = 1 ]; then
    c=${spin:i%4:1}; i=$((i+1))
    printf "\r%s[build] makeblastdb… %ds %s%s" "$C_INFO" "$((SECONDS-start))" "$c" "$C_RST" >&2
  fi
  sleep 1
done
if wait "$mbpid"; then
  [ "$TTY" = 1 ] && printf "\r\033[K" >&2
  echo "${C_OK}[build] makeblastdb done ($((SECONDS-start))s)${C_RST}"
  rm -f "$mblog"
else
  [ "$TTY" = 1 ] && printf "\r\033[K" >&2
  echo "${C_WARN}[build] makeblastdb FAILED:${C_RST}" >&2; cat "$mblog" >&2; rm -f "$mblog"; exit 1
fi
fi

echo "$REL" > "$DBDIR/GTDB_VERSION.txt"
echo "${C_OK}[build] done: $(grep -c '^>' "$FASTA") sequences, GTDB r${REL}${C_RST}"
