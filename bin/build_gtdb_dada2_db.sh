#!/usr/bin/env bash
# Convert HiFiTaxa's BLCA-format GTDB reference into the TWO DADA2 references the
# Naive-Bayes path needs (canonical DADA2 two-step design):
#
#   1. GENUS reference  : headers are the 6-rank lineage Kingdom..Genus
#        >Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia
#        used by dada2::assignTaxonomy() (bootstrap, per-rank confidence).
#   2. SPECIES reference : headers are ">accession Genus species"
#        >RS_GCF_000123.1 Escherichia coli
#        used by dada2::addSpecies() (exact 100% match, genus-consistency checked).
#
# WHY TWO FILES (not species-in-headers like PacBio): assignTaxonomy with the full
# species lineage in headers does NOT scale to full GTDB SSU r232 (~82k unique
# species labels) — the bootstrap collapses to Kingdom-only assignments. Genus-level
# training + exact-match species is the standard DADA2 design (as used by the
# SILVA/RDP training sets) and runs against the same full GTDB release as BLCA and
# Emu. See scripts/dada2_assign_tax.R for the matching classifier logic.
#
# Both references contain the SAME sequences (and the SAME BLCA-parsed GTDB that
# BLCA + Emu use), so all three classifiers stay anchored on one reference.
#
# Usage:
#   bash bin/build_gtdb_dada2_db.sh <blca_fasta> <blca_tax> <out_genus_fa_gz> <out_species_fa_gz>
#
# Example:
#   bash bin/build_gtdb_dada2_db.sh \
#       db/gtdb_ssu_BLCAparsed.fasta \
#       db/gtdb_ssu_BLCAparsed.taxonomy \
#       db_nb/gtdb_ssu_dada2_genus.fa.gz \
#       db_nb/gtdb_ssu_dada2_species.fa.gz
#
# The BLCA taxonomy file is `accession<TAB>superkingdom:X;phylum:Y;...;species:Z;`.

set -eo pipefail

BLCA_FASTA="${1:?missing arg: BLCA-parsed FASTA}"
BLCA_TAX="${2:?missing arg: BLCA taxonomy TSV}"
OUT_GENUS="${3:?missing arg: output genus dada2 fasta.gz}"
OUT_SPECIES="${4:?missing arg: output species dada2 fasta.gz}"

[ -f "$BLCA_FASTA" ] || { echo "ERROR: $BLCA_FASTA not found" >&2; exit 1; }
[ -f "$BLCA_TAX"   ] || { echo "ERROR: $BLCA_TAX not found"   >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }

mkdir -p "$(dirname "$OUT_GENUS")" "$(dirname "$OUT_SPECIES")"
echo "[dada2-db] reformatting BLCA reference -> genus + species DADA2 references"

python3 - "$BLCA_FASTA" "$BLCA_TAX" "$OUT_GENUS" "$OUT_SPECIES" <<'PY'
import sys, gzip
fasta_in, tax_in, genus_out, species_out = sys.argv[1:5]

# Load accession -> parsed ranks dict
acc2ranks = {}
with open(tax_in) as fh:
    for line in fh:
        line = line.rstrip("\n").rstrip("\r")
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) < 2:
            continue
        acc, blca_tax = parts
        ranks = {}
        for tok in blca_tax.split(";"):
            tok = tok.strip()
            if ":" in tok:
                rk, name = tok.split(":", 1)
                ranks[rk.strip()] = name.strip()
        acc2ranks[acc] = ranks

print(f"  loaded {len(acc2ranks):,} taxonomy entries", file=sys.stderr)

# BLCA uses "superkingdom" for Kingdom.
GENUS_ORDER = ["superkingdom", "phylum", "class", "order", "family", "genus"]

def genus_header(ranks):
    """6-rank lineage Kingdom..Genus; trailing empty ranks trimmed. None if empty."""
    lineage = [ranks.get(r, "") for r in GENUS_ORDER]
    while lineage and lineage[-1] == "":
        lineage.pop()
    if not lineage:
        return None
    return ";".join(lineage)

def species_header(ranks):
    """'Genus species' for addSpecies. None unless we have genus + a binomial."""
    g  = ranks.get("genus", "").strip()
    sp = ranks.get("species", "").strip()
    if not g or not sp:
        return None
    # GTDB species is normally already the binomial ("Escherichia coli"); if it is
    # epithet-only, prepend the genus so token2=genus, token3=epithet for addSpecies.
    if not sp.startswith(g + " "):
        sp = g + " " + sp
    if len(sp.split()) < 2:
        return None
    return sp

written_g = written_s = skipped = 0
with open(fasta_in) as fh_in, \
     gzip.open(genus_out, "wt") as fh_g, \
     gzip.open(species_out, "wt") as fh_s:
    cur_acc = None
    cur_seq = []

    def flush():
        global written_g, written_s, skipped
        if cur_acc is None:
            return
        ranks = acc2ranks.get(cur_acc)
        if ranks is None:
            skipped += 1
            return
        seq = "".join(cur_seq)
        gh = genus_header(ranks)
        if gh is not None:
            fh_g.write(">" + gh + "\n" + seq + "\n")
            written_g += 1
        sh = species_header(ranks)
        if sh is not None:
            fh_s.write(">" + cur_acc + " " + sh + "\n" + seq + "\n")
            written_s += 1

    for line in fh_in:
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            flush()
            cur_acc = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)
    flush()

print(f"  genus ref  : wrote {written_g:,} seqs -> {genus_out}",   file=sys.stderr)
print(f"  species ref: wrote {written_s:,} seqs -> {species_out}", file=sys.stderr)
print(f"  skipped {skipped:,} (no taxonomy entry for accession)",  file=sys.stderr)
PY

echo "[dada2-db] done."
ls -lh "$OUT_GENUS" "$OUT_SPECIES"
echo "[dada2-db] use: Rscript scripts/dada2_assign_tax.R asvs.fa <threads> $OUT_GENUS $OUT_SPECIES 80"
