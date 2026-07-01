#!/usr/bin/env bash
# Fungal/ITS analogue of build_gtdb_dada2_db.sh, but for the SINGLE-STEP
# Naive-Bayes design that ITS requires. Use this when params.marker == 'ITS'.
#
# Convert the BLCA-format UNITE reference into the ONE DADA2 reference the
# single-step ITS Naive-Bayes path needs:
#
#   unite_full_singlestep_ref.fa.gz : headers are the full 7-rank lineage
#       Kingdom..Species (de-underscored binomial species), terminated by ';'
#         >Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Aspergillus;Aspergillus niger;
#       consumed by dada2::assignTaxonomy() with
#         taxLevels = Kingdom,Phylum,Class,Order,Family,Genus,Species
#       (NO addSpecies). See scripts/dada2_assign_tax_singlestep.R for the
#       matching classifier logic, and process nb_classify_singlestep in
#       modules/taxonomy_nb.nf for the wiring.
#
# WHY SINGLE-STEP (not the 16S two-step genus + addSpecies): on ITS the
# two-step addSpecies exact-match path collapses to ~0 species (ITS is variable
# and the UNITE species labels do not align with addSpecies' exact-100%-match
# assumption). The single-step 7-rank assignTaxonomy with the species lineage in
# the headers is the ITS-appropriate NB design and is what the fungal benchmark
# uses. ITS reference sets (UNITE ~100k records) are small enough that the
# bootstrap does not collapse the way full GTDB SSU r232 (~82k species) does for
# 16S, so the species rank can live directly in the training headers here.
#
# The reference contains the SAME sequences as the BLCA-parsed UNITE that BLCA
# and EMITS use, so all three classifiers stay anchored on one UNITE release.
#
# Usage:
#   bash bin/build_unite_dada2_db.sh <blca_fasta> <blca_tax> <out_singlestep_fa_gz>
#
# Example:
#   bash bin/build_unite_dada2_db.sh \
#       db_unite/unite_BLCAparsed.fasta \
#       db_unite/unite_BLCAparsed.taxonomy \
#       db_unite/unite_full_singlestep_ref.fa.gz
#
# Paths are relative to the HiFiTaxa projectDir; the pipeline expects the output
# at ${projectDir}/db_unite/unite_full_singlestep_ref.fa.gz
# (params.unite_dada2_singlestep_db).
#
# The BLCA taxonomy file is `accession<TAB>superkingdom:X;phylum:Y;...;species:Z;`.

set -eo pipefail

BLCA_FASTA="${1:?missing arg: BLCA-parsed FASTA}"
BLCA_TAX="${2:?missing arg: BLCA taxonomy TSV}"
OUT_SINGLESTEP="${3:?missing arg: output single-step dada2 fasta.gz}"

[ -f "$BLCA_FASTA" ] || { echo "ERROR: $BLCA_FASTA not found" >&2; exit 1; }
[ -f "$BLCA_TAX"   ] || { echo "ERROR: $BLCA_TAX not found"   >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }

mkdir -p "$(dirname "$OUT_SINGLESTEP")"
echo "[dada2-db] reformatting BLCA reference -> single-step 7-rank DADA2 reference"

python3 - "$BLCA_FASTA" "$BLCA_TAX" "$OUT_SINGLESTEP" <<'PY'
import sys, gzip
fasta_in, tax_in, out = sys.argv[1:4]

# Full 7-rank order. BLCA stores Kingdom under the "superkingdom" key (UNITE k__
# value "Fungi"); we render it as the leading taxLevel for assignTaxonomy.
ORD = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

# Load accession -> ';'-joined 7-rank header (trailing empty ranks trimmed).
acc2h = {}
with open(tax_in) as fh:
    for line in fh:
        line = line.rstrip("\n").rstrip("\r")
        if "\t" not in line:
            continue
        acc, lin = line.split("\t", 1)
        d = {}
        for tok in lin.split(";"):
            tok = tok.strip()
            if ":" in tok:
                k, v = tok.split(":", 1)
                d[k.strip()] = v.strip()
        vals = [d.get(r, "") for r in ORD]
        while vals and vals[-1] == "":
            vals.pop()
        if vals:
            acc2h[acc] = ";".join(vals) + ";"

print(f"  loaded {len(acc2h):,} taxonomy entries", file=sys.stderr)

written = skipped = 0
with open(fasta_in) as fh_in, gzip.open(out, "wt") as fh_out:
    cur_acc = None
    cur_seq = []

    def flush():
        global written, skipped
        if cur_acc is None:
            return
        header = acc2h.get(cur_acc)
        if header is None:
            skipped += 1
            return
        fh_out.write(">" + header + "\n" + "".join(cur_seq) + "\n")
        written += 1

    for line in fh_in:
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            flush()
            cur_acc = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)
    flush()

print(f"  single-step ref: wrote {written:,} seqs -> {out}", file=sys.stderr)
print(f"  skipped {skipped:,} (no taxonomy entry for accession)", file=sys.stderr)
PY

echo "[dada2-db] done."
ls -lh "$OUT_SINGLESTEP"
echo "[dada2-db] use: Rscript scripts/dada2_assign_tax_singlestep.R asvs.fa <threads> $OUT_SINGLESTEP 80"
