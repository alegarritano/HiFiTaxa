#!/usr/bin/env python3
"""Merge a DADA2 ASV frequency table (samples x ASVs) with a classifier's
taxonomy output (BLCA, NB best_taxonomy, Emu-on-ASVs, etc.) into a single
wide TSV:

    Feature ID  <sample1>  <sample2>  ...  <sampleN>  Domain  Phylum  ...  Species  (Confidence)

Input — frequency table (TSV; QIIME2 biom export OR `qiime feature-table
transpose` export OR a plain TSV). Two header conventions accepted:
    a) #OTU ID  sample1  sample2  ...      (biom-style)
    b) Feature ID  sample1  sample2  ...   (QIIME2 export)
The first column is the ASV ID; the rest are sample frequencies.

Input — taxonomy. Two formats accepted; auto-detected:
    1) Wide CSV/TSV with one column per rank (HiFiTaxa's blca_taxonomy_table.csv,
       nb_taxonomy_table.csv): Feature ID, Domain, Phylum, Class, Order, Family,
       Genus, Species, [Confidence]
    2) QIIME2 TSVTaxonomyFormat (best_taxonomy.tsv): Feature ID  Taxon  Confidence
       where Taxon is `d__X; p__Y; ...; s__Z`.

Usage:
    merge_asv_freq_taxonomy.py \\
        --asv-freq results/dada2/feature-table.tsv  \\
        --taxonomy results/taxonomy_blca/blca_taxonomy_table.csv \\
        --out      results/blca_asv_freq_taxonomy.tsv

Optional:
    --lineage-col    write the taxonomy as ONE column ('Taxonomy') in the
                     d__X;p__Y;...;s__Z style instead of 7 separate columns
    --keep-confidence  include a Confidence column at the end (only when the
                     taxonomy file has one; ignored otherwise)
"""
import argparse, csv, sys, os, re

RANKS = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
GTDB_RANK_LETTERS = "dpcofgs"


def _sniff_delim(path):
    """TSV vs CSV based on the first non-empty line."""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            return "," if line.count(",") > line.count("\t") else "\t"
    return "\t"


def load_freq(path):
    """Returns (sample_names: list, freq: {feature_id: {sample: count}})."""
    delim = _sniff_delim(path)
    samples = None
    freq = {}
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter=delim)
        for row in rdr:
            if not row:
                continue
            first = row[0].strip()
            if first.startswith("# Constructed from") or first == "":
                continue
            if first in ("#OTU ID", "Feature ID", "OTU ID", "#OTU_ID"):
                samples = [c.strip() for c in row[1:]]
                continue
            if samples is None:
                # No explicit header — assume first non-comment row IS the header
                samples = [c.strip() for c in row[1:]]
                continue
            fid = first
            vals = {}
            for s, v in zip(samples, row[1:]):
                try:
                    fv = float(v)
                    vals[s] = int(fv) if fv.is_integer() else fv
                except (ValueError, AttributeError):
                    vals[s] = 0
            freq[fid] = vals
    if samples is None:
        sys.exit(f"[merge] couldn't find a header row in {path}")
    return samples, freq


def _parse_gtdb_taxon(taxon):
    """`d__Bacteria; p__X; ... ; s__Z` → 7-element list of names (empty for missing)."""
    out = [""] * 7
    if not taxon:
        return out
    for tok in taxon.split(";"):
        tok = tok.strip()
        if len(tok) >= 3 and tok[1:3] == "__":
            idx = GTDB_RANK_LETTERS.find(tok[0])
            if idx >= 0:
                out[idx] = tok[3:].strip()
    return out


def load_taxonomy(path):
    """Returns (feat_id → (lineage_list[7], confidence_str_or_empty)).

    Handles both HiFiTaxa's wide CSV (Feature ID,Domain,Phylum,...) and the
    QIIME2 TSVTaxonomyFormat (Feature ID, Taxon, Confidence) interchangeably.
    """
    delim = _sniff_delim(path)
    tax = {}
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter=delim)
        header = next(rdr, None)
        if header is None:
            return tax
        # Skip QIIME2's optional metadata-types row
        header = [h.strip() for h in header]
        cols = {h.lower(): i for i, h in enumerate(header)}
        # Decide format
        wide_form = all(r.lower() in cols for r in ("domain", "phylum", "species"))
        is_taxon = "taxon" in cols
        # Skip an optional `#q2:types` row
        peek = None
        try:
            peek = next(rdr)
            if peek and peek[0].startswith("#"):
                peek = None
        except StopIteration:
            pass

        def emit(row):
            if not row:
                return
            fid = row[0].strip()
            if not fid or fid.startswith("#"):
                return
            if wide_form:
                lin = [row[cols[r.lower()]].strip() if cols[r.lower()] < len(row) else ""
                       for r in RANKS]
                conf = row[cols["confidence"]].strip() if "confidence" in cols and cols["confidence"] < len(row) else ""
                tax[fid] = (lin, conf)
            elif is_taxon:
                taxon = row[cols["taxon"]].strip() if cols["taxon"] < len(row) else ""
                lin = _parse_gtdb_taxon(taxon)
                conf = row[cols["confidence"]].strip() if "confidence" in cols and cols["confidence"] < len(row) else ""
                tax[fid] = (lin, conf)

        if peek is not None:
            emit(peek)
        for row in rdr:
            emit(row)
    return tax


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--asv-freq", required=True,
                    help="TSV of ASV freq per sample (biom export or feature-table.tsv)")
    ap.add_argument("--taxonomy", required=True,
                    help="taxonomy table (wide CSV or QIIME2 TSVTaxonomyFormat)")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--lineage-col", action="store_true",
                    help="write taxonomy as a single 'Taxonomy' column "
                         "(d__X;p__Y;...) instead of 7 rank columns")
    ap.add_argument("--keep-confidence", action="store_true",
                    help="include a Confidence column if the taxonomy file has one")
    args = ap.parse_args()

    print(f"[merge] freq    : {args.asv_freq}", file=sys.stderr)
    print(f"[merge] taxonomy: {args.taxonomy}", file=sys.stderr)

    samples, freq = load_freq(args.asv_freq)
    print(f"[merge]   {len(samples)} samples, {len(freq):,} ASVs in freq table",
          file=sys.stderr)

    tax = load_taxonomy(args.taxonomy)
    print(f"[merge]   {len(tax):,} ASVs in taxonomy", file=sys.stderr)

    # Order ASVs: same order as the freq table (preserve DADA2 ranking)
    asv_order = list(freq.keys())
    overlap = sum(1 for a in asv_order if a in tax)
    print(f"[merge]   {overlap:,} / {len(asv_order):,} ASVs have a taxonomy entry",
          file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        # header
        header = ["Feature ID"] + samples
        if args.lineage_col:
            header.append("Taxonomy")
        else:
            header.extend(RANKS)
        if args.keep_confidence:
            header.append("Confidence")
        w.writerow(header)

        for fid in asv_order:
            row = [fid]
            for s in samples:
                row.append(freq[fid].get(s, 0))
            lin, conf = tax.get(fid, ([""] * 7, ""))
            if args.lineage_col:
                # GTDB-prefixed single string; drop empty trailing ranks
                parts = [f"{GTDB_RANK_LETTERS[i]}__{name}" for i, name in enumerate(lin) if name]
                row.append(";".join(parts))
            else:
                row.extend(lin)
            if args.keep_confidence:
                row.append(conf)
            w.writerow(row)

    print(f"[merge] wrote {args.out}  ({overlap:,} ASVs with taxonomy)", file=sys.stderr)


if __name__ == "__main__":
    main()
