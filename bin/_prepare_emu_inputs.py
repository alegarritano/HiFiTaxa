#!/usr/bin/env python3
"""Convert the GTDB BLCA-parsed (FASTA + taxonomy TSV) into Emu-format inputs.

Emu's `build-database` (custom mode) expects:
  sequences.fasta : FASTA where each record id (first whitespace-delimited
                    token of the header) is a seq id present in seq2tax.map.
                    The headers do NOT need to embed the tax id; Emu builds its
                    internal `species_taxid.fasta` itself, keyed via seq2tax.
  seq2tax.map     : two-column TSV  (seq_id <tab> tax_id)
  taxonomy.tsv    : header row + one row per tax_id; first column = tax_id,
                    remaining columns are the rank labels (any rank scheme).

This script:
  1) assigns a unique synthetic tax_id to each unique GTDB species lineage,
  2) writes seq2tax.map and taxonomy.tsv,
  3) filters the input FASTA to records whose accession has a known taxonomy
     (passing the original headers through unchanged).

It is called by bin/build_gtdb_emu_db.sh and is not intended to be run alone.
"""
import argparse
import os
import sys


RANK_PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")


def strip_prefix(s):
    s = s.strip()
    for p in RANK_PREFIXES:
        if s.startswith(p):
            return s[len(p):]
    return s


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--blca-fasta", required=True,
                    help="BLCA-parsed GTDB SSU FASTA (e.g. db/gtdb_ssu_BLCAparsed.fasta)")
    ap.add_argument("--blca-tax", required=True,
                    help="BLCA-parsed GTDB taxonomy TSV (accession<tab>lineage)")
    ap.add_argument("--out-dir", required=True,
                    help="Destination dir for Emu-format inputs")
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Read accession -> lineage; assign one tax_id per unique lineage.
    acc2tax = {}
    lin2id = {}
    next_id = 1
    with open(args.blca_tax) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            acc, lin = parts[0], parts[1]
            tid = lin2id.get(lin)
            if tid is None:
                tid = next_id
                lin2id[lin] = tid
                next_id += 1
            acc2tax[acc] = tid

    print(f"[emu-prep] {len(acc2tax)} sequences mapped to {len(lin2id)} unique lineages",
          file=sys.stderr)

    # 2. Write seq2tax.map
    seq2tax_path = os.path.join(args.out_dir, "seq2tax.map")
    with open(seq2tax_path, "w") as fh:
        for acc, tid in acc2tax.items():
            fh.write(f"{acc}\t{tid}\n")

    # 3. Write taxonomy.tsv (tax_id + 7 stripped ranks)
    tax_path = os.path.join(args.out_dir, "taxonomy.tsv")
    with open(tax_path, "w") as fh:
        fh.write("tax_id\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
        for lin, tid in lin2id.items():
            ranks = [strip_prefix(x) for x in lin.split(";")]
            # pad or truncate to exactly 7 ranks
            ranks = (ranks + [""] * 7)[:7]
            fh.write(f"{tid}\t" + "\t".join(ranks) + "\n")

    # 4. Filter FASTA: keep only records whose accession is in seq2tax.
    #    Pass through original headers (Emu keys seq2tax_dict by record.id, i.e.
    #    the first whitespace token of the header).
    fa_path = os.path.join(args.out_dir, "sequences.fasta")
    n_kept = 0
    n_skipped_records = 0
    skip_block = False
    with open(args.blca_fasta) as src, open(fa_path, "w") as dst:
        for line in src:
            if line.startswith(">"):
                hdr = line[1:].rstrip("\n")
                acc = hdr.split()[0] if hdr else ""
                if acc not in acc2tax:
                    n_skipped_records += 1
                    skip_block = True
                    continue
                skip_block = False
                dst.write(line)
                n_kept += 1
            elif not skip_block:
                dst.write(line)

    print(f"[emu-prep] FASTA: kept {n_kept} sequences, "
          f"skipped {n_skipped_records} (no taxonomy mapping)", file=sys.stderr)
    print(f"[emu-prep] wrote {seq2tax_path}, {tax_path}, {fa_path}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main() or 0)
