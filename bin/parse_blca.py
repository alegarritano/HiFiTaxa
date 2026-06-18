#!/usr/bin/env python3
"""Parse a BLCA .blca.out file into tidy taxonomy tables.

BLCA line format:
  <ASVID>\tsuperkingdom:Name;conf;phylum:Name;conf;...;species:Name;conf;
  <ASVID>\tUnclassified
  <ASVID>\tSkipped

Outputs:
  blca_taxonomy_table.csv       Feature ID + 7 ranks (names only; for phyloseq)
  blca_taxonomy_confidence.csv  same + per-rank bootstrap confidence
"""
import csv, sys

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
COLS  = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

def parse_line(rest):
    """Return (names dict, conf dict) keyed by rank."""
    names = {r: "NA" for r in RANKS}
    confs = {r: "" for r in RANKS}
    if rest in ("Unclassified", "Skipped", ""):
        return names, confs, rest
    toks = [t for t in rest.split(";") if t != ""]
    # tokens alternate: 'rank:name', 'conf', 'rank:name', 'conf', ...
    i = 0
    while i < len(toks):
        if ":" in toks[i]:
            rank, name = toks[i].split(":", 1)
            conf = toks[i + 1] if (i + 1) < len(toks) and ":" not in toks[i + 1] else ""
            if rank in names:
                names[rank] = name
                confs[rank] = conf
            i += 2
        else:
            i += 1
    return names, confs, "OK"

def main(infile):
    rows_tax, rows_conf = [], []
    n_ok = n_unc = 0
    with open(infile) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            asv = parts[0]
            rest = parts[1] if len(parts) > 1 else ""
            names, confs, status = parse_line(rest)
            if status == "OK":
                n_ok += 1
            else:
                n_unc += 1
            rows_tax.append([asv] + [names[r] for r in RANKS])
            rows_conf.append([asv] + [v for r in RANKS for v in (names[r], confs[r])])

    with open("blca_taxonomy_table.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["Feature ID"] + COLS); w.writerows(rows_tax)
    conf_header = ["Feature ID"] + [c for col in COLS for c in (col, col + "_conf")]
    with open("blca_taxonomy_confidence.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(conf_header); w.writerows(rows_conf)

    print(f"parsed {n_ok + n_unc} ASVs: {n_ok} classified, {n_unc} unclassified/skipped")
    print("-> blca_taxonomy_table.csv, blca_taxonomy_confidence.csv")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "dada2_ASV.blca.out")
