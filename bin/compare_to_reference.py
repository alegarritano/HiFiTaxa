#!/usr/bin/env python3
"""Compare a local pipeline result against the bundled reference results and
print a GREEN (similar) / YELLOW (close) / RED (too different) verdict.

Always prints BOTH our reference and your local numbers.

usage: compare_to_reference.py <local_outdir> [reference_dir]
  <local_outdir>   a pipeline --outdir (contains dada2/dada2_ASV.fasta and
                   taxonomy_blca/blca_taxonomy_table.csv)
  [reference_dir]  default: <repo>/example/reference

Suggested thresholds (tune to taste):
  PASS  : ASV Jaccard >= 0.95  AND taxonomy agreement >= 0.95
          AND species Jaccard >= 0.90  (AND Bray-Curtis <= 0.10 if available)
  FAIL  : ASV Jaccard <  0.80  OR  taxonomy agreement <  0.85
          OR species Jaccard <  0.75  OR  Bray-Curtis > 0.30
  WARN  : anything in between
For the same data + same pinned GTDB release you should see ~1.0 everywhere
(denoise is deterministic; BLCA assignments are identical, only bootstrap
confidence jitters), so deviations point to a different DB/params/platform.
"""
import csv, sys
from pathlib import Path

G = "\033[1;32m"; Y = "\033[1;33m"; R = "\033[1;31m"; B = "\033[1m"; X = "\033[0m"

def find(d, *names):
    d = Path(d)
    for n in names:
        if (d / n).is_file():
            return d / n
        hit = next((p for p in d.rglob(n) if p.is_file()), None)
        if hit:
            return hit
    return None

def asv_ids(fa):
    if not (fa and Path(fa).is_file()):
        return set()
    return {l[1:].split()[0].strip() for l in open(fa) if l.startswith(">")}

def tax_map(csvf):
    m = {}
    if not csvf:
        return m
    with open(csvf) as f:
        for row in csv.DictReader(f):
            fid = row.get("Feature ID") or row.get("id")
            m[fid] = (row.get("Species") or row.get("species") or "").strip()
    return m

def jaccard(a, b):
    u = a | b
    return len(a & b) / len(u) if u else 1.0

def load_ft(tsv):
    d = {}
    if not tsv:
        return d
    for ln in open(tsv):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        try:
            d[p[0]] = sum(float(x) for x in p[1:])
        except ValueError:
            pass
    return d

def species_abund(ft, tax):
    out = {}
    for a, r in ft.items():
        out[tax.get(a, "NA")] = out.get(tax.get(a, "NA"), 0.0) + r
    tot = sum(out.values()) or 1.0
    return {k: v / tot for k, v in out.items()}

def bray_curtis(a, b):
    keys = set(a) | set(b)
    num = sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys)
    den = sum(a.values()) + sum(b.values())
    return num / den if den else 0.0

def main():
    if len(sys.argv) < 2:
        print("usage: compare_to_reference.py <local_outdir> [reference_dir]"); return 2
    local = sys.argv[1]
    ref = sys.argv[2] if len(sys.argv) > 2 else str(Path(__file__).resolve().parent.parent / "example" / "reference")

    L_fa, R_fa = find(local, "dada2_ASV.fasta"), find(ref, "dada2_ASV.fasta")
    L_csv, R_csv = find(local, "blca_taxonomy_table.csv"), find(ref, "blca_taxonomy_table.csv")
    L_ft, R_ft = find(local, "feature_table.tsv", "feature-table.tsv"), find(ref, "feature_table.tsv")
    if not (L_fa and L_csv):
        print(f"{R}ERROR: could not find dada2_ASV.fasta + blca_taxonomy_table.csv under {local}{X}"); return 2

    la, ra = asv_ids(L_fa), asv_ids(R_fa)
    lt, rt = tax_map(L_csv), tax_map(R_csv)
    shared = la & ra
    asv_j = jaccard(la, ra)
    tax_agree = (sum(1 for a in shared if lt.get(a) == rt.get(a)) / len(shared)) if shared else 0.0
    lsp = {v for v in lt.values() if v and v != "NA"}
    rsp = {v for v in rt.values() if v and v != "NA"}
    sp_j = jaccard(lsp, rsp)
    bc = None
    if L_ft and R_ft:
        bc = bray_curtis(species_abund(load_ft(L_ft), lt), species_abund(load_ft(R_ft), rt))

    print(f"\n{B}=== HiFiTaxa : local vs reference ==={X}")
    print(f"{'':26}{'REFERENCE (ours)':>18}{'LOCAL (yours)':>18}")
    print(f"{'ASVs':26}{len(ra):>18}{len(la):>18}")
    print(f"{'named species':26}{len(rsp):>18}{len(lsp):>18}")
    print(f"{'shared ASVs':26}{len(shared):>18}{'':>18}")
    print(f"\n{B}similarity metrics{X}")
    print(f"  ASV-set Jaccard            : {asv_j:.3f}")
    print(f"  taxonomy agreement (shared): {tax_agree:.3f}")
    print(f"  species-set Jaccard        : {sp_j:.3f}")
    if bc is not None:
        print(f"  species Bray-Curtis        : {bc:.3f}   (0 = identical)")
    else:
        print(f"  species Bray-Curtis        : n/a (no feature_table.tsv; presence-only)")

    is_pass = asv_j >= 0.95 and tax_agree >= 0.95 and sp_j >= 0.90 and (bc is None or bc <= 0.10)
    is_fail = asv_j < 0.80 or tax_agree < 0.85 or sp_j < 0.75 or (bc is not None and bc > 0.30)
    print()
    if is_pass:
        print(f"{G}PASS  Your results match the reference within tolerance.{X}\n"); return 0
    if is_fail:
        print(f"{R}FAIL  Your results differ substantially from the reference.{X}")
        print(f"{R}      Check the GTDB release/build, params, and platform.{X}\n"); return 1
    print(f"{Y}WARN  Close to the reference but outside the strict tolerance — review the metrics above.{X}\n")
    return 0

if __name__ == "__main__":
    sys.exit(main())
