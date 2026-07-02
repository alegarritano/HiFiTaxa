#!/usr/bin/env python3
"""Scan PacBio HiFi reads for the amplicon primers actually present.

Given a FASTQ(.gz), test a panel of candidate primers against the 5' and 3' ends
of a sample of reads (both read orientations, IUPAC-aware, up to a mismatch
budget) and report which pair best explains the data and where it sits. Use this
to pick --forward_primer / --reverse_primer for a HiFiTaxa run when the primers
used to generate a read set are unknown or undocumented.

No dependencies beyond the Python 3 standard library.

Usage:
    python scan_primers.py reads.fastq.gz
    python scan_primers.py reads.fastq.gz --sample 5000 --window 60 --max-mismatch 3

The primer panel defaults to HiFiTaxa's fungal/eukaryote candidates (see
docs/PRIMERS.md) plus the 16S 27F/1492R pair; edit PRIMERS below or pass
--primers name=SEQ,name=SEQ to supply your own.
"""
import argparse
import gzip
import io
import re
import sys

# name -> (sequence 5'->3'). Forward primers anchor the read 5' end; reverse
# primers appear as their reverse complement at the read 3' end (and vice versa
# on the opposite strand, which --both-orientations covers).
PRIMERS = {
    # fungal / eukaryote long-read (docs/PRIMERS.md)
    "ITS1catta":    "ACCWGCGGARGGATCATTA",   # fwd, 3' 18S upstream of ITS1
    "LR5_TW14ngs":  "TCCTGAGGGAAACTTCG",     # rev, 28S/LSU
    "ITS4ngsUni":   "CCTSCSCTTANTDATATGC",   # rev, end of ITS2
    "ITS1F":        "CTTGGTCATTTAGAGGAAGTAA", # fwd, fungal ITS1
    "ITS4":         "TCCTCCGCTTATTGATATGC",   # rev, fungal ITS4
    # HiFiTaxa ITS config default pair
    "1391F":        "GTACACACCGCCCGTC",       # fwd (config forward_primer ITS)
    "ITS4ngsUni_rc":"GCATATHANTAAGSGSAGG",    # rev, config reverse_primer ITS (= revcomp ITS4ngsUni)
    # 16S (for ATCC / bacterial read sets)
    "27F":          "AGRGTTYGATYMTGGCTCAG",
    "1492R":        "AAGTCGTAACAAGGTARCY",
}

# known role of each panel primer, used only to label the verdict fwd/rev.
ROLE = {
    "ITS1catta": "fwd", "1391F": "fwd", "ITS1F": "fwd", "27F": "fwd",
    "LR5_TW14ngs": "rev", "ITS4ngsUni": "rev", "ITS4": "rev",
    "ITS4ngsUni_rc": "rev", "1492R": "rev",
}

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}
_COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                      "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def iupac_regex(primer):
    return re.compile("".join(IUPAC.get(b, re.escape(b)) for b in primer.upper()))


def fuzzy_find(pattern_seq, text, max_mismatch):
    """Return (start, mismatches) of the best IUPAC match of pattern_seq in text
    within the mismatch budget, else None. Simple sliding window — primers are
    short and the search window is small, so this is fast enough."""
    plen = len(pattern_seq)
    classes = [IUPAC.get(b, b) for b in pattern_seq.upper()]
    best = None
    for i in range(0, len(text) - plen + 1):
        mism = 0
        for j, cls in enumerate(classes):
            ch = text[i + j]
            ok = (ch in cls) if len(cls) > 1 else (ch == cls)
            if not ok:
                mism += 1
                if mism > max_mismatch:
                    break
        else:
            if best is None or mism < best[1]:
                best = (i, mism)
                if mism == 0:
                    break
    return best


def open_maybe_gz(path):
    if path == "-":
        return io.TextIOWrapper(sys.stdin.buffer)
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path)


def iter_reads(path, limit):
    n = 0
    with open_maybe_gz(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip().upper()
            fh.readline()  # +
            fh.readline()  # qual
            if not seq:
                break
            yield seq
            n += 1
            if limit and n >= limit:
                break


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq", help="reads FASTQ or FASTQ.gz ('-' for stdin)")
    ap.add_argument("--sample", type=int, default=2000,
                    help="reads to scan (default 2000; 0 = all)")
    ap.add_argument("--window", type=int, default=50,
                    help="bp at each read end to search (default 50)")
    ap.add_argument("--max-mismatch", type=int, default=3,
                    help="mismatch budget per primer match (default 3)")
    ap.add_argument("--primers",
                    help="override panel: name=SEQ,name=SEQ,...")
    args = ap.parse_args()

    primers = dict(PRIMERS)
    if args.primers:
        primers = {}
        for kv in args.primers.split(","):
            name, seq = kv.split("=", 1)
            primers[name.strip()] = seq.strip().upper()

    # counters: for each primer, how often it hits the 5' end (as given) and the
    # 3' end (as reverse complement), summed over both read orientations.
    hits5 = {p: 0 for p in primers}
    hits3 = {p: 0 for p in primers}
    pos5 = {p: [] for p in primers}
    pos3 = {p: [] for p in primers}
    n = 0

    for seq in iter_reads(args.fastq, args.sample):
        n += 1
        for orient in (seq, revcomp(seq)):
            head = orient[:args.window]
            tail = orient[-args.window:]
            for name, pseq in primers.items():
                m = fuzzy_find(pseq, head, args.max_mismatch)
                if m:
                    hits5[name] += 1
                    pos5[name].append(m[0])
                rc = revcomp(pseq)
                m = fuzzy_find(rc, tail, args.max_mismatch)
                if m:
                    hits3[name] += 1
                    pos3[name].append(m[0] - args.window)  # negative = from 3' end

    if n == 0:
        sys.exit("no reads read")

    def pct(x):
        return 100.0 * x / n

    print(f"# scanned {n} reads | window {args.window} bp each end | "
          f"max_mismatch {args.max_mismatch}")
    print(f"# a read is counted in both orientations, so >100% is possible for a "
          f"primer that sits at BOTH ends\n")
    print(f"{'primer':14} {'5prime_end':>11} {'3prime_end(rc)':>15}   {'median_5pos':>11}")
    rows = sorted(primers, key=lambda p: -(hits5[p] + hits3[p]))
    for p in rows:
        med5 = (sorted(pos5[p])[len(pos5[p]) // 2] if pos5[p] else "-")
        print(f"{p:14} {pct(hits5[p]):9.1f}% {pct(hits3[p]):13.1f}%   {str(med5):>11}")

    # Verdict. Reads are scanned in both orientations, so a true forward primer
    # shows near-equal 5' and 3' rates (its rc appears at the 3' end of the
    # flipped copy) and likewise for the reverse primer. The amplicon is therefore
    # bounded by the two DISTINCT top-scoring primers, not one primer at both ends.
    # Cross-reactive near-duplicates (e.g. ITS4 vs ITS4ngsUni) are collapsed by
    # preferring, within a role, the highest scorer.
    print()
    total = {p: hits5[p] + hits3[p] for p in primers}
    ranked = [p for p in sorted(primers, key=lambda p: -total[p]) if pct(total[p]) >= 40]
    if not ranked:
        print("# no candidate primer explains >=40% of read ends — reads may already "
              "be primer-trimmed, or the primers are not in this panel "
              "(pass --primers name=SEQ,...).")
        return
    fwd = next((p for p in ranked if ROLE.get(p) == "fwd"), None)
    rev = next((p for p in ranked if ROLE.get(p) == "rev" and p != fwd), None)
    if fwd:
        print(f"# forward primer: {fwd} = {primers[fwd]} ({pct(hits5[fwd]):.0f}% of reads)")
    if rev:
        print(f"# reverse primer: {rev} = {primers[rev]} ({pct(hits3[rev]):.0f}% of reads, as reverse complement)")
    if fwd and rev:
        print(f"#   -> HiFiTaxa: --forward_primer {primers[fwd]} --reverse_primer {revcomp(primers[rev])}")
        print(f"#      (reverse_primer is given to cutadapt as the reverse complement of {rev})")
    else:
        print(f"# top primers present: {', '.join(ranked)} — assign fwd/rev by hand.")


if __name__ == "__main__":
    main()
