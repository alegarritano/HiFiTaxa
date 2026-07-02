#!/usr/bin/env python3
"""Remove a set of species from a reference FASTA — the "clade exclusion" step of
the mock-vs-depleted-DB benchmark. Given a list of species known to be in a mock
community, drop every reference sequence assigned to those species so the mock
reads must be classified WITHOUT their true species in the database (testing
graceful genus-level fallback).

Handles the two reference formats HiFiTaxa uses:

  1. BLCA-parsed (GTDB or UNITE): pass --in-taxonomy. Species come from the
     `species:` field of the taxonomy TSV, matched to the FASTA by accession
     (first whitespace token of the header). Writes a depleted FASTA and, if
     --out-taxonomy is given, the matching depleted taxonomy. Feed the outputs
     to makeblastdb (BLCA) and to build_gtdb_dada2_db.sh / build_unite_dada2_db.sh
     (NB) / build_gtdb_emu_db.sh (Emu).

  2. Lineage-in-header (the EMITS minimap2 target unite.fasta): omit
     --in-taxonomy. Species = the last ';'-delimited field of the header
     (e.g. '>Fungi;...;Abrothallus subhalei;'). Writes a depleted FASTA.

Species matching is case-insensitive and whitespace-trimmed. --exclude-species
is a text file, one species name per line ('#' comments and blanks ignored); or
pass --species-col to pull the names from a column of a TSV truth table.

Usage:
  # BLCA-parsed (GTDB) — for BLCA/NB/Emu:
  python deplete_reference.py --in-fasta db/gtdb_ssu_BLCAparsed.fasta \\
      --in-taxonomy db/gtdb_ssu_BLCAparsed.taxonomy \\
      --exclude-species atcc_species.txt \\
      --out-fasta DEP/gtdb_ssu_BLCAparsed.fasta \\
      --out-taxonomy DEP/gtdb_ssu_BLCAparsed.taxonomy

  # UNITE EMITS target (lineage in header):
  python deplete_reference.py --in-fasta db_unite/unite.fasta \\
      --exclude-species fungi_species.txt --out-fasta DEP/unite.fasta

  # pull species from a truth-table column instead of a plain list:
  python deplete_reference.py ... --exclude-species truth.tsv --species-col species
"""
import argparse
import sys
from pathlib import Path


def norm(s):
    return s.strip().lower()


def read_fasta_records(path):
    rec_id, header, lines = None, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if rec_id is not None:
                    yield rec_id, header, lines
                header = line[1:].rstrip("\n")
                rec_id = header.strip().split()[0] if header.strip() else ""
                lines = []
            else:
                lines.append(line.rstrip("\n"))
        if rec_id is not None:
            yield rec_id, header, lines


def species_from_taxonomy_field(lineage):
    for part in lineage.rstrip(";").split(";"):
        if part.startswith("species:"):
            return part.split(":", 1)[1].strip()
    return ""


def species_from_header(header):
    # UNITE-style: '>Fungi;...;Genus species;' -> last non-empty ';' field
    parts = [p for p in header.rstrip(";").split(";") if p.strip()]
    return parts[-1].strip() if parts else ""


def load_exclude(path, species_col):
    names = set()
    with open(path) as fh:
        if species_col is not None:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                ci = header.index(species_col)
            except ValueError:
                sys.exit(f"--species-col '{species_col}' not in header: {header}")
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) > ci and cols[ci].strip():
                    names.add(norm(cols[ci]))
        else:
            for line in fh:
                s = line.strip()
                if s and not s.startswith("#"):
                    names.add(norm(s))
    if not names:
        sys.exit(f"no species names read from {path}")
    return names


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-fasta", required=True)
    ap.add_argument("--in-taxonomy", help="BLCA-parsed taxonomy TSV; omit for lineage-in-header FASTAs")
    ap.add_argument("--exclude-species", required=True,
                    help="species to remove: one per line, or a TSV column via --species-col")
    ap.add_argument("--species-col", help="if --exclude-species is a TSV, the species column name")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-taxonomy", help="write depleted taxonomy (BLCA-parsed mode only)")
    args = ap.parse_args()

    exclude = load_exclude(args.exclude_species, args.species_col)
    print(f"[deplete] excluding {len(exclude)} species")

    tax = {}
    if args.in_taxonomy:
        with open(args.in_taxonomy) as fh:
            for line in fh:
                if "\t" not in line:
                    continue
                acc, lineage = line.rstrip("\n").split("\t", 1)
                tax[acc] = lineage

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    n_total = n_kept = n_dropped = 0
    matched = set()
    kept_acc = set()

    with open(args.out_fasta, "w") as out:
        for acc, header, lines in read_fasta_records(args.in_fasta):
            n_total += 1
            if args.in_taxonomy:
                sp = species_from_taxonomy_field(tax.get(acc, ""))
            else:
                sp = species_from_header(header)
            if norm(sp) in exclude:
                n_dropped += 1
                matched.add(norm(sp))
                continue
            out.write(f">{header}\n")
            for s in lines:
                out.write(s + "\n")
            n_kept += 1
            kept_acc.add(acc)

    if args.out_taxonomy and args.in_taxonomy:
        Path(args.out_taxonomy).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_taxonomy, "w") as out:
            for acc, lineage in tax.items():
                if acc in kept_acc:
                    out.write(f"{acc}\t{lineage}\n")

    print(f"[deplete] {n_total:,} records -> kept {n_kept:,}, dropped {n_dropped:,}")
    print(f"[deplete] {len(matched)}/{len(exclude)} excluded species matched >=1 reference record")
    missing = sorted(exclude - matched)
    if missing:
        print(f"[deplete] WARNING: {len(missing)} excluded species not found in the reference "
              f"(name mismatch or already absent), e.g.: {missing[:5]}")
    print(f"[deplete] wrote {args.out_fasta}"
          + (f" + {args.out_taxonomy}" if (args.out_taxonomy and args.in_taxonomy) else ""))


if __name__ == "__main__":
    main()
