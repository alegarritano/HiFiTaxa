#!/usr/bin/env python3
"""Convert a UNITE general-release FASTA into the BLCA-parsed reference format
that the whole HiFiTaxa benchmark consumes (FASTA + taxonomy TSV).

This is the fungal/ITS analogue of bin/build_gtdb_blca_db.sh's parse step. The
downstream tooling (build_unite_dada2_db.sh, build_emits_db.sh, blca_main.py,
the scorers) was written for GTDB, where the top rank key is the literal string
"superkingdom". BLCA in particular HARD-CODES its seven accepted rank keys as
    superkingdom, phylum, class, order, family, genus, species
(blca_main.py: `levels = [...]`). So to reuse every downstream script verbatim
we map UNITE's kingdom (k__Fungi) onto the "superkingdom" key. It still holds the
fungal kingdom value ("Fungi"); only the *display* label changes to "Kingdom" in
the fungal scorers. Nothing else needs the rank renamed.

UNITE header format (general release "dynamic"):
    >SpeciesName|INSDC_acc|SH_id|type|k__Fungi;p__Ascomycota;c__..;o__..;f__..;g__..;s__Genus_epithet
We build:
    * a BLAST-safe, unique record id   = the SH id with '.' -> '_'  (e.g. SH1227328_10FU)
    * a BLCA-format lineage            = superkingdom:Fungi;phylum:..;..;species:Genus epithet;
      - the per-rank 'x__' prefix is stripped
      - the SPECIES value is de-underscored ("Abrothallus_subhalei" -> "Abrothallus subhalei")
        so it is a real binomial (genus + epithet), which is what dada2::assignTaxonomy and the
        GTDB-derived species-reference builder both expect.

Outputs (into --out-dir):
    unite_BLCAparsed.fasta       bare-id headers (+ a BLAST index built by the wrapper)
    unite_BLCAparsed.taxonomy    id <TAB> superkingdom:..;...;species:..;
    UNITE_PARSE_STATS.json       counts + diagnostics (incl. placeholder-species fraction)

Placeholder species: ~62% of UNITE species labels are not named to species
("Tomentella_sp", "Ascomycota_sp", "unidentified", "..._Incertae_sedis"). By
default we KEEP them verbatim (faithful to the database, mirroring how the GTDB
benchmark kept GTDB's placeholder genera). They are *flagged* in the stats and
written to placeholder_species_ids.txt so the scorer can optionally exclude them
from species-rank metrics (see --exclude-placeholder-species there). Use
--drop-placeholder-species here to instead blank the species rank at build time.
"""
import argparse
import json
import re
from pathlib import Path


# UNITE 'x__' rank prefix -> the BLCA/GTDB rank key. NOTE the kingdom (k__) maps
# to "superkingdom" on purpose (see module docstring / blca_main.py levels).
PREFIX2KEY = {
    "k__": "superkingdom",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}
RANK_ORDER = ["superkingdom", "phylum", "class", "order",
              "family", "genus", "species"]

# Epithet tokens that mark "not actually identified to species".
_PLACEHOLDER_EPITHETS = {"sp", "sp.", "spp", "spp.", "unidentified",
                         "uncultured", "incertae", "indet", "indet."}


def is_placeholder_species(species_value, genus_value=""):
    """True if a (de-underscored, binomial) species string is a UNITE
    placeholder rather than a named species. Examples that are placeholders:
    'Tomentella sp', 'Ascomycota sp', 'Fungi sp', 'unidentified',
    'Helotiales Incertae sedis'. 'Abrothallus subhalei' is NOT a placeholder."""
    s = species_value.strip()
    if not s:
        return True
    low = s.lower()
    if "incertae" in low or "unidentified" in low or "uncultured" in low:
        return True
    toks = s.split()
    if len(toks) < 2:           # single word -> not a binomial
        return True
    epithet = toks[1].lower().rstrip(".")
    if epithet in {"sp", "spp", "indet"}:
        return True
    # genus repeated as epithet placeholder e.g. "Russula sp" already caught;
    # also catch "<rank-name> sp" where token0 is itself a higher-rank label.
    return False


def parse_header(header):
    """Return (record_id, {rank_key: value}) for one UNITE header line (no '>').
    record_id is BLAST-safe and unique. Missing/empty ranks are simply absent."""
    # Split the pipe-delimited prefix from the trailing 'k__..;..;s__..' lineage.
    # The lineage always starts at the first 'k__' token.
    m = re.search(r"k__", header)
    if m is None:
        return None, None
    meta = header[:m.start()].rstrip("|").split("|")
    lineage_str = header[m.start():].strip()

    # record id: prefer the SH id (meta[2]); it is unique across the release.
    # Fall back to the INSDC accession, then the whole meta joined.
    sh_id = meta[2] if len(meta) >= 3 else ""
    insdc = meta[1] if len(meta) >= 2 else ""
    raw_id = sh_id or insdc or "|".join(meta)
    # BLAST-safe: keep alnum, replace everything else (incl '.') with '_'.
    rec_id = re.sub(r"[^A-Za-z0-9]+", "_", raw_id).strip("_")

    ranks = {}
    for tok in lineage_str.split(";"):
        tok = tok.strip()
        if len(tok) < 3 or tok[:3] not in PREFIX2KEY:
            continue
        key = PREFIX2KEY[tok[:3]]
        val = tok[3:].strip()
        if not val:
            continue
        if key == "species":
            # UNITE species = "Genus_epithet" -> "Genus epithet" (real binomial).
            val = val.replace("_", " ").strip()
        ranks[key] = val
    return rec_id, ranks


def lineage_to_blca(ranks, drop_placeholder_species=False):
    """Render the ordered 'superkingdom:..;...;species:..;' string. Trailing
    empty ranks are dropped. Returns ('', is_placeholder_flag)."""
    out = []
    is_ph = False
    for key in RANK_ORDER:
        val = ranks.get(key, "")
        if key == "species" and val:
            is_ph = is_placeholder_species(val, ranks.get("genus", ""))
            if is_ph and drop_placeholder_species:
                val = ""
        if val:
            out.append(f"{key}:{val}")
    return ";".join(out) + (";" if out else ""), is_ph


def iter_fasta(path):
    rec_id, lines = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if rec_id is not None:
                    yield rec_id, lines
                rec_id = line[1:].rstrip("\n")
                lines = []
            else:
                lines.append(line.strip())
        if rec_id is not None:
            yield rec_id, lines


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-fasta", required=True,
                    help="UNITE general-release dynamic FASTA")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out-prefix", default="unite_BLCAparsed",
                    help="basename for the .fasta / .taxonomy outputs")
    ap.add_argument("--min-length", type=int, default=0,
                    help="drop reference sequences shorter than this (default 0 "
                         "= keep all; ITS is short, unlike 16S, so do NOT use the "
                         "1000 bp 16S threshold here)")
    ap.add_argument("--drop-placeholder-species", action="store_true",
                    help="blank the species rank for UNITE placeholders "
                         "('X sp', 'unidentified', 'Incertae sedis') at build time")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fa_out = out_dir / f"{args.out_prefix}.fasta"
    tx_out = out_dir / f"{args.out_prefix}.taxonomy"
    ph_out = out_dir / "placeholder_species_ids.txt"

    n_total = n_kept = n_no_lineage = n_dup = 0
    n_placeholder = 0
    seen = set()
    kingdoms = {}
    placeholder_ids = []

    with open(fa_out, "w") as fa, open(tx_out, "w") as tx:
        for header, seq_lines in iter_fasta(args.in_fasta):
            n_total += 1
            seq = "".join(seq_lines)
            if args.min_length > 0 and len(seq) < args.min_length:
                continue
            rec_id, ranks = parse_header(header)
            if rec_id is None or not ranks:
                n_no_lineage += 1
                continue
            if rec_id in seen:
                n_dup += 1
                # disambiguate rather than silently overwrite
                rec_id = f"{rec_id}_{n_dup}"
            seen.add(rec_id)

            kingdoms[ranks.get("superkingdom", "")] = \
                kingdoms.get(ranks.get("superkingdom", ""), 0) + 1

            lineage, is_ph = lineage_to_blca(
                ranks, drop_placeholder_species=args.drop_placeholder_species)
            if is_ph:
                n_placeholder += 1
                placeholder_ids.append(rec_id)

            fa.write(f">{rec_id}\n{seq}\n")
            tx.write(f"{rec_id}\t{lineage}\n")
            n_kept += 1

    ph_out.write_text("\n".join(placeholder_ids) + ("\n" if placeholder_ids else ""))

    stats = {
        "in_fasta": args.in_fasta,
        "min_length": args.min_length,
        "drop_placeholder_species": args.drop_placeholder_species,
        "n_total_records": n_total,
        "n_written": n_kept,
        "n_skipped_no_lineage": n_no_lineage,
        "n_duplicate_ids_disambiguated": n_dup,
        "kingdom_values": kingdoms,
        "n_placeholder_species": n_placeholder,
        "placeholder_species_fraction": round(n_placeholder / max(1, n_kept), 4),
        "rank_key_note": ("kingdom (k__) is stored under the 'superkingdom' key so "
                          "BLCA's hard-coded rank list accepts it; display label is "
                          "'Kingdom' in the fungal scorers."),
        "placeholder_note": ("placeholder_species_ids.txt lists records whose species "
                             "is a UNITE placeholder (e.g. 'Tomentella sp'); the scorer "
                             "can optionally exclude these from species-rank metrics."),
    }
    (out_dir / "UNITE_PARSE_STATS.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2))
    print(f"\n[parse] wrote {fa_out}")
    print(f"[parse] wrote {tx_out}")
    print(f"[parse] wrote {ph_out}  ({n_placeholder:,} placeholder-species records)")


if __name__ == "__main__":
    main()
