#!/usr/bin/env python3
"""Compare BLCA vs Emu species composition as side-by-side stacked bars.

Inputs (assumes the standard HiFiTaxa <outdir>/ layout):
  <outdir>/dada2/dada2-ccs_table_filtered.qza   (BLCA feature table)
  <outdir>/taxonomy_blca/blca_taxonomy_table.csv (BLCA ASV -> taxonomy)
  <outdir>/taxonomy_emu/emu_species_table.tsv    (Emu species rel. abundance)

Outputs:
  <out_dir>/blca_vs_emu_species_taxaplot.png     (the plot)
  <out_dir>/blca_vs_emu_species_taxaplot.pdf     (vector copy)
  <out_dir>/blca_vs_emu_species_relabund.tsv     (joined table, for the record)

Run with the qiime2 conda env (gives us biom + pandas + matplotlib + h5py).

Usage:
    python bin/compare_blca_emu_taxaplot.py <pipeline_outdir> <out_dir> [--top N]
"""
import argparse
import os
import re
import sys
import tempfile
import zipfile

import biom
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ATCC MSA-2006-style 20-strain mock community (expected proportions, fractions
# summing to 1.0). Triggered by --truth atcc on the CLI.
ATCC_TRUTH = {
    "Acinetobacter baumannii":      0.0018,
    "Bacillus pacificus":           0.0180,
    "Phocaeicola vulgatus":         0.0002,
    "Bifidobacterium adolescentis": 0.0002,
    "Clostridium beijerinckii":     0.0180,
    "Cutibacterium acnes":          0.0018,
    "Deinococcus radiodurans":      0.0002,
    "Enterococcus faecalis":        0.0002,
    "Escherichia coli":             0.1800,
    "Helicobacter pylori":          0.0018,
    "Lactobacillus gasseri":        0.0018,
    "Neisseria meningitidis":       0.0018,
    "Porphyromonas gingivalis":     0.1800,
    "Pseudomonas paraeruginosa":    0.0180,
    "Cereibacter sphaeroides":      0.1800,
    "Schaalia odontolytica":        0.0002,
    "Staphylococcus aureus":        0.0180,
    "Staphylococcus epidermidis":   0.1800,
    "Streptococcus agalactiae":     0.0180,
    "Streptococcus mutans":         0.1800,
}


def normalise_species(name):
    """Reduce GTDB-style 'Genus_A species' / 'Genus_A' to 'Genus species' /
    'Genus' for matching. Handles both standalone genus labels (e.g. when
    aggregating BLCA tax to the Genus rank) and full 'Genus species' labels.
    Leaves the original label intact if no '_X' suffix is present."""
    if not isinstance(name, str):
        return name
    return re.sub(r"^([A-Z][a-zA-Z]+)_[A-Z](\s|$)", r"\1\2", name)


# GTDB-vs-NCBI species name swaps for the SAME organism (by ATCC reference).
# Two tiers:
#
# ALWAYS_APPLY_ALIASES — silently applied in every figure and every
# comparison. These map GTDB **placeholder / unnamed-MAG accessions** onto
# their NCBI species. Placeholders carry no biological information that a
# reader could decode from the label itself ('G047199095 sp047199095' is just
# a UUID), so showing them in a published figure helps no one. They are NOT
# silently applied in the underlying data — only at presentation time.
#
# GTDB_TO_NCBI_ALIAS — opt-in (--apply-aliases). These map GTDB **named**
# species onto their NCBI synonyms. They are off by default because the
# difference between, say, 'Escherichia boydii' (GTDB) and 'Escherichia coli'
# (NCBI) is itself a real and reportable finding for any closed-reference
# benchmark against an NCBI-defined truth.
#
# Direction: `GTDB label` -> `NCBI / ATCC label`.

ALWAYS_APPLY_ALIASES = {
    # GTDB placeholder genera + species for genomes that NCBI calls Escherichia
    # coli but GTDB's 95% ANI species threshold splits out into "not ideal"
    # nameless clusters (acknowledged by the GTDB developers themselves,
    # https://forum.gtdb.ecogenomic.org/t/newly-formed-ecma0423-out-of-escherichia-coli-and-representative-choice/827).
    #
    # - G047199095 sp047199095 : an Escherichia MAG (e.g. GCA_000226585.1,
    #   "E. coli XH140A" by NCBI) that falls outside the named Escherichia
    #   clade at GTDB's ANI threshold.
    # - ECMA0423 sp047199055   : the ~13,000-genome split-out from the
    #   GTDB r232 Escherichia coli cluster (reference GCA_047199055.1 hits
    #   94.9% ANI to the E. coli type strain, just under the 95% cutoff).
    #
    # Species-level aliases (full Genus + species labels):
    "G047199095 sp047199095": "Escherichia coli",
    "ECMA0423 sp047199055":   "Escherichia coli",
    # Genus-level aliases (used when BLCA / NB aggregate to the Genus rank,
    # at which point the placeholder genus appears on its own):
    "G047199095":             "Escherichia",
    "ECMA0423":               "Escherichia",
}

GTDB_TO_NCBI_ALIAS = {
    # Escherichia coli aliases (NCBI: E. coli):
    #   - 'Escherichia boydii' — GTDB merges Shigella into Escherichia and
    #     re-splits the E. coli tree by ANI; many NCBI E. coli strains
    #     (incl. the ATCC E. coli reference) sit under E. boydii in GTDB.
    "Escherichia boydii":         "Escherichia coli",
    # Pseudomonas paraeruginosa aliases:
    #   - 'Pseudomonas aeruginosa' — GTDB / NCBI rename for the ATCC 9027 reference.
    "Pseudomonas aeruginosa":     "Pseudomonas paraeruginosa",
}


def apply_aliases(name, use_aliases=False):
    """Apply normalise_species + ALWAYS_APPLY_ALIASES (placeholder MAGs);
    additionally apply the opt-in GTDB↔NCBI named-species synonym map when
    use_aliases=True."""
    n = normalise_species(name)
    n = ALWAYS_APPLY_ALIASES.get(n, n)
    if use_aliases:
        n = GTDB_TO_NCBI_ALIAS.get(n, n)
    return n


# ---------- I/O helpers --------------------------------------------------

def extract_biom_from_qza(qza_path):
    """Pull feature-table.biom out of the .qza (which is just a zip)."""
    tmpdir = tempfile.mkdtemp(prefix="hifiblca_biom_")
    with zipfile.ZipFile(qza_path) as zf:
        member = next((n for n in zf.namelist() if n.endswith("feature-table.biom")),
                      None)
        if not member:
            sys.exit(f"ERROR: no feature-table.biom inside {qza_path}")
        zf.extract(member, tmpdir)
    return os.path.join(tmpdir, member)


def load_blca_relabund_by_species(qza_path, tax_csv):
    """Join the per-ASV BLCA taxonomy with the feature table, aggregate counts
    by species, then divide by per-sample totals to get relative abundance."""
    biom_path = extract_biom_from_qza(qza_path)
    table = biom.load_table(biom_path)
    df = table.to_dataframe(dense=True)         # rows = ASV ids, cols = samples
    df.index.name = "Feature ID"

    tax = pd.read_csv(tax_csv)
    tax = tax[["Feature ID", "Genus", "Species"]].fillna("Unassigned")
    merged = df.reset_index().merge(tax, on="Feature ID", how="left")
    merged["Species"] = merged["Species"].fillna("Unassigned")
    # Some Species entries are empty strings -> tag them too
    merged.loc[merged["Species"].str.strip() == "", "Species"] = "Unassigned"

    sample_cols = [c for c in df.columns]
    species_counts = merged.groupby("Species")[sample_cols].sum()
    totals = species_counts.sum(axis=0).replace(0, np.nan)
    species_relabund = species_counts.div(totals, axis=1).fillna(0.0)
    return species_relabund


def load_emu_relabund(tsv_path):
    """Emu species table: first 7 cols are taxonomy ranks, remaining cols are
    per-sample relative abundance. Strip the `species:` prefix."""
    df = pd.read_csv(tsv_path, sep="\t")
    tax_cols = ["species", "genus", "family", "order", "class", "phylum",
                "superkingdom"]
    sample_cols = [c for c in df.columns if c not in tax_cols]
    out = df.set_index("species")[sample_cols].copy()
    out.index = out.index.str.replace(r"^species:", "", regex=True)
    out.index.name = "Species"
    # Drop rows where all values are 0 (Emu sometimes emits empty trailing rows)
    out = out.loc[out.sum(axis=1) > 0]
    return out


# ---------- Plot ---------------------------------------------------------

def build_species_order(blca, emu, truth):
    """Union of all species across BLCA, Emu, and (optionally) the truth.

    Sort: ATCC species first in descending expected abundance, then any
    extra species detected by BLCA or Emu in descending max mean abundance.
    """
    union = set(blca.index) | set(emu.index)
    if truth is not None:
        union |= set(truth.index)

    truth_order = (truth.sort_values(ascending=False).index.tolist()
                   if truth is not None else [])
    truth_set = set(truth_order)

    extras = sorted(union - truth_set,
                    key=lambda s: -max(blca.get(s, pd.Series([0])).mean()
                                       if s in blca.index else 0,
                                       emu.get(s, pd.Series([0])).mean()
                                       if s in emu.index else 0))
    return [s for s in truth_order if s in union] + extras


def build_palette(n):
    """Distinct, high-contrast palette for n species. Combines tab20 + tab20b +
    tab20c (60 colours total) and falls back to hsv beyond that."""
    pools = []
    for name in ("tab20", "tab20b", "tab20c"):
        pools += list(plt.get_cmap(name).colors)
    if n <= len(pools):
        return pools[:n]
    hsv = plt.get_cmap("hsv")
    extra = [hsv(i / max(1, n - len(pools))) for i in range(n - len(pools))]
    return pools + extra


def _add_method_bar(ax, x_pos, df, sample, species_order, colours, width,
                    method_label):
    """Stack one bar (one method × one sample) at x_pos."""
    bottom = 0.0
    for sp, col in zip(species_order, colours):
        v = float(df.at[sp, sample]) if sp in df.index and sample in df.columns else 0.0
        if v == 0.0:
            continue
        ax.bar(x_pos, v, width=width, bottom=bottom,
               color=col, edgecolor="white", linewidth=0.3)
        bottom += v
    ax.text(x_pos, -0.022, method_label, ha="center", va="top",
            fontsize=7, color="#444", rotation=90)


def make_plot(blca, emu, truth, out_png, out_pdf, top=None, use_aliases=False):
    """Plot per-sample {BLCA, Emu[, Expected]} stacked bars.

    `top`: if a positive int, restrict to top-N species and pool the rest as
    'Other'. If None or 0, show every detected species (no 'Other').
    `truth`: pandas Series of expected proportions (per species). When given,
    one extra 'Expected' bar is added at the right with the truth composition.
    `use_aliases`: when True, collapse known GTDB↔NCBI synonyms (E. boydii ->
    E. coli, etc.) before plotting. Off by default so the discrepancy stays
    visible in the figure (you can describe it in the manuscript text).
    """
    common_samples = [s for s in emu.columns if s in blca.columns]
    if not common_samples:
        sys.exit("ERROR: no overlapping samples between BLCA and Emu outputs")

    # Always normalise GTDB '_A' suffixes (Genus_A species -> Genus species)
    # so the plain GTDB ↔ NCBI base names line up. Optional GTDB↔NCBI synonym
    # collapsing only happens when use_aliases=True.
    blca = blca.copy(); blca.index = [apply_aliases(s, use_aliases) for s in blca.index]
    emu  = emu.copy();  emu.index  = [apply_aliases(s, use_aliases) for s in emu.index]
    # Merge any rows that collapsed onto the same name after aliasing.
    blca = blca.groupby(level=0).sum()
    emu  = emu.groupby(level=0).sum()

    species_order_full = build_species_order(blca, emu, truth)

    if top:  # legacy top-N + Other mode
        species_order = species_order_full[:top] + ["Other"]
        def _restrict(df):
            kept = df.reindex(species_order_full[:top]).fillna(0.0)
            other = df.loc[~df.index.isin(species_order_full[:top])].sum(axis=0).to_frame("Other").T
            return pd.concat([kept, other])
        blca_p = _restrict(blca[common_samples])
        emu_p  = _restrict(emu[common_samples])
        if truth is not None:
            t_kept = truth.reindex(species_order_full[:top]).fillna(0.0)
            t_other = truth.loc[~truth.index.isin(species_order_full[:top])].sum()
            truth_p = pd.concat([t_kept, pd.Series({"Other": t_other})])
        else:
            truth_p = None
    else:    # show ALL species
        species_order = species_order_full
        blca_p = blca[common_samples].reindex(species_order).fillna(0.0)
        emu_p  = emu[common_samples].reindex(species_order).fillna(0.0)
        truth_p = (truth.reindex(species_order).fillna(0.0)
                   if truth is not None else None)

    colours = build_palette(len(species_order))

    n_methods = 3 if truth_p is not None else 2
    n_samples = len(common_samples)
    bar_w = 0.26 if n_methods == 3 else 0.36
    # one cluster per sample; methods sit inside the cluster
    cluster_x = np.arange(n_samples) * (bar_w * (n_methods + 0.6))
    offsets = {"BLCA": -bar_w, "Emu": 0.0, "Expected": +bar_w} if n_methods == 3 \
              else {"BLCA": -bar_w / 2 - 0.02, "Emu": +bar_w / 2 + 0.02}

    fig_w = max(11, 0.9 * n_samples * n_methods)
    fig_h = max(7, 0.18 * len(species_order))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for i, sample in enumerate(common_samples):
        _add_method_bar(ax, cluster_x[i] + offsets["BLCA"],
                        blca_p, sample, species_order, colours, bar_w, "BLCA")
        _add_method_bar(ax, cluster_x[i] + offsets["Emu"],
                        emu_p,  sample, species_order, colours, bar_w, "Emu")
        if truth_p is not None:
            # Treat the Expected truth as identical for every sample.
            truth_df = pd.DataFrame({sample: truth_p})
            _add_method_bar(ax, cluster_x[i] + offsets["Expected"],
                            truth_df, sample, species_order, colours, bar_w,
                            "Expected")

    ax.set_xticks(cluster_x)
    ax.set_xticklabels(common_samples, rotation=30, ha="right")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Relative abundance")
    title_suffix = (" (vs ATCC expected)" if truth_p is not None else "")
    if top:
        ax.set_title(f"BLCA vs Emu species composition — top {top} + Other{title_suffix}")
    else:
        ax.set_title(f"BLCA vs Emu species composition — all detected{title_suffix}")
    ax.spines[["top", "right"]].set_visible(False)

    # Legend with all species. Use 1–2 columns depending on count.
    legend_handles = [matplotlib.patches.Patch(color=c, label=s)
                      for s, c in zip(species_order, colours)]
    ncol = 1 if len(species_order) <= 18 else 2
    ax.legend(handles=legend_handles,
              loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7,
              title="Species", frameon=False, ncol=ncol)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {out_png}")
    print(f"[plot] wrote {out_pdf}")

    # joined TSV: BLCA + Emu + (optionally) Expected, all species
    blca_p_n = blca_p.copy(); blca_p_n.columns = [f"BLCA::{s}" for s in blca_p_n.columns]
    emu_p_n  = emu_p.copy();  emu_p_n.columns  = [f"Emu::{s}"  for s in emu_p_n.columns]
    pieces = [blca_p_n, emu_p_n]
    if truth_p is not None:
        pieces.append(truth_p.to_frame("Expected"))
    joined = pd.concat(pieces, axis=1)
    return joined


# ---------- Main ---------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir",  help="pipeline output dir (must contain dada2/, taxonomy_blca/, taxonomy_emu/)")
    ap.add_argument("dest",    help="where to write the plot + joined TSV")
    ap.add_argument("--top", type=int, default=0,
                    help="top N species to show; rest pooled as 'Other'. "
                         "0 (default) shows every detected species.")
    ap.add_argument("--truth", choices=("atcc",), default=None,
                    help="overlay the ground-truth composition. 'atcc' uses the bundled "
                         "20-strain ATCC MSA-2006-style mock community.")
    ap.add_argument("--apply-aliases", action="store_true",
                    help="collapse known GTDB↔NCBI species-name swaps (E. boydii / "
                         "G047199095 sp047199095 -> E. coli; P. aeruginosa -> "
                         "P. paraeruginosa). OFF by default — for mock-community work "
                         "the discrepancy is a finding to report, not to hide.")
    args = ap.parse_args()

    outdir = os.path.abspath(args.outdir)
    dest   = os.path.abspath(args.dest)
    os.makedirs(dest, exist_ok=True)

    blca = load_blca_relabund_by_species(
        os.path.join(outdir, "dada2", "dada2-ccs_table_filtered.qza"),
        os.path.join(outdir, "taxonomy_blca", "blca_taxonomy_table.csv"))
    emu  = load_emu_relabund(
        os.path.join(outdir, "taxonomy_emu", "emu_species_table.tsv"))

    truth = None
    if args.truth == "atcc":
        truth = pd.Series(ATCC_TRUTH, name="Expected").sort_values(ascending=False)
        # truth uses NCBI labels already; only normalise '_A' suffixes (no-op here).
        truth.index = [normalise_species(s) for s in truth.index]

    print(f"[load] BLCA: {blca.shape[0]} species across {blca.shape[1]} samples")
    print(f"[load] Emu : {emu.shape[0]} species across {emu.shape[1]} samples")
    if truth is not None:
        print(f"[load] ATCC truth: {truth.shape[0]} expected species")
    if args.apply_aliases:
        print(f"[load] applying GTDB↔NCBI alias map ({len(GTDB_TO_NCBI_ALIAS)} entries)")

    suffix = "all" if not args.top else f"top{args.top}"
    if truth is not None:
        suffix += "_vs_atcc"
    if args.apply_aliases:
        suffix += "_aliased"
    out_png = os.path.join(dest, f"blca_vs_emu_species_taxaplot_{suffix}.png")
    out_pdf = os.path.join(dest, f"blca_vs_emu_species_taxaplot_{suffix}.pdf")
    joined = make_plot(blca, emu, truth, out_png, out_pdf,
                       top=args.top, use_aliases=args.apply_aliases)

    out_tsv = os.path.join(dest, f"blca_vs_emu_species_relabund_{suffix}.tsv")
    joined.to_csv(out_tsv, sep="\t")
    print(f"[plot] wrote {out_tsv}")


if __name__ == "__main__":
    main()
