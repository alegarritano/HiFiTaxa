# Benchmark scoring and figure scripts

The R scripts in this folder score the four benchmark runs of the HiFiTaxa manuscript
and draw its figures: the 16S rRNA gene mock (ATCC MSA-1003, eight Kinnex Revio
libraries), the fungal ITS mock (103 species; Tedersoo et al. 2026, Mol Ecol Resour 26,
e70189), and the GTDB r232 and UNITE leave 10% out holdouts. `identity.pbs` and
`identity_resume.pbs` are the Gadi PBS jobs that measured the identity of every holdout
test sequence to its closest training sequence (Supplementary Table S12).

## Requirements

R 4.5 with ggplot2, dplyr, tidyr, stringr, readr, patchwork and scales (loaded by
`00_common.R`); `13_agreement_sankey.R` also needs ggalluvial, and `09_figure2.R`,
`08_holdout_retrieval.R` and `11_holdout_decomposition.R` use cowplot, ggplotify or
ggbreak. The published tables were produced with R 4.5.2, ggplot2 4.0.2, dplyr 1.2.0,
tidyr 1.3.2, stringr 1.6.0, readr 2.2.0, patchwork 1.3.2 and scales 1.4.0.

## Inputs

The scripts read the output directories of the four pipeline runs. They are not part of
the repository; `HIFITAXA_DATA` must point at a directory that holds them under these
names:

| directory | run | produced by |
|---|---|---|
| `results_atcc_16s` | 16S rRNA gene mock, HiFiTaxa v1.2.3 (commit 9d6dcf3), 5 to 7 October 2026 | `benchmark/pbs/run_atcc_16s.pbs` |
| `results_fungi_its` | fungal ITS mock, 17 July 2026 | `benchmark/pbs/run_fungi_its.pbs` |
| `benchmark_holdout_16s_seed42` | GTDB r232 holdout, 3 to 5 July 2026 | `benchmark/pbs/holdout.pbs` with `MARKER=16S` |
| `benchmark_holdout_its_seed42` | UNITE holdout, 17 July 2026 | `benchmark/pbs/holdout.pbs` with `MARKER=ITS` |

The three July runs were made from the development checkout between versions 1.1.0 and
1.2.0; the classifier code is the same as in v1.2.3 apart from commit 49aceef (the BLCA
identity threshold is passed to BLAST), which was in place for the two ITS runs. Each
mock directory must hold `dada2/feature_table.tsv`, `taxonomy_blca/`, `taxonomy_nb/` and
`taxonomy_emu/` (16S) or `taxonomy_emits/per_sample/` (ITS); each holdout directory
`test_10.taxonomy`, `reference_90.taxonomy`, `split_stats.json` and
`results/taxonomy_{blca,nb}/`. `16_holdout_identity_bins.R` also needs
`identity_16s.tsv` and `identity_its.tsv`, written by the identity job.

## Environment variables

| variable | meaning | default |
|---|---|---|
| `HIFITAXA_DATA` | directory holding the four run directories above | required |
| `HIFITAXA_OUT` | where `tables/` and `figures/` are written | `benchmark/analysis/regenerated/` (ignored by git) |
| `HIFITAXA_REPO` | a HiFiTaxa checkout, for `benchmark/mock_exclude/fungi_103_species.txt` | found automatically when run from `benchmark/analysis/` of a checkout; otherwise required |
| `IDENT_DIR` | folder holding `identity_16s.tsv` and `identity_its.tsv` (script 16 only) | this folder |
| `FP_THRESH` | false positive threshold as a fraction of the reads classified at that rank | `0.001` (0.1%) |

Every script locates `00_common.R` next to itself, so it can be run from any directory:

```
HIFITAXA_DATA=/path/to/run_dirs HIFITAXA_OUT=/path/to/out Rscript benchmark/analysis/01_holdout_accuracy.R
```

## Scripts, outputs and where they appear in the manuscript

Run the table scripts first (01, 02, 04, 05, 07, 14, 16, 17), then the figure scripts. All
outputs go to `HIFITAXA_OUT/tables/` and `HIFITAXA_OUT/figures/`.

| script | tables | figures | manuscript |
|---|---|---|---|
| `01_holdout_accuracy.R` | `holdout_accuracy.csv`, `holdout_orphan_backoff.csv`, `holdout_blca80_coverage_precision.csv`, `holdout_its_placeholder_denominators.csv` | `fig_holdout_accuracy` | Supplementary Tables S4 and S6; Table S5 is derived from `holdout_accuracy.csv` (coverage = 100 minus pct_unassigned; precision = correct / (correct + incorrect)) |
| `02_mock_detection.R` | `mock_detection.csv` | `fig_mock_detection` | Supplementary Table S1 |
| `03_clade_exclusion.R` | `clade_exclusion.csv` | `fig_clade_exclusion` | not in the manuscript (needs the clade exclusion runs of `benchmark/pbs/clade_exclusion.pbs`) |
| `04_mock_composition.R` | `atcc_composition_genus.csv` | `fig_atcc_deviation_heatmap` | Supplementary Table S3 (its composition error row is the sum of the absolute deviations in that file) and Supplementary Fig. S2 |
| `05_mock_readweighted.R` | `mock_readweighted_accuracy.csv` | | Supplementary Table S7 |
| `06_mock_taxaplots.R` | | `fig_mock_taxaplot_16s`, `fig_mock_taxaplot_its` | Supplementary Fig. S1 (16S) |
| `07_mock_retrieval.R` | `mock_retrieval.csv`; with `FP_THRESH=0` and `FP_THRESH=0.0001`, `mock_retrieval_threshold_any_read.csv` and `mock_retrieval_threshold_0.01pct.csv` | `fig_mock_retrieval_16s`, `fig_mock_retrieval_its` | Supplementary Tables S2 and S11 |
| `08_holdout_retrieval.R` | reads `holdout_accuracy.csv` and `holdout_orphan_backoff.csv` written by 01 | `fig_holdout_retrieval_16s`, `fig_holdout_retrieval_its` | working figures |
| `09_figure2.R` | `mock_detection_unassigned.csv` | `fig2_main` (drawn at 180 mm, 300 dpi; the PDF needs `pdftocairo` or Ghostscript to embed fonts) | Figure 2 |
| `10_holdout_taxaplot.R` | | `fig_holdout_taxaplot_16s` | Supplementary Fig. S3 |
| `11_holdout_decomposition.R` | reads `holdout_accuracy.csv` written by 01 | `fig_holdout_decomposition` | Supplementary Fig. S4 |
| `12_runtimes.R` | `classifier_runtimes.csv` | `fig_classifier_runtimes` | not in the manuscript (parses `pipeline_trace.txt` of each run) |
| `13_agreement_sankey.R` | | `fig_agreement_sankey` | working figure |
| `14_agreement_taxa.R` | `agreement_genus_16s_in_reference.csv`, `agreement_genus_16s_orphan.csv` | `fig_agreement_taxa_16s_in_reference`, `fig_agreement_taxa_16s_orphan` | Supplementary Tables S8 and S9; Supplementary Figs. S5 and S6 |
| `15_holdout_16s_genus_weighted.R` | `holdout_16s_genus_weighted.csv`, `holdout_16s_genus_weighted_ecoli.csv` | | not in the manuscript (a table of an earlier draft) |
| `16_holdout_identity_bins.R` | `holdout_accuracy_by_identity.csv` | | Supplementary Table S12 |
| `17_mock_per_library.R` | `atcc_composition_genus_per_library.csv` (needs `atcc_composition_genus.csv` from 04) | | Supplementary Table S13 |

Supplementary Table S10 (read tracking) is compiled from the pipeline's own read
tracking tables (`stats.tsv`, `primer_removal_stats.tsv`) and not from these scripts.
The scoring rules (name folding, the rule for crediting species calls on the 16S mock,
placeholder exclusion on the ITS holdout, the BLCA@80 mask) are in `00_common.R` and
described in the Supplementary Methods.

## Identity of the holdout test sequences to the training fraction

`identity.pbs` takes `test_ids_16s.txt` and `test_ids_its.txt` (the ids of the test
sequences, one per line, the first column of each `test_10.taxonomy`) and the
BLCA parsed reference FASTA files of GTDB r232 and UNITE. It first finds the test
sequences that have an identical training sequence by exact comparison, then runs
megablast for the rest (at most five hits per query, minimum identity 80%) and keeps the
best hit whose alignment covers at least 80% of the query length. `identity_resume.pbs`
does the same with the remaining queries split over 48 single thread blastn processes,
and resumes after a walltime kill. Both write `identity_16s.tsv` and `identity_its.tsv`
(`id`, `best_hit`, `pident`, `alnlen`, `qlen`), which `16_holdout_identity_bins.R` reads
from `IDENT_DIR`. The PBS headers and reference paths are those of the Gadi install used
for the manuscript (project is54); change them for another system.
