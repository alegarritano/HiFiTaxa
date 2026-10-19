# HiFiTaxa Gadi test harness

Three benchmarks, run as non-interactive Gadi PBS jobs. Assumes the databases,
container images and conda envs are already staged on Gadi as in
[../docs/gadi.md](../docs/gadi.md) (steps 1–5b), and that the reads (ATCC 16S,
103-species fungal ITS) live under `/scratch/<proj>/$USER/...`.

| # | Benchmark | Query | Reference | Classifiers | Score at |
|---|---|---|---|---|---|
| 1 | Mock vs full DB | mock reads | full DB | BLCA, NB, Emu/EMITS | species (recovery) |
| 2 | Mock vs clade-excluded DB | mock reads | DB with mock species removed | BLCA, NB, Emu/EMITS | genus (fallback) |
| 3 | Leave-10%-out holdout | removed 10% of reference (as ASVs) | 90% train | BLCA, NB only | species/genus |

Experiment 3 is BLCA+NB only: its queries are reference *sequences*, not reads, so
the read-level profilers (Emu/EMITS) don't apply. Experiments 1–2 cover all four.

Every script has a `CONFIG (edit me)` block at the top. Set `PROJ` to your NCI
project code and the read/repo paths; each var is also overridable at submit time
(`qsub -v READS_DIR=/path,PROJ=ab12 ...`).

## Data provenance

- **ATCC MSA-1003 16S mock** — the eight staggered mix libraries of PacBio's public
  Kinnex 16S Revio data set
  <https://downloads.pacbcloud.com/public/dataset/Kinnex-16S/DATA-ZymoATCCMix-32plex-Revio/>
  (`2-Demux/fastq_files`; forward barcodes F01 to F08 paired with reverse barcodes R37
  to R44, listed as the ATCC staggered mix in `sample_primer_pairs_key.txt`). All reads
  of each library (2.9 to 4.0 million) were used.
- **103-species fungal mock** — PacBio Revio HiFi, full-ITS (ITS9MUNngs / ITS4ngsUni),
  ENA **PRJEB108994**. Reference:
  > Tedersoo L, *et al.* (2026). *Benchmarking full-length ITS metabarcoding across
  > Illumina 2 × 500, PacBio, and Oxford Nanopore sequencing using mock and soil
  > communities.* Molecular Ecology Resources 26, e70189. doi:10.1111/1755-0998.70189
  > Code: <https://github.com/Mycology-Microbiology-Center/fullITS-multiplatform-eval>

## 1. Find the primers used (`scan_primers.py`)

The holdout benchmark uses no reads or PCR, but the real-read runs need the
amplicon primers. Scan one read file to detect them (IUPAC-aware, both
orientations, stdlib-only):

```bash
python benchmark/scan_primers.py /scratch/<proj>/$USER/reads_fungi_its/<one>.fastq.gz
```

It prints a per-primer hit table and a suggested `--forward_primer` /
`--reverse_primer`. On the bundled ITS mock it reports **1391F**
(`GTACACACCGCCCGTC`) + **ITS4ngsUni** (rev-comp `GCATATHANTAAGSGSAGG`) — the
pipeline's ITS defaults. If your Gadi fungal reads are the ~1.8 kb ITS–28S
amplicon they will instead show **ITS1catta / LR5_TW14ngs** (see
[../docs/PRIMERS.md](../docs/PRIMERS.md)); set those in `run_fungi_its.pbs`.

## Experiment 1 — mock vs full DB (`run_*.pbs`)

```bash
qsub -v READS_DIR=/scratch/<proj>/$USER/ATCC_bacteria benchmark/pbs/run_atcc_16s.pbs
qsub -v READS_DIR=/scratch/<proj>/$USER/fungi_103      benchmark/pbs/run_fungi_its.pbs
```

Each builds `samples.tsv` + `metadata.tsv` from its `READS_DIR`, then runs
`bin/run_pipeline.py --classifier all` offline against the full database. The 16S job
disables both ASV filters (`--min_asv_totalfreq 0 --min_asv_sample 0`) so that the
0.02% tier of the staggered mock is kept, and uses `--emu_threads 10`; the ITS job uses
`--min_asv_totalfreq 1` with the default `--min_asv_sample 1`. These are the settings of
the manuscript runs.

## Experiment 2 — mock vs clade-excluded DB (`clade_exclusion.pbs`)

Removes the mock's known species from the reference, **rebuilds every classifier
DB from the depleted source**, then classifies the mock reads. All classifiers see
the same depleted DB (so they're comparable); score at **genus** since species was
removed on purpose. Needs a species list / truth table per mock.

Curated exclusion lists are committed under `benchmark/mock_exclude/` (derived
from the ATCC MSA-1003 composition and Table S1, with the manuscript's synonym
folding baked in — old+new names, and the Cortinarius segregates that UNITE uses):

```bash
qsub -v MARKER=16S,EXCLUDE_LIST=benchmark/mock_exclude/atcc_msa1003_species.txt,READS_DIR=/scratch/<proj>/$USER/ATCC_bacteria benchmark/pbs/clade_exclusion.pbs
qsub -v MARKER=ITS,EXCLUDE_LIST=benchmark/mock_exclude/fungi_103_species.txt,READS_DIR=/scratch/<proj>/$USER/fungi_103    benchmark/pbs/clade_exclusion.pbs
```

`deplete_reference.py` handles both the BLCA-parsed references (BLCA/NB/Emu) and
the lineage-in-header EMITS target, and folds GTDB `_<UPPER>` suffixes so
`Escherichia coli` also matches `Escherichia coli_F`, `Cereibacter sphaeroides`
matches `Cereibacter_A sphaeroides`, etc. `EXCLUDE_LIST` can instead be your own
plain list (one species/line) or a TSV truth table (add `SPECIES_COL=<column>`).
Check the `N/M excluded species matched` line in the job log — unmatched names
mean a nomenclature mismatch to reconcile.

## Experiment 3 — leave-10%-out holdout, BLCA + NB (`holdout.pbs`)

The two splitters are the exact scripts that produced the published GTDB/UNITE
splits. Standalone:

```bash
python benchmark/split_gtdb_holdout.py  --in-fasta db/gtdb_ssu_BLCAparsed.fasta \
  --in-taxonomy db/gtdb_ssu_BLCAparsed.taxonomy --out-dir holdout_16s --seed 42
python benchmark/split_unite_holdout.py --in-fasta db_unite/unite_BLCAparsed.fasta \
  --in-taxonomy db_unite/unite_BLCAparsed.taxonomy --out-dir holdout_its --seed 42
```

Or split + rebuild BLCA/NB DBs from the 90% + classify the held-out 10% in one job:

```bash
qsub -v MARKER=16S,TEST_SUBSAMPLE=0 benchmark/pbs/holdout.pbs   # NB = two-step (GTDB)
qsub -v MARKER=ITS                  benchmark/pbs/holdout.pbs   # NB = single-step (UNITE)
```

NB handling matches the pipeline design: **two-step** (genus + addSpecies) for
16S/GTDB via `taxonomy_only`; **single-step** 7-rank `assignTaxonomy` for ITS/UNITE
(run directly, since `taxonomy_only`'s NB path is two-step only). Each split writes:

| file | contents |
|---|---|
| `reference_90.fasta` / `.taxonomy` | 90% training reference — the DBs are rebuilt from this |
| `test_10.fasta` / `.taxonomy` | 10% held-out queries + ground-truth lineages |
| `orphan_test_accessions.txt` | test species absent from train (species-level accuracy ceiling) |
| `split_stats.json` | counts + species-coverage diagnostics |

`TEST_SUBSAMPLE=N` caps how many held-out queries are classified; `0` (the default,
and the manuscript setting) classifies the whole held-out set (95,731 GTDB r232 and
10,213 UNITE sequences). Predictions land in
`benchmark_holdout_<marker>_seed42/results/taxonomy_{blca,nb}/`; score them with the
scripts in `benchmark/analysis/` (see below).

> **Cross-group reproducibility.** `random.Random(seed).shuffle` is stable across
> Python 3 versions and machines, so the split is fully determined by two things:
> the **input reference file** and the **`--seed`**. Two groups get an identical
> split only if they feed in the *same* reference (same GTDB release / same UNITE
> DOI, built the same way — record order matters) with the same seed. To guarantee
> it, either pin and share the exact `*_BLCAparsed.{fasta,taxonomy}` used, or
> distribute the resulting `reference_90`/`test_10` files (and `split_stats.json`)
> directly. Python is used rather than R precisely because R's `sample()` RNG
> changed between R versions (the `sample.kind` change in R 3.6).

## Scoring and figures (`analysis/`)

The R scripts that turn the run directories of experiments 1 and 3 into the
manuscript's supplementary tables, Figure 2 and the supplementary figures are in
[`analysis/`](analysis/README.md), together with the Gadi jobs that measured the
identity of every holdout test sequence to its closest training sequence. They read the
run directories from `HIFITAXA_DATA` and write to `HIFITAXA_OUT`.
