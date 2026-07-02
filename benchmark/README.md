# HiFiTaxa Gadi test harness

Scripts to (a) run the pipeline on real reads as non-interactive PBS jobs, and
(b) run the leave-10%-out reference holdout benchmark. Assumes the databases,
container images and conda envs are already staged on Gadi as in
[../docs/gadi.md](../docs/gadi.md) (steps 1–5b), and that the reads (ATCC 16S,
103-species fungal ITS) live under `/scratch/<proj>/$USER/...`.

Every script has a `CONFIG (edit me)` block at the top. Set `PROJ` to your NCI
project code and the read/repo paths; each var is also overridable at submit time
(`qsub -v READS_DIR=/path,PROJ=ab12 ...`).

## Data provenance

- **ATCC 16S mock** — ATCC MSA 16S, 192-plex PacBio HiFi. Downloaded from
  <http://downloads.pacbcloud.com/public/dataset/atcc_msa/16S_192plex_HiFi.fastq.tar.gz>
- **103-species fungal mock** — PacBio Revio HiFi, full-ITS (ITS9mun / ITS4ngsUni),
  ENA **PRJEB108994**. Reference:
  > Tedersoo L, *et al.* (2026). *Benchmarking full-length ITS metabarcoding across
  > Illumina 2×500, PacBio and Nanopore using mock and soil communities.* bioRxiv.
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

## 2. Run the pipeline on real reads (non-interactive)

```bash
qsub benchmark/pbs/run_atcc_16s.pbs      # ATCC 16S mock  -> results_atcc_16s/
qsub benchmark/pbs/run_fungi_its.pbs     # fungal ITS     -> results_fungi_its/
```

Each script builds `samples.tsv` + `metadata.tsv` from its `READS_DIR`, then runs
`bin/run_pipeline.py --classifier all` offline (`--skip-gtdb-check`,
`-profile singularity`, `--publish_dir_mode copy`).

## 3. Leave-10%-out holdout benchmark

One splitter per holdout (both Python, both standalone). They are the exact
scripts that produced the published GTDB/UNITE splits:

```bash
# 16S / GTDB (default --min-length 1000):
python benchmark/split_gtdb_holdout.py \
  --in-fasta   db/gtdb_ssu_BLCAparsed.fasta \
  --in-taxonomy db/gtdb_ssu_BLCAparsed.taxonomy \
  --out-dir holdout_16s --seed 42

# ITS / UNITE (default --min-length 0; ITS is short):
python benchmark/split_unite_holdout.py \
  --in-fasta   db_unite/unite_BLCAparsed.fasta \
  --in-taxonomy db_unite/unite_BLCAparsed.taxonomy \
  --out-dir holdout_its --seed 42

# or split + makeblastdb + classify held-out queries in one PBS job:
qsub -v MARKER=16S benchmark/pbs/holdout_blca.pbs    # uses split_gtdb_holdout.py
qsub -v MARKER=ITS benchmark/pbs/holdout_blca.pbs    # uses split_unite_holdout.py
```

Each writes into `--out-dir`:

| file | contents |
|---|---|
| `reference_90.fasta` / `.taxonomy` | 90% training reference — rebuild the classifier DB from this |
| `test_10.fasta` / `.taxonomy` | 10% held-out queries + ground-truth lineages |
| `orphan_test_accessions.txt` | test species absent from train (species-level accuracy ceiling) |
| `split_stats.json` | counts + species-coverage diagnostics |

The split is a pure seeded random split (not stratified by species), so the
orphan rate is measured, not engineered away. `holdout_blca.pbs` then
`makeblastdb`s `reference_90` and classifies `test_10` with BLCA via the
pipeline's `taxonomy_only` entry. `TEST_SUBSAMPLE=N` caps the number of queries
(BLCA is per-query expensive; a full GTDB 10% is ~95k seqs). Score the
predictions in `blca_results/taxonomy_blca/` against `test_10.taxonomy`.

> **Cross-group reproducibility.** `random.Random(seed).shuffle` is stable across
> Python 3 versions and machines, so the split is fully determined by two things:
> the **input reference file** and the **`--seed`**. Two groups get an identical
> split only if they feed in the *same* reference (same GTDB release / same UNITE
> DOI, built the same way — record order matters) with the same seed. To guarantee
> it, either pin and share the exact `*_BLCAparsed.{fasta,taxonomy}` used, or
> distribute the resulting `reference_90`/`test_10` files (and `split_stats.json`)
> directly. Python is used rather than R precisely because R's `sample()` RNG
> changed between R versions (the `sample.kind` change in R 3.6).
