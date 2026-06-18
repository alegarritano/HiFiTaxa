# HiFiTaxa

HiFiTaxa is a Nextflow pipeline for PacBio HiFi full length 16S rRNA data. It
denoises reads with QIIME2 and DADA2, then assigns taxonomy against GTDB with up
to three classifiers you can run in one invocation and compare directly: BLCA,
Emu, and DADA2 Naive Bayes.

Pipeline: read QC, primer trimming (cutadapt), DADA2 denoising, ASV filtering,
then taxonomy with any combination of BLCA, Emu, and the DADA2 two step Naive
Bayes (genus `assignTaxonomy` plus exact match `addSpecies`). All three
classifiers are anchored on the latest GTDB SSU release (r232).

## Requirements

- Nextflow 24.04 to 24.x (the launcher pins 24.10.9; Nextflow 25/26 will not parse the config)
- Singularity/Apptainer (on an HPC) or Docker (on a laptop), or conda/mamba
- About 3 GB of free disk for the GTDB database, plus image cache, and internet for the first run

## Install

```
git clone https://github.com/alegarritano/HiFiTaxa.git
cd HiFiTaxa
mamba env create -f environment.yml
conda activate HiFiTaxa
source set_apptainer_cache.sh     # optional on HPC: the launcher sets repo-local image caches itself; source only to override
```

`environment.yml` installs a small driver environment: Nextflow, the launcher's
Python dependencies, and makeblastdb for building the database. The actual
analysis tools come from containers (`-profile singularity` or `-profile docker`)
or from conda environments you point at (`-profile standard`).

### Zero-config setup (no bioinformatics experience required)

The launcher auto-installs every analysis env it needs on first use. You do not
have to know any tool names or conda commands. The first `python bin/run_pipeline.py`
call splits into three phases:

**Phase 1 (about 1 min, interactive).** You answer 5 quick questions:

  1. Confirm the GTDB SSU reference release.
  2. Drop reference sequences shorter than 1000 bp? (recommended for full length 16S)
  3. Add the **Emu** classifier? (EM based species profiling on raw reads)
  4. Add the **NB** classifier (`nb`)? DADA2 ASVs go to DADA2 Naive Bayes (genus
     `assignTaxonomy`) plus exact match `addSpecies`.
  5. Run the bundled 8 sample ATCC mock to validate the install when ready?

**Phase 2, install (about 20 to 60 min, unattended).** Based on your answers the launcher:

  1. Bootstraps `mamba` into the conda base env if missing (about 3 to 5 min).
  2. Downloads and parses the GTDB SSU database (about 10 min, about 3 GB).
  3. If you said yes to Emu: auto creates `hifitax_emu` (emu, minimap2, samtools
     from bioconda, about 3 to 5 min) and builds the Emu DB (about 5 min).
  4. If you said yes to NB: builds the two DADA2 references (genus and species)
     from the GTDB SSU database in about 1 to 2 minutes. No training step.
     On `-profile singularity/docker` the classifier runs inside the QIIME2
     amplicon image you already cache for denoising, so no extra env is needed.
     On `-profile standard` the launcher auto creates a `qiime2-amplicon` conda
     env once (downloads QIIME2's own install YAML, about 10 to 15 min).
  5. Pre pulls the Singularity/Docker images for QC, denoising, and BLCA.

**Phase 3, validate and run (only if you said yes to Q5).** Runs the bundled
8 sample ATCC mock end to end with every classifier you enabled, compares the
result to the reference shipped in `example/`, prints PASS/WARN/FAIL, then
proceeds to your own data.

If an env is already on your system (any name, any prefix), the launcher
auto detects it via `conda env list` and skips the install. Override the
detected env any time with `--emu_env <path>`, `--qiime2_env <path>`, or
`--blca_env <path>`.

## For subsequent runs

After the first setup, every new shell or HPC job needs the environment prepared
again before you launch the pipeline:

```
cd HiFiTaxa
conda activate HiFiTaxa
source set_apptainer_cache.sh     # optional on HPC: the launcher sets repo-local image caches itself; source only to override
```

## Quick test

Run the bundled 8 sample ATCC mock and compare it against the reference results
shipped in `example/`:

```
python bin/run_pipeline.py --test --profile singularity
```

(On a laptop with Docker instead of Singularity, use `--profile docker`.) Add
`--classifier all` to exercise all three classifiers in the smoke test.

The first run builds the GTDB r232 database (about 10 minutes, once), runs the
example, and prints PASS/WARN/FAIL comparing your output to `example/reference/`.

## Using conda instead of containers

If you cannot use containers, run with `-profile standard`. The launcher will
auto detect (or auto install) the QIIME2 amplicon env, the Emu env, and the
BLCA env, so you do not need to supply any flags:

```
python bin/run_pipeline.py --test --profile standard
```

If you already have envs you would rather reuse (for example on a shared HPC),
point the launcher at them explicitly:

```
python bin/run_pipeline.py --test --profile standard \
  --qiime2_env /path/to/qiime2-amplicon-env \
  --emu_env    /path/to/hifitax_emu \
  --blca_env   /path/to/blca-env
```

## Run your own data

Run it interactively (it builds the sample sheet and asks you things) or
non interactively for batch/HPC jobs (everything via flags). Reads may be
`.fastq.gz` or plain `.fastq`.

### Interactive

1. Put **all** your reads (`.fastq.gz` or `.fastq`) in a folder. The default
   location the launcher looks at is `00_Reads/` inside the project, but you can
   keep them anywhere and pass `--reads_dir /path/to/reads`.
2. Start the pipeline (no `--input` needed):

   ```
   # reads in 00_Reads/ (default)
   python bin/run_pipeline.py --profile singularity

   # reads anywhere else
   python bin/run_pipeline.py --profile singularity --reads_dir /scratch/me/my_run/fastqs
   ```

   On the first run it offers the quick test above; once that finishes it sets
   up your run.
3. The launcher scans `00_Reads` and writes two files for you:
   - `samples.tsv`, one row per read file (sample id plus absolute path)
   - `metadata.tsv`, one row per sample, `sample_name` filled in, `condition` left blank.
4. It pauses and asks you to fill the `condition` column in `metadata.tsv`, save,
   and confirm. It will not proceed until every sample has a condition.
5. It then asks you to confirm the primers (full length 16S 27F / 1492R by
   default, or type your own forward then reverse sequence).
6. The analysis runs: QC, primer trim, DADA2, ASV filter, then your chosen classifier(s).

### Batch / HPC scheduler jobs (non interactive)

A scheduler job has no terminal to answer prompts, so make `samples.tsv` and
`metadata.tsv` yourself and pass everything as flags. `samples.tsv` is
`sample-id <tab> absolute-filepath`; `metadata.tsv` is `sample_name <tab> condition`.

```
python bin/run_pipeline.py \
  --input samples.tsv --metadata metadata.tsv \
  --classifier all \
  --forward_primer AGRGTTYGATYMTGGCTCAG --reverse_primer AAGTCGTAACAAGGTARCY \
  --skip-gtdb-check \
  --outdir results --profile singularity
```

Flags for batch runs:

| flag | effect |
|---|---|
| `--input` / `--metadata` | sample sheet plus metadata (required; avoids all prompts) |
| `--classifier` | which classifier(s) to run (see below); default `blca` |
| `--forward_primer` / `--reverse_primer` | primers (omit to use the 27F / 1492R defaults) |
| `--skip_primer_trim` | skip cutadapt (for example reads already trimmed) |
| `--skip-gtdb-check` | use the existing DB as is (no network check, no build) |
| `--assume-yes` | build/rebuild the DB if missing or outdated |
| `--assume-no` | never build; proceed with the existing DB |
| `--skip-test` | never run the first run example test |
| `--publish_dir_mode copy` | write real files into `results/` instead of symlinks into `work/` (recommended for archiving) |
| `-resume` | resume a previous run (passed through to Nextflow) |

With no terminal the launcher never prompts; if you pass no database flag it
builds the DB only when it is missing. Example `sbatch` wrapper:

```
#!/bin/bash
#SBATCH -c 16 --mem 64G -t 12:00:00
source set_apptainer_cache.sh        # keep image caches off $HOME (HPC)
python bin/run_pipeline.py \
  --input samples.tsv --metadata metadata.tsv \
  --classifier all \
  --skip-gtdb-check --outdir results --profile singularity
```

### Taxonomy only (from an existing ASV FASTA)

Skip QC and denoising and classify an existing ASV/sequence FASTA directly. This
entry runs **BLCA, Emu, and/or NB**, chosen with `--classifier`:

```
# BLCA (default)
python bin/run_pipeline.py --asv_fasta my_ASVs.fasta --outdir tax_out --profile singularity

# all three
python bin/run_pipeline.py --asv_fasta my_ASVs.fasta --classifier all --outdir tax_out --profile singularity

# just NB (genus assignTaxonomy + exact-match addSpecies)
python bin/run_pipeline.py --asv_fasta my_ASVs.fasta --classifier nb --outdir tax_out --profile singularity
```

Each classifier writes a per-sequence taxonomy. NB here runs `assignTaxonomy`
plus `addSpecies` and writes `best_taxonomy.tsv` and `gtdb_nb.tsv`; it skips the
ASV-by-sample frequency merge, because a standalone FASTA carries no frequency
table. Emu gives a per-sequence assignment whose abundance numbers are not
quantitative on dereplicated input.

## Classifier choice (BLCA, Emu, DADA2 NB)

The pipeline can classify ASVs or reads with any combination of three
classifiers, all anchored on the same GTDB release:

- **BLCA** (default). Bayesian LCA over BLAST hits with **per rank bootstrap
  confidence** for every ASV. Runs on DADA2 denoised ASVs (Gao *et al.* 2017,
  *mSystems*).
- **Emu**. EM based species level abundance estimator. Runs on the trimmed raw
  reads (no DADA2), with strong species level resolution and a low false
  positive rate (Curry *et al.* 2022, *Nat. Methods*).
- **NB**. DADA2 Naive Bayes (the RDP algorithm of Wang *et al.* 2007, *AEM*,
  [10.1128/aem.00062-07](https://journals.asm.org/doi/10.1128/aem.00062-07)),
  run as the canonical DADA2 two step: `assignTaxonomy()` against a genus level
  GTDB reference (per rank bootstrap, `minBoot=80`) followed by exact match
  `addSpecies()` for species. A genus level reference is required because
  `assignTaxonomy()` with full species lineages does not scale to the roughly
  82k species labels in full GTDB SSU r232 (the bootstrap dilutes and collapses
  to Kingdom). Species recovered by `addSpecies` are high precision but lower
  recall on full length reads, since the ASV must exactly match a reference 16S.

Pick with `--classifier`: a single name, a comma list, or a shorthand. `all`
expands to `blca,emu,nb`; `both` expands to `blca,emu` (kept for back compat):

```
python bin/run_pipeline.py --input samples.tsv --metadata metadata.tsv \
  --classifier blca            # default: DADA2 then BLCA / GTDB
# or
  --classifier emu             # trimmed reads then Emu / GTDB
# or
  --classifier nb              # DADA2 ASVs then DADA2 NB genus + addSpecies / GTDB  (aliases: dada2_nb, qiime2_nb)
# or
  --classifier blca,emu        # both BLCA and Emu (legacy 'both' still works)
# or
  --classifier all             # run all three branches
```

The three classifiers complement each other: BLCA gives per rank bootstrap
confidence on denoised ASVs; Emu gives species level abundance estimates
directly from reads with a low false positive rate; NB is the widely used
DADA2/RDP Naive Bayes approach, useful as a reference comparison. Running
several branches keeps one reference (GTDB) and compares them on the same
samples in a single invocation.

The `-entry taxonomy_only` mode (triggered by `--asv_fasta`) supports BLCA, Emu,
and NB on a standalone FASTA.

## GTDB database

### BLCA database

Lives in `--gtdb_db_dir` (default `db/`): `gtdb_ssu_BLCAparsed.fasta` with its
BLAST index, `gtdb_ssu_BLCAparsed.taxonomy`, and `GTDB_VERSION.txt`. Default
release is r232. The Emu and NB references are both derived from these same
files, so all three classifiers share one reference build.

On the first run the preflight builds it. With `--gtdb-release latest` it checks
the GTDB site for a newer release and offers to rebuild. If the site is
unreachable it uses whatever is already there. During a build it asks whether to
drop reference sequences shorter than 1000 bp; set this non interactively with
`--min-ref-len` (0 keeps everything).

Build by hand:
```
bash bin/build_gtdb_blca_db.sh 232 db 1000     # 1000 drops shorter refs, 0 keeps all
```

### Emu database (only when `--classifier` includes `emu`)

Lives in `--emu_db_dir` (default `db_emu/`). It is built from the same
BLCA parsed GTDB files so both classifiers see the same reference release, by
remapping each GTDB lineage to a synthetic tax id and calling `emu
build-database`. On first use with an Emu enabled run the launcher builds it;
`--assume-yes` builds it non interactively, `--skip-emu-db-check` uses whatever
is already there.

On `-profile singularity/docker` the build runs **inside the Emu container** (it
ships `emu` + `python3`), so no conda env is created. This avoids the brittle
emu/samtools/pysam conda solve that fails on many HPCs. On `-profile standard`
the launcher creates a small conda env from conda-forge + bioconda instead.

Build by hand (conda env with `emu` active, or inside the Emu image):
```
bash bin/build_gtdb_emu_db.sh \
  db/gtdb_ssu_BLCAparsed.fasta \
  db/gtdb_ssu_BLCAparsed.taxonomy \
  db_emu

# or, with no conda, straight in the Emu container:
singularity exec --bind "$PWD" \
  docker://quay.io/biocontainers/emu@sha256:61ea3336f12d41930d73e57ce1b041bce48d66b4011a165bf1f0efce9d684777 \
  bash bin/build_gtdb_emu_db.sh \
    db/gtdb_ssu_BLCAparsed.fasta db/gtdb_ssu_BLCAparsed.taxonomy db_emu
```

### NB references (only when `--classifier` includes `nb`)

The DADA2 NB path needs **no training**. It uses two gzipped FASTA references,
built from the same BLCA parsed GTDB that BLCA and Emu use, under `--nb_db_dir`
(default `db_nb/`):

- `gtdb_ssu_dada2_genus.fa.gz` (6 rank Kingdom to Genus, for `assignTaxonomy`)
- `gtdb_ssu_dada2_species.fa.gz` (`>acc Genus species`, for `addSpecies`)

Building is a fast reformat (about 1 to 2 minutes, low memory), not a training
step. On the first NB enabled run the launcher builds both automatically; pass
`--skip-nb-classifier-check` to assume they are already in place.

Build by hand:
```
bash bin/build_gtdb_dada2_db.sh \
  db/gtdb_ssu_BLCAparsed.fasta \
  db/gtdb_ssu_BLCAparsed.taxonomy \
  db_nb/gtdb_ssu_dada2_genus.fa.gz \
  db_nb/gtdb_ssu_dada2_species.fa.gz
```

This needs only `python3` (no `qiime`, no container). The genus level reference
is used for `assignTaxonomy` because the full species lineage does not scale to
GTDB r232 (about 82k species labels dilute the bootstrap to Kingdom only);
species is then recovered by exact match `addSpecies`.

## Outputs (`<outdir>/`)

Always:
- `stats.tsv`, read tracking through the pipeline.

When `--classifier` includes `blca`:
- `dada2/dada2_ASV.fasta`, filtered ASV sequences
- `dada2/dada2-ccs_table_filtered.qza`, ASV table
- `taxonomy_blca/blca_taxonomy_table.csv`, 7 rank taxonomy, one row per ASV
- `taxonomy_blca/blca_taxonomy_confidence.csv`, same, with per rank bootstrap confidence
- `taxonomy_blca/ASV_blca.out`, raw BLCA output

When `--classifier` includes `emu`:
- `taxonomy_emu/emu_species_table.tsv`, species level relative abundance matrix (samples by species)
- `taxonomy_emu/emu_genus_table.tsv`, genus level relative abundance matrix
- `taxonomy_emu/per_sample/<sample>_rel-abundance.tsv`, raw Emu per sample output

When `--classifier` includes `nb`:
- `taxonomy_nb/best_taxonomy.tsv`, per ASV lineage (`d__..s__`) and confidence (QIIME2 TSVTaxonomyFormat)
- `taxonomy_nb/best_taxonomy_withDB.tsv`, same, plus the source database column
- `taxonomy_nb/gtdb_nb.tsv`, per rank detail: Kingdom to Genus with bootstrap support, plus exact match Species
- `taxonomy_nb/best_tax_merged_freq_tax.tsv`, ASV by sample frequency table merged with the taxonomy
- `taxonomy_nb/best_tax.qza`, taxonomy as a QIIME2 artifact

By default these are symlinks into `work/`. Pass `--publish_dir_mode copy` to
write standalone files (safe to keep after `work/` is cleaned).

## Parameters

| param | default | meaning |
|---|---|---|
| `--input` / `--metadata` | (none) | samples TSV / metadata TSV |
| `--asv_fasta` | (none) | taxonomy only mode (`-entry taxonomy_only`; BLCA, Emu, and/or NB via `--classifier`) |
| `--classifier` | `blca` | classifier branch(es): `blca`, `emu`, `nb` (aliases: `dada2_nb`, `qiime2_nb`), a comma list, or `all` / `both`. The `nb` branch denoises with DADA2 and classifies ASVs with DADA2's `assignTaxonomy` against a genus level GTDB reference, then recovers species by exact match `addSpecies`. On `-profile singularity/docker` the tools are already in the QIIME2 image, so no extra env install is needed. |
| `--forward_primer` / `--reverse_primer` | 27F / 1492R | primer sequences |
| `--min_len` / `--max_len` | 1000 / 1600 | DADA2 length filter (BLCA and NB branches) |
| `--pooling_method` | pseudo | DADA2 pooling (BLCA and NB branches) |
| `--min_asv_totalfreq` / `--min_asv_sample` | 5 / 1 | ASV level noise filter applied after DADA2 (BLCA and NB branches). `min_asv_totalfreq` drops ASVs whose **summed read count across all samples** is below the threshold (default 5, likely sequencing noise). `min_asv_sample` drops ASVs present in **fewer than this many samples** (default 1, keep anything in at least one sample). Both auto drop to 0 for single sample runs so nothing gets filtered. |
| `--blca_chunk_size` | auto | ASVs per BLCA chunk; auto splits into about (available cores minus 2) chunks |
| `--max_cpus` | available cores minus 2 | usable cores for auto chunking the BLCA step (leaves headroom) |
| `--blca_minid` | 90 | BLCA minimum percent identity |
| `--emu_db_dir` | `db_emu/` | prebuilt Emu DB directory (auto built from the GTDB BLCA DB on first Emu enabled run) |
| `--emu_type` | `map-hifi` | minimap2 preset Emu uses. `map-hifi` is the PacBio HiFi preset (minimap2 2.19 or newer) and the default here. Other choices: `map-ont` (Oxford Nanopore), `map-pb` (legacy PacBio CLR), `sr` (short reads), `lr:hq` (high quality long reads). |
| `--emu_threads` | 4 | CPU threads per `emu abundance` call |
| `--emu_min_abundance` | Emu default | minimum relative abundance threshold (for example `1e-5`) for Emu's noise trim |
| `--skip-emu-db-check` | off | use the existing Emu DB as is; never build |
| `--nb_db_dir` | `db_nb/` | directory for the two DADA2 NB references (`gtdb_ssu_dada2_genus.fa.gz` and `gtdb_ssu_dada2_species.fa.gz`); both auto built from the GTDB BLCA DB on first NB enabled run (fast reformat, about 1 to 2 min, low memory, no training) |
| `--nb_min_bootstrap` | 80 | DADA2 `minBoot` for the genus level `assignTaxonomy` step |
| `--skip-nb-classifier-check` | off | use the existing NB references as is; never build |
| `--gtdb-release` | `latest` | GTDB release to build (`latest` queries the GTDB site for the newest available release). Pass an explicit number (for example `232`) to pin a version for reproducibility. |
| `--min-ref-len` | 1000 | drop reference seqs shorter than this (0 keeps all) |
| `--qiime2_env` / `--blca_env` / `--emu_env` | auto | conda envs for `-profile standard`. All auto detected and auto installed if missing; supply explicit paths only to reuse pre existing envs. `--qiime2_env` is also used for the QIIME2 export/tabulate steps in the NB branch. |
| `--blca_container` / `--emu_container` | image refs | per branch container images for `-profile docker` / `singularity` |
| `--publish_dir_mode` | symlink | how outputs land in `<outdir>/`: `symlink` (default) or `copy` (standalone files) |
| `--progress-bar` | off | show a tqdm bar instead of Nextflow's native display |

## Notes

- Conda works on Linux, but QIIME2 and BLAST have no Apple Silicon builds. On
  macOS, reuse an existing QIIME2 amplicon env (add seqkit and csvtk to it) and a
  clustalo/blast/biopython env, and pass them with `--qiime2_env` / `--blca_env`.
- BLCA scores each ASV independently, so chunk size only affects speed, not the
  taxonomy.
- BLCA runs with `-p 1` (the conda clustalo build has no OpenMP); parallelism
  comes from running chunks at the same time.
- NB classifies each ASV independently against the references, so the result
  does not depend on how many ASVs are in the run.

## References

- Bolyen E *et al.* *Reproducible, interactive, scalable and extensible
  microbiome data science using QIIME 2.* Nature Biotechnology 37, 852 to 857
  (2019). The QIIME 2 framework. For plugin citations, run
  `qiime <plugin-name> --citations` after activating the QIIME 2 environment.
- Callahan BJ *et al.* *DADA2: High resolution sample inference from Illumina
  amplicon data.* Nature Methods 13, 581 to 583 (2016). The DADA2 denoiser, and
  the `assignTaxonomy` plus exact match `addSpecies` used in the NB branch.
- Wang Q *et al.* *Naive Bayesian classifier for rapid assignment of rRNA
  sequences into the new bacterial taxonomy.* Applied and Environmental
  Microbiology 2007, [10.1128/aem.00062-07](https://journals.asm.org/doi/10.1128/aem.00062-07).
  The RDP Naive Bayes 16S algorithm that DADA2's `assignTaxonomy` implements.
- Gao X, Lin H, Dong Q. *A Bayesian taxonomic classification method for 16S
  rRNA gene sequences with improved species level accuracy.* mSystems 2017.
  The BLCA classifier.
- Curry KD *et al.* *Emu: species level microbial community profiling of full
  length 16S rRNA Oxford Nanopore sequencing data.* Nature Methods 2022. The
  Emu classifier.
- Parks DH *et al.* *GTDB: a complete and systematic taxonomy.* Nucleic Acids
  Research 2025. The reference database.

## License and citation

The pipeline code is MIT (see `LICENSE`). Bundled tools keep their own licenses:
BLCA (`bin/blca_main.py`, GPL) and the denoise modules adapted from PacBio's
HiFi-16S-workflow (BSD-3-Clause-Clear). If you use it, please cite this repo
(`CITATION.cff`) and the tools it wraps: QIIME2, DADA2, cutadapt, BLCA, Emu,
GTDB, and Nextflow.
