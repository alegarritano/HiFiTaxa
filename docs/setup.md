# First-run setup (zero-config)

[← back to README](../README.md)

The launcher auto-installs every analysis env it needs on first use. You do not
have to know any tool names or conda commands. The first `python bin/run_pipeline.py`
call splits into three phases.

## Phase 1 (about 1 min, interactive)

You answer 5 quick questions:

1. Confirm the GTDB SSU reference release.
2. Drop reference sequences shorter than 1000 bp? (recommended for full length 16S)
3. Add the **Emu** classifier? (EM based species profiling on raw reads)
4. Add the **NB** classifier (`nb`)? DADA2 ASVs go to DADA2 Naive Bayes (genus
   `assignTaxonomy`) plus exact match `addSpecies`.
5. Run the bundled 8 sample ATCC mock to validate the install when ready?

## Phase 2 — install (about 20 to 60 min, unattended)

Based on your answers the launcher:

1. Bootstraps `mamba` into the conda base env if missing (about 3 to 5 min).
2. Downloads and parses the GTDB SSU database (about 10 min, about 3 GB).
3. If you said yes to Emu: auto creates `hifitax_emu` (emu, minimap2, samtools
   from bioconda) and builds the Emu DB (about 5 min).
4. If you said yes to NB: builds the two DADA2 references (genus and species)
   from the GTDB SSU database in about 1 to 2 minutes. No training step. On
   `-profile singularity/docker` the classifier runs inside the QIIME2 amplicon
   image you already cache for denoising, so no extra env is needed. On
   `-profile standard` the launcher auto creates a `qiime2-amplicon` conda env
   once (downloads QIIME2's own install YAML, about 10 to 15 min).
5. Pre pulls the Singularity/Docker images for QC, denoising, and BLCA.

## Phase 3 — validate and run (only if you said yes to Q5)

Runs the bundled 8 sample ATCC mock end to end with every classifier you
enabled, compares the result to the reference shipped in `example/`, prints
PASS/WARN/FAIL, then proceeds to your own data.

## Reusing existing environments

If an env is already on your system (any name, any prefix), the launcher
auto-detects it via `conda env list` and skips the install. Override the
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
