# HiFiTaxa

HiFiTaxa is a Nextflow pipeline for PacBio HiFi full-length amplicon data. It
handles two markers, selected with `--marker`:

- **16S** (default): QC, primer trimming (cutadapt), DADA2 denoising, ASV
  filtering, then taxonomy against GTDB (r232) with any combination of BLCA, Emu,
  and the DADA2 two-step Naive Bayes (genus `assignTaxonomy` plus exact-match
  `addSpecies`) — run in one invocation and compared on the same samples.
- **ITS** (fungal): QC, then itsxrust pulls out the ITS region, and taxonomy is
  assigned against UNITE. EMITS (a read-level EM profiler, the fungal analogue of
  Emu) is the default; BLCA and single-step Naive Bayes run on the DADA2 ASVs if
  you ask for them (DADA2 runs only then).

See **[docs/MARKER_AWARE.md](docs/MARKER_AWARE.md)** for the full marker design.

## Requirements

- Nextflow 24.04 to 24.x (the launcher pins 24.10.9; Nextflow 25/26 will not parse the config)
- Singularity/Apptainer (on an HPC) or Docker (on a laptop), or conda/mamba
- About 3 GB of free disk for the GTDB database (16S) or ~0.5 GB for UNITE (ITS), plus image cache, and internet for the first run

## Install

```
git clone https://github.com/alegarritano/HiFiTaxa.git
cd HiFiTaxa
mamba env create -f environment.yml
conda activate HiFiTaxa
source set_apptainer_cache.sh     # HPC: keep the image + conda caches inside the repo (off your $HOME quota)
```

`environment.yml` installs a small driver environment: Nextflow, the launcher's
Python dependencies, and makeblastdb for building the database. The actual
analysis tools come from containers (`-profile singularity` or `-profile docker`)
or from conda environments you point at (`-profile standard`).

The launcher auto-installs every analysis env and reference it needs on first
use — you do not have to know any tool names or conda commands. For the full
first-run walkthrough (the interactive questions, what gets built, reusing
existing envs) see **[docs/setup.md](docs/setup.md)**. To run without containers,
see **[docs/profiles.md](docs/profiles.md)**.

## Run your own data

Run interactively (it builds the sample sheet and prompts you) or non
interactively for batch/HPC jobs (everything via flags). Reads may be
`.fastq.gz` or plain `.fastq`.

### Interactive

1. Put **all** your reads in a folder (default `00_Reads/`, or pass `--reads_dir /path`).
2. Start the pipeline:

   ```
   python bin/run_pipeline.py --profile singularity
   # fungal ITS instead of 16S:
   python bin/run_pipeline.py --profile singularity --marker ITS
   # reads elsewhere:
   python bin/run_pipeline.py --profile singularity --reads_dir /scratch/me/fastqs
   ```
3. The launcher scans the reads and writes `samples.tsv` (one row per read file)
   and `metadata.tsv` (one row per sample). It pauses for you to fill the
   `condition` column, then asks you to confirm the primers (full length 16S
   27F / 1492R by default). It then runs QC → primer trim → DADA2 → ASV filter →
   your chosen classifier(s).

### Batch / HPC (non interactive)

A scheduler job has no terminal to answer prompts, so make `samples.tsv`
(`sample-id <tab> absolute-filepath`) and `metadata.tsv`
(`sample_name <tab> condition`) yourself and pass everything as flags:

```
python bin/run_pipeline.py \
  --input samples.tsv --metadata metadata.tsv \
  --classifier all \
  --forward_primer AGRGTTYGATYMTGGCTCAG --reverse_primer AAGTCGTAACAAGGTARCY \
  --skip-gtdb-check \
  --outdir results --profile singularity --publish_dir_mode copy
```

The full flag and parameter reference is in **[docs/parameters.md](docs/parameters.md)**.

### Taxonomy only (from an existing ASV FASTA)

Skip QC and denoising and classify an existing ASV/sequence FASTA directly with
`--asv_fasta` (runs BLCA, Emu, and/or NB, chosen with `--classifier`):

```
python bin/run_pipeline.py --asv_fasta my_ASVs.fasta --classifier all --outdir tax_out --profile singularity
```

Each classifier writes a per-sequence taxonomy; see
**[docs/parameters.md](docs/parameters.md)** for the per-branch outputs and caveats.

## Classifier choice (BLCA, Emu, DADA2 NB)

The pipeline can classify ASVs or reads with any combination of three
classifiers, all anchored on the same GTDB release:

- **BLCA** (default). Bayesian LCA over BLAST hits with **per rank bootstrap
  confidence** for every ASV. Runs on DADA2 denoised ASVs (Gao *et al.* 2017,
  *mSystems*).
- **Emu**. EM based species level abundance estimator. Runs on the trimmed raw
  reads (no DADA2), with strong species level resolution and a low false
  positive rate (Curry *et al.* 2022, *Nat. Methods*).
- **NB**. DADA2 Naive Bayes (the RDP algorithm of Wang *et al.* 2007, *AEM*),
  run as the canonical DADA2 two step: `assignTaxonomy()` against a genus level
  GTDB reference (per rank bootstrap, `minBoot=80`) followed by exact match
  `addSpecies()` for species. A genus level reference is required because
  `assignTaxonomy()` with full species lineages does not scale to the roughly
  82k species labels in full GTDB SSU r232 (the bootstrap dilutes and collapses
  to Kingdom). Species from `addSpecies` are high precision but lower recall on
  full length reads, since the ASV must exactly match a reference 16S.

For **fungal ITS** (`--marker ITS`) the read-level classifier is **EMITS** (the
fungal analogue of Emu), and it is the default. BLCA and single-step NB run on the
itsxrust-extracted DADA2 ASVs if you add them; Emu and the 16S two-step NB do not
apply to ITS.

Pick with `--classifier`: a single name, a comma list, or a shorthand. For 16S,
`all` expands to `blca,emu,nb`; `both` expands to `blca,emu` (kept for back compat):

```
python bin/run_pipeline.py --input samples.tsv --metadata metadata.tsv \
  --classifier blca            # default: DADA2 then BLCA / GTDB
# or  --classifier emu         # trimmed reads then Emu / GTDB
# or  --classifier nb          # DADA2 ASVs then DADA2 NB genus + addSpecies (aliases: dada2_nb, qiime2_nb)
# or  --classifier blca,emu    # both (legacy 'both' still works)
# or  --classifier all         # run all three branches
```

The three complement each other: BLCA gives per rank bootstrap confidence on
denoised ASVs; Emu gives species level abundance estimates directly from reads
with a low false positive rate; NB is the widely used DADA2/RDP Naive Bayes
approach, useful as a reference comparison. Running several branches keeps one
reference (GTDB) and compares them on the same samples in a single invocation.

## Databases

For 16S, all classifiers share one GTDB reference build (default r232), downloaded
and built automatically on the first run. For **ITS**, download a UNITE release
from <https://unite.ut.ee/repository.php> and the pipeline builds the BLCA / NB /
EMITS references from it. To pin a release, rebuild, or build any reference by
hand, see **[docs/databases.md](docs/databases.md)**.

## Offline / air-gapped HPC

On compute nodes with no internet, pre-stage the container images on a login
node once, then run offline. See **[docs/offline.md](docs/offline.md)** (and
**[docs/gadi.md](docs/gadi.md)** for the Gadi-specific version).

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
  microbiome data science using QIIME 2.* Nature Biotechnology 37, 852–857 (2019).
- Callahan BJ *et al.* *DADA2: High resolution sample inference from Illumina
  amplicon data.* Nature Methods 13, 581–583 (2016). The denoiser, and the
  `assignTaxonomy` + `addSpecies` used in the NB branch.
- Wang Q *et al.* *Naive Bayesian classifier for rapid assignment of rRNA
  sequences into the new bacterial taxonomy.* Applied and Environmental
  Microbiology 2007, [10.1128/aem.00062-07](https://journals.asm.org/doi/10.1128/aem.00062-07).
- Gao X, Lin H, Dong Q. *A Bayesian taxonomic classification method for 16S rRNA
  gene sequences with improved species level accuracy.* mSystems 2017. The BLCA classifier.
- Curry KD *et al.* *Emu: species level microbial community profiling of full
  length 16S rRNA Oxford Nanopore sequencing data.* Nature Methods 2022. The Emu classifier.
- Parks DH *et al.* *GTDB: a complete and systematic taxonomy.* Nucleic Acids
  Research 2025. The reference database.

## License and citation

The pipeline code is MIT (see `LICENSE`). Bundled tools keep their own licenses:
BLCA (`bin/blca_main.py`, GPL) and the denoise modules adapted from PacBio's
HiFi-16S-workflow (BSD-3-Clause-Clear). If you use it, please cite this repo
(`CITATION.cff`) and the tools it wraps: QIIME2, DADA2, cutadapt, BLCA, Emu,
GTDB, and Nextflow.
