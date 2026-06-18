# Parameters and flags

[← back to README](../README.md)

## Batch-run flags

For non-interactive scheduler jobs, pass everything as flags (no terminal to
answer prompts). `samples.tsv` is `sample-id <tab> absolute-filepath`;
`metadata.tsv` is `sample_name <tab> condition`.

| flag | effect |
|---|---|
| `--input` / `--metadata` | sample sheet plus metadata (required; avoids all prompts) |
| `--classifier` | which classifier(s) to run (see table below); default `blca` |
| `--forward_primer` / `--reverse_primer` | primers (omit to use the 27F / 1492R defaults) |
| `--skip_primer_trim` | skip cutadapt (for example reads already trimmed) |
| `--skip-gtdb-check` | use the existing DB as is (no network check, no build) |
| `--assume-yes` | build/rebuild the DB if missing or outdated |
| `--assume-no` | never build; proceed with the existing DB |
| `--skip-test` | never run the first-run example test |
| `--publish_dir_mode copy` | write real files into `results/` instead of symlinks into `work/` (recommended for archiving) |
| `-resume` | resume a previous run (passed through to Nextflow) |

With no terminal the launcher never prompts; if you pass no database flag it
builds the DB only when it is missing.

## Full parameter reference

| param | default | meaning |
|---|---|---|
| `--input` / `--metadata` | (none) | samples TSV / metadata TSV |
| `--asv_fasta` | (none) | taxonomy-only mode (`-entry taxonomy_only`; BLCA, Emu, and/or NB via `--classifier`) |
| `--classifier` | `blca` | classifier branch(es): `blca`, `emu`, `nb` (aliases: `dada2_nb`, `qiime2_nb`), a comma list, or `all` / `both`. The `nb` branch denoises with DADA2 and classifies ASVs with DADA2's `assignTaxonomy` against a genus-level GTDB reference, then recovers species by exact match `addSpecies`. On `-profile singularity/docker` the tools are already in the QIIME2 image, so no extra env install is needed. |
| `--forward_primer` / `--reverse_primer` | 27F / 1492R | primer sequences |
| `--min_len` / `--max_len` | 1000 / 1600 | DADA2 length filter (BLCA and NB branches) |
| `--pooling_method` | pseudo | DADA2 pooling (BLCA and NB branches) |
| `--min_asv_totalfreq` / `--min_asv_sample` | 5 / 1 | ASV-level noise filter after DADA2 (BLCA and NB branches). `min_asv_totalfreq` drops ASVs whose summed read count across all samples is below the threshold; `min_asv_sample` drops ASVs present in fewer than this many samples. Both auto-drop to 0 for single-sample runs. |
| `--blca_chunk_size` | auto | ASVs per BLCA chunk; auto splits into about (available cores − 2) chunks |
| `--max_cpus` | cores − 2 | usable cores for auto-chunking the BLCA step (leaves headroom) |
| `--blca_minid` | 90 | BLCA minimum percent identity |
| `--emu_db_dir` | `db_emu/` | prebuilt Emu DB directory (auto-built from the GTDB BLCA DB on first Emu-enabled run) |
| `--emu_type` | `map-hifi` | minimap2 preset Emu uses. `map-hifi` is the PacBio HiFi preset and the default. Others: `map-ont`, `map-pb`, `sr`, `lr:hq`. |
| `--emu_threads` | 4 | CPU threads per `emu abundance` call |
| `--emu_min_abundance` | Emu default | minimum relative abundance threshold (e.g. `1e-5`) for Emu's noise trim |
| `--skip-emu-db-check` | off | use the existing Emu DB as is; never build |
| `--nb_db_dir` | `db_nb/` | directory for the two DADA2 NB references; both auto-built from the GTDB BLCA DB on first NB-enabled run (fast reformat, no training) |
| `--nb_min_bootstrap` | 80 | DADA2 `minBoot` for the genus-level `assignTaxonomy` step |
| `--skip-nb-classifier-check` | off | use the existing NB references as is; never build |
| `--gtdb-release` | `latest` | GTDB release to build (`latest` queries the GTDB site). Pass an explicit number (e.g. `232`) to pin for reproducibility. |
| `--min-ref-len` | 1000 | drop reference seqs shorter than this (0 keeps all) |
| `--qiime2_env` / `--blca_env` / `--emu_env` | auto | conda envs for `-profile standard`; auto-detected/installed. Supply paths only to reuse existing envs. `--qiime2_env` is also used for the QIIME2 export/tabulate steps in the NB branch. |
| `--blca_container` / `--emu_container` | image refs | per-branch container images for `-profile docker` / `singularity` |
| `--publish_dir_mode` | symlink | how outputs land in `<outdir>/`: `symlink` (default) or `copy` (standalone files) |
| `--progress-bar` | off | show a tqdm bar instead of Nextflow's native display |

## Taxonomy-only outputs

Each classifier writes a per-sequence taxonomy. NB runs `assignTaxonomy` plus
`addSpecies` and writes `best_taxonomy.tsv` and `gtdb_nb.tsv`; it skips the
ASV-by-sample frequency merge, because a standalone FASTA carries no frequency
table. Emu gives a per-sequence assignment whose abundance numbers are not
quantitative on dereplicated input.

## Notes

- BLCA scores each ASV independently, so chunk size only affects speed, not the taxonomy.
- BLCA runs with `-p 1` (the conda clustalo build has no OpenMP); parallelism comes from running chunks at the same time.
- NB classifies each ASV independently, so the result does not depend on how many ASVs are in the run.
