# Databases (GTDB/Emu/NB for 16S, UNITE for ITS)

[← back to README](../README.md)

All three 16S classifiers share one GTDB reference build (default release r232).
On the first run the preflight builds it automatically; the sections below cover
pinning a release and building each reference by hand. For the ITS marker, UNITE
replaces GTDB (see [UNITE references](#unite-references-its-marker) at the end).

## Build everything interactively

`bash bin/build_gtdb_db.sh [release] [db-dir]` downloads and parses GTDB **once**,
then asks four questions — drop sequences <1000 bp, and whether to format for
**BLCA**, **NB**, and **Emu** — and builds only the formats you choose:

```
bash bin/build_gtdb_db.sh 232 db
```

Handy for preparing the database by hand (e.g. on an HPC login node before an
offline run). BLCA and NB build with the driver env (`makeblastdb`, `python`);
the Emu build needs `emu`, so it runs inside the Emu container if `emu` is not on
PATH. The per-format scripts below are the non-interactive primitives this
builder (and the launcher) call.

**Timing (~5–10 min total):** GTDB download + BLCA parse/index ~2–3 min (~76 MB
download, ~957k sequences kept after the 1000 bp filter), NB references ~1–2 min,
Emu database ~2–3 min. (Database builds are quick — the slow first-run step is
converting the container images to SIFs; see [offline.md](offline.md).)

## BLCA database

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

## Emu database (only when `--classifier` includes `emu`)

Lives in `--emu_db_dir` (default `db_emu/`). It is built from the same
BLCA-parsed GTDB files so both classifiers see the same reference release, by
remapping each GTDB lineage to a synthetic tax id and calling `emu
build-database`. On first use with an Emu-enabled run the launcher builds it;
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

## NB references (only when `--classifier` includes `nb`)

The DADA2 NB path needs **no training**. It uses two gzipped FASTA references,
built from the same BLCA-parsed GTDB that BLCA and Emu use, under `--nb_db_dir`
(default `db_nb/`):

- `gtdb_ssu_dada2_genus.fa.gz` (6 rank Kingdom to Genus, for `assignTaxonomy`)
- `gtdb_ssu_dada2_species.fa.gz` (`>acc Genus species`, for `addSpecies`)

Building is a fast reformat (about 1 to 2 minutes, low memory), not a training
step. On the first NB-enabled run the launcher builds both automatically; pass
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

## UNITE references (ITS marker)

With `--marker ITS` the pipeline classifies fungal ITS against UNITE instead of
GTDB. Three references are built into `--unite_db_dir` (default `db_unite/`) from
one UNITE general FASTA release, so BLCA, the single-step NB, and EMITS all anchor
on the same sequences:

- `unite_BLCAparsed.fasta` (+ BLAST index) and `unite_BLCAparsed.taxonomy`, for BLCA
- `unite_full_singlestep_ref.fa.gz` (full 7-rank lineage headers), for the single-step ITS Naive-Bayes
- `unite.fasta` (the raw UNITE release, staged verbatim), as the EMITS minimap2 target
- `UNITE_VERSION.txt`, recording the release tag

Unlike GTDB, UNITE has no machine-readable release API, so there is no
auto-download: you fetch the FASTA yourself. Download the **general FASTA release
for Fungi** from the UNITE repository page and pick the `Current` row (version
10.0, 2025-02-19 at time of writing):

<https://unite.ut.ee/repository.php>

On the first ITS run the preflight prints that link, asks for the path to your
downloaded `sh_general_release_dynamic_*.fasta`, and builds the three references
from it. Non interactively (batch/HPC), pass the path with `--unite-fasta <path>`;
`--skip-gtdb-check` uses whatever is already in `db_unite/` as-is.

Build by hand (the three primitives the preflight calls):

```
FASTA=~/Downloads/sh_general_release_dynamic_19.02.2025.fasta
bash bin/build_unite_blca_db.sh  "$FASTA" db_unite 0    # 0 keeps all refs; a number drops shorter ones
bash bin/build_unite_dada2_db.sh db_unite/unite_BLCAparsed.fasta db_unite/unite_BLCAparsed.taxonomy db_unite/unite_full_singlestep_ref.fa.gz
bash bin/build_emits_db.sh       "$FASTA" db_unite
```

These need only the driver env (`makeblastdb`, `python3`); no container, no
`qiime`. ITS uses a single-step Naive-Bayes (one `assignTaxonomy` straight to
species against the 7-rank reference) because UNITE species are largely
sequence-disjoint, so the two-step `addSpecies` recovery that 16S uses does not
apply.

### On an HPC with offline compute nodes (e.g. Gadi)

Prepare `db_unite/` on a login node before the job, the same as GTDB. Either build
it there from a FASTA you have copied across (the block above), or transfer a
`db_unite/` you built elsewhere:

```
# copy a prebuilt db_unite/ (-L dereferences the unite.fasta symlink)
rsync -avL db_unite/ user@gadi-dm.nci.org.au:/scratch/<proj>/$USER/HiFiTaxa/db_unite/
```

See [gadi.md](gadi.md) for the full offline workflow.
