# Running HiFiTaxa on Gadi (NCI)

[← back to README](../README.md) · [generic offline guide](offline.md)

Gadi's compute and interactive nodes have **no outbound internet**; only login
(head) nodes do. So the workflow is: **do all the downloading on a login node,
then run the analysis in an interactive job.** Home quota is small — put conda,
the repo, the database, and the image cache on `/scratch/<proj>` or
`/g/data/<proj>`. Replace `<proj>` with your NCI project code throughout.

## On the login node (has internet)

### 1. Install conda + mamba

```bash
cd /scratch/<proj>/$USER
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /scratch/<proj>/$USER/miniforge3
source /scratch/<proj>/$USER/miniforge3/etc/profile.d/conda.sh
# Miniforge ships mamba. (If you used plain Miniconda: conda install -n base -c conda-forge mamba)
```

### 2. Load Singularity

```bash
module load singularity
```

### 3. Install HiFiTaxa, its driver env, and the Nextflow runtime

```bash
git clone https://github.com/alegarritano/HiFiTaxa.git
cd HiFiTaxa
mamba env create -f environment.yml
conda activate HiFiTaxa
source set_apptainer_cache.sh         # repo-local image + temp caches (keeps them off $HOME)

# Cache the pinned Nextflow runtime + any plugins here, so the offline compute
# node has them. NXF_HOME on /scratch is shared with the compute node.
export NXF_HOME="$PWD/.nextflow"
NXF_VER=24.10.9 nextflow info         # downloads the Nextflow runtime into NXF_HOME
```

### 4. Compile all container images

Compute nodes can't pull images, so cache them here **first** — the Emu database
build in the next step reuses the Emu image. Pull each into the repo-local cache
(full list and exact filenames in [offline.md](offline.md)).

> **Budget ~1–1.5 h on the first build** (one-time; cached after). The downloads
> are fast — the slow part is Apptainer converting each image to a `.sif` on
> `/scratch`: roughly QIIME2 amplicon **~45 min**, pb-16s-nf-tools **~15 min**,
> BLCA **~5 min**, Emu **~5 min**.

```bash
cd "$NXF_SINGULARITY_CACHEDIR"
singularity pull kpinpb-pb-16s-nf-tools-latest.img docker://kpinpb/pb-16s-nf-tools:latest
singularity pull 'quay.io-qiime2-amplicon@sha256-4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2.img' \
  docker://quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2
singularity pull ghcr.io-alegarritano-hifitax-1.0.0.img docker://ghcr.io/alegarritano/hifitax:1.0.0
# only if you'll run Emu:
singularity pull 'quay.io-biocontainers-emu@sha256-61ea3336f12d41930d73e57ce1b041bce48d66b4011a165bf1f0efce9d684777.img' \
  docker://quay.io/biocontainers/emu@sha256:61ea3336f12d41930d73e57ce1b041bce48d66b4011a165bf1f0efce9d684777
cd -
```

### 5. Download and build the GTDB database

```bash
bash bin/build_gtdb_db.sh 232 db
```

This downloads and parses GTDB **once**, then asks whether to drop sequences
shorter than 1000 bp and which formats to build — **BLCA**, **NB**, **Emu**.
Build the ones matching the `--classifier` you'll run. BLCA and NB build with the
driver env; the Emu database build reuses the Emu image you cached in step 4.
All of this needs internet, so do it here on the login node.

Quick (~5–10 min total): the GTDB download + BLCA parse/index is ~2–3 min, the
NB references ~1–2 min, and the Emu database ~2–3 min.

### 5b. Download and build the UNITE database (ITS marker only)

For fungal ITS runs (`--marker ITS`), download UNITE instead of GTDB. UNITE has no
`latest` API, so pick the release you want from the repository page, resolve its
DOI to the file URL through the PlutoF API, and pull it straight down:

```
# Browse https://unite.ut.ee/repository.php and take the DOI of the 'Current'
# general FASTA release for Fungi (v10.0 / 2025-02-19 -> 10.15156/BIO/3301229 here).
DOI=10.15156/BIO/3301229
URL=$(curl -s "https://api.plutof.ut.ee/v1/public/dois/?format=vnd.api%2Bjson&identifier=$DOI" \
      | python3 -c "import sys,json;print(json.load(sys.stdin)['data'][0]['attributes']['media'][-1]['url'])")
wget -O unite.tgz "$URL"
mkdir -p unite_src && tar xzf unite.tgz -C unite_src
FASTA=$(find unite_src -name 'sh_general_release_dynamic_*.fasta' | head -1)

bash bin/build_unite_blca_db.sh  "$FASTA" db_unite 0
bash bin/build_unite_dada2_db.sh db_unite/unite_BLCAparsed.fasta db_unite/unite_BLCAparsed.taxonomy db_unite/unite_full_singlestep_ref.fa.gz
bash bin/build_emits_db.sh       "$FASTA" db_unite
```

This builds all three UNITE references (BLCA + single-step NB + EMITS) into
`db_unite/`, all with the driver env (`makeblastdb`, `python3`), no container
needed. It needs internet, so do it here on the login node. (The itsxrust and
EMITS steps also need their images cached in step 4: pull
`ghcr.io/ayobi/itsxrust:latest` and, for EMITS, `ghcr.io/ayobi/emits:latest`.)

Quick (~2–5 min): the ~29 MB download plus the BLCA parse/index and the two
reformatted references.

### 5c. Warm the ITS conda envs (ITS marker only)

The fungal itsxrust and EMITS steps run via small conda envs that bundle nhmmer +
minimap2, not the images (the upstream images lack `ps` and nhmmer). Nextflow builds
them from `envs/itsxrust.yml` / `envs/emits.yml` on first use, so build them here on
the login node (internet) and the offline job reuses them:

```
source set_apptainer_cache.sh                 # sets NXF_CONDA_CACHEDIR (repo-local, shared with the job)
nextflow run . -profile test_its,singularity  # first run builds the two envs + validates the ITS path
```

No `module load hmmer`, no extra `-c` config: nhmmer ships inside the conda env.

## In an interactive job (offline)

### 6. Request the job

```bash
qsub -I -P <proj> -q normal \
  -l walltime=06:00:00,ncpus=16,mem=64GB,jobfs=100GB,storage=scratch/<proj>+gdata/<proj>
```

### 7. Re-prepare the shell and run

Each fresh job starts clean, so reload everything (no internet needed now):

```bash
source /scratch/<proj>/$USER/miniforge3/etc/profile.d/conda.sh
module load singularity
cd /scratch/<proj>/$USER/HiFiTaxa
conda activate HiFiTaxa
source set_apptainer_cache.sh
export NXF_HOME="$PWD/.nextflow"      # the Nextflow runtime you cached on the login node
export NXF_OFFLINE=true               # don't reach the network for Nextflow plugins/updates

python bin/run_pipeline.py \
  --input samples.tsv --metadata metadata.tsv \
  --classifier all \
  --skip-gtdb-check \
  --profile singularity --publish_dir_mode copy \
  --outdir results
```

`--skip-gtdb-check` (database already built), the pre-cached images, and the
cached `NXF_HOME` + `NXF_OFFLINE=true` mean the run never needs the network. Run
from `/scratch` so Nextflow's `work/` lands there, not on your home quota.
