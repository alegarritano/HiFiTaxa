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

### 3. Install HiFiTaxa and its driver env

```bash
git clone https://github.com/alegarritano/HiFiTaxa.git
cd HiFiTaxa
mamba env create -f environment.yml
conda activate HiFiTaxa
source set_apptainer_cache.sh    # repo-local image + temp caches (keeps them off $HOME)
```

### 4. Compile all container images

Compute nodes can't pull images, so cache them here **first** — the Emu database
build in the next step reuses the Emu image. Pull each into the repo-local cache
(full list and exact filenames in [offline.md](offline.md)):

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

python bin/run_pipeline.py \
  --input samples.tsv --metadata metadata.tsv \
  --classifier all \
  --skip-gtdb-check \
  --profile singularity --publish_dir_mode copy \
  --outdir results
```

`--skip-gtdb-check` (database already built) + the pre-cached images mean the run
never needs the network. Run from `/scratch` so Nextflow's `work/` lands there,
not on your home quota.
