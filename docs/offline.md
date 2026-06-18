# Offline / air-gapped nodes (pre-staging the container images)

[← back to README](../README.md) · [Gadi-specific version](gadi.md)

Many HPC interactive and compute nodes have **no outbound internet**, so the
images cannot be pulled at run time. The fix is to download them **once on a
login node that does have internet**, into a persistent cache, then point every
later run at that same cache. The pipeline uses these images:

| Image | Used by | Pulled when |
|-------|---------|-------------|
| `kpinpb/pb-16s-nf-tools:latest` | QC / import | always |
| `quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2` | denoising, NB | always |
| `ghcr.io/alegarritano/hifitax:1.0.0` | BLCA steps (`params.blca_container`) | `--classifier` includes `blca` (default) |
| `quay.io/biocontainers/emu@sha256:61ea3336f12d41930d73e57ce1b041bce48d66b4011a165bf1f0efce9d684777` | Emu steps (`params.emu_container`) | `--classifier` includes `emu` |

## Easiest: let the first run pull them on a login node

`NXF_SINGULARITY_CACHEDIR` is where Nextflow keeps the converted `.img` files
(default: `singularity_cache/` inside the repo; `set_apptainer_cache.sh` sets it
and the layer/temp caches to repo-local paths). On a **login node with
internet**, populate that cache once — the launcher pre-pulls every image:

```bash
source set_apptainer_cache.sh                 # persistent, repo-local caches
python bin/run_pipeline.py --profile singularity --reads_dir /path/to/reads
```

Then on the **offline interactive node**, source the same file (or export the
same `NXF_SINGULARITY_CACHEDIR`) and run normally. Nextflow finds the cached
`.img` files and never touches the network:

```bash
source set_apptainer_cache.sh                 # same cache dir as above
python bin/run_pipeline.py --input samples.tsv --metadata metadata.tsv \
  --skip-gtdb-check --profile singularity
```

To put the cache on shared storage instead of the repo, export a fixed path on
both nodes before sourcing/running, e.g.
`export NXF_SINGULARITY_CACHEDIR=/scratch/$USER/hifitaxa_sif`.

## Manual: pull each image yourself

If you'd rather stage them by hand, pull (on the internet-connected node) into
`$NXF_SINGULARITY_CACHEDIR` using Nextflow's exact cache filenames (`/` and `:`
become `-`, extension `.img`), so Nextflow reuses them instead of re-pulling:

```bash
source set_apptainer_cache.sh
cd "$NXF_SINGULARITY_CACHEDIR"
singularity pull kpinpb-pb-16s-nf-tools-latest.img \
  docker://kpinpb/pb-16s-nf-tools:latest
singularity pull 'quay.io-qiime2-amplicon@sha256-4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2.img' \
  docker://quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2
singularity pull ghcr.io-alegarritano-hifitax-1.0.0.img \
  docker://ghcr.io/alegarritano/hifitax:1.0.0
# only if you run Emu:
singularity pull 'quay.io-biocontainers-emu@sha256-61ea3336f12d41930d73e57ce1b041bce48d66b4011a165bf1f0efce9d684777.img' \
  docker://quay.io/biocontainers/emu@sha256:61ea3336f12d41930d73e57ce1b041bce48d66b4011a165bf1f0efce9d684777
```

After staging, run once on the connected node with `-profile singularity` to
confirm nothing re-downloads; if a file name is off, Nextflow re-pulls it and
caches it under the correct name.
