# Execution profiles (containers vs conda)

[← back to README](../README.md)

Pick the runtime with `-profile`:

- `-profile singularity` / `-profile docker` — analysis tools come from pinned
  container images (recommended on HPC / laptops).
- `-profile standard` — analysis tools come from conda environments the launcher
  auto-detects or auto-installs.

## Using conda instead of containers

If you cannot use containers, run with `-profile standard`. The launcher will
auto-detect (or auto-install) the QIIME2 amplicon env, the Emu env, and the BLCA
env, so you do not need to supply any flags:

```
python bin/run_pipeline.py --profile standard
```

If you already have envs you would rather reuse (for example on a shared HPC),
point the launcher at them explicitly:

```
python bin/run_pipeline.py --profile standard \
  --qiime2_env /path/to/qiime2-amplicon-env \
  --emu_env    /path/to/hifitax_emu \
  --blca_env   /path/to/blca-env
```

## macOS note

Conda works on Linux, but QIIME2 and BLAST have no Apple Silicon builds. On
macOS, reuse an existing QIIME2 amplicon env (add seqkit and csvtk to it) and a
clustalo/blast/biopython env, and pass them with `--qiime2_env` / `--blca_env`.
