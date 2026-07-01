# Keep ALL container-image caches inside this cloned repo, not your $HOME.
#
# By default Apptainer/Singularity download image layers into
# $HOME/.apptainer/cache (or $HOME/.singularity/cache) and build SIFs in a temp
# dir under $HOME or /tmp -- which quickly fills a home quota on an HPC.
# SOURCE this file before running the pipeline to redirect every cache here:
#
#     source set_apptainer_cache.sh
#     python bin/run_pipeline.py --input samples.tsv --metadata metadata.tsv -profile singularity
#
# Make it permanent by adding the `source ...` line to ~/.bashrc or your job script.
# (Nothing to run/install -- this only sets environment variables; it must be
#  SOURCED, not executed, so the variables persist in your shell.)

# Directory containing this file = the cloned repo root. Works whether sourced
# from bash (BASH_SOURCE) or zsh ($0 is set to the sourced path).
_self="${BASH_SOURCE[0]:-$0}"
_here="$( cd "$( dirname "$_self" )" >/dev/null 2>&1 && pwd )"

# Nextflow's converted-SIF cache (matches nextflow.config's default location)
export NXF_SINGULARITY_CACHEDIR="$_here/singularity_cache"
# Apptainer/Singularity's own layer-download cache (both names, for either tool)
export APPTAINER_CACHEDIR="$_here/apptainer_cache"
export SINGULARITY_CACHEDIR="$_here/apptainer_cache"
# Temp space used while building/converting images (keeps it off $HOME and /tmp)
export APPTAINER_TMPDIR="$_here/apptainer_tmp"
export SINGULARITY_TMPDIR="$_here/apptainer_tmp"

# Repo-local conda cache for the itsxrust/EMITS envs (ITS marker runs them via
# conda, not the images). Keeping it here means an offline compute node reuses the
# envs the internet-connected login node built.
export NXF_CONDA_CACHEDIR="$_here/conda_cache"

mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$NXF_CONDA_CACHEDIR"

echo "[cache] container caches now live under: $_here"
echo "[cache]   NXF_SINGULARITY_CACHEDIR = $NXF_SINGULARITY_CACHEDIR"
echo "[cache]   APPTAINER_CACHEDIR       = $APPTAINER_CACHEDIR   (= SINGULARITY_CACHEDIR)"
echo "[cache]   APPTAINER_TMPDIR         = $APPTAINER_TMPDIR   (= SINGULARITY_TMPDIR)"
echo "[cache]   NXF_CONDA_CACHEDIR       = $NXF_CONDA_CACHEDIR"

unset _self _here
