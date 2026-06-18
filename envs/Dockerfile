# BLCA runtime image for HiFiTaxa (clustalo + BLAST + BioPython).
# The BLCA/parse scripts are NOT baked in — Nextflow stages them from the
# pipeline's bin/ onto the task PATH, so this image only carries the tools.
# Built linux/amd64 and pushed to GHCR by .github/workflows/build-blca-image.yml.
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="hifitax" \
      org.opencontainers.image.description="BLCA taxonomy tools (clustalo, blast, biopython) for HiFiTaxa"

USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.10 \
        biopython=1.84 \
        clustalo=1.2.4 \
        "blast>=2.14" \
        perl \
        procps-ng && \
    micromamba clean --all --yes
# procps-ng provides `ps`, which Nextflow needs to collect task metrics; without
# it the task wrapper aborts before the command runs.

# make the env tools available regardless of how the container is invoked
ENV PATH=/opt/conda/bin:$PATH

# sanity check at build time
RUN clustalo --version && blastn -version && ps --version && python -c "import Bio; print('biopython', Bio.__version__)"
