# paralog2viewer: open-source tools only. ANNOVAR, InterVar, MViewer.CLI,
# HGMD, OMIM, Taiwan Biobank and all patient data are mounted at run time
# (/opt/annovar, /opt/intervar, and the directory passed with --data).
#
# Base image pinned by version tag. The digest is recorded after the first CI
# build (job summary of the `image` job); then replace the tag with
# tag@sha256:<digest>.
ARG BASE_IMAGE=mambaorg/micromamba:2.9.0-debian13-slim

FROM ${BASE_IMAGE} AS runtime
USER root
SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
ENV MAMBA_ROOT_PREFIX=/opt/conda \
    MPLCONFIGDIR=/opt/mplconfig

# Liberation Sans (Arial metrics) for figures; procps provides ps for process
# listings. Debian package versions follow the pinned base release.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        fonts-liberation fontconfig procps \
    && fc-cache -f \
    && rm -rf /var/lib/apt/lists/*

# Caller and analysis tools, exact versions.
RUN micromamba create -y -p /opt/recall -c conda-forge -c bioconda \
        --override-channels \
        clair3=2.0.3 sawfish=2.2.1 minimap2=2.31 samtools=1.24 \
        bcftools=1.24 pysam=0.24.1 matplotlib-base=3.10 python=3.11.16 \
        pytorch=2.13.0 longphase=2.0.2 numpy=2.4.6 blast=2.17.0 \
        pip=26.2.2 curl ca-certificates perl \
    && micromamba clean -a -y

# VEP in its own environment (Perl dependencies conflict with the caller env).
RUN micromamba create -y -p /opt/vep -c conda-forge -c bioconda \
        --override-channels ensembl-vep=116.1 perl-compress-raw-zlib=2.214 \
    && micromamba clean -a -y \
    && printf '#!/bin/sh\nPATH=/opt/vep/bin:$PATH exec /opt/vep/bin/vep "$@"\n' \
        > /usr/local/bin/vep \
    && chmod +x /usr/local/bin/vep

# Clair3 ONT model r1041_e82_400bps_sup_v520 (HKU). sha256 computed on
# 2026-11-13 from a download of the same URLs; sizes 8304746 and 11983414 bytes.
RUN model_url=https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/r1041_e82_400bps_sup_v520 \
    && mkdir -p /opt/models/recall \
    && /opt/recall/bin/curl -fsSL --retry 3 "${model_url}/pileup.pt" \
        -o /opt/models/recall/pileup.pt \
    && /opt/recall/bin/curl -fsSL --retry 3 "${model_url}/full_alignment.pt" \
        -o /opt/models/recall/full_alignment.pt \
    && printf '%s  %s\n' \
        ec73e7fd3d955d8d64efbab29ba81f45c6129774a6c18beecc67c5829611fb4b /opt/models/recall/pileup.pt \
        3d45d34401613df5c4dcd31f3c444e13338a5cb0c8beb395d4f6ab71682f25d3 /opt/models/recall/full_alignment.pt \
        > /opt/models/recall/SHA256SUMS \
    && sha256sum -c /opt/models/recall/SHA256SUMS

ENV PATH=/opt/annovar:/opt/recall/bin:/usr/local/bin:/usr/bin:/bin \
    CONDA_PREFIX=/opt/recall \
    BCFTOOLS_PLUGINS=/opt/recall/libexec/bcftools \
    MPLBACKEND=Agg

# The package goes last so code changes reuse the tool layers above. Runtime
# dependencies (pysam, numpy, matplotlib) are already satisfied by conda.
COPY pyproject.toml VERSION README.md LICENSE /src/
COPY src/ /src/src/
RUN pip install --no-cache-dir --no-deps /src \
    && rm -rf /src \
    && python -c "import matplotlib.font_manager as fm; fm._load_fontmanager(try_read_cache=False); fm.findfont('Liberation Sans', fallback_to_default=False)" \
    && chmod -R a+rX /opt/mplconfig \
    && paralog2viewer --help >/dev/null \
    && blastn -version >/dev/null \
    && run_clair3.sh --help >/dev/null \
    && sawfish --version >/dev/null \
    && vep --help >/dev/null

ENTRYPOINT []
CMD ["paralog2viewer", "--help"]

# Test stage: dev tools and the tests, run inside the image by CI.
FROM runtime AS test
WORKDIR /repo
RUN pip install --no-cache-dir pytest==9.2.0 ruff==0.17.5
COPY pyproject.toml VERSION README.md LICENSE ./
COPY src/ src/
COPY tests/ tests/
RUN pip install --no-cache-dir --no-deps -e .
