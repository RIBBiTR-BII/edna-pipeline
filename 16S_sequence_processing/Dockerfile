# =============================================================================
# Dockerfile: QIIME2 + R pipeline environment
# Base: official QIIME2 amplicon image (includes conda + all QIIME2 deps)
# Adds: R, tidyverse, yaml, phylotools, openxlsx, yq
#
# Build:
#   docker build -t edna-pipeline .
#
# Run (interactive shell, with current dir mounted to /data):
#   docker run -it --rm -v $(pwd):/data my-pipeline bash
#
# NOTE: Pin the QIIME2 epoch below to match the version you use in production.
# =============================================================================

FROM quay.io/qiime2/amplicon:2025.4

# -----------------------------------------------------------------------------
# System labels
# -----------------------------------------------------------------------------
LABEL description="QIIME2 + R pipeline environment"
LABEL qiime2.version="2026.1"

# -----------------------------------------------------------------------------
# Environment: suppress interactive prompts during apt installs
# -----------------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# -----------------------------------------------------------------------------
# System dependencies required by R packages
#   - r-base: R itself (via apt, available in the base image's Ubuntu/Debian)
#   - libcurl4-openssl-dev: needed by curl / httr (tidyverse dependency)
#   - libssl-dev: needed by openssl (tidyverse dependency)
#   - libxml2-dev: needed by xml2 (tidyverse dependency)
#   - libfontconfig1-dev: needed by systemfonts (ragg/ggplot2 dependency)
#   - libharfbuzz-dev, libfribidi-dev: needed by textshaping (ggplot2)
#   - libfreetype6-dev, libpng-dev, libtiff5-dev, libjpeg-dev: image libs
#   - libzip-dev: needed by openxlsx
#   - pandoc: optional but useful if you render Rmarkdown
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libzip-dev \
    pandoc \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# yq — YAML processor for bash (mikefarah/yq, Go-based)
#
# The QIIME2 base image ships a Python-based yq (kislyuk/yq) via conda,
# which uses different syntax and will break scripts written for mikefarah/yq.
# We: (1) remove the conda yq, (2) install the correct static binary,
# (3) place it in /usr/bin so it wins regardless of how conda reorders PATH.
# -----------------------------------------------------------------------------
RUN conda run pip uninstall -y yq 2>/dev/null || true \
    && conda remove -n base --force -y yq 2>/dev/null || true
RUN wget -qO /usr/local/bin/yq \
    https://github.com/mikefarah/yq/releases/download/v4.44.2/yq_linux_amd64 \
    && chmod +x /usr/local/bin/yq \
    && ln -sf /usr/local/bin/yq /usr/bin/yq

# -----------------------------------------------------------------------------
# R packages
# - Installed in a single RUN layer to keep the image cache efficient.
# - phylotools is not on CRAN; it is installed from GitHub via remotes.
# - Set Ncpus to use parallel compilation where available.
# -----------------------------------------------------------------------------
RUN Rscript -e " \
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c('tidyverse', 'yaml', 'openxlsx', 'remotes'), \
                     Ncpus = parallel::detectCores()); \
    remotes::install_github('helixcn/phylotools'); \
"

# -----------------------------------------------------------------------------
# Working directory — mount your project here at runtime with -v $(pwd):/data
# -----------------------------------------------------------------------------
WORKDIR /data

# -----------------------------------------------------------------------------
# Default entrypoint: drop into bash so the conda environment is active
# (QIIME2 base image sets up conda init; bash will source it automatically)
# -----------------------------------------------------------------------------
CMD ["/bin/bash"]
