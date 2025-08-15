# Dockerfile.r
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Speed & stability helpers
ENV CRAN_REPO=https://cloud.r-project.org
WORKDIR /app

# Ensure Git for remotes and common dev libs
RUN apt-get update && apt-get install -y --no-install-recommends \
      libgit2-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
      build-essential \
    && rm -rf /var/lib/apt/lists/*

# Force modern C++ so packages using BH (Boost) and others compile cleanly
RUN mkdir -p /root/.R && printf '%s\n' \
  'CXX11STD=-std=gnu++14' \
  'CXX14STD=-std=gnu++14' \
  'CXX17STD=-std=gnu++17' \
  > /root/.R/Makevars

# Install CRAN + Bioc packages, then GitHub-only AnnotationGx, and verify
RUN R -q -e "options(timeout=1200, repos=c(CRAN='${CRAN_REPO}')); \
             install.packages(c('optparse','dplyr','remotes','BH','piano'))" \
 && R -q -e "if (!requireNamespace('BiocManager', quietly=TRUE)) \
               install.packages('BiocManager', repos='${CRAN_REPO}'); \
             BiocManager::install(c('Biobase','biomaRt','CoreGx','PharmacoGx','fgsea','Xeva'), \
                                  ask=FALSE, update=FALSE, \
                                  Ncpus=parallel::detectCores(), \
                                  build_vignettes=FALSE)" \
 && R -q -e "remotes::install_github('bhklab/AnnotationGx', upgrade='never', dependencies=TRUE)" \
 && R -q -e "pkgs <- c('optparse','dplyr','Biobase','biomaRt','CoreGx','PharmacoGx','fgsea','Xeva','AnnotationGx','piano'); \
             ok <- vapply(pkgs, requireNamespace, logical(1), quietly=TRUE); \
             if (!all(ok)) { cat('Missing after install: ', paste(pkgs[!ok], collapse=', '), '\n'); quit(status=1) }"

CMD ["bash"]
