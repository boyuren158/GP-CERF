FROM rocker/verse:4.2.3

LABEL org.opencontainers.image.authors="nkhoshnevis@g.harvard.edu"

# Linux internal commands
RUN apt-get update
RUN apt-get install -y --no-install-recommends vim
RUN apt-get install -y --no-install-recommends zsh
RUN apt-get install -y --no-install-recommends libgsl-dev
RUN apt-get clean

# R package installation from CRAN
RUN R -e "install.packages(c( \
    'covr', \
    'roxygen2', \
    'devtools'), repos='https://cloud.r-project.org')"

# R package installation from CRAN
RUN R -e "install.packages(c( \
    'parallel', \
    'xgboost', \
    'stats', \
    'MASS', \
    'spatstat.geom', \
    'logger', \
    'Rcpp', \
    'RcppArmadillo', \
    'ggplot2', \
    'cowplot', \
    'rlang', \
    'Rfast', \
    'SuperLearner', \
    'wCorr'), repos='https://cloud.r-project.org')"

ENV RENV_VERSION 0.15.1
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

