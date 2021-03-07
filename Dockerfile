
FROM rocker/rstudio-stable:3.4.4
MAINTAINER Andrew Oliver, UCI
LABEL version="1.0.0"

ENV RENV_VERSION 0.9.3-106
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /project
COPY renv.lock renv.lock
RUN R -e 'renv::restore()'



