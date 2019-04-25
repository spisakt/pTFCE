FROM rocker/tidyverse
LABEL maintainer="Tamas Spisak"

RUN R -e "install.packages(c('devtools', 'roxygen2', 'mmand', 'oro.nifti'), repos='http://cran.us.r-project.org')"

# run tests

RUN git clone https://github.com/spisakt/pTFCE.git /home/tester/src/pTFCE
WORKDIR /home/tester/src/pTFCE
RUN R -e "devtools::test()"

# test if it can be installed and loaded
RUN R -e "devtools::install_github('spisakt/pTFCE@v0.1.3'); library(pTFCE))"

