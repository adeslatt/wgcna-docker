# Full contents of Dockerfile
FROM rocker/tidyverse
LABEL description="Base docker image with tidyverse and hence R and util libraries"
ARG ENV_NAME="wgcna"

# Install wgcna
# Since our base image is an R docker base we will use cran for instal

RUN apt-get update && \ 
    Rscript -e "install.packages('WGCNA')" && \
    apt-get clean -y


#
# Jupytext can be used to convert this to a a notebook
# as a script it can be used in a workflow
#
# to run as a script you can use Rscript wgcna.Rscript
ADD ./wgcna.R /usr/local/bin/

RUN chmod +x /usr/local/bin/wgcna.R


