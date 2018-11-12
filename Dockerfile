FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y \
    r-base \
    r-cran-devtools \
    r-cran-plyr \
    pandoc

RUN R -e "devtools::install_github('fbreitwieser/sankeyD3')"
RUN R -e "install.packages('networkD3', repos = 'http://cran.us.r-project.org')"

COPY krakenSankey.R /

ARG UNAME=testuser
ARG UID=1000
ARG GID=1000
RUN groupadd -g $GID -o $UNAME
RUN useradd -m -u $UID -g $GID -o -s /bin/bash $UNAME
USER $UNAME

WORKDIR /data
ENTRYPOINT ["Rscript", "/krakenSankey.R"]
