FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y \
    r-base \
    pandoc \
    wget \
    r-cran-devtools \
    r-cran-htmltools \
    r-cran-htmlwidgets \
    r-cran-yaml \
    r-cran-igraph \
    r-cran-magrittr \
    r-cran-plyr

RUN wget https://github.com/fbreitwieser/sankeyD3/archive/v0.2.tar.gz
RUN wget https://cran.r-project.org/src/contrib/networkD3_0.4.tar.gz

RUN R CMD INSTALL v0.2.tar.gz
RUN R CMD INSTALL networkD3_0.4.tar.gz

COPY krakenSankey.R /
WORKDIR /data
ENTRYPOINT ["Rscript", "/krakenSankey.R"]
