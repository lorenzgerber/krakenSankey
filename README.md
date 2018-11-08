# krakenSankey

### Commandline
This repository contains an R file that can be run with `Rscript` from the cli
to sankey html plots from kraken report files:
```
Rscript krakenSankey.R kraken.report
```
Functions in `krakenSankey.R` require the 
[sankeyD3](https://github.com/fbreitwieser/sankeyD3) package which itself
has a number of dependencies. 

### Docker
The repository also contains a `Dockerfile` to build a Docker image that already
contains all required R packages. The Docker image can also be obtained from 
the author's [DockerHub](https://hub.docker.com/r/lorenzgerber/krakensankey/).

The syntax for running dockerized:
```
docker run -v /host/path/to/krakenreport/:/data/ lorenzgerber/krakensankey kraken.report
```
The result file, `kraken.html` is written in the mounted directory.

### References / Acknowledgments
Most of the code in this repo is 1:1 extracted from Florian Breitwieser's repo
[pavian](https://github.com/fbreitwieser/pavian).


