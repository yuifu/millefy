# <img src="img/logo.svg" width="30%" align="right" /> Millefy


Millefy is a tool for visualizing __read coverage of single-cell RNA sequensing (scRNA-seq) datasets__ in genomic contexts. By dynamically and automatically reorder single cells based on 'locus-specific' pseudotime, Millefy highlights cell-to-cell heterogeneity in read covreage of scRNA-seq data.

Millefy is available as an R package and [a Docker image with JupyterLab](https://github.com/yuifu/datascience-notebook-millefy).

## Example of millefy plot

<img src="img/Pou5f1.png" width="100%" />


## Getting started
### Installation

```
install.packages(c("data.table", "dtplyr", "proxy", "viridisLite"), dependencies=TRUE)

source("https://bioconductor.org/biocLite.R"); biocLite()
biocLite(c('Rsamtools', 'GenomicRanges',  'rtracklayer', 'destiny'))

devtools::install_github("yuifu/millefy")
```

#### Requirements

- R (version 3.2.2 or higher)


### Quick example

To quickly learn to use Millefy, try out [Quick example](tutorial/Quick_example.md) with the included example dataset.


## Learning more
See [Tutorial](tutorial/Tutorial.md) for details.

