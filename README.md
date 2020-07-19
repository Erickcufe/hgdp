# hgdp
Package with data cleaned from the Human Genome Diversity Panel (HGDP) dataset. 

# Source
Datasets are from Stanford U, contains ~ 660,918 tag SNPs (Illumina HuHap 650k), in autosomes, chromosome X and Y, the pseudoautosomal region and mitochondrial DNA, typed across 1043 individuals from all panel populations (Li JZ et al. Science 319: 1100-4, 2008). 

## Installation

```
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github("Erickcufe/hgdp", ref = "erick")
```
