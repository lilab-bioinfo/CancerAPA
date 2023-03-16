# CancerAPA
[![Github Release](https://img.shields.io/badge/release-v1.0-brightgreen)](https://github.com/3UTR/CancerAPA)
[![python Release](https://img.shields.io/badge/python-2.7.14-brightgreen)](https://www.python.org/downloads/)
[![R Release](https://img.shields.io/badge/R-3.6.2-brightgreen)](https://cran.r-project.org/)
[![bedtools Release](https://img.shields.io/badge/bedtools-v2.25.0-brightgreen)](https://github.com/arq5x/bedtools2)
[![Samtools Release](https://img.shields.io/badge/samtools-v1.9-brightgreen)](http://www.htslib.org/)

In this study, we performed the first large-scale and systematic analysis assessing the genetic effects of APA on 25 cancer types in 49 human tissues from the Genotype-Tissue Expression(GTEx).

This repository contains all source code for the analyses in manuscript ["Pan-cancer GWAS analysis".](https://medrxiv.org/cgi/content/short/2023.02.28.23286554v1)

## System requirements
### Hardware and operating system requirements
A workstation or computer cluster running a POXIS system (Unix, Linux, or macOS) is required (we used the Linux distribution Ubuntu 18.04.6 LTS). The minimum requirements for the test dataset are an 8-core processor, 16 GB of RAM, and 1 TB of hard disk space.

### Software and dependencies
R and dependent packages

* `R (>= v3.6.2)`
* `PEER (v1.3)`
* `data.table`
* `coloc(v5.2.1)`
* `glmnet`

Python and dependent packages
* `Python (version >= v2.7.14)`
* `numpy`
* `scipy`

External software
* `tabix (v1.7)`
* `bedtools (>= v2.30.0)`
* `samtools (>= v1.10)`
* `plink 1.9 beta`
* `CAVIAR (v2.2)`

## Installation guide (Time required: 30 min)
We highly recommend using conda to setup and manage the software environment. The instructions below introduce how to setup the running environment with conda.

1. Install python and conda. Download installer script “Anaconda2-2019.10-Linux-x86_64.sh” from [Anaconda Repository](https://repo.anaconda.com/archive/). Python and conda will be installed in $HOME/anaconda2/bin by default.
```
> bash Anaconda2-2019.10-Linux-x86_64.sh
```
2. Install other software and dependencies with conda
```
> conda install -c conda-forge r-base
> conda install -c conda-forge r-dplyr
> conda install -c bioconda r-optparse
> conda install -c bioconda Bioconductor-impute
> conda install -c bioconda r-peer
> conda install -c bioconda r-matrixeqtl
> conda install -c bioconda bedtools
> conda install -c bioconda plink=1.90
> conda install -c bioconda samtools
> conda install -c bioconda vcftools
> conda install -c bioconda tabix
```
## Demo data for running the codes

## Instructions for use


1.GWAS clump and fine-mapping
```
> bash gwas_clump.sh
> bash gwas_finemap.sh
```

2.Estimate heritability and genetic correlation
```
> bash gwas_heritability.sh
> bash gwas_genetic_cor_ldsc.sh
```

3. APA quantification and 3'aQTL mapping
We have developed a pipeline to analyze APA and call 3'aQTL before. Please refer to our [3aQTL-pipe]() for detailed instructions.

4. Enrichment of 3'aQTL in cancer GWAS signals
```
bash enrich_aQTL_in_GWAS.sh
```

  
5. Build and run 3aTWAS model
```
> bash aTWAS_Model_fusion.sh
```
  




## Authors

Hui Chen, Zeyang Wang, Jia Wang, Wenyan Chen, Xuelian Ma, Xudong Zou, Mireya Plass, Cheng Lian, Ting Ni, Gong-Hong Wei,  Wei Li, Lin Deng, Lei Li

Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

## Citation
* Code and Execution:





## Contact
For any issues, please create a GitHub Issue.

## Funding
This work was supported by National Natural Science Foundation of China (no. 32100533) and startup funds from Shenzhen Bay Laboratory to L.L.
