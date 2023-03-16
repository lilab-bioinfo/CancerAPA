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
* `tabix (v1.7)
* `bedtools (v2.30.0)`
* `samtools (v1.10)`
* `plink 1.9 beta`

## Installation guide
We highly recommend using conda to setup and manage the software environment. The instructions below introduce how to setup the running environment with conda.



## Outline of the analyses

### 1.GWAS clump and fine-mapping
* GWAS clump using [plink](https://www.cog-genomics.org/plink/) and fine-mapping using [CAUSALdb](https://github.com/mulinlab/CAUSALdb-finemapping-pip)

### 2.Estimate heritability and genetic correlation
* GWAS heritabiltiy estimates and genetic correlaiton using [LDSC](https://github.com/bulik/ldsc)

### 3. APA quantification
* APA quantification using [DarPars2]()

### 4. 3'aQTL mapping
* 3aQTL calling and fine-mapping using [CAVIAR](https://github.com/fhormoz/caviar)

### Analysis of 3'aQTL enrichment in cancer GWAS signals
* Integreation of cancer GWAS and 3aQTLs to assess the 3aQTL for their contributions to disease susceptibility
  * Partitioned heritabiltiy estimates using [LDSC](https://github.com/bulik/ldsc) and heritablity enrichment esitimates using [fgwas](https://github.com/joepickrell/fgwas)
  
### Build and run 3aTWAS model
  * Transcriptome-wide association analysis using [FUSION](http://gusevlab.org/projects/fusion/)
  




## Authors

Hui Chen, Zeyang Wang, Jia Wang, Wenyan Chen, Xuelian Ma, Xudong Zou, Mireya Plass, Cheng Lian, Ting Ni, Gong-Hong Wei,  Wei Li, Lin Deng, Lei Li

Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

## Citation
* Code and Execution:





## Contact
For any issues, please create a GitHub Issue.

## Funding
This work was supported by National Natural Science Foundation of China (no. 32100533) and startup funds from Shenzhen Bay Laboratory to L.L.
