# CancerAPA
[![Github Release](https://img.shields.io/badge/release-v1.0-brightgreen)](https://github.com/3UTR/CancerAPA)
[![python Release](https://img.shields.io/badge/python-2.7.14-brightgreen)](https://www.python.org/downloads/)
[![R Release](https://img.shields.io/badge/R-3.6.2-brightgreen)](https://cran.r-project.org/)
[![bedtools Release](https://img.shields.io/badge/bedtools-v2.25.0-brightgreen)](https://github.com/arq5x/bedtools2)
[![Samtools Release](https://img.shields.io/badge/samtools-v1.9-brightgreen)](http://www.htslib.org/)

In this study, we performed the first large-scale and systematic analysis assessing the genetic effects of APA on 25 cancer types in 49 human tissues from the Genotype-Tissue Expression(GTEx).

This repository contains all source code for the analyses in manuscript "[A distinct class of pan-cancer susceptibility genes revealed by alternative polyadenylation transcriptome-wide association study](https://medrxiv.org/cgi/content/short/2023.02.28.23286554v1) ".

## System requirements
### Hardware and operating system requirements
A workstation or computer cluster running a POXIS system (Unix, Linux, or macOS) is required (we used the Linux distribution Ubuntu 18.04.6 LTS). The minimum requirements for the test dataset are an 8-core processor, 16 GB of RAM, and 1 TB of hard disk space.

### Software and dependencies
R and dependent packages

* `R (>= v3.6.1)`
* `PEER (v1.3)`
* `data.table (1.14.2)`
* `coloc (v5.1.1)`
* `glmnet (4.1-3)`
* `methods`
* `MatrixEQTL(v2.1.0)`
* `RColorBrewer`
* `optparse (1.7.1)`
* `plink2R (1.1)`

Python and dependent packages
* `Python (version >= v2.7.14)`
* `numpy`
* `scipy`

External software
* `LDSC (v1.0.1)`
* `plink 1.9 beta`
* `CAUSALdb`
* `tabix (v1.7)`
* `bedtools (>= v2.30.0)`
* `samtools (>= v1.10)`
* `DaPars (v2.1)`
* `CAVIAR (v2.2)`
* `fgwas(0.3.6)`
* `GCTA (1.93.3beta2)`
* `GEMMA (0.98.5)`
* `FUSION`

## Installation guide (Time required: 30 min)
We highly recommend using conda to setup and manage the software environment. The instructions below introduce how to setup the running environment with conda.

* Install python and conda. Download installer script “Anaconda2-2019.10-Linux-x86_64.sh” from [Anaconda Repository](https://repo.anaconda.com/archive/). Python and conda will be installed in $HOME/anaconda2/bin by default.
```
> bash Anaconda2-2019.10-Linux-x86_64.sh
```
* Install other software and dependencies with conda
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
* Install CASUALdb
```
# Clone the repo:
> git clone https://github.com/mulinlab/CAUSALdb-finemapping-pip.git

# set up conda environment and download fine-mapping tools:
> cd bin
> bash 00_set_up.sh
> cd ..

# build reference panel:
> cd ref
> python 01_prepare_reference.py
> cd ..
```
* Install FUSION software package from github:
```
> wget -c  https://github.com/gusevlab/fusion_twas/archive/master.zip

```
## Demo data for running the codes

## Instructions for use


1. GWAS clump and fine-mapping
```
> bash gwas_clump.sh
> bash gwas_finemap.sh
```

2. Estimate heritability and genetic correlation
```
> bash gwas_heritability.sh
> bash gwas_genetic_cor_ldsc.sh
```

3. APA quantification and 3'aQTL mapping
We have developed a pipeline to analyze APA and call 3'aQTL before. Please refer to our [3aQTL-pipe](https://github.com/3UTR/3aQTL-pipe) for detailed instructions.

4. Enrichment of 3'aQTL in cancer GWAS signals
```
>bash enrich_aQTL_in_GWAS.sh
```

5. Build and run 3aTWAS model
```
> bash aTWAS_Model_fusion.sh
```

## Authors

Hui Chen, Wenyan Chen, Xudong Zou, Lei Li (Institute of Systems and Physical Biology, Shenzhen Bay Laboratory)

## Citation

* Hui Chen^, Zeyang Wang^, Jia Wang^, Wenyan Chen, Xuelian Ma, Xudong Zou, Mireya Plass, Cheng Lian, Ting Ni, Gong-Hong Wei, Wei Li#, Lin Deng#, Lei Li# 
**A distinct class of pan-cancer susceptibility genes revealed by alternative polyadenylation transcriptome-wide association study** **(2023)** medrxiv 10.1101/2023.02.28.23286554


* Xudong Zou, Ruofan Ding, Wenyan Chen, Gao Wang, Shumin Cheng, Qin Wang, Wei Li, Lei Li. **Using population-scale transcriptomic and genomic data to map 3' UTR alternative polyadenylation quantitative trait loci** ***STAR Protocols***,3(3):101566 **(2022)**.
DOI: https://doi.org/10.1016/j.xpro.2022.101566


* Lei Li, Kai-Lieh Huang, Yipeng Gao, Ya Cui, Gao Wang, Nathan D. Elrod, Yumei Li, Yiling Elaine Chen, Ping Ji, Fanglue Peng, William K. Russell, Eric J. Wagner & Wei Li.**An atlas of alternative polyadenylation quantitative trait loci contributing to complex trait and disease heritability** ***Nature Genetics***,53,994-1005 **(2021)**. DOI:https://doi.org/10.1038/s41588-021-00864-5


## Contact
For any issues, please create a GitHub Issue.

## Funding
This work was supported by National Natural Science Foundation of China (no. 32100533) and startup funds from Shenzhen Bay Laboratory to L.L.
