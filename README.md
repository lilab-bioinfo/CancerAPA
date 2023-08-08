# CancerAPA
[![Github Release](https://img.shields.io/badge/release-v1.0-brightgreen)](https://github.com/3UTR/CancerAPA)
[![python Release](https://img.shields.io/badge/python-2.7.14-brightgreen)](https://www.python.org/downloads/)
[![R Release](https://img.shields.io/badge/R-3.6.2-brightgreen)](https://cran.r-project.org/)
[![bedtools Release](https://img.shields.io/badge/bedtools-v2.25.0-brightgreen)](https://github.com/arq5x/bedtools2)
[![Samtools Release](https://img.shields.io/badge/samtools-v1.9-brightgreen)](http://www.htslib.org/)

In this study, we performed the first large-scale and systematic analysis assessing the genetic effects of APA on 22 cancer types in 49 human tissues from the Genotype-Tissue Expression(GTEx v8).

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
* `MatrixEQTL(v2.1.0)`
* `RColorBrewer`
* `optparse (1.7.1)`
* `plink2R (1.1)`

Python and dependent packages
* `Python (version >= v2.7.14)`
* `numpy`
* `scipy`
* `argparse`

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

## Installation guide
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
> conda install -c bioconda r-coloc
> conda install -c bioconda bedtools
> conda install -c bioconda plink=1.90
> conda install -c bioconda samtools
> conda install -c bioconda vcftools
> conda install -c bioconda tabix
```
* Install CASUALdb for GWAS summary statistics based fine-mapping

  - Clone the repo:
  ```
  > git clone https://github.com/mulinlab/CAUSALdb-finemapping-pip.git
  ```
  - set up conda environment and download fine-mapping tools:
  ```
  > cd bin
  > bash 00_set_up.sh
  > cd ..
  ```
  - build reference panel:
  ```
  > cd ref
  > python 01_prepare_reference.py
  > cd ..
  ```

* Install LDSC for estimating heritability and genetic correlation from GWAS summary statisctics.

  - Clone the repo:
  ```
  > git clone https://github.com/bulik/ldsc.git
  > cd ldsc
  ```
  
  - Create an Anaconda environment with LDSC's dependencies:
  ```
  > conda env create --file environment.yml
  > source activate ldsc
  ```

* Install FUSION software package from github:
  - Download and unpack the FUSION software package from github:
  ```
  > wget -c  https://github.com/gusevlab/fusion_twas/archive/master.zip
  > unzip master.zip
  > cd fusion_twas-master
  ```
  - Download and unpack the 1000Genome LD reference data:
  ```
  > wget -c https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
  > tar xjvf LDREF.tar.bz2  
  ```
  - Download and unpack the plink2R library：
  ```
  > wget https://github.com/gabraham/plink2R/archive/master.zip
  > unzip master.zip
  ```
  - Launch R and install required libraries:
  ```
  install.packages(c('optparse','RColorBrewer'))
  install.packages('plink2R-master/plink2R/',repos=NULL)
  ```
  - the following steps are required for computing weights
    * Add the bundled `GCTA` binary to path
    * Add `plink` to path
    * Launch R and install the required libraries：
    ```
    install.packages(c('glmnet','methods'))
    ```
  - If using BSLMM, download and install `GEMMA` software, add to path. Generate a symbolic link to the output by calling ln -s ./ output in the directory where you will run FUSION.weights.R (this is a workaround because GEMMA requires results to go into an `output` subdirectory).

## Demo data for running the codes

Each analysis has a demo input data under the src/ folder.

## Instructions for use
**1. Definition of lead SNPs and trait-assocaited loci using `plink` clumping:**
```
> bash run_plink_clump.sh
```

**2. Fine mapping of trait-associated loci to defined credible SNPs:**

```
> bash run_finemapping.sh
```

**3. Estimate SNP heritability and genetic correlation:**

```
> bash run_ldsc_heritability.sh
> bash run_ldsc_genetic_correlation.sh
```

**4. APA quantification and 3'aQTL mapping:**

  We have developed a pipeline to analyze APA and call 3'aQTL before. Please refer to our [3aQTL-pipe](https://github.com/3UTR/3aQTL-pipe) for detailed instructions.

**5. Enrichment of 3'aQTL in cancer GWAS signals**
```
> bash LDSC-Partitioned-Heritability.sh
> bash run-fgwas-cancer.sh  
```
  
**6. Colocalization of trait-associated loci and 3'aQTL:**

```
> bash run_aQTL_colocalization.sh
```

**7. Build 3'TWAS model and run transcriptome-wide association analysis:**

* Compute GTEx v8 APA predictive models:
    ```
    Rscript FUSION_compute_weights.R \
    --bfile ${tissName}.${GNAME} \
    --tmp ${tissName}.${GNAME}.tmp \
    --out ${tissName}.${GNAME} \
    --PATH_plink $PLINK \
    --PATH_gcta $GCTA \
    --PATH_gemma $GEMMA \
    --models top1,blup,lasso,enet \
    --covar ${COVAR} \
    --hsq_p 0.05 \
    ```
    
* APA-based transcriptome wide association analysis
  - Input: GWAS summary statistics
        
    The input of GWAS summary statistics is in LD-score format, which mainly contain the following columns:
    1. `SNP` : SNP identifier(rsID)
    2. `A1` : first allele (effect allele)
    3. `A2` : second allel (other allele)
    4. `Z` : Z-score,sign with respect to `A1`

  - Input: APA weights
          
    The reference data are loaded using the `./WEIGHTS/${tissueName}.pos` which points to the individual `*RDat` weights files, their gene identifiers, physical positions. Only weights in the file will be evaluated.
      
   - Performing the APA level-disease association
        
        ```
        Rscript FUSION.assoc_test.R \
        --sumstats $WORK_DIR/GWAS/*.sumstats
        --weights $WORK_DIR/output/WEIGHTS/'${tissueName}'.pos \
        --weights_dir ./WEIGHTS/ \
        --ref_ld_chr ./LDREF/1000G.EUR. \
        --chr ${chr} \
        --out $WORK_DIR/results/'${GWAS}'/'${tissueName}'/'${GWAS}'.'${tissueName}'.chr${chr}.dat 
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
If you have any comments, suggestions, questions, etc, please feel free to create a GitHub Issue.

## Funding
This work was supported by National Natural Science Foundation of China (no. 32100533) and startup funds from Shenzhen Bay Laboratory to L.L.
