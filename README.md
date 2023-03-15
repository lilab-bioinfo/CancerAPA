# Pan-cancer GWAS APA analysis
Analysis shared by ["A distinct class of pan-cancer susceptibility genes revealed by alternative polyadenylation transcriptome-wide association study"](https://medrxiv.org/cgi/content/short/2023.02.28.23286554v1)
# Introduction
In this study, we performed the first large-scale and systematic analysis assessing the genetic effects of APA on 25 cancer types in 49 human tissues from the Genotype-Tissue Expression(GTEx) (v8).
# Outline
* GWAS clump using [plink](https://www.cog-genomics.org/plink/) and fine-mapping using [CAUSALdb](https://github.com/mulinlab/CAUSALdb-finemapping-pip)
* GWAS heritabiltiy estimates and genetic correlaiton using [LDSC](https://github.com/bulik/ldsc)
* APA quantification using [DaPars2](https://github.com/3UTR/DaPars2), aQTL calling using [matrixEQTL](https://github.com/andreyshabalin/MatrixEQTL) and fine-mapping on 3aQTL using [CAVIAR](https://github.com/fhormoz/caviar)
* Integreation of cancer GWAS and 3aQTLs to assess the 3aQTL for their contributions to disease susceptibility
  * Partitioned heritabiltiy estimates using [LDSC](https://github.com/bulik/ldsc) and heritablity enrichment esitimates using [fgwas](https://github.com/joepickrell/fgwas)
  * Transcriptome-wide association analysis using [FUSION](http://gusevlab.org/projects/fusion/)
  
## Prerequisites
1.Download and install python3
```
wget -c https://repo.anaconda.com/archive/Anaconda3-5.3.0-Linux-x86_64.sh
bash Anaconda3-5.3.0-Linux-x86_64.sh
```
2.Download and install other neccessary tools with conda
```
conda install -c r r-base
conda install -c conda-forge r-dplyr
conda install -c bioconda r-peer
conda install -c bioconda Bioconductor-impute
conda install -c bioconda matrixeqtl
conda install -c bioconda bedtools
conda install -c bioconda vcftools
conda install -c bioconda samtools
conda install -c bioconda plink=1.90
```
