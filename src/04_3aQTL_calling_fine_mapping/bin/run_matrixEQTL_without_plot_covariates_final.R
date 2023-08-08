#!/usr/bin/env Rscript

# Be sure to use an up to date version of R and Matrix eQTL.


# print usage
usage <- function() {
    cat(
          'usage: run_matrixEQTL.R <tissue>
	run_matrixEQTL.R
    tissue       tissue
')
}
lp <- "~/R/x86_64-conda_cos6-linux-gnu-library/3.6"
.libPaths(lp)
# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
CIS_DISTANCE <- 1e6
TISSUE <- args[1]
#TISSUE <- 'GTEx_EBV_transformed'
DIR <- '~/aQTL_pipeline'
OUT_PREFIX <- TISSUE

# Check input args
if (is.na(CIS_DISTANCE)) {
  usage()
  quit(save='no', status=1)
}
if (is.na(TISSUE)) {
  usage()
  quit(save='no', status=1)
}
if (is.na(DIR)) {
  usage()
  quit(save='no', status=1)
}
if (is.na(OUT_PREFIX)) {
  usage()
  quit(save='no', status=1)
}


# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
# install.packages("MatrixEQTL")
library('MatrixEQTL')

## Location of the package with the data files.
base.dir = DIR;
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste0(base.dir, "/output_filter/", TISSUE, "_SNPs.vcf_genotype.vcf");
snps_location_file_name = paste0(base.dir, "/Input_SNP_locs/", "GTEx_snp_loc.txt");

# Gene expression file name
expression_file_name = paste0(base.dir, "/peer_pdui/",TISSUE, ".pdui.impute.txt.v2");
gene_location_file_name = paste0(base.dir, "/Input_gene_locs/", "gene_3UTR_loc_hg38.txt");

# Covariates file name
# Set to character() for no covariates
#covariates_file_name = character();
covariates_file_name = paste0(base.dir, '/peer_pdui/', TISSUE, ".pdui.peer.covariates_line6.txt");

# Output file name
output_file_name_cis = paste0(base.dir, '/output_QTL/', TISSUE, ".cis_eqtl_all_approachb.txt");
output_file_name_tra = paste0(base.dir, '/output_QTL/', TISSUE, ".trans_eqtl_approachb.txt");
output_figure_name_cis = paste0(base.dir, '/output_QTL/', TISSUE, ".cis_eqtl_genotype_info_approachb.pdf");
#pdf(output_figure_name_cis)


# Only associations significant at this level will be saved
#pvOutputThreshold_cis = 1e-2;
#pvOutputThreshold_tra = 1e-5;
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = CIS_DISTANCE;

## Load gene expression data
expr = read.csv(expression_file_name,header=TRUE,sep="\t",row.names=1,check.names = F)  # GxN where G is number of genes and N is number of samples
#expr = expr[, -1]


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 500000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

# extract the samples with both RNAseq and genotyped, and reorder both snps and expression
snps$ColumnSubsample(match(intersect(colnames(snps), colnames(expr)), colnames(snps)))
expr=expr[,intersect(colnames(snps), colnames(expr))]

message("# filter out SNPs with MAF<=0.01 ...");
# Note: here Minor allele is the minor one for the samples, not necessary the same one as the population)
maf.list = vector('list', length(snps))
na.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
  na.list[[sl]] = is.na(rowMeans(slice));
}
maf = unlist(maf.list)
na  = unlist(na.list)
cat('SNPs before filtering:',nrow(snps), "\n")
#snps$RowReorder(maf>0.01);  # remove
snps$RowReorder(maf>0.1 & dim(slice)[2]*maf>10)
cat('SNPs after filtering:',nrow(snps), "\n")

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

message("# loading SNP and gene position files...")
######################
snpspos = read.csv(snps_location_file_name, header = TRUE, sep="\t");
snpspos=snpspos[,1:3]
genepos = read.csv(gene_location_file_name, header = TRUE, sep="\t");

expr_subset = as.matrix(expr)
gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1      # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$CreateFromMatrix(expr_subset)

me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);


