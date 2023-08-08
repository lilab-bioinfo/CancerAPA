#/opt/apps/R/3.6.2/bin/R

# get_coloc_summaries.R
# Caleb Matthew Radens
# cradens@mail.med.upenn.edu
# Last Update: 2016_1_5

system("echo ===========================")
R_ver <- substr(version$version.string,1,15) # Get R version
system("echo inside get_coloc_summaries")
system(paste("echo",R_ver))
system("echo ===========================")

# Choose a path to load R packages from:
lp <- "~/R/x86_64-conda_cos6-linux-gnu-library/3.6"

# Add your library path to the current session of R's library path variable
.libPaths(lp)

base <- getwd()
setwd(base)

# Get functions from import.R
source("import.R")
source("coloc_analysis.R")
require(coloc,lib=lp)

system("echo wrapper package and script dependencies loaded",wait=FALSE)

# Extract cmd argument
args <- commandArgs(trailingOnly = TRUE)
# See function description in import.R
args <- get_command_args(args)

merged_file <- args[["file"]]

merged_table<-read.table(merged_file,header=TRUE,stringsAsFactors = FALSE)

if(length(merged_table[,1])==0){
  write.table("hello world",file=paste(merged_file,"_completed_but_was_empty",sep=""))
  # Delete file once it has been imported
  #system(paste("rm",merged_file),wait=FALSE)
  stop(paste("Table was empty:",merged_file))
}

# Delete file once it has been imported
#system(paste("rm",merged_file),wait=FALSE)
merged_table<-merged_table[order(merged_table[,"PV_eQTL"]),]#edit by HYM
#remove duplicate snps
merged_table<-merged_table[!duplicated(merged_table[,"chr_pos"]),]#edit by HYM

cat(length(merged_table[,1])," shared SNPs to be analyzed for colocalisation ","\n",sep='')

# Build lists for coloc.abf() input
#   SNPs will be named according to their hg19 position

snp <- as.character(merged_table[,'chr_pos'])

PV_GWAS <- merged_table[,'PV_GWAS'] #edit by zhec
names(PV_GWAS) <- snp

PV_eQTL <- merged_table[,'PV_eQTL']  #edit by zhec
names(PV_eQTL) <- snp

#varbeta_eQTL <- merged_table[,'varbeta_eQTL']
#names(varbeta_eQTL) <- snp                    #edit by zhec

#varbeta_GWAS <- merged_table[,'varbeta_GWAS']
#names(varbeta_GWAS) <- snp                    #edit by zhec

merged_table$MAF[merged_table$MAF<0.001] <- 0.001
MAF <- merged_table[,'MAF']

names(MAF) <- snp

N_eQTL <- merged_table[,'N_eQTL']
names(N_eQTL) <- snp

N_GWAS <- merged_table[,'N_GWAS']
names(N_GWAS) <- snp

eQTL_list <- list(pvalues = PV_eQTL, type = "quant", N = N_eQTL, snp = snp)
GWAS_list <- list(pvalues = PV_GWAS, type = "quant", N = N_GWAS , snp = snp)

result <- coloc.abf(dataset1 = eQTL_list, dataset2 = GWAS_list, MAF = MAF)

summary<- result$summary

summary<-data.frame(nsnps=summary[1],
                    hyp0=summary[2],
                    hyp1=summary[3],
                    hyp2=summary[4],
                    hyp3=summary[5],
                    hyp4=summary[6],
                    stringsAsFactors = FALSE)
write.table(summary,paste(merged_file,"_summary",sep=""),row.names=FALSE)

write.table("hello world",file=paste(merged_file,"_completed",sep=""))

split_underscores <- unlist(strsplit(merged_file,split="_"))
chr <- matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-2] #edit by zhec
if(length(split_underscores)==8){ #edit by HYM
  NM<-matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-7] #edit by HYM
  NMloc<-matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-6] #edit by HYM
  gene <- paste(NM,NMloc,sep = "_") #edit by HYM
  trait <- matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-5] #edit by zhec
  }
if(length(split_underscores)==7){ #edit by HYM
  NM<-matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-6] #edit by HYM
  NMloc<-matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-5] #edit by HYM
  gene <- paste(NM,NMloc,sep = "_") #edit by HYM
  trait <- matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-4] #edit by zhec
  }
if(length(split_underscores)==6){ #edit by HYM
    gene <- matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-5] #edit by zhec
    trait <- matrix(split_underscores,ncol=length(split_underscores))[,length(split_underscores)-4] #edit by zhec
  }#edit by HYM

gene_trait <- paste(gene,trait,sep="_")


bedified <- bedify_coloc(result, chr, gene_trait)

write.table(bedified,file=paste(gene_trait,"_BEDIFIED",sep=""),row.names=FALSE,col.names=FALSE)
