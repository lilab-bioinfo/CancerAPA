#/opt/apps/R/3.6.2/bin/R

# intersect_analysis.R
# Caleb Matthew Radens
# cradens@mail.med.upenn.edu
# 2016_1_6

# This script plots the top results from wrapper.R (identifies high hypothesis 4
#  genes from summary_table.txt, and then runs coloc.abf on them, then plots them.)

# I tested this script using the PMACS R module 3.1.1
system("echo ===========================",wait=FALSE)
system("echo inside make_plots.R",wait=FALSE)
R_ver <- substr(version$version.string,1,15) # Get R version
system(paste("echo",R_ver),wait=FALSE)
system("echo ===========================",wait=FALSE)

# Choose a path to load R packages from:
lp <- "~/R/x86_64-conda_cos6-linux-gnu-library/3.6"

# Add your library path to the current session of R's library path variable
.libPaths(lp)

source("import.R")
require("data.table",lib=lp)

system("echo wrapper package and script dependencies loaded and checked",wait=FALSE)

if (length(list.files(pattern = "coloc_hepG2_chromHMM.wo"))==0){
  stop("Please make sure coloc_hepG2_chromHMM.wo is in the same directory as intersect_analysis.R or adjust script accordingly")
}

# Note: the snp_hyp4 column may have least 1 very small number that R converts to 0. That's OK in this case.
#  (that's why you get a big error about a very small number and using rmpfr etc...)
table<- fread(input="coloc_hepG2_chromHMM.wo",stringsAsFactors = FALSE)
setnames(table,1:24, c("chr", "start", "end", "gene_full", "trait", "hyp3", "hyp4", "snp_hyp4","eQTL_beta", "eQTL_lABF", "eQTL_credible",
                  "GWAS_beta", "GWAS_lABF","GWAS_credible", "both_credible",
                  "chr2", "start_chromHMM", "end_chromHMM", "state", "idunno", "idunno2", "overlap1", "overlap2", "idunno3"))

Cutoff_hyp4_1 <- 0.75
Cutoff_hyp4_2 <- 0.3
States<-c("1_Active_Promoter","2_Weak_Promoter", "3_Poised_Promoter",
 "4_Strong_Enhancer",       "5_Strong_Enhancer", "6_Weak_Enhancer",
 "7_Weak_Enhancer",       "8_Insulator",   "9_Txn_Transition",
 "10_Txn_Elongation", "11_Weak_Txn",  "12_Repressed",
 "13_Heterochrom/lo", "14_Repetitive/CNV", "15_Repetitive/CNV")
both_credibles <- table[which(table$both_credible),]
snp_hyp4_75 <- table[which(table$snp_hyp4>Cutoff_hyp4_2),]
snp_hyp4_75_gene_25 <- snp_hyp4_75[which(snp_hyp4_75$hyp4>Cutoff_hyp4_1),]


all_state_counts <- list()
both_cred_state_counts <- list()
snp_hyp4_75_state_counts <- list()
snp_hyp4_75_gene_25_state_counts <- list()
for (state in States){
  all_state_counts[[state]] <- length(which(table$state == state))
  both_cred_state_counts[[state]] <- length(which(both_credibles$state == state))
  snp_hyp4_75_state_counts[[state]] <- length(which(snp_hyp4_75$state == state))
  snp_hyp4_75_gene_25_state_counts[[state]] <- length(which(snp_hyp4_75_gene_25$state == state))
}

# Get the ratio of each state's count over the sum of all state counts
ratiolize <- function(List){
  names <- names(List)
  total <- 0
  for (name in names){
    total <- total + List[[name]]
  }
  for (name in names){
    List[[name]] <- round(100*(List[[name]] / total),1)
  }
  return(List)
}
all_state_counts <- ratiolize(all_state_counts)
both_cred_state_counts <- ratiolize(both_cred_state_counts)
snp_hyp4_75_state_counts <- ratiolize(snp_hyp4_75_state_counts)
snp_hyp4_75_gene_25_state_counts <- ratiolize(snp_hyp4_75_gene_25_state_counts)


png("State histograms.png", width=10, height=7.5, units="in", res=200)
layout(
  mat=matrix(c(1,2,3,4,5), 5, 1, byrow = TRUE)
)
par(mar=c(0,2,1,0))
boxplot(all_state_counts,xaxt="n",ylim=c(0,50))
legend("topright",legend = "All",bty="n")
par(mar=c(0,2,1,0))
boxplot(both_cred_state_counts,xaxt="n",ylim=c(0,50))
legend("topright",legend = "Both Credibles",bty="n")
par(mar=c(0,2,1,0))
boxplot(snp_hyp4_75_state_counts,xaxt="n",ylim=c(0,50))
legend("topright",legend = paste("SNP Hyp4 Above",100*Cutoff_hyp4_2),bty="n")
par(mar=c(0,2,1,0))
boxplot(snp_hyp4_75_gene_25_state_counts,las=2,ylim=c(0,50))
legend("topright",legend = paste("SNP Hyp4 Above",100*Cutoff_hyp4_2,"and Gene Hyp4 Above",100*Cutoff_hyp4_1),bty="n")
par(mar=c(0,0,0,0))
dev.off()

HDL <- table[which(table$trait=="HDL"),]
HDL_range <- c(range(abs(HDL$eQTL_lABF)),range(abs(HDL$GWAS_lABF)))

LDL <- table[which(table$trait=="LDL"),]
LDL_range <- c(range(abs(LDL$eQTL_lABF)),range(abs(LDL$GWAS_lABF)))

TC <- table[which(table$trait=="TC"),]
TC_range <- c(range(abs(TC$eQTL_lABF)),range(abs(TC$GWAS_lABF)))

TG <- table[which(table$trait=="TG"),]
TG_range <- c(range(abs(TG$eQTL_lABF)),range(abs(TG$GWAS_lABF)))

HDL_sig <- HDL[which(HDL$snp_hyp4>0.3),]
HDL_sig <- HDL_sig[which(HDL_sig$hyp4>0.75),]

LDL_sig <- LDL[which(LDL$snp_hyp4>0.3),]
LDL_sig <- LDL_sig[which(LDL_sig$hyp4>0.75),]

TC_sig <- TC[which(TC$snp_hyp4>0.3),]
TC_sig <- TC_sig[which(TC_sig$hyp4>0.75),]

TG_sig <- TG[which(TG$snp_hyp4>0.3),]
TG_sig <- TG_sig[which(TG_sig$hyp4>0.75),]

library(grid)

library(hexbin)

png("HDL_lABF_comparison_top_snps.png", width=10, height=7.5, units="in", res=200)
p<-plot(hexbin(abs(HDL$eQTL_lABF),abs(HDL$GWAS_lABF)),xlab="eQTL lABF",ylab="HDL lABF",main="HDL with Top SNPs",colramp= BTY)
pushHexport(p$plot.vp)
grid.points(abs(HDL_sig$eQTL_lABF), abs(HDL_sig$GWAS_lABF), pch=23, gp=gpar(col="skyblue",fill="red"))
dev.off()

png("LDL_lABF_comparison_top_snps.png", width=10, height=7.5, units="in", res=200)
p<-plot(hexbin(abs(LDL$eQTL_lABF),abs(LDL$GWAS_lABF)),xlab="eQTL lABF",ylab="LDL lABF",main="LDL with Top SNPs",colramp= BTY)
pushHexport(p$plot.vp)
grid.points(abs(LDL_sig$eQTL_lABF), abs(LDL_sig$GWAS_lABF), pch=23, gp=gpar(col="skyblue",fill="red"))
dev.off()

png("TC_lABF_comparison_top_snps.png", width=10, height=7.5, units="in", res=200)
p<-plot(hexbin(abs(TC$eQTL_lABF),abs(TC$GWAS_lABF)),xlab="eQTL lABF",ylab="TC lABF",main="TC with Top SNPs",colramp= BTY)
pushHexport(p$plot.vp)
grid.points(abs(TC_sig$eQTL_lABF), abs(TC_sig$GWAS_lABF), pch=23, gp=gpar(col="skyblue",fill="red"))
dev.off()

png("TG_lABF_comparison_top_snps.png", width=10, height=7.5, units="in", res=200)
p<-plot(hexbin(abs(TG$eQTL_lABF),abs(TG$GWAS_lABF)),xlab="eQTL lABF",ylab="TG lABF",main="TG with Top SNPs",colramp= BTY)
pushHexport(p$plot.vp)
grid.points(abs(TG_sig$eQTL_lABF), abs(TG_sig$GWAS_lABF), pch=23, gp=gpar(col="skyblue",fill="red"))
dev.off()


#=======================



png("HDL_beta_comparison_top_snps.png", width=10, height=7.5, units="in", res=200)
p<-plot(hexbin(abs(HDL$eQTL_beta),abs(HDL$GWAS_beta)),xlab="eQTL beta",ylab="HDL beta",main="HDL with Top SNPs",colramp= BTY)
pushHexport(p$plot.vp)
grid.points(abs(HDL_sig$eQTL_beta), abs(HDL_sig$GWAS_beta), pch=23, gp=gpar(col="skyblue",fill="red"))
dev.off()

png("LDL_beta_comparison_top_snps.png", width=10, height=7.5, units="in", res=200)
p<-plot(hexbin(abs(LDL$eQTL_beta),abs(LDL$GWAS_beta)),xlab="eQTL beta",ylab="LDL beta",main="LDL with Top SNPs",colramp= BTY)
pushHexport(p$plot.vp)
grid.points(abs(LDL_sig$eQTL_beta), abs(LDL_sig$GWAS_beta), pch=23, gp=gpar(col="skyblue",fill="red"))
dev.off()

png("TC_beta_comparison_top_snps.png", width=10, height=7.5, units="in", res=200)
p<-plot(hexbin(abs(TC$eQTL_beta),abs(TC$GWAS_beta)),xlab="eQTL beta",ylab="TC beta",main="TC with Top SNPs",colramp= BTY)
pushHexport(p$plot.vp)
grid.points(abs(TC_sig$eQTL_beta), abs(TC_sig$GWAS_beta), pch=23, gp=gpar(col="skyblue",fill="red"))
dev.off()

png("TG_beta_comparison_top_snps.png", width=10, height=7.5, units="in", res=200)
p<-plot(hexbin(abs(TG$eQTL_beta),abs(TG$GWAS_beta)),xlab="eQTL beta",ylab="TG beta",main="TG with Top SNPs",colramp= BTY)
pushHexport(p$plot.vp)
grid.points(abs(TG_sig$eQTL_beta), abs(TG_sig$GWAS_beta), pch=23, gp=gpar(col="skyblue",fill="red"))
dev.off()

# This is a really ugly way to get all of the above plots on one image...
# library(gridExtra)
# png("lABF_comparison.png", width=10, height=7.5, units="in", res=200)
# data_list <-list(HDL,LDL,TC,TG)
# legendlist <- c("HDL","LDL","TC","TG")
# plotList <- lapply(1:4, function(i) {
#   eQTL_lABF = abs(data_list[[i]]$eQTL_lABF)
#   GWAS_lABF = abs(data_list[[i]]$GWAS_lABF)
#   hexbinplot(eQTL_lABF ~ GWAS_lABF)
# })
# do.call(grid.arrange, c(plotList, ncol=2))
# dev.off()
