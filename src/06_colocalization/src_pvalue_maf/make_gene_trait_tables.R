#/opt/apps/R/3.6.2/bin/R

# make_gene_trait_tables.R
# Caleb Matthew Radens
# cradens@mail.med.upenn.edu
# Last Update: 2016_1_26
# modified by YoSon

# This is a wrapper script that runs imports a GWAS trait table and all eQTL tables that
#   fall within 1e6 bp of sentinal SNP regions in the GWAS table, then exports tables
#   that are a merge of the GWAS and each eQTL table.

# I tested this script using the PMACS R module 3.1.2
system("echo ===========================",wait=FALSE)
system("echo inside make_gene_trait_tables.R",wait=FALSE)
R_ver <- substr(version$version.string,1,15) # Get R version
system(paste("echo",R_ver),wait=FALSE)
system("echo ===========================",wait=FALSE)

# Loads the required packages
require(data.table)

base <- getwd()  #edit by zhec: change path
setwd(base)
source("import.R")

system("echo wrapper package and script dependencies loaded",wait=FALSE)

# Extract cmd arguments (see get_command_args in import.R for details)
# Too Long, Didn't Read: you can add trailing arguments to this function if
#  called in PMACS.. those arguments are captured by get_command_args

args <- commandArgs(trailingOnly = TRUE)
# See function description in import.R
args <- get_command_args(args)

GWAS_file <- args[["file"]]
Trait <- args[["trait"]]
eQTL_directory <- args[["eQTL_dir"]]
sentinal_file <- args[["sentinal_file"]]


# Retrieve list of genes matching the pattern below and assign filepath to each gene
genes <- get_gene_names(directory= eQTL_directory, Pattern=".txt")

# Import the GWAS table (this will be slow)
start_time_GWAS<-proc.time()[3]
trait_table <- read_GWAS(File_path=GWAS_file,
                          Columns=9,        #edit by zhec: change columns=11 to 9 and change the col number of each parameters
                          Skip=1,
                          Sep="\t",
                          Chr_col=1,
                          Chr_pos_col=3,
                          Rsid_col=4,
                          Pos_col=2,
                          MAF_col=9,
                          N_col=7,
                          PV_col=8,
                          Beta_col=5,
                          Varbeta_col=6,
                          Var_is_SE = TRUE)
#
 end_time_GWAS<-proc.time()[3]
 time_GWAS<- end_time_GWAS-start_time_GWAS
 print("time_GWAS")
 print(time_GWAS)
# Identify the sentinal SNPs (those SNPs that are smaller than some cutoff
#   and not within range of each other)
# By default, look for SNPs that are at least range away from each other and
# are less than or equal to 5e-8
cutoff_PV <- 5e-8
range <- 1e6
range_gene <-1e5 ## modified from 1e6

trait_table_PVs <- read_GWAS(File_path=sentinal_file,
                             Columns=9,        #edit by zhec: change columns=11 to 9 and change the col number of each parameters
                             Skip=1,
                             Sep="\t",
                             Chr_col=1,
                             Chr_pos_col=3,
                             Rsid_col=4,
                             Pos_col=2,
                             MAF_col=9,
                             N_col=7,
                             PV_col=8,
                             Beta_col=5,
                             Varbeta_col=6,
                             Var_is_SE = TRUE)
trait_table_PVs<-subset(trait_table_PVs,MAF<0.99 & MAF>0.01)
#print(trait_table_PVs[1:3,])

# Import the microArray alignment table associated with the eQTL
microArray_alnTable_file <- "~/2021-10-31-cancer-GWAS/aQTL_coloc/output/microArray/final/allTissue_hg19_microArray.txt"    #edit by Lei: change path
microArray_alnTable<-get_microArray_table(microArray_alnTable_file)
#print (microArray_alnTable)

# Extract the TESs from microalignment table
# g <- microArray_alnTable$gene
# unique_gene_indeces <- order(g)[!duplicated(sort(g))]
# microArray_alnTable <- microArray_alnTable[unique_gene_indeces,]
# t <- microArray_alnTable$TEScolClasses[Pos_col]
# unique_TES_indeces <- order(t)[!duplicated(sort(t))]
# microArray_alnTable <- microArray_alnTable[unique_TES_indeces,]

genes <- merge(genes, microArray_alnTable, by=c("gene"), all = FALSE)

unique_genes <- unique(genes$gene)
#print(unique_genes)

start_time_findgene<-proc.time()[3]
# Determine which genes are within range bp of the sentinal SNPs
sig_genes <- c()
tss <- c()
sig_trait_table_PVs <- c()
for (Gene in unique_genes){
  gene_data <- genes[which(genes$gene==Gene),]
  g_chr <- unique(gene_data$chr)
  if (length(g_chr)>1){
    #stop(paste("More than 1 chr associated with:",Gene))   #edit by zhec, add #
    next
  }
  # Use the mode of the gene TSS from the table (if equal numbers of different TSSs,
  #   it picks one TSS at randon, I think)
  u_gene_tss <- unique(gene_data$TSS)
  g_tss <- u_gene_tss[which.max(tabulate(match(gene_data$TSS, u_gene_tss)))]
  # Calculate the range to search for P-values
  bottom<-g_tss-range_gene
  top<-g_tss+range_gene

  # If chromosome of gene and the top GWAS PV_SNP match:
  if (TRUE%in%grepl(g_chr,trait_table_PVs$chr)){
    # Get the loci from trait_table_PV that match the chromosome of the gene
    t_pos <- trait_table_PVs$chr_pos[which(grepl(g_chr,trait_table_PVs$chr))]
    t_pos <- as.integer(matrix(unlist(strsplit(t_pos,split=":")),ncol=2,byrow=TRUE)[,2])
    # If the g_tss is within 1000kb of the genome-wide significant PV:
    if (TRUE%in%(t_pos >= bottom & t_pos <= top)){
      # Add the gene to the sig_gene list if any of the trait PVs are close to the gene
      sig_genes <- c(sig_genes,Gene)
      tss <- c(tss, g_tss)
    }
  }
}
end_time_findgene<-proc.time()[3]
time_findgene<-end_time_findgene-start_time_findgene
print("time_findgene")
print(time_findgene)
# Make sure there was at least 1 gene added to the sig_gene list
if (length(sig_genes) == 0){
  # If not, stop the script and write the _completed file so wrapper.R may move on2
  write.table("hello world",file=paste(Trait,"_completed",sep=""))
  stop(paste(Trait,"has no significant PVs close to the genes tested."))
}

genes_to_analyze <- data.frame(gene=sig_genes,
                               TSS=tss,
                               stringsAsFactors=FALSE)
genes_to_analyze_ori <- merge(genes_to_analyze, genes, by=c("gene","TSS"))
genes_to_analyze <-genes_to_analyze_ori[!duplicated(genes_to_analyze_ori$gene),]
#print(genes_to_analyze)
gene_files_to_analyze <- genes_to_analyze$filepath
chromosomes_to_analyze <- genes_to_analyze$chr
n_genes_to_analyze <- length(genes_to_analyze$filepath)
#print("n_genes_to_analyze")
#print(n_genes_to_analyze)
gene_trait_tss <- paste(genes_to_analyze$gene,Trait,genes_to_analyze$TSS,sep="_")

start_time_merge<-proc.time()[3]
countSubNum <- 1
for (gene_row in seq(from=1,to=n_genes_to_analyze)){
  gene_table <- read_eQTL(File_path=gene_files_to_analyze[gene_row],
                           Columns=5,            #edit by zhec: change columns to 5 and change the col number of each parameters
                           Skip=1,
                           Sep="\t",
                           Gene=genes_to_analyze$gene[gene_row],
                           Chr_col=1,
                           Rsid_col=5,
                           Pos_col=2,
                           N_col=4 ,
                           PV_col=3,
                           Var_is_SE = FALSE)
  # # Merge GWAS and eQTL table columns:
  # #   only keep SNPs that have a hg19 position between GWAS and eQTL.
  #
  merged_table <- merge(gene_table,
                         trait_table,
                         by = c("chr_pos"),
                         all = FALSE)
  print(gene_row)
  print(gene_files_to_analyze[gene_row])
  #sys_arg <- paste("python /dfs4/weil21-lab1/leil22/Projects/2018-12-31-3QTL-immune/2020-08-17-Blueprint-coloc/merge.py ",
   #                gene_files_to_analyze[gene_row], " ",
   #               GWAS_file, " ",
   #                gene_trait_tss[gene_row],"_",chromosomes_to_analyze[gene_row],"_analyze_me"," ",
   #                sep="")
  #system(sys_arg)
  #if(countSubNum %% 6 !=0 && countSubNum != n_genes_to_analyze){  #add by zhec for controlling to wait every 50 scripts submitted
  #  system(sys_arg,wait=FALSE)
  #}else{
  #  system(sys_arg,wait=TRUE)  #add by zhec for controlling to wait every 50 scripts submitted
  #}
 # print(countSubNum)
 # print(genes_to_analyze$gene[gene_row])
  countSubNum<-countSubNum+1

  write.table(merged_table, file = paste(gene_trait_tss[gene_row],"_",chromosomes_to_analyze[gene_row],"_analyze_me",sep=""), quote=FALSE)
}
#while(length(list.files(pattern="_analyze_me"))<n_genes_to_analyze){
 # print(length(list.files(pattern="_analyze_me")))
#  Sys.sleep(5)
#}
#print("out")
write.table("hello world",file=paste(Trait,"_completed",sep=""))
end_time_merge<-proc.time()[3]
time_merge<-end_time_merge-start_time_merge
print("time_merge")
print(time_merge)
