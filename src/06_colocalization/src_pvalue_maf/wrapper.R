#/opt/apps/R/3.6.2/bin/R

# wrapper.R
# Caleb Matthew Radens
# cradens@mail.med.upenn.edu
# Last Update: 2016_1_28
# Modified by YoSon

# This is a wrapper script that runs coloc.abf() on eQTL<->GWAS, genome wide
#  The output is
# chr|start|end|gene|trait|gene_trait_hyp0|gene_trait_hyp1|gene_trait_hyp2|gene_trait_hyp3|
#   gene_trait_hyp4|snp_hyp4|eQTL_lABF|eQTL_credible|GWAS_lABF|GWAS_credible|both_credible
#
# This script depends on import.R, coloc_analysis, get_coloc_summaries.R, and make_gene_trait_tables.R
#   (see those scripts for details)
# All three scripts (those mentioned ^^^ and make_gene_trait_tables.R) should be located in
#   the same directory.

# example usage:
# Rscript wrapper.R

# I tested this script using the PMACS R module 3.1.2
system("echo ===========================",wait=FALSE)
R_ver <- substr(version$version.string,1,15) # Get R version
system(paste("echo",R_ver),wait=FALSE)
system("echo ===========================",wait=FALSE)


# Start a stopwatch!
start_time <- proc.time()[3]

# Loads the required packages
require(data.table)
require(hash)
require(coloc)

args <- commandArgs(TRUE)

if (length(args)==0) {
  stop("input files are not specified", call.=FALSE)
}

base <- args[1]
setwd(base)

# Get functions from import.R
source("import.R")

system("echo wrapper package and script dependencies loaded and checked",wait=FALSE)


GWAS_directory <- args[3]

# Retreive list of traits from GWAS_directory and assign full filepaths to them
traits <- get_GWAS_traits(directory = GWAS_directory, Pattern = "GLGC")

eQTL_directory <- args[2]

sentinal_file <-args[4]



# For each trait, launch make_gene_trait_tables.R script
for (row in seq(from=1,to=length(traits$filepath))){
  file<-traits$filepath[row]
  trait<-traits$trait[row]
  sys_arg <- paste("Rscript make_gene_trait_tables.R ",
                   "--","file ",file, " ",
                   "--","trait ",trait, " ",
                   "--","eQTL_dir ",eQTL_directory, " ",
                   "--","sentinal_file ",sentinal_file," ",
                   sep="")
  # NOTE: wait=TRUE means R will submit a job, wait for it to finish, and then
  # continue. Setting it to FALSE utilizes the parallel computing capabilities of
  # the OS because it will execute the command and continue running.
  # (thus submitting as many system commands as there are traits, in this case)
  system(sys_arg,wait=FALSE)
}
#print("back to wrapper")
# Wait for the trait scripts to finish writing gene_trait tables to file
while(length(list.files(pattern="_completed"))<length(traits$trait)){
  Sys.sleep(5)
}

# Time check
elapsed_gt = proc.time()[3] - start_time

# Remove the files that have _complete in their names
#  (these files were used to determine with make_gene_trait_tables was finished)
for(completed in list.files(pattern="_completed")){
  system(paste("rm",completed))
}

# Get a list of the file names for the mergred gene_trait files
merged_table_files <- list.files(pattern="_analyze_me")

# Launch get_coloc_summaries for each merged gene_trait file
countSubNum <- 1      #add by zhec for controlling to wait every 50 scripts submitted
for (file in merged_table_files){

  sys_arg <- paste("Rscript get_coloc_summaries.R ",
                   "--","file ",file, " ",
                   sep="")
  # NOTE: wait=TRUE means R will submit a job, wait for it to finish, and then
  # continue. Setting it to FALSE utilizes the parallel computing capabilities of
  # the OS because it will execute the command and continue running.
  # (thus submitting as many system commands as there are traits, in this case)
  system(sys_arg)
  #if(countSubNum %% 50 !=0){  #add by zhec for controlling to wait every 50 scripts submitted
  #  system(sys_arg,wait=FALSE)
  #}else{
   # system(sys_arg,wait=TRUE)  #add by zhec for controlling to wait every 50 scripts submitted
  #}
  countSubNum<-countSubNum+1
}

# Wait for all of the get_coloc_summaries scripts to finish
#while ((length(list.files(pattern="_completed"))+length(list.files(pattern = "_completed_but_was_empty")))<length(merged_table_files)){
#  Sys.sleep(5)
#}

# Time check
elapsed_coloc = proc.time()[3] - start_time

# If, after attempting to merge an eQTL and GWAS table, the merge is empty,
#  make_gene_trait_tables writes a file that indicates this:
empty_files = list.files(pattern="_completed_but_was_empty")

# Delete files that have _completed in their name
#   (these files were used to determine when all the merged gene_traits had been analyzed)
for(completed in list.files(pattern="_completed")){
  system(paste("rm",completed), wait=TRUE)
}

pattern <- "_summary"
catch <- 0
# Initiate the final summary table
summary_table <- data.frame()

# The number of summary tables to add to the final table is equal to the number of
#  gene_traits analyzed minus the number of merged tables that ended up being empty:
n_summary_tables <- length(merged_table_files) - length(empty_files)

# Build the final summary table file.
# Note: previous version of summary_table.txt is overwritten!
file_name <- paste(base,"summary_table.txt",sep="/")
write.table(NULL,file_name,row.names=FALSE,col.names = FALSE)
while (catch < n_summary_tables){
  summary_files <- list.files(path=base, pattern = pattern)
  n_files_to_analyze <- length(summary_files)
  if (n_files_to_analyze==0){
    system("echo Waiting for summary files.", wait=FALSE)
    Sys.sleep(1)
    next
  }

  catch <- catch + n_files_to_analyze
  percent_done <- as.character(100*round(catch/n_summary_tables,digits=2))
  system(paste("echo percent done: ",percent_done,"%",sep=""),wait=FALSE)

  summary_split <- matrix(unlist(strsplit(summary_files,split="_")),ncol=7,byrow=TRUE) #edit by zhec
  for (index in seq(from=1,to=n_files_to_analyze)){
    gene_trait_tss <- paste(summary_split[index,1],summary_split[index,2],summary_split[index,3],sep="_")    #edit by zhec:tes->tss, name show be tss,it's mean it
    table <- read.table(summary_files[index],stringsAsFactors = FALSE, header = TRUE)
    #system(paste("rm",summary_files[index]))
    row.names(table) <- gene_trait_tss        #edit by zhec:tes->tss, name show be tss,it's mean it
    write.table(table,file=file_name,sep="\t",append=TRUE,col.names = FALSE,quote = FALSE)

  }
}

pattern <- "BEDIFIED"
catch <- 0

# The number of bed tables to add to the final table is equal to the number of
#  gene_traits analyzed minus the number of merged tables that ended up being empty:
n_bed_tables <- length(merged_table_files) - length(empty_files)

# Build the final bedfile.
# Note: previous version of coloc_bed_table.BED is overwritten!
file_name <- paste(base,"coloc_bed_table.BED",sep="/")
write.table(NULL,file=file_name,row.names=FALSE,col.names=FALSE)

while (catch < n_bed_tables){
  bed_files <- list.files(path=base, pattern = pattern)
  n_files_to_analyze <- length(bed_files)
  if (n_files_to_analyze==0){
    system("echo Waiting for bed files.", wait=FALSE)
    Sys.sleep(1)
    next
  }

  catch <- catch + n_files_to_analyze
  percent_done <- as.character(100*round(catch/n_bed_tables,digits=2))
  system(paste("echo percent done: ",percent_done,"%",sep=""),wait=FALSE)

  for (index in seq(from=1,to=n_files_to_analyze)){
    table <- read.table(bed_files[index],stringsAsFactors = FALSE,skip = 0)
    # Warning: don't set wait=False, or you'll potentially add a bed file 2+ times..
    #system(paste("rm",bed_files[index]))
    write.table(table,file="coloc_bed_table.BED",
                sep="\t",
                append=TRUE,
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
  }
}

# Stop a stopwatch!
end_time <- proc.time()[3]
elapsed <- end_time - start_time

system(paste("echo","Writing gene_tables elapsed time:",elapsed_gt),wait=FALSE)
system(paste("echo","Running coloc elapsed time:",elapsed_coloc),wait=FALSE)
system(paste("echo",n_summary_tables,"gene-by-traits were analyzed for colocalization. Total elapsed time:",elapsed),wait=FALSE)
