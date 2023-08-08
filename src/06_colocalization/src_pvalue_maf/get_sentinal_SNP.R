#/opt/apps/R/3.6.2/bin/R
base <- getwd()
setwd(base)
start_time <- proc.time()[3]
source("import.R")
args <- commandArgs(trailingOnly = TRUE)
args <- get_command_args(args)
GWAS_file <- args[["file"]]
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

cutoff_PV <- 5e-8
range <- 1e6

trait_table_PVs <- get_sentinal_snp_regions(trait_table, Range = range, Cutoff_PV = cutoff_PV)
write.table(trait_table_PVs,"sentinalSNP.txt",sep="\t",row.names=FALSE)
end_time <- proc.time()[3]
elapsed <- end_time - start_time
system(paste("echo","Find sentinal SNP elapsed time:",elapsed),wait=FALSE)
