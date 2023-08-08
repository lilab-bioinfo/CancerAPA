#/opt/apps/R/3.6.2/bin/R

# coloc_analysis.R
# Caleb Matthew Radens
# cradens@mail.med.upenn.edu
# 2016_1_8

# coloc_analysis.R contains functions that, given outputs from coloc_prep.py,
#   will use the R coloc package to identify potentially colocalizing SNPs
#   from eQTL and GWAS experiments.

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

require(coloc, lib=lp)
require(hash, lib=lp)
source("import.R")

system("echo wrapper package and script dependencies loaded",wait=FALSE)

batch_run_coloc_abf <- function(eQTL_directory,
                                GWAS_directory,
                                genes,
                                traits,
                                significant_only = TRUE,
                                summary_only = TRUE){
	# Run coloc_abf() on all combinations of genes and traits.
  #
  # If significant_only = TRUE, only run coloc.abf() if there is a genome-wide significant pvalue in the GWAS
  #   that is within 1e6 of the gene. See run_coloc_abf() for details
	#
	#	If summary_only = FALSE, the output will include both coloc.abf() summaries and results.
	#		By default, this function sets summary_only to TRUE. The results tables are much larger
	#		than a summary. See run_coloc_abf() for details
	#
	#	Arguments:
	#	eQTL_directory: points to CLI eQTL snptest .gz files
	#		typeof: character
	#
	#	GWAS_directory: points to GLGC Joint GWAS summary stat files
	#		typeof: character
	#
	#	genes: list of genes (CAPITALIZED!)
	#		typeof: character vector
	#		format: c("SORT1", "GENE2", ...)
	#		length: 1 or more
	#
	#	traits: list of traits
	#		typeof: character vector
	#		format: c("HDL", "LDL", ...)
	#		length: 1 or more
  #
  # significant_only: Boolean
	#
	#	summary_only: Boolean
	#
	# Return a dictionary whereby the keys are gene_trait and values are the coloc results
	#
	# NOTE:
	# Do not use this function for more than 100 genes or 4 traits when summary_only = FALSE!
	# 	Otherwise, the output will start to get really really big.
	# 	(100 genes/4 traits was an arbitrarily chosen max: if you realllly want to increase these
  #   defaults, just change them below:)

  MAX_genes <- 300
  MAX__traits <- 4

	print("Beginning batch import")

	# Prevent function from being run with too many gene<->trait combos when summary_only=FALSE
	if (!summary_only){
		if (length(genes) > MAX_genes){
			stop(paste("The number of genes is ",length(genes),", which is too many ","Do not run this function with too many genes; it will use too much memory!"))
		}
		if (length(traits) > MAX__traits){
			stop(paste("The number of traits is ",length(traits),", which is too many ","Do not run this function with too many traits; it will use too much memory!"))
		}
	}
	# Start the clock!
	batch_time <- proc.time()

	result_dict <- hash()

	# Extract the eQTL gene and GWAS trait tables to be run in coloc_abf
	# Start clock!
	gene_extract_time <- proc.time()

	gene_tables <- import_genes(genes, eQTL_directory)

	print("Gene extraction time elapsed:")
	# Lap update!
	print(proc.time() - gene_extract_time)

	# Start clock!
	GWAS_extract_time <- proc.time()

	trait_tables <- import_GWAS_tables(traits, GWAS_directory) #edit by zhec

	print("GWAS extraction time elapsed:")
	# Lap update!
	print(proc.time() - GWAS_extract_time)

	# Start clock!
	dictionary_build_time <- proc.time()

	for (trait in names(trait_tables)){
  	for (gene in names(gene_tables)){
  	  print(paste("Coloc analyzing:",gene,"and",trait,sep=" "))

  	  summary_and_results <- run_coloc_abf(gene_tables[[gene]], trait_tables[[trait]], significant_only = significant_only)
  	  gene_trait <- paste(gene,trait,sep="_")
  	  # Expected output of run_coloc_abf is list of length 2 [$summary and $results]
  	  if (length(summary_and_results)==2){
  	    if (summary_only){
  	      result_dict[gene_trait] <- summary_and_results$summary
  	    } else{
  	      result_dict[gene_trait] <- summary_and_results
  	    }
  	    # If significant_only = TRUE, the summary table will indicate gene-by-trait was not signiciant:
  	  } else if (summary_and_results == "NO_SIGNIFICANCE"){
  	      result_dict[gene_trait] <- "NO_SIGNIFICANCE"
  	  }
  	}
	}
	print("dictionary building time elapsed:")
	print(proc.time() - dictionary_build_time)

	print("Batch run time elapsed:")
	# Lap update!
	print(proc.time() - batch_time)
	return(result_dict)
}

run_coloc_abf <- function(eQTL_table, GWAS_table, significant_only = TRUE, GENOME_WIDE_SIG = 5e-8) {
  # Given eQTL and GWAS data with beta/varbeta/N and pvalues/N/MAF, respectively,
  #   utilize coloc.abf() to calculate the probability of hypothesis 4:
  #   (the eQTL and GWAS share a causal variant)
  #
  # Assumptions:
  #   pvalues is a character vector that needs to be converted to a double vector
  #   0 < pvalues <= 1
  #   0 < MAF <= 0.5
  #
  # Notes:
  #   MAF is derived from the GWAS because the GWAS has a much greater
  #     sample size than the eQTL and is a more accurate estimate of MAF.
  #
  # Returns coloc.abf() output

  # Merge GWAS and eQTL table columns:
  #   only keep SNPs that have a shared hg19 position between GWAS and eQTL.
  merged_table <- merge(eQTL_table, GWAS_table, by = c("chr_pos"), all = FALSE)  #edit by zhec

  cat(length(merged_table[,1])," shared SNPs to be analyzed for colocalisation ","\n",sep='')

  # Build lists for coloc.abf() input
  #   SNPs will be named according to their hg19 position
  snp <- as.character(merged_table[,'chr_pos'])  #edit by zhec: position->chr_pos

 # beta <- merged_table[,'beta']
  #names(beta) <- snp                    #edit by zhec: add #

  PV_eQTL <- merged_table[,'PV_eQTL']    #edit by zhec: add PV_eQTL
  names(PV_eQTL) <- snp

  pvalues <- merged_table[,'PV_GWAS']    #edit by zhec: pvalue->PV_GWAS
  # Due to buggy nature of fread(), importing tables with very small numbers crashes R.
  #   to get around this, pvalues was coerced into a character during import.

  # Convert pvalues from character to double.
  #   Note: R cannot handle numbers smaller than ~1e-300, so this function
  #   converts these numbers to 1e-300 instead of letting R turn them into 0s.
  pvalues<-not_too_tiny(pvalues)

  # If coloc should only run on genes for which there is a genome-wide significant p-value
  # in the same region in the GWAS:
  if (significant_only){
    if (min(pvalues)>GENOME_WIDE_SIG){
      # No genome-wide significant SNPs are close to this gene for this trait
      #   Return null
      cat("NO_SIGNIFICANCE","\n")
      return("NO_SIGNIFICANCE")
    }
  }

  names(pvalues) <- snp

  MAF <- merged_table[,'MAF']
  names(MAF) <- snp

  N_eQTL <- merged_table[,'N_eQTL']
  names(N_eQTL) <- snp

  N_GWAS <- merged_table[,'N_GWAS']
  names(N_GWAS) <- snp

  eQTL_list <- list(pvalues = PV_eQTL, N = N_eQTL,  type = "quant", snp = snp)    #edit by zhec: delete beta and varbeta, add pvalues.
  GWAS_list <- list(pvalues = pvalues, N = N_GWAS,  type = "quant", snp = snp)

  result <- coloc.abf(dataset1 = eQTL_list, dataset2 = GWAS_list, MAF = MAF)

  return(result)
}

coloc_abf_top_snps <- function(coloc_result, cutoff_PP = 0.5){
  # Given one output from coloc.abf(), returns a table of SNPs with
  #   hypothesis 4 posterior probabilities (PP). Only returns SNPs
  #   whose PP are greater than cutoff_PP, which is default set to 0.5.
  #
  #	Arguments:
  #		coloc_result: one output from coloc.abf() with
  #		  typeof: list()
  #			elements:
  #				$summary: coloc.abf() summary
  #				$results: coloc.abf() results
  #
  # Returns:
  #		Table of SNPs with a hyp4_PP > cutoff_PP
  #		snp | SNP.PP.H4

  # Extract hypothesis 4 SNP PPs
  hyp_4_PP <- coloc_result$results$SNP.PP.H4

  # Rank the PPs
  pp_rank <- cbind(hyp_4_PP,as.numeric(factor(hyp_4_PP)))

  # The row index from the result$PP corresponds to the row index of result$snp..
  result_index <- seq(from=1, to=length(hyp_4_PP), by=1)

  # Extract PP
  hyp_4_PP <- pp_rank[,1]

  # Extract the PP rank (higher is better)
  rank <- pp_rank[,2]

  # Extract the snp names
  snps <- coloc_result$results$snp

  # Building table with snp names, hyp_4_PP, and rank
  combined_results  <- data.frame(snps,hyp_4_PP, stringsAsFactors = FALSE)

  # Rank the results base don their hyp_4_PP (higher is better)
  ranked_combined_results <- combined_results[order(hyp_4_PP,decreasing=TRUE),]

  # Remove the row names
  rownames(ranked_combined_results)=NULL

  # Find indeces of SNPs that have hyp4 PP > cutoff_PP
  cutoff_SNPs <- which(ranked_combined_results$hyp_4_PP>cutoff_PP)
  if(length(cutoff_SNPs)==0){
    return(list(snps="No SNPs"))
  }

  # Update results with cutoff SNPs
  ranked_combined_results<-ranked_combined_results[cutoff_SNPs,]
  return(ranked_combined_results)
}

batch_credibility_set <- function(Result_dict, Cutoff = 0.95){
  # Cycles through the result_dictionary generated from batch_run_coloc_abf
  # and runs credibility_set on each result.
  #
  # Suggested usage:
  # > output_from_batch_run_coloc_abf <- batch_credibility_set(output_from_batch_run_coloc_abf)
  #     Using it in this way overwrites the output with the updated results
  coloc_results <- names(Result_dict)
  for (result in coloc_results){
    Result_dict[[result]] <- credibility_set(Result_dict[[result]], Cutoff)
  }
  return(Result_dict)
}

credibility_set <- function(Coloc_result, Cutoff = 0.95) {
  # Given the results from a coloc analysis, return the results inluding
  #  a sublist of of SNPs that  are, as a set, Cutoff percent likely to
  #  contain the causal variant.
  #  Based on Peter Donnelly's 2012 Nature paper (doi:10.1038).
  if(!(length(Coloc_result)==2)
     && !names(Coloc_result)[1]=="summary"
     && !names(Coloc_result)=="Result"){
    stop("Expected a default coloc.abf() output, but that's not what was given.")
  }
  copy <- Coloc_result

  # get number of SNPs in coloc result
  n_snps <- length(Coloc_result$results$lABF.df1)

  # Get list of snps
  if (FALSE %in% grepl(pattern=":",x=Coloc_result$results$snp)){
    stop("Expected coloc_result$results$snp to be formatted as such: 'chr#:####'")
  }

  split  <- matrix(unlist(strsplit(x=as.character(Coloc_result$results$snp), split=":")),ncol=2,byrow=TRUE) #edit by zhec
  all_snps <- split[,2]

  # Get log(ABF)s, take antilog
  abf_df1 <- data.frame(abf=not_too_big(10^Coloc_result$results$lABF.df1),index=seq(from=1,to=n_snps))

  attach(abf_df1)
  # Sort by ABF
  abf_df1 <- abf_df1[order(-abf),]
  detach(abf_df1)

  # Get log(ABF)s, take antilog
  abf_df2 <- data.frame(abf=not_too_big(10^Coloc_result$results$lABF.df2),index=seq(from=1,to=n_snps))

  attach(abf_df2)
  # Sort by ABF
  abf_df2 <- abf_df2[order(-abf),]
  detach(abf_df2)

  sum_abf_df1 <- sum(abf_df1$abf)
  sum_abf_df2 <- sum(abf_df2$abf)

  probability_1 = 0
  set_1 = c()

  probability_2 = 0
  set_2 = c()

  i <- 1
  # Add ABFs to the credibile set until the set includes Cutoff percent of the total ABF
  while (probability_1 < Cutoff && i<n_snps){
    # Update percent of ABF accounted for so far in credibility set 1
    probability_1 <- probability_1 + abf_df1$abf[i]/sum_abf_df1

    # Add next highest SNP to credible set 1
    set_1 <- c(set_1, all_snps[abf_df1$index[i]])
    i <- i+1
  }

  i <- 1
  # Add ABFs to the credibile set until the set includes Cutoff percent of the total ABF
  while (probability_2 < Cutoff && i<n_snps){
    # Update percent of ABF accounted for so far in credibility set 2
    probability_2 <- probability_2 + abf_df2$abf[i]/sum_abf_df2

    # Add next highest SNP to credible set 2
    set_2 <- c(set_2, all_snps[abf_df2$index[i]])
    i <- i+1
  }

  overlap <- duplicated(c(set_1,set_2))

  overlap <- c(set_1,set_2)[overlap]

  # Add a sublist to the coloc result list (it already has $summary and $result)
  copy$credible_snps <- list(credible_1=list(probability = probability_1, snps = set_1),
                                     credible_2=list(probability = probability_2, snps = set_2),
                                     in_both=overlap)
  return(copy)
}

batch_bedify_coloc <- function(Result_dict, Cutoff_PP = 0.75, Write_to_file = FALSE, File_out = "coloc_bed_table.BED"){
  # Given the output from batch_run_coloc_abf, generate gene_trait bed tables and
  #  returns a list of the bedified results.
  #
  # If Write_to_file = TRUE:
  #   This funciton writes all of the bedtables to a single file.
  #   Note: there aren't column headers in the bed file. See bedify_coloc for column
  #     name information.

  gene_traits <- names(Result_dict)

  # Identify which gene_traits are significant
  significant_gene_traits <- c()
  for (gene_trait in gene_traits){
    if (Result_dict[[gene_trait]]$summary[6]>=Cutoff_PP){
      significant_gene_traits <- c(significant_gene_traits, gene_trait)
    }
  }

  # String split all of the "gene.#_trait" to grab the gene names
  genes <- unique(matrix(unlist(strsplit(significant_gene_traits,split="\\.")),ncol=2,byrow=TRUE)[,1])

  gene_dict <- get_gene_TssTes(Genes=genes,only_TssTes=FALSE)

  ngenes <- length(genes)
  # Extract chromosome info from gene dictionary
  chrs <- sapply(1:ngenes, function(index, dictionary) dictionary[[genes[index]]]$chr,dictionary=gene_dict)

  bedified_results <- list()
  chromosomes <- c()
  n_siggenes <- length(significant_gene_traits)
  if(Write_to_file){
    if (length(list.files(pattern=File_out))>0){
      stop(paste(File_out,"already exists. Please move, rename, or delete it first."))
    }
  }
  for (sig_index in 1:n_siggenes){
    sig_g_trait <- significant_gene_traits[sig_index]
    matched <- sapply(1:ngenes, function(index, to_match) grepl(genes[index], to_match), to_match=sig_g_trait)
    matched<-genes[matched]
    matched_chr <- gene_dict[[matched]]$chr
    chromosomes <- c(chromosomes, matched_chr)
    bedified <- bedify_coloc(Result_dict[[sig_g_trait]], matched_chr, Gene_trait = sig_g_trait)
    if (!Write_to_file){
      bedified_results[[sig_g_trait]] <- bedified
    } else{
      write.table(bedified,file=File_out,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t",quote = FALSE)
    }
  }

  if (!Write_to_file){
    return(bedified_results)
  }
}

bedify_coloc <- function(Coloc_result, Chromosome, Gene_trait){
  # Given an output from coloc.abf() and the chromosome,
  # build a bed table that includes individual SNP level data for:
  #       chr: chromosome
  #       start: loci of SNP
  #       end: start+1
  #       gene
  #       trait
  #       gene_trait_hyp3: overall hypothesis 3 PP
  #       gene_trait_hyp4: overall hypothesis 4 PP
  #       Hyp4: Each SNP's individual hypothesis 4 PP
  #       eQTL_beta: beta from the eQTL
  #       eQTL_lABF: Each SNP's log ABF for the eQTL
  #       eQTL_credible: boolean: is the SNP in the eQTL 95% credibility set?
  #       GWAS_beta: beta from the GWAS
  #       GWAS_lABF: Each SNP's log ABF for the GWAS
  #       GWAS_credible: boolean: is the SNP in the GWAS 95% credibility set?
  #       both_credible: boolean: is the SNP in both 95% credibility sets?
  #

  Coloc_result <- credibility_set(Coloc_result, .95)
  #print(Coloc_result)

  split_gt <- matrix(unlist(strsplit(Gene_trait,split="_")),ncol=2)
  gene <- split_gt[,1]
  trait <- split_gt[,2]

  # Get the number of SNPs in the result
  nsnps <- Coloc_result$summary[1]

  # Make sure snp is formatted as expected
  if (FALSE %in% grepl(pattern=":",x=Coloc_result$results$snp)){
    stop("Expected coloc_result$results$snp to be formatted as such: 'chr#:####'")
  }
  # Split chr#:#### into a matrix
  split <- matrix(unlist(strsplit(x=as.character(Coloc_result$results$snp),split=":")),byrow = TRUE, ncol=2) #edit by zhec
  start <- as.integer(split[,2])
  end <- start+1
  snp_hyp4 <- Coloc_result$results$SNP.PP.H4
  eQTL_beta <- (Coloc_result$results$z.df1 * (sqrt(Coloc_result$results$V.df1)))
  eQTL_lABF <- Coloc_result$results$lABF.df1
  eQTL_credible <- as.integer(Coloc_result$credible_snps$credible_1$snps)
  eQTL_credible <- sapply(1:nsnps, function(snp_index, to_match) TRUE%in%grepl(start[snp_index],to_match),to_match=eQTL_credible)
  GWAS_beta <- (Coloc_result$results$z.df2 * (sqrt(Coloc_result$results$V.df2)))
  GWAS_lABF <- Coloc_result$results$lABF.df2
  GWAS_credible <- as.integer(Coloc_result$credible_snps$credible_2$snps)
  GWAS_credible <- sapply(1:nsnps, function(snp_index, to_match) TRUE%in%grepl(start[snp_index],to_match),to_match=GWAS_credible)
  both_credible <- as.integer(Coloc_result$credible_snps$in_both)
  both_credible <- sapply(1:nsnps, function(snp_index, to_match) TRUE%in%grepl(start[snp_index],to_match),to_match=both_credible)

  # Remove column names to avoid error message when making data frame below
  names(Coloc_result$summary) <- NULL

  bed_table <- data.frame(chr=Chromosome,
                          start=start,
                          end=end,
                          gene=gene,
                          trait=trait,
                          gene_trait_hyp3=Coloc_result$summary[5],
                          gene_trait_hyp4=Coloc_result$summary[6],
                          snp_hyp4=snp_hyp4,
                          eQTL_beta=eQTL_beta,
                          eQTL_lABF=eQTL_lABF,
                          eQTL_credible=eQTL_credible,
                          GWAS_beta=GWAS_beta,
                          GWAS_lABF=GWAS_lABF,
                          GWAS_credible=GWAS_credible,
                          both_credible=both_credible,
                          stringsAsFactors = FALSE)
  return(bed_table)
}
