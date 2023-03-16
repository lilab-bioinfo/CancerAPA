#/opt/apps/R/3.6.2/bin/R

# plot_coloc.R
# Caleb Matthew Radens
# 2016_1_5
# modified by YoSon

# plot_coloc contains functions for running coloc.abf() and visualizing
#   the input and output data and results.

# Given genes, traits, directories, run coloc.abf() and output plots

# I tested this script using the PMACS R module 3.1.2
system("echo ===========================",wait=FALSE)
system("echo inside plot_coloc.R",wait=FALSE)
R_ver <- substr(version$version.string,1,15) # Get R version
system(paste("echo",R_ver),wait=FALSE)
system("echo ===========================",wait=FALSE)

source('coloc_analysis.R')
require('hash')
require('scales')
require('data.table')
library(Cairo)    #added by zhec

system("echo wrapper package and script dependencies loaded and checked",wait=FALSE)

plot_coloc_results <- function(Coloc_result,Range=2e5,Cutoff_PP=0.3,X_axis_ticks=10){
  #   Given coloc.abf() outputs, plot the gene vs trait ABFs and hyp4 PPs
  #     - For each gene, a series of gene vs trait ABF and hyp4 PP plots will
  #         be generated and saved as a .png
  #     - This is a wrapper function that calls on get_plotting_data, get_genes,
  #         and plot_genes. See each of those functions for further details.
  #     - Note: this function will generate .pngs in the current working directory.
  #
  #   Arguments:
  #     Coloc_result: output from batch_run_coloc_abf
  # 	Range: how far from the top SNP should the plots go?
  #		Cutoff_PP: the Posterior Probability over which a SNP is defined as a "top SNP"
  # 	X_axis_ticks: the number of ticks to be shown on the x-axes of the plots
  if (length(results[[names(results)[1]]]) < 2){
    stop("Coloc_result needs to have summary *AND* results: set summary_only=FALSE in batch_run_coloc_abf()")
  }

  # Retreive 95% credibility sets
  Coloc_result <- batch_credibility_set(Coloc_result, Cutoff=0.95)

  plotting_data <- get_plotting_data(Coloc_result= Coloc_result,
                                     Range= Range,
                                     Cutoff_PP= Cutoff_PP,
                                     X_axis_ticks= X_axis_ticks)
  genes <- get_genes(Plotting_data= plotting_data)
  plot_genes(Plotting_data= plotting_data,
             Genes= genes)
}

get_plotting_data <- function(Coloc_result,Range=2e5,Cutoff_PP=0.3,X_axis_ticks=10){
  # Generate lists of data needed to plot coloc.abf() results.
  #
  # Arguments
  #   Coloc_result: output from batch_run_coloc_abf
  # 	Range: how far from the top SNP should the plots go?
  #		Cutoff_PP: the Posterior Probability over which a SNP is defined as a "top SNP"
  # 	X_axis_ticks: the number of ticks to be shown on the x-axes of the plots
  #
  # Returns a dictionary of lists with the following elements:
  #     Note: each dictionary key is a gene_trait
  # 	$"eQTL_lABF"        - Log(ABF) calculated in coloc.abf()
  #		$"GWAS_lABF"        - Log(ABF) calculated in coloc.abf()
  #		$"hyp4"             - Hypothesis 4 posterior probability (PP) of shared variant
  #   $"chr"              - Chrom qosome of region
  #		$"loci"             - Loci info for the SNPs
  #   $"TSS"              - All TSSs from V19 genome build matching gene name
  #   $"TES"              - All TESs from V19 genome build matching gene name
  #		$"top_snps"         - Top SNPs (with individual hyp4 PP > cutoff PP)
  #   $"eQTL_credible_set"- 95% credibility set of eQTL SNPs
  #   $"GWAS_credible_set"- 95% credibility set of GWAS SNPs
  #   $"credible_both"    - SNPs in both the eQTL and GWAS credible sets
  #		$"cutoff_PP"        - PP above which a SNP is considered a "top SNP"
  #		$"eQTL"             - Name of the gene
  #		$"trait"            - Name of the trait
  #		$"maj_axis_ticks"   - Loci of the major x-axis ticks
  #		$"mai_axis_ticks"   - Loci of the minor x-axis ticks
  #   $"left_legend_pos"  - x,y coordinates where left legend should be
  #   $"right_legend_pos" - x,y coordinates where right legend should be
  #		$"summary"          - Coloc.abf() summary (# of SNPs and hyp0-4 PPs)

  # Initialize dictionary of gene_traits
  all_lABFs<-hash()

  microArrayTable <- import_microarray_table()

  # Get the gene_traits
  gene_traits<-names(Coloc_result)
  gene_names <- matrix(unlist(strsplit(gene_traits,split="\\.")),ncol=2,byrow=TRUE)[,1]
  gene_TssTes_dict <- get_gene_TssTes(Genes<-gene_names,only_TssTes = FALSE)
  for (key in gene_traits){
    data <- Coloc_result[[key]]
    positions <- as.integer(c(data$results[,"snp"]))
    bottom <- round(min(positions),0)
    top <- round(max(positions),0)
    middle <- round(mean(c(top,bottom)),0)
    top_snps <- coloc_abf_top_snps(Coloc_result[[key]],cutoff_PP = Cutoff_PP)
    gene_name <- matrix(unlist(strsplit(key,split="\\.")),ncol=2,byrow=TRUE)[,1]
    TssTes<-gene_TssTes_dict[[gene_name]]
    if(length(TssTes)>1){
      TSS<- unique(gene_TssTes_dict[[gene_name]]$TSS)
      TES<- unique(gene_TssTes_dict[[gene_name]]$TES)
      chr<- gene_TssTes_dict[[gene_name]]$chr
    } else{
      TSS<-"Not Found"
      TES<-"Not Found"      #edit be zhec: TSS->TES
      chr<-microArrayTable[[gene_name]]$chr
    }
    loci_range <- c(middle-Range, middle+Range)
    if(loci_range[1]<bottom || loci_range[2]>top){
      print(loci_range)
      print(bottom)
      print(top)
      print(paste("User-defined range is wider than ",gene_name,"'s available positions.",sep=""))
    }
    range_indeces <- which(positions >= loci_range[1] & positions <= loci_range[2])
    split_name <- matrix(unlist(strsplit(key,split="_")),nrow=1,ncol=2)
    ticks<- get_axis_labels(Range=Range, N_ticks=X_axis_ticks, Center=middle)
    maj_ticks<- ticks$"major"
    min_ticks<- ticks$"minor"
    all_lABFs[key] <- list("eQTL_lABF"=data$results[,"lABF.df1"][range_indeces],
                           "GWAS_lABF"=data$results[,"lABF.df2"][range_indeces],
                           "hyp4"=data$results[,"SNP.PP.H4"][range_indeces],
                           "chr"=chr,
                           "loci"=positions[range_indeces],
                           "TSS"=TSS,
                           "TES"=TES,
                           "top_snps"=as.character(top_snps$snps),
                           "eQTL_credible_set"=data$credible_snps$credible_1$snps,
                           "GWAS_credible_set"=data$credible_snps$credible_2$snps,
                           "credible_both"=data$credible_snps$in_both,
                           "cutoff_PP"=Cutoff_PP,
                           "eQTL"=split_name[1],
                           "trait"=split_name[2],
                           "maj_axis_ticks"=maj_ticks,
                           "min_axis_ticks"=min_ticks,
                           "left_legend_pos"=get_legend_pos(as.integer(positions[range_indeces]),"left"),
                           "right_legend_pos"=get_legend_pos(as.integer(positions[range_indeces]),"right"),
                           "summary"=data$summary)
  }
  all_gene_traits <- matrix(unlist(strsplit(gene_traits,split="_")),
                            nrow=length(gene_traits),ncol=2,byrow=TRUE)
  all_lABFs["genes"]<-unique(all_gene_traits[,1])
  return(all_lABFs)
}

get_axis_labels <- function(Range, N_ticks, Center){
  # Given a center point, a range from the center,
  #   and a number of ticks, return tick locations.
  N_ticks<- N_ticks-1
  spacing<-(Range*2)/N_ticks
  minimum <- Center-Range
  maximum <- Center+Range
  maj<-seq(from=minimum,to=maximum,by=spacing)
  min<-c((maj[1]-0.5*spacing),(maj+(0.5*spacing)))
  ticks<-list("major"=maj,"minor"=min)

  return(ticks)
}

get_legend_pos <- function (X_axis_data, Side="left", P=0.2){
  # Given the X-variable data and a side ('left' or 'right'),
  #   Return a postion that is P% away from the left or right edge.
  #
  # Arguments:
  #   X_axis_data: the data to be plotted on the X-axis
  #     typeof: integer
  #   Side: which side of the plot to calculate a position from
    left_most <- min(X_axis_data)
    right_most <- max(X_axis_data)
    range <- right_most-left_most
  if (Side == "left") {
    pos<-left_most+round(P*range)
  } else if (Side == "right"){
    pos<-right_most-round(P*range)
  } else {
    stop(paste(Side,"is not an acceptable argument. Only 'left' or 'right' accepted."))
  }
  return(pos)
}

get_genes <- function(Plotting_data){
  # Build dictionary of each gene pointing at all the coloc results
  #  run with that gene.
  #
  # Argument
  #   Plotting_data: output from batch_run_coloc_abf:
  #     Dictionary of coloc.abf() results with gene_traits as keys
  #
  # Returns a dictionary of character->character; gene->gene_traits
  #   where the all gene_traits with a matching gene name are
  #   assigned to the gene.
  genes <- hash()
  gene_traits <- sort(names(Plotting_data))
  gene_traits <- gene_traits[gene_traits != "genes"]
  for (gene in Plotting_data$"genes"){
    genes[gene] <- gene_traits[grep(gene,gene_traits)]
  }
  return(genes)
}

plot_genes <- function (Plotting_data, Genes){
  for (gene in names(Genes)){
    plot_name <- paste("Colocalization Analysis of ",
                       gene,
                       ".png",
                       sep="")
    n_gene_traits <- length(Genes[[gene]])
    CairoPNG(file=plot_name, width=10, height=7.5, units="in", res=200) #edit by zhec: png->CairoPNG
    plot_layout_coords<-c()
    bottom_label_ratio<-200-floor(185/(n_gene_traits*2))*n_gene_traits*2
    for(i in 1:(n_gene_traits*2)){
      plot_layout_coords<-c(plot_layout_coords,rep(i,floor(185/(n_gene_traits*2))))
    }
    plot_layout_coords<-c(plot_layout_coords,rep(((n_gene_traits*2)+1),bottom_label_ratio))
    layout(
      mat=matrix(plot_layout_coords, 200, 1, byrow = TRUE)
    )
    counter <- 1
    for (gene_trait in Genes[[gene]]){
      data <- Plotting_data[[gene_trait]]
      credible_loci <- data$loci[unique(grep(
                              paste(data$eQTL_credible_set,collapse="|"),
                              data$loci,
                              value=FALSE))]
      credible_loci <- paste(credible_loci,collapse="|")
      par(mar=c(0,3,0.5,3))
      plot(data$loci, data$eQTL_lABF,
           yaxt="n",
           # col=ifelse(grepl(credible_loci,data$loci), "deepskyblue4", "deepskyblue"),
           col="deepskyblue",
           pch=ifelse(grepl(credible_loci,data$loci), 0, 8),
           cex=ifelse(grepl(credible_loci,data$loci), 1.2, .3),
           ylab="",
           xaxt="n",
           xlab=""
      )
      axis(side=2,col.axis="deepskyblue4",cex.axis=1)
      # Plot all known TSS and TES info, if the info was found in the transcripts table
      tss<-data$TSS
      if(!"Not Found"%in%tss){
        abline(v=tss,lty=2,col="dark green")
      }
      tes<-data$TES
      if(!"Not Found"%in%tes){
        abline(v=tes,lty=3)
      }
      #To make it say gene name on the side:
#       mtext(paste(data$eQTL,"log(ABF)"),
#             side=2,line=-1.2,cex=0.9,
#             col = "blue")
      # To make it say "eQTL log(ABF)" on the side
      mtext("eQTL log(ABF)",
            side=2,line=-1.2,cex=0.65,
            col = "deepskyblue")

      credible_loci <- data$loci[unique(grep(
                                  paste(data$GWAS_credible_set,collapse="|"),
                                  data$loci,
                                  value=FALSE))]
      credible_loci <- paste(credible_loci,collapse="|")
      par(new=T)
      par(mar=c(0,3,0.5,3))
      plot(data$loci,data$GWAS_lABF,
           yaxt="n",
           # col=ifelse(grepl(credible_loci,data$loci), "firebrick1", "firebrick4"),
           col="firebrick4",
           pch=ifelse(grepl(credible_loci,data$loci), 1, 8),
           cex=ifelse(grepl(credible_loci,data$loci), 1.2, .3),
           ylab="",
           xaxt="n"
      )
      axis(side=4, col.axis="firebrick4",line=0,cex.axis=1)
      mtext(paste(data$trait,"log(ABF)"),
            side=4,line=-1.2,cex=.75,
            col = "firebrick4")
      legend(x=data$"right_legend_pos",
             y=max(data$GWAS_lABF),
             legend=c(paste("Hyp4 PP: ",round(100*data$summary[6],1),"%"),
                      paste("Hyp3 PP: ",round(100*data$summary[5],1),"%")),
             cex=1.3,
             text.col="grey60",
             text.font=2,
             bty="n")
      legend(x=data$"right_legend_pos",
             y=max(data$GWAS_lABF),
             legend=c(paste("Hyp4 PP: ",round(100*data$summary[6],1),"%"),
                      paste("Hyp3 PP: ",round(100*data$summary[5],1),"%")),
             cex=1.28,
             text.col="black",
             text.font=2,
             bty="n")
      if(counter==1){
        legend(x="topleft",
               legend=c(gene),
               cex=1.5,
               text.col="grey60",
               text.font=2,
               bty="n")
        legend(x="topleft",
               legend=c(gene),
               cex=1.48,
               text.col="black",
               text.font=2,
               bty="n")
      }
      par(new=F)
#       if (counter == n_gene_traits){
#         par(mar=c(.5,3,.5,3))
#       } else{
#         par(mar=c(1,3,.5,3))
#       }
      par(mar=c(1,3,.5,3))
      plot(data$loci,100*data$hyp4,pch=4,
           xaxt="n",
           yaxt="n",
           xlab="",
           ylab="",
           ylim=c(0,100))
      percentages<-c(0,50,100)
      axis(side=2,cex.axis=1,
           at=percentages,
           labels= percentages)
      mtext("Hyp4 PP",side=2,line=-1.2,cex=0.65)
#       legend(x=data$"right_legend_pos",
#              y=max(data$hyp4),
#              legend=data$top_snps,
#              title=paste("SNP(s) with PP >",data$"cutoff_PP"))
      if (counter == n_gene_traits){
        par(mar=c(0,0,0,0))
        # axis(side=1,at=X,labels=X)
        # mtext("hello",side=1,line=2)
        mtext(paste(data$chr, "position"),side=1,line=2,cex=1.5)
        axis(side=1,at=data$"min_axis_ticks",tick=TRUE,labels=FALSE,line=-1)
        axis(side=1,at=data$"maj_axis_ticks",tick=TRUE,labels=TRUE,line=-1)
      }

      counter = counter + 1
    }
    dev.off()
  }
}
