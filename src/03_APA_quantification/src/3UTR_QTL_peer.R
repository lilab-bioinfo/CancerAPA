#!/usr/bin/env Rscript

# 2017-12-09

# print usage
usage <- function() {
  cat(
    'usage: peer.R <tissue>
peer.R

author: Lei Li
  tissue             tissue name
')
}

# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
TISSUE <- args[1]
print(TISSUE)
# Check input args
if (is.na(TISSUE)) {
  usage()
  quit(save='no', status=1)
}

library(gplots)
library(peer)

# simulate the colors of Matlab
myCols=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

# set your display range here
pairs.breaks <- seq(-3, 3, by=.05)
## pairs.breaks
## length(pairs.breaks)

OUT.PREFIX <- 'peer_pdui'

# ----------------------------------

PATH.DIR <- '/dfs5/weil21-lab2/leil22/GTEx/'
EXPR.DIR <- '/dfs5/weil21-lab2/leil22/GTEx/PDUIs'

mat <- read.table(paste0(EXPR.DIR, '/', TISSUE, '_combined_All_PDUIs_clean.txt'), stringsAsFactors=FALSE, header=FALSE)
mat <- mat[,-2]

colnames(mat) <- mat[1,]
rownames(mat) <- mat[,1]
mat <- as.matrix(mat[-1,-1])

mat <- mat[, colMeans(is.na(mat)) <= 0.8]
mat <-  mat[rowMeans(is.na(mat)) < 0.5,]

class(mat) <- 'numeric'

print("Loading finished")


## prior to PEER correction
#png(paste0(PATH.DIR, '/', OUT.PREFIX, '/', TISSUE, '.expr.clust.png'), width=8, height=8, res=150, units='in')
#heatmap.2(as.matrix(mat), distfun=function(x) dist(x,method='euclidian'), hclustfun=function(x) hclust(x,method='ward.D2'),
#                    trace='none', dendrogram='both', Rowv=TRUE, Colv=TRUE, breaks=pairs.breaks, col=colorRampPalette(myCols), scale='none', symkey=T, na.color='grey', density.info='histogram', cexRow=0.2, cexCol=0.5, main=paste0(TISSUE, '\nexpr clustering'))
#dev.off()

## ----------------------------------------
# peer correct

model <- PEER()

# set covariates
COV.DIR <- '/dfs5/weil21-lab2/leil22/GTEx/Covariates_by_tissue/'

covs = read.table(paste0(COV.DIR, '/', TISSUE,'..v8.covariates.v2.txt'), header=TRUE, check.names = F,sep="\t")
rownames(covs) <- covs[,1]
covs <- covs[,-1]
covs <- covs[,colnames(covs) %in% colnames(mat),]
covs_se <- t(covs)
colnames(covs_se) <- NULL

PEER_setCovariates(model, as.matrix(covs_se))

dim(PEER_getCovariates(model))
#impute missing values in 3'UTR expression
library(impute)
mat.ds <- mat[,colnames(mat) %in% rownames(covs_se)]
#t(as.matrix(mat)) <-  impute.knn(t(as.matrix(mat)))
mat_impute <- impute.knn(mat.ds)
PEER_setPhenoMean(model, t(as.matrix(mat_impute$data)))

dim(PEER_getPhenoMean(model))

# set number of peer factors
if (ncol(mat.ds) < 150) { 
	numcov <- 15
} else if (ncol(mat.ds) < 250) {
	numcov <- 30
} else if (ncol(mat.ds) >= 250) {
	numcov <- 35 
}

PEER_setNk(model, numcov)
PEER_getNk(model)

PEER_update(model)

# diag
#pdf(paste0(PATH.DIR, '/',OUT.PREFIX, '/', TISSUE, '.peer.diag.pdf'), width=6, height=8)
#PEER_plotModel(model)
#dev.off()
## cor(PEER_getCovariates(model)[,1], PEER_getX(model)[,2])

factors = t(PEER_getX(model))
weights = PEER_getW(model)
precision = PEER_getAlpha(model)

residuals = t(PEER_getResiduals(model))
rownames(residuals) <- rownames(mat.ds)
colnames(residuals) <- colnames(mat.ds)

rownames(factors) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "pcr", "platform", "sex",paste0("PEER_top_factor_",1:numcov))
colnames(factors) <- colnames(mat.ds)

#residuals.ds <- residuals[seq(1, nrow(residuals), ceiling(nrow(residuals)/num.ds)),]
residuals.ds <- residuals

#png(paste0(PATH.DIR, '/',OUT.PREFIX, '/', TISSUE, '.expr.peer.clust.png'), width=8, height=8, res=150, units='in')
#heatmap.2(as.matrix(residuals.ds), distfun=function(x) dist(x,method='euclidian'), hclustfun=function(x) hclust(x,method='ward.D2'),
#          trace='none', dendrogram='both', Rowv=TRUE, Colv=TRUE, breaks=pairs.breaks, col=colorRampPalette(myCols), scale='none', symkey=T, na.color='grey', density.info='histogram', cexRow=0.2, cexCol=0.5, main=paste0(TISSUE, '\nexpr clustering'))
#dev.off()

covariate_file <- paste0(PATH.DIR, '/', OUT.PREFIX, '/', TISSUE, ".pdui.peer.covariates.txt")
write.table(cbind(rownames(factors), factors), file=covariate_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

gz1 <- paste0(PATH.DIR, '/', OUT.PREFIX, '/', TISSUE, ".expr.peer.txt")
write.table(cbind(rownames(residuals), residuals), file=gz1, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
