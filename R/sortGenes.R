#' sortGenes
#'
#' sortGenes is the main function of the genesorteR package. It takes a gene
#' expression matrix and cell cluster assignments. It binarizes the expression
#' matrix and calculates empirical statistics on gene expression in each
#' cluster including a specificity score that can be used to rank genes in cell
#' clusters.
#'
#' When \code{binarizeMethod} is "median", expression matrix binarization is
#' done by estimating a cutoff value uniformly on all values in the matrix.
#' This is equal to the median of all non-zero entries in the matrix and is
#' returned in \code{cutoff}. When binarizeMethod is "naive", all non-zero
#' entries are kept and the minimum value of non-zero entries is returned in
#' \code{cutoff}. When the input matrix \code{x} has already been binarized, 
#' set \code{binarizeMethod} to "naive". You can set a specific cutoff value 
#' for binarization, by setting \code{binarizeMethod} to a numeric value >= 0.
#'
#' The specificity scores balance the posterior probability of observing a cell
#'  cluster given the gene (gene-cluster specificity) with its conditional 
#' probability given the cluster (a measure of gene expression). This ensures 
#' that highly specific genes are also highly expressed. The \code{specScore} 
#' matrix is considered the main output of this function, and on which many of the
#' remaining calculations by other functions in the \code{genesorteR} package 
#' are performed. The values in this matrix can be used to rank features 
#' (genes in scRNA-Seq) in clusters.
#' 
#' Note that if returnInput is set to FALSE (input expression matrix will no be
#' returned in the output), many of the other functions that accept the output of
#' \code{sortGenes} will break. 
#'
#' \code{sortGenes} can in principle be applied to both a raw count matrix or a
#' normalized log-count expression matrix.
#' @param x A numeric (sparse) matrix. It will be coerced to a dgCMatrix sparse
#'   matrix. Rows represent genes, columns represent cells.
#' @param classLabels A numeric or character vector or a factor of the same
#'   length as ncol(x) that represents cell cluster assignments. It will be
#'   coerced to a factor whose levels are the cell cluster names.
#' @param binarizeMethod Either "median" (default) or "naive" or a numeric cutoff. 
#' See Details.
#' @param returnInput Return the input matrix and cell classes? \code{TRUE} by 
#' default. See Details.
#' @param cores An integer greater than zero (1 by default) that indicates how
#'   many cores to use for parallelization using mclapply.
#' @return \code{sortGenes} returns a list with the following components:
#'   \item{binary}{The binarized gene expression matrix. A sparse matrix of
#'   class dgCMatrix.} \item{cutoff}{The cutoff value
#'   used to binarize the gene expression matrix. Anything lower than this value
#'   was set to zero.} \item{removed}{A numeric vector containing the row
#'   indeces of genes that were removed because they were not expressed in any
#'   cells after binarization. If none were removed, this will be \code{NULL}.} 
#'   \item{geneProb}{A numeric vector whose length is equal to nrow(binary) that 
#'   lists the fraction of cells in which a gene was detected.} \item{condGeneProb}{A
#'   sparse matrix of class dgCMatrix with as many rows as nrow(binary) and as
#'   many columns as the number of cell clusters. It includes the conditional
#'   probability of observing a gene in a cluster.} \item{postClustProb}{A
#'   sparse matrix of class dgCMatrix with the same size as \code{condGeneProb},
#'   containing the posterior probability that a cell belongs to a certain
#'   cluster or type given that the gene was observed.} \item{specScore}{A
#'   sparse matrix of class dgCMatrix with the same size as \code{condGeneProb},
#'   containing a specificity score for each gene in each cell cluster. See
#'   Details.} \item{classProb}{A numeric vector whose length is equal to the
#'   number of cell clusters, containing the fraction of cells belonging to each
#'   cluster.} \item{inputMat}{the input \code{x} matrix, after being coerced
#'   to a sparse matrix of class dgCMatrix.} \item{inputClass}{the input
#'   \code{classLabels} after being coerced to a factor.}
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @export
#' @examples
#' data(kidneyTabulaMuris)
#' #basic functionality
#' gs = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#'
#' #the top 10 genes for each cluster by specificity scores
#' top_genes = apply(gs$specScore, 2, function(x) names(head(sort(x, decreasing = TRUE), n = 10)))
#' 
#' #the same top 10 genes but using the plotTopMarkerHeat function
#' plotTopMarkerHeat(gs, top_n = 10, outs = TRUE, plotheat = FALSE)
#'
#' #naive binarization keeps any non-zero input in the input matrix
#' gs_naive = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType, binarizeMethod = "naive")
#'
#' #different genes?
#' plotTopMarkerHeat(gs_naive, top_n = 10, outs=TRUE, plotheat=FALSE)
sortGenes = function(x, classLabels, binarizeMethod = "median", returnInput = TRUE, cores = 1) {
	
	classLabels = as.factor(classLabels)
	if (length(setdiff(levels(classLabels), classLabels)) > 0) {
		classLabels = classLabels[,drop = TRUE] #drop nonoccuring levels
		warning("A Friendly Warning: The cell type factor had some empty levels. You probably don't need to do anything. Cell types were likely filtered beforehand.")
	}

	if (length(levels(classLabels)) < 2) {
		stop("Sorry but that's an error. classLabels specifies only one cell cluster!")
	}	

	ww = which(as.vector(table(classLabels)) == 1)	
	if (length(ww) > 0) {
		stop("Sorry but that's an error. sortGenes() won't continue. There were some cell types comprised of only one cell. Although it's technically possible to have this, something is likely off with your cell clustering.")	
	}
	
	classLabelsNum = as.integer(classLabels)
	classLabelsLab = levels(classLabels)

	x = as(x, "dgCMatrix")

	xbin = binarize(x, method = binarizeMethod)

	rem = which((Matrix::rowSums(xbin$mat)) == 0)
	if (length(rem) > 0) {
		xbin$mat = xbin$mat[-rem,]
		warning("A Friendly Warning: Some genes were removed because they were zeros in all cells after binarization. You probably don't need to do anything but you might want to look into this. Maybe you forgot to pre-filter the genes? You can also use a different binarization method. Excluded genes are available in the output under '$removed'.")
	}

	classProb = getClassProb(classLabelsNum)
	names(classProb) = classLabelsLab

	geneProb = getGeneProb(xbin$mat)

	condGeneCluster = getGeneConditionalCluster_stats(xbin$mat, classProb, classLabels, cores = cores)
	clusterPostGene = getClusterPostGene(condGeneCluster, geneProb, classProb)
	specScore = getSpecScore(clusterPostGene, condGeneCluster)

	if (returnInput) {
		return(list(binary = xbin$mat, cutoff = xbin$cutoff, removed = rem, geneProb = geneProb, condGeneProb = condGeneCluster, postClustProb = clusterPostGene, specScore = specScore, classProb = classProb, inputMat = x, inputClass = classLabels))
	} else {
		return(list(binary = xbin$mat, cutoff = xbin$cutoff, removed = rem, geneProb = geneProb, condGeneProb = condGeneCluster, postClustProb = clusterPostGene, specScore = specScore, classProb = classProb))
	}
}
