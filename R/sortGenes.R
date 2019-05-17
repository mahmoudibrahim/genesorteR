#' sortGenes
#'
#' sortGenes is the main function of the genesorteR package. It takes a gene
#' expression matrix and cell cluster assignments. It binarizes the expression
#' matrix and calculates empirical statistics on gene expression in each
#' cluster.
#'
#' When \code{binarizeMethod} is "median", expression matrix binarization is
#' done by estimating a cutoff value uniformally on all values in the matrix.
#' This is equal to the median of all non-zero entries in the matrix and is
#' returned in \code{cutoff}. When binarizeMethod is "naive", all non-zero
#' entries are kept and the minimum value of non-zero entries is returned in
#' \code{cutoff}. When the input matrix \code{x} has already been binarized, set
#' \code{binarizeMethod} to "naive".
#'
#' The specificity score matrix \code{specScore} is calculated as
#' \code{condGeneProb * postClustProb}. Therefore, it simply balances the
#' posterior probability of observing a cell cluster given the gene
#' (gene-cluster specificity) with its conditional probability given the cluster
#' (a measure of gene expression). This ensures that highly specific genes are
#' also highly expressed. The \code{specScore} matrix is considered the main
#' output of this function, and on which many of the remaining calculations by
#' other functions in the \code{genesorteR} package are performed. The values in
#' this matrix can be used to rank features (genes in scRNA-Seq) in clusters.
#'
#' \code{sortGenes} can in principle be applied to both a raw count matrix or a
#' normalized log-count expression matrix.
#' @param x A numeric (sparse) matrix. It will be coerced to a dgCMatrix sparse
#'   matrix. Rows represent genes, columns represent cells.
#' @param classLabels A numeric or character vector or a factor of the same
#'   length as ncol(x) that represents cell cluster assignments. It will be
#'   coerced to a factor whose levels are the cell cluster names.
#' @param binarizeMethod Either "median" (default) or "naive". See Details.
#' @param cores A number greater than zero (1 by default) that indicates how
#'   many cores to use for parallelization using mclapply.
#' @return \code{sortGenes} returns a list with the following components:
#'   \item{binary}{The binarized gene expression matrix. A sparse matrix of
#'   class dgCMatrix of the same size as "x".} \item{cutoff}{The cutoff value
#'   used to binarize the gene expression matrix. Anything lower than this value
#'   was set to zero.} \item{removed}{A numeric vector containing the row
#'   indeces of genes that were removed because they were not expressed in any
#'   cells. If none were removed, this will be \code{NULL}.} \item{geneProb}{A
#'   numeric vector whose length is equal to nrow(binary) that lists the
#'   fraction of cells in which a gene was detected.} \item{condGeneProb}{A
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
#' #basic functionality
#' set.seed(1234)
#' exp = matrix(sample(0:20,1000,replace=TRUE), ncol = 20)
#' rownames(exp) = sapply(1:50, function(x) paste0("g", x))
#' classLab = sample(c("cell type 1","cell type 2"),20,replace=TRUE)
#' sg = sortGenes(exp, classLab)
#'
#' #naive binarization keeps any non-zero input in the input matrix
#' sg_naive = sortGenes(exp, classLab, binarizeMethod = "naive")
sortGenes = function(x, classLabels, binarizeMethod = "median", cores = 1) {
	
	classLabels = as.factor(classLabels)
	classLabels = classLabels[,drop = TRUE] #drop nonoccuring levels (issue a warning?)
	
	classLabelsNum = as.integer(classLabels)
	classLabelsLab = levels(classLabels)

	x = as(x, "dgCMatrix")

	xbin = binarize(x, method = binarizeMethod)

	rem = which((Matrix::rowSums(xbin$mat)) == 0)
	if (length(rem) > 0) {
		xbin$mat = xbin$mat[-rem,]
		##issue a warning here##
	}

	classProb = getClassProb(classLabelsNum)
	names(classProb) = classLabelsLab

	geneProb = getGeneProb(xbin$mat)

	condGeneCluster = getGeneConditionalCluster_stats(xbin$mat, classProb, classLabels, cores = cores)
	clusterPostGene = getClusterPostGene(condGeneCluster, geneProb, classProb)
	specScore = getSpecScore(clusterPostGene, condGeneCluster)

	return(list(binary = xbin$mat, cutoff = xbin$cutoff, removed = rem, geneProb = geneProb, condGeneProb = condGeneCluster, postClustProb = clusterPostGene, specScore = specScore, classProb = classProb, inputMat = x, inputClass = classLabels))

}
