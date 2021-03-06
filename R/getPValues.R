#' getPValues
#'
#' getPValues performs a permutation test on the gene-cell type specificity
#' score to obtain a p-value value on each gene-cell type specificity score. The
#' permutation keeps everything the same except that cell type assignments are
#' permuted between the cells (but note that cell type proportions are also kept
#' the same). This function is a way to select all differentially expressed
#' genes between all cell classes globally in a dataset in one go or to define 
#' differentially expressed genes in a specific cell cluster.
#'
#' @param gs The output of \code{sortGenes()}.
#' @param numPerm The number of permutations to do. The default value 5 is to
#'   ensure quick running time for very large datasets but can be increased to
#'   increase the power of the permutation test, see Details.
#' @param correctMethod The method used to correct p-values for multiple
#'   hypothesis testing. Any valid input to "method" in \code{p.adjust} is
#'   allowed. p-value correction is done on a gene-by-gene basis.
#' @param testGenes A character vector of gene names to restrict the output
#'   p-values to.
#' @param subsetCells A numeric vector of cell indeces to restrict the
#'   permutation to. Note that the selected cells should still contain at least
#'   one cell from each of the cell clusters contained in \code{$specScore}.
#'   Note this option is still experimental and might not give consistent
#'   results.
#' @param cores An integer greater than zero (1 by default) that indicates how
#'   many cores to use for parallelization using mclapply.
#' @param seed The seed for random permutations.
#'
#' @return \code{getPValues} returns a list with the following components:
#'   \item{permuteVal}{A sparse matrix with as many rows as genes and as many
#'   columns as \code{ncol($specScore) * numPerm}, containing the specificity
#'   score matrices for all permutations concatenated after each other
#'   column-wise.} \item{startIndeces}{A numeric vector of length \code{numPerm}
#'   that indicates the starting column for each performed permutation in
#'   \code{permuteVal}} \item{pval}{A matrix with as many rows as genes and
#'   columns as \code{ncol($specScore)}, containing the p-values for the null
#'   hypothesis that the gene is not highly specific to a cell cluster is true.}
#'   \item{adjpval}{A matrix with the same size as \code{pval}, containing the
#'   corrected p-values using the method specified in \code{correctMethod}}
#' @export
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @examples
#' data(kidneyTabulaMuris)
#' gs = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' pp = getPValues(gs)
#' #obtain genes that are "differentially expressed" in at least one cluster
#' markers = names(which(apply(pp$adjpval, 1, function(x) any(x < 0.01))))
#'
#' @seealso
#' getMarkers
getPValues = function(gs, numPerm = 5, correctMethod = "BH", testGenes = NULL, subsetCells = NULL, cores = 1, seed = 111) {

	#recompute without non-binary genes
	takeout = which(Matrix::rowSums(gs$binary) == ncol(gs$binary))
	if (length(takeout) > 0) {
		gs = sortGenes(gs$binary[-takeout,], gs$inputClass, binarizeMethod = "naive", cores = cores)
	}

	ngroup = ncol(gs$specScore)
	ngene = nrow(gs$specScore)
	mm = Matrix::sparseMatrix(1, 1, dims = c(ngene, (ngroup * numPerm)))
	mm = mm * 1
	dd1L = rep(0, length = numPerm)

	#permutation
	set.seed(seed)
	for (i in 1:numPerm) {
		dd1 = (ngroup * i) - (ngroup - 1)
		dd2 = (ngroup * i)
		mm[,dd1:dd2] = getPerma(gs, subsetCells = subsetCells, cores = cores)
		dd1L[i] = dd1
	}
	rownames(mm) = rownames(gs$specScore)
	ee = ecdf(as.vector(mm@x))
	nn = length(as.vector(mm@x))

	#getting pvalues
	if (is.null(testGenes)) {
		pval = matrix(ee(as.matrix(gs$specScore)), ncol = ncol(gs$specScore), byrow = FALSE) * nn
		pval = nn - pval
		pval = (pval + 1) / (nn + 1)
	
		if (cores == 1) {
			padj = do.call(rbind, lapply(1:ngene, function(k) p.adjust(pval[k,], method = correctMethod)))
		} else {
			padj = do.call(rbind, mclapply(1:ngene, function(k) p.adjust(pval[k,], method = correctMethod), mc.cores = cores))
		}
		rownames(padj) = rownames(gs$specScore)
		rownames(pval) = rownames(gs$specScore)

	} else {
		wgene = which(rownames(gs$specScore) %in% testGenes)
		ngeneT = length(wgene)
		if (ngeneT != length(testGenes)) {
			warning("Some of the genes you specified with 'testGenes' were not found. Were they filtered out by 'sortGenes()?'")
		}
		namesT = rownames(gs$specScore)[wgene]

		pval = matrix(ee(as.matrix(gs$specScore[wgene,])), ncol = ncol(gs$specScore), byrow = FALSE) * nn
		pval = nn - pval
		pval = (pval + 1) / (nn + 1)


		if (cores == 1) {
			padj = do.call(rbind, lapply(1:ngeneT, function(k) p.adjust(pval[k,], method = correctMethod)))
		} else {
		padj = do.call(rbind, mclapply(1:ngeneT, function(k) p.adjust(pval[k,], method = correctMethod), mc.cores = cores))
		}
		rownames(padj) = namesT
		rownames(pval) = namesT
	}

	colnames(padj) = colnames(gs$specScore)
	colnames(pval) = colnames(gs$specScore)
	colnames(mm) = rep(colnames(gs$specScore), numPerm)
	return(list(permuteVal = mm, startIndeces = dd1L, pval = pval, adjpval = padj))

}
