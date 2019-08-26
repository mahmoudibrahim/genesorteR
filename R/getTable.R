#' getTable
#'
#' Summarize \code{sortGenes} and \code{getPValues} results in one data frame.
#'
#' \code{getTable} combines the results of sortGenes and getPValues in one table
#' and also calculates the average fold-change of expression, allowing to select 
#' a fold-change cutoff and a p-value cutoff to determine variable genes.
#'
#'
#' @param gs The output of \code{sortGenes}.
#' @param pp The output of \code{getPValues}.
#' @param fc_cutoff Default is 0, which means only genes that have an average 
#' fold-change higher than 0 for a given cluster are reported. Positive and negative 
#' numbers are allowed. Se to FALSE to switch off filtering on the average fold-change
#' value.
#' @param adjpval_cutoff A numeric adjusted p-value cutoff value. Default is 0.05.
#' @param islog A logical value. TRUE (default) means the expression matrix supplied 
#' previously to \code{sortGenes} is in log space.
#' @param pseudocount A numeric value. A pseudocount to add to the expression matrix before
#' taking the log if \code{islog} is FALSE. 
#'
#' @return \code{getTable} returns a data frame containing the average log fold change, the adjusted 
#' p-value and the specificity score for all variable genes in each cluster. The data frame is sorted 
#' by the specificity score. Note that a gene may appear multiple times for multiple clusters.
#' 
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @export
#' @examples
#' data(kidneyTabulaMuris)
#' gs = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' pp = getPValues(gs)
#' tab = getTable(gs, pp)
#' 
#' #A quick diagnostic plot (fold change correlates with specificity score)
#' plot(tab$Average.Log.Fold.Change, tab$Specificity.Score, col = as.factor(tab$Cluster), pch = 20)
#'
#' #all variable genes
#' unique(rownames(tab))
#'
#' #To get all genes without any cutoffs, set adjpval_cutoff to 1 and fc_cutoff to FALSE
#' tab = getTable(gs, pp, fc_cutoff = FALSE, adjpval_cutoff = 1)
getTable = function(gs, pp, fc_cutoff = 0, adjpval_cutoff = 0.05, islog = TRUE, pseudocount = 1) {

	cl = getClassIndeces(gs$inputClass)
	
	if (islog) {
		mat = gs$inputMat
	} else {
		mat = log(gs$inputMat + pseudocount)
	}
	
	fc = do.call(cbind, lapply(cl, function(clx) Matrix::rowMeans(mat[,clx]) - Matrix::rowMeans(mat[,-clx])))


	nm = intersect(rownames(fc), rownames(pp$adjpval))
	fc = fc[which(rownames(fc) %in% nm),]
	pp$adjpval = pp$adjpval[which(rownames(pp$adjpval) %in% nm),]

	if (isFALSE(fc_cutoff)) {
		machine = lapply(1:ncol(fc), function(i) names(pp$adjpval[which(pp$adjpval[,i] < adjpval_cutoff),i]))
	} else {
		machine = lapply(1:ncol(fc), function(i) intersect(names(pp$adjpval[which(pp$adjpval[,i] < adjpval_cutoff),i]), names(fc[which(fc[,i] > fc_cutoff),i])))
	}
	
	tab = do.call(rbind, lapply(1:ncol(fc), function(i) data.frame(fc[which(rownames(fc) %in% machine[[i]]),i], pp$adjpval[which(rownames(pp$adjpval) %in% machine[[i]]),i], gs$specScore[which(rownames(gs$specScore) %in% machine[[i]]),i], colnames(gs$specScore)[i])))
	colnames(tab) = c("Average.Log.Fold.Change","Adjusted.pvalue","Specificity.Score","Cluster")
	
	tab = tab[order(tab$Specificity.Score, decreasing = TRUE),]
	
	
	return(tab)
}
