#' getTable
#'
#' Summarize \code{sortGenes} and \code{getPValues} results in one table.
#'
#' \code{getTable} combines the results of sortGenes and getPValues in one table
#' and also calculates the average fold-change of expression, allowing to select 
#' a fold-change cutoff and a p-value cutoffto determine differential expression.
#'
#'
#' @param sg The output of \code{sortGenes}.
#' @param mm The output of \code{getPValues}.
#' @param fc_cutoff Default is FALSE which means no fold change cut off is used to 
#' filter the results. If it is a numeric value, it means only genes that have a 
#' an average fold-change higher than this number for a given cluster are reported.
#'
#' @param adjpval_cutoff A numeric adjusted p-value cutoff value. Default is 0.05.
#' @param log A logical value. TRUE (default) means the expression matrix supplied 
#' previously to \code{sortGenes} is in log space.
#' @param pseudocount A numeric value. A pseudocount to add to the expression matrix before
#' taking the log if \code{log} is FALSE. 
#'
#' @return \code{getTable} returns a data frame containing the average log fold change, the adjusted 
#' p-value and the specificity score for all variable genes in all clusters. 
#' 
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @export
#' @examples
#' data(kidneyTabulaMuris)
#' gs = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' pp = getPValues(gs)
#' tab = getTable(gs, pp)
getTable = function(gs, pp, fc_cutoff = FALSE, adjpval_cutoff = 0.05) {

	cl = getClassIndeces(gs$inputClass)
	
	if (log) {
		mat = gs$inputMat
	} else {
		mat = log(gs$inputMat + pseudocount)
	}
	
	fc = do.call(cbind, lapply(cl, function(clx) Matrix::rowMeans(mat[,clx]) - Matrix::rowMeans(mat[,-clx])))

	if (!fc_cutoff) {
		machine = lapply(1:ncol(fc), function(i) names(pp$adjpval[which(pp$adjpval[,i] < adjpval_cutoff),i]))
	} else {
		machine = lapply(1:ncol(fc), function(i) intersect(names(pp$adjpval[which(pp$adjpval[,i] < adjpval_cutoff),i]), names(fc[which(fc[,i] < fc_cutoff),i])))
	}
	
	tab = do.call(rbind, lapply(1:ncol(fc), function(i) data.frame(fc[which(rownames(fc) %in% machine[[i]]),i], pp$adjpval[which(rownames(pp$adjpval) %in% machine[[i]]),i], gs$specScore[which(rownames(gs$specScore) %in% machine[[i]]),i], colnames(gs$specScore)[i])))
	
	colnames(tab) = c("Average Log Fold Change","Adjusted p-value","Specificity Score","Cluster")
	
	return(tab)
}
