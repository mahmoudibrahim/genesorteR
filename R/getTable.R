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
#' numbers are allowed. Set to FALSE to switch off filtering on the average fold-change
#' value.
#' @param adjpval_cutoff A numeric adjusted p-value cutoff value. Default is 0.05.
#' @param islog A logical value. TRUE (default) means the expression matrix supplied 
#' previously to \code{sortGenes} is in log space.
#' @param pseudocount A numeric value. A pseudocount to add to the expression matrix before
#' taking the log if \code{islog} is FALSE. 
#'
#' @return \code{getTable} returns a data frame containing the gene name, the average log fold 
#' change, the adjusted p-value and the specificity score for all variable genes in each cluster.
#' The data frame is sorted by the specificity score. Note that a gene may appear multiple times 
#' for multiple clusters.
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
#' unique(tab$Gene.Name)
#'
#' #To get all genes without any cutoffs, set adjpval_cutoff to >1 and fc_cutoff to FALSE
#' tab = getTable(gs, pp, fc_cutoff = FALSE, adjpval_cutoff = 1.1)
getTable = function(gs, pp, fc_cutoff = 0, adjpval_cutoff = 0.05, islog = TRUE, pseudocount = 1) {

	#pval
	tooearly = apply(pp$adjpval, 2, function(x) length(x[x < adjpval_cutoff]))
	toolate = which(tooearly == 0)

	if (length(toolate) > 0) {
		cl = getClassIndeces(gs$inputClass[-(which(gs$inputClass %in% names(tooearly[toolate])))])
		nomen = names(gs$classProb[-toolate])
		sp = gs$specScore[,-toolate]
		
		tooearly = names(which(apply(pp$adjpval[,-toolate], 1, function(x) any(x < adjpval_cutoff))))
	} else {
		cl = getClassIndeces(gs$inputClass)
		nomen = names(gs$classProb)
		sp = gs$specScore
		
		tooearly = names(which(apply(pp$adjpval, 1, function(x) any(x < adjpval_cutoff))))
	}


	if (length(tooearly) > 0) {
		pv = pp$adjpval[which(rownames(pp$adjpval) %in% tooearly),]
	} else {
		stop("No genes passing the p-value cutoff were found!")
	}
	
	

	#fc
	toolate = which(rownames(gs$inputMat) %in% tooearly)
	if (!islog) {
		mat = log(gs$inputMat[toolate,] + pseudocount)
	} else {
		mat = gs$inputMat[toolate,]
	}
	
	fc = do.call(cbind, lapply(cl, function(clx) (Matrix::rowMeans(mat[,clx])) - (Matrix::rowMeans(mat[,-clx]))))
	colnames(fc) = nomen
	if (is.numeric(fc_cutoff)) {
		tooearly = names(which(apply(fc, 1, function(x) any(x > fc_cutoff))))
		toolate = which(rownames(fc) %in% tooearly)
		if (length(toolate) > 0) {
			fc = fc[toolate,]
		} else {
			stop("No genes passing the p-value and fold-change cutoffs were found!")
		}
	}

	#table
	pv = pv[,which(colnames(pv) %in% colnames(fc))]
	nm = intersect(rownames(fc), rownames(pv))
	fc = fc[which(rownames(fc) %in% nm),]
	pv = pv[which(rownames(pv) %in% nm),]
	sp = sp[which(rownames(sp) %in% nm),]
	fc = fc[order(rownames(fc)),]
	pv = pv[order(rownames(pv)),]
	sp = sp[order(rownames(sp)),]
	
	tab = data.frame(rep(rownames(fc), ncol(fc)), c(fc), c(pv), c(as.matrix(sp)), rep(colnames(fc), each = nrow(fc)))
	colnames(tab) = c("Gene.Name", "Average.Log.Fold.Change","Adjusted.pvalue","Specificity.Score","Cluster")
	
	if (is.numeric(fc_cutoff)) {
		toolate = which((tab[,2] > fc_cutoff) & (tab[,3] < adjpval_cutoff))
		tab = tab[toolate,]
	} else {
		toolate = which(tab[,3] < adjpval_cutoff)
		tab = tab[toolate,]
	}
	tab = tab[order(tab$Specificity.Score, decreasing = TRUE),]
	
	
	return(tab)
}
