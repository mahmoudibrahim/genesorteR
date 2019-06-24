#' getClassAUC
#'
#' getClassAUC implements one way to investigate clustering quality. It processes the
#' output of \code{sortGenes} to obtain a curve for each cell cluster for all
#' gene specificity scores against their ranking in the cluster. The Area Under
#' the Curve (AUC) can be used as a measure of clustering quality in terms of the
#' possibility to identify cell clusters using a few marker genes. See Details.
#'
#' Given the specificity score for all genes in a certain cell cluster, we can
#' assume that a well-separated easily-identified cell cluster will have a
#' relatively small number of genes that have a very high specificity score. Top
#' marker genes for a  cluster that is poorly separated from other cell
#' clusters will have average or low specificity scores. Sorting the genes for
#' each cell cluster by their specificity scores and plotting the scaled scores
#' in order creates a curve that should be far from the diagonal for
#' well-separated clusters but close to the diagonal for poorly-separated
#' clusters. The AUC of this curve can be used to quantify this intuition and
#' estimate a clustering quality metric.
#'
#'
#' @param gs A list containing \code{$specScore} sparse matrix. Typically the
#'   output of \code{sortGenes()}.
#' @param markers A character vector of gene names to restrict this analysis to.
#'   See Details.
#' @param plotCurves Should a plot be drawn? default value is TRUE.
#' @param colors Color palette for the plot.
#'
#' @seealso \code{getMarkers} returns a cell cluster Shannon index that tends to
#' correlate well with the AUC metric returned by \code{getClassAUC}.
#' @return \code{getClassAUC} returns a numeric vector of length
#'   \code{ncol($specScore)} that contains the AUC for each cell cluster.
#' @export
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @examples
#' #randomly generated expression matrix and cell clusters
#' set.seed(1234)
#' exp = matrix(sample(0:20,1000,replace=TRUE), ncol = 20)
#' rownames(exp) = sapply(1:50, function(x) paste0("g", x))
#' cellType = sample(c("cell type 1","cell type 2"),20,replace=TRUE)
#' sg = sortGenes(exp, cellType)
#' classAUC = getClassAUC(sg)
#'
#' #"reasonably" separated clusters
#' data(sim)
#' sg = sortGenes(sim$exp, sim$cellType)
#' classAUC = getClassAUC(sg)
#'
#' #real data with three well separated clusters
#' data(kidneyTabulaMuris)
#' sg = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' classAUC = getClassAUC(sg)
getClassAUC = function(gs, markers = NULL, plotCurves = TRUE, colors = NULL) {


	if (is.null(colors)) {
		colors = c("#d64e3c7F","#7dd5547F","#8641c67F","#cfca457F","#7778cb7F","#59803d7F","#d04d9c7F","#73d6a87F","#492f607F","#ccc4977F","#7f343b7F","#72acc07F","#b97d407F","#c796b57F","#45483a7F", "#A020F07F","#00FF007F","#FFFF007F")
	}
	
	auc = rep(0, length = ncol(gs$specScore))
	if (is.null(markers)) {
		for (i in 1:ncol(gs$specScore)) {
			auc[i] = fastAUC((score(gs$specScore[order(gs$specScore[,i], decreasing=F),i])),(score(1:nrow(gs$specScore))))
		}

		if (plotCurves == TRUE) {
			plot((score(gs$specScore[order(gs$specScore[,1], decreasing=F),1])), (score(1:nrow(gs$specScore))), type = "l", col = colors[1], lwd = 4, ylab = "Genes (Increasing Specificity Score)", xlab = "Scaled Specificity Score")
			lines(seq(0,1,length.out = 100), (seq(0,1,length.out=100)), lty = 2, col = "#99999940", lwd = 3)
			for (i in 2:ncol(gs$specScore)) {
				lines((score(gs$specScore[order(gs$specScore[,i], decreasing=F),i])), (score(1:nrow(gs$specScore))), type = "l", col = colors[i], lwd = 4)
			}
			legend("bottomright", legend = colnames(gs$specScore), pch=23, col=colors, pt.bg=colors, pt.cex=2, cex=1, bty="n")
		}

	} else {
		temp = gs$specScore[which(rownames(gs$specScore) %in% markers),]
		for (i in 1:ncol(gs$specScore)) {
			auc[i] = fastAUC((score(temp[order(temp[,i], decreasing=F),i])), score(1:nrow(temp)))
		}

		if (plotCurves == TRUE) {
			plot((score(temp[order(temp[,1], decreasing=F),1])), score(1:nrow(temp)), type = "l", col = colors[1], lwd = 4, ylab = "Genes (Increasing Specificity Score)", xlab = "Scaled Specificity Score")
			lines(seq(0,1,length.out = 100), (seq(0,1,length.out=100)), lty = 2, col = "#99999940", lwd = 3)
			for (i in 2:ncol(gs$specScore)) {
				lines((score(temp[order(temp[,i], decreasing=F),i])), score(1:nrow(temp)), type = "l", col = colors[i], lwd = 4)
			}	
			legend("bottomright", legend = colnames(gs$specScore), pch=23, col=colors, pt.bg=colors, pt.cex=2, cex=1, bty="n")
		}
	}

	names(auc) = colnames(gs$specScore)
	return(auc)

}
