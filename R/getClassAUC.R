#' getClassAUC
#'
#' getClassAUC implements one way to investigate clustering quality. It processes the
#' output of \code{sortGenes} to obtain a curve for each cell cluster of all
#' gene specificity scores against their ranking in the cluster. The Area Under
#' the Curve (AUC) can be used as a measure of clustering quality in terms of the
#' possibility to identify cell clusters using a few marker genes. See Details.
#'
#' Given the specificity score for all genes in a certain cell cluster, we can
#' assume that a well-separated easily-identified cell cluster will have a
#' relatively small number of genes that have a very high specificity score. Top
#' marker genes for a cell cluster that is poorly separated from other cell
#' clusters will have average or low specificity scores. Sorting the genes for
#' each cell cluster by their specificity score and plotting the scaled scores
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
#' @author Mahmoud M Ibrahim <m3i@@selfishscience.net>
#' @examples
#' #randomly generated cell clusters
#' set.seed(1234)
#' exp = matrix(sample(0:20,1000,replace=TRUE), ncol = 20)
#' rownames(exp) = sapply(1:50, function(x) paste0("g", x))
#' cellType = sample(c("cell type 1","cell type 2"),20,replace=TRUE)
#' sg = sortGenes(exp, cellType)
#' classAUC = getClassAUC(sg)
#' mean(classAUC) #overall clustering quality
#'
#' #"reasonably" separated clusters, with a few clear markers
#' data(sim)
#' sg = sortGenes(sim$exp, sim$cellType)
#' classAUC = getClassAUC(sg)
#' mean(classAUC)
#'
#' #real data with three well separated clusters
#' data(kidneyTabulaMuris)
#' sg = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' classAUC = getClassAUC(sg)
#' mean(classAUC)
getClassAUC = function(gs, markers = NULL, plotCurves = TRUE, colors = c("#d64e3c7F","#7dd5547F","#8641c67F","#cfca457F","#7778cb7F","#59803d7F","#d04d9c7F","#73d6a87F","#492f607F","#ccc4977F","#7f343b7F","#72acc07F","#b97d407F","#c796b57F","#45483a7F", "#A020F07F","#00FF007F","#FFFF007F")) {

  auc = rep(0, length = ncol(gs$specScore))
  if (is.null(markers)) {
    for (i in 1:ncol(gs$specScore)) {
      auc[i] = getAUC(score(1:nrow(gs$specScore)), score(gs$specScore[order(gs$specScore[,i], decreasing=F),i]))
    }

    if (plotCurves == TRUE) {
      plot(score(1:nrow(gs$specScore)), score(gs$specScore[order(gs$specScore[,1], decreasing=F),1]), ylim = c(1,0), type = "l", col = colors[1], lwd = 4, xlab = "Genes (Increasing Specificity Score)", ylab = "Scaled Specificity Score")
      for (i in 2:ncol(gs$specScore)) {
        lines(score(1:nrow(gs$specScore)), score(gs$specScore[order(gs$specScore[,i], decreasing=F),i]), ylim = c(1,0), type = "l", col = colors[i], lwd = 4)
      }
      legend("bottomleft", legend = colnames(gs$specScore), pch=23, col=colors, pt.bg=colors, pt.cex=2, cex=1, bty="n")
    }

  } else {
    temp = gs$specScore[which(rownames(gs$specScore) %in% markers),]
    for (i in 1:ncol(gs$specScore)) {
      auc[i] = getAUC(score(1:nrow(temp)), score(temp[order(temp[,i], decreasing=F),i]))
    }

    if (plotCurves == TRUE) {
      plot(score(1:nrow(temp)), score(temp[order(temp[,1], decreasing=F),1]), ylim = c(1,0), type = "l", col = colors[1], lwd = 4, xlab = "Genes (Increasing Specificity Score)", ylab = "Scaled Specificity Score")
      for (i in 2:ncol(gs$specScore)) {
        lines(score(1:nrow(temp)), score(temp[order(temp[,i], decreasing=F),i]), ylim = c(1,0), type = "l", col = colors[i], lwd = 4)
      }
      legend("bottomleft", legend = colnames(gs$specScore), pch=23, col=colors, pt.bg=colors, pt.cex=2, cex=1, bty="n")
    }
  }

  names(auc) = colnames(gs$specScore)
  return(auc)

}
