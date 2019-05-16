#' plotMarkerScores
#'
#' plotMarkerScores plots scatter plots of the specificity score and the gene
#' specificity Shannon index to investigate \code{getMarkers} results.
#' @param mm A list. The output of \code{getMarkers}
#' @param gs A list. The output of \code{sortGenes}
#' @param markers A character vector containing the genes to highlight in the
#'   scatter plot.
#' @param colors A color palette to use for the scatter plot.
#' @author Mahmoud M Ibrahim <m3i@@selfishscience.net>
#'
#' @export
#' @examples
#' data(kidneyTabulaMuris)
#' sg = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' mm = getMarkers(sg,quant=0.99)
#' plotMarkerScores(mm, sg)
plotMarkerScores = function(mm, gs, markers = mm$markers, colors = c("#99999920","orangered4")) {
  coloring = rep(colors[1],length(mm$gene_entropy))
  coloring[which(names(mm$gene_entropy) %in% markers)] = colors[2]

  par(mfrow = c(1,2))
  plot(mm$gene_entropy, mm$maxScaledSpecScore, col = coloring, pch = 19, xlab = "Specificity Shannon Index", ylab = "Max. Scaled Specificity Score")
  plot(mm$gene_entropy, apply(gs$specScore, 1, max), col = coloring, pch = 19, xlab = "Specificity Shannon Index", ylab = "Max. Specificity Score")
}
