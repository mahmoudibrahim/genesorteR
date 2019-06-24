#' plotCorrelationHeat
#'
#' \code{plotCorrelationHeat} uses the specificity scores returned by \code{sortGenes} to 
#' correlate cell clusters with each other and plot a heatmap.
#'
#' @param gs The output of \code{sortGenes}.
#' @param markers Restrict correlation analysis to those genes. A character
#'   vector.
#' @param corMethod Correlation method, will passed to \code{method} in the
#'   function \code{cor}.
#' @param colors Color palette for drawing the heatmap
#' @param outs Should the \code{pheatmap} object and correlation matrix be 
#' returned? FALSE by default.
#' @param displayNumbers Should correlation values be displayed on the heatmap?
#'   TRUE by default.
#' @return If \code{outs} is TRUE, the pheatmap object and the correlation matrix will 
#' be returned.
#' @export
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @examples
#' data(kidneyTabulaMuris)
#' gs = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' plotCorrelationHeat(gs)
#'
#' #user only marker genes and spearman correlation
#' mm = getMarkers(gs, quant = 0.95)
#' plotCorrelationHeat(gs, markers = mm$markers, corMethod = "spearman")
#'
#' #do not write correlation values, useful if there are many cell clusters
#' plotCorrelationHeat(gs, markers = mm$markers, corMethod = "spearman", displayNumbers = FALSE)
plotCorrelationHeat = function(gs, markers = NULL, corMethod = "pearson", colors = colorRampPalette(rev(c("orangered4","orangered","gray90","dodgerblue","dodgerblue4")))(n=100), outs = FALSE, displayNumbers = TRUE) {

  if (is.null(markers)) {
    cc = cor(as.matrix(gs$specScore), method = corMethod)
  } else {
    cc = cor(as.matrix(gs$specScore[which(rownames(gs$specScore) %in% markers),]), method = corMethod)
  }

  rownames(cc) = colnames(gs$specScore)
  colnames(cc) = colnames(gs$specScore)

  p=pheatmap(cc, breaks = seq(-1,1,length.out = 99), show_rownames=T, show_colnames=T, display_numbers=displayNumbers, color = colors)

  if (outs) {
    return(list(pheat = p, corr = cc))
  }

}
