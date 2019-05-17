#' plotMarkerHeat
#'
#' Plot a heatmap of expression values for a select set of genes. Optionally
#' also cluster those genes.
#' @param exp A matrix of expression values. Typically the one supplied to
#'   \code{sortGenes}.
#' @param classes A vector or factor of cell classes, whose length is equal to
#'   \code{ncol(exp)}.
#' @param markers A character vector of gene names to plot in the heatmap.
#' @param colors Color palette used for the heatmap.
#' @param newOrder Reorder the clusters in the heatmap? See Examples.
#' @param clusterGenes Cluster genes before plotting?
#' @param clusterGenesK How many clusters should genes by clustered into? See
#'   Details.
#' @param outs Should gene cluster output and pheatmap object be returned? FALSE
#'   by default.
#' @param plotheat Should the heatmap be drawn? TRUE by default.
#' @param seed Randomization seed used for gene clustering initialization.
#' @return If \code{outs} is TRUE, \code{plotMarkerHeat} returns a list
#'   containing: \item{p}{The pheatmap object corresponding to the plot.}
#'   \item{gene_class_info}{Gene cluster assignments if gene clustering was
#'   requested. See Details.} \item{new_class_info}{The new order of cell
#'   clusters. See Examples.}
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @export
#' @examples
#' data(kidneyTabulaMuris)
#' sg = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' mm = getMarkers(sg, quant = 0.999)
#'
#' #this plots a heatmap without reordering genes
#' plotMarkerHeat(sg$inputMat, sg$inputClass, mm$markers)
#'
#' #so now cluster genes and return the clustering: a better looking plot
#' pp = plotMarkerHeat(sg$inputMat, sg$inputClass, mm$markers, clusterGenes=TRUE, outs=TRUE)
#' pp$gene_class_info #cell clusters
#'
#' #only cluster genes, do not make plots
#' pp = plotMarkerHeat(sg$inputMat, sg$inputClass, mm$markers, clusterGenes=TRUE, 
#' outs=TRUE, plotheat=FALSE)
#'
#'
#' #reorder cell clusters in the heatmap
#' pp = plotMarkerHeat(sg$inputMat, sg$inputClass, mm$markers, clusterGenes=TRUE, newOrder = c(2,1,3))
plotMarkerHeat = function(exp, classes, markers, colors = colorRampPalette(rev(c("orangered4","orangered","gray90","dodgerblue","dodgerblue4")))(n=100), newOrder = 1:length(unique(classes)), clusterGenes = FALSE, clusterGenesK = length(unique(classes)), outs = FALSE, plotheat = TRUE, seed = 10) {

  if (is.character(classes) || is.factor(classes)) {
    theoldorder = levels(as.factor(classes))
    classes = as.vector(classes)

    for (i in 1:length(newOrder)) {
      classes[which(classes == theoldorder[i])] = newOrder[i]
    }
    classes = as.numeric(classes)
  } else if (is.numeric(classes)) {
    map = data.frame(oldO = 1:length(unique(classes)), newO = newOrder)
    classes = map$oldO[match(classes, map$newO)]
  }


  temp = as.matrix(exp[which(rownames(exp) %in% markers),order(classes, decreasing = FALSE)])
  temp = temp[order(match(rownames(temp), markers)),]

  if (clusterGenes) {
    go = t(apply(temp, 1, function(x) x / (sqrt(sum(x^2)))))
    set.seed(seed)
    kk = kmeans(go, clusterGenesK, nstart = 10)$cluster
    temp = temp[order(kk, decreasing = FALSE),]
  }

  temp = t(apply(temp, 1, scale))
  cut = quantile(as.vector(abs(temp)), probs = 0.95)
  temp[temp > cut] = cut
  temp[temp < -cut] = -cut

  colnames(temp) = (1:ncol(temp))
  anno = data.frame(cellType = t((classes)))
  colnames(anno) = colnames(temp)


  if (clusterGenes) {
    if (plotheat) {
      p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, gaps_col =  cumsum(table(sort(classes))), gaps_row =  cumsum(table(sort(kk))), display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(-cut,cut,length.out = 101))
    } else {p = NULL}
    if (outs) {
      return(list(pheat = p, new_class_info = classes, gene_class_info = kk))
    }
  } else {
    if (plotheat) {
      p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, gaps_col =  cumsum(table(sort(classes))), display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(-cut,cut,length.out = 101))
    } else {p = NULL}
    if (outs) {
      return(list(pheat = p, new_class_info = classes))
    }
  }

}
