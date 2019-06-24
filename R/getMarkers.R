#' getMarkers
#'
#' getMarkers processes the output of \code{sortGenes} to select a relatively
#' small set of marker genes. 
#'
#'
#' getMarkers relies on calculating an entropy-like metric (called here the Gene 
#' Specificity Shannon Index) calculated on the scaled gene-cluster specificity 
#' score. The intuition is that genes with low Shannon Index are either lowly 
#' expressed or highly specific to one or few cell clusters. Therefore, we select
#' the top n\% of genes according to the scaled specificity score and then cluster 
#' those genes based on their Shannon Index into genes with high Shannon Index 
#' and genes with low Shannon Index. \code{n} is controlled by \code{quant}, 
#' the default value 0.99 means top 1 \code{%} of genes. Using lower values up to 
#' 0.95 may be ok but lower values than 0.95 are not recommended nor necessary. It 
#' can be increased to 0.999 for example to obtain a smaller set of genes, especially
#' for data with many cell clusters. Scaling the specificity score for each cluster
#' is done to try to guarantee that markers for each cluster will eventually be 
#' selected even if the cluster is not absolutely well-separated.
#'
#' It can also return the mutual information between each gene and cell clusters, 
#' as well as a Cluster Shannon Index, indicating how well separated cell clusters are 
#' from each other. Cluster Shannon Index is calculated on the scaled specificity score 
#' of selected marker genes for each cluster. The intuition is that, when restricted to
#' top marker genes, well-defined clusters will have few high scoring genes 
#' (low Cluster Shannon Index).
#'
#' @param gs The output of \code{sortGenes()}.
#' @param quant A number greater than zero and smaller than one. 0.99 by
#'   default. See Details.
#' @param mutualInfo Logical. If \code{TRUE}, the mutual information between
#'   gene expression and cell clustering will be returned for each gene. \code{FALSE}
#'   by default.
#' @param classEnt Logical. If \code{TRUE}, a "Cluster Shannon Index" is
#'   returned for each cluster. See Details. \code{FALSE} by default.
#' @return \code{getMarkers} returns a list with the following components:
#'   \item{markers}{A character vector containing the names of selected marker
#'   features.} \item{maxScaledSpecScore}{A numeric vector of length
#'   \code{nrow(gs$specScore)}, including the maximum scaled specificity score for
#'   each gene across all cell clusters.} \item{gene_shannon_index}{A numeric vector
#'   of the same length as \code{maxScaledSpecScore}, including the Gene Shannon
#'   Index for each gene. See Details.} \item{mutInfo}{A numeric vector of the
#'   same length as \code{maxScaledSpecScore}, including the mutual information
#'   between the expression of each gene and the cell clustering.}
#'   \item{classEntropy}{A numeric vector of the same length as
#'   \code{ncol(specScore)}, including the Class Shannon Index for each cluster
#'   based on the selected marker genes.}
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @export
#' @examples
#' #randomly generated data and cell clusters, almost no markers are found
#' set.seed(1234)
#' exp = matrix(sample(0:20,1000,replace=TRUE), ncol = 20)
#' rownames(exp) = sapply(1:50, function(x) paste0("g", x))
#' cellType = sample(c("cell type 1","cell type 2"),20,replace=TRUE)
#' sg = sortGenes(exp, cellType)
#' mm = getMarkers(sg,quant=0.95)
#' length(mm$markers) #only one marker gene was found
#'
#'
#' #"reasonably" separated clusters, with a few clear markers
#' data(sim)
#' gs = sortGenes(sim$exp, sim$cellType)
#' mm = getMarkers(gs,quant=0.95)
#' length(mm$markers)
#'
#' #real data with three well separated clusters
#' data(kidneyTabulaMuris)
#' gs = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' mm = getMarkers(gs, quant = 0.99)
#' length(mm$markers) #we found 109 candidate markers
#' #we want to get a more focused list:
#' mm = getMarkers(gs, quant = 0.999)
#' length(mm$markers) #11 genes that can alone descriminate between the cell types
#' @seealso
#' getPValues
getMarkers = function(gs, quant = 0.99, mutualInfo = FALSE, classEnt = FALSE) {

  scored = apply(gs$specScore, 2, score)
  ent = apply(scored, 1, getEntropy)
  maxi = apply(scored, 1, max)
  ww = which(maxi > quantile(maxi, probs=quant))
  m = Mclust(ent[ww], 2, modelNames="E", verbose = FALSE)
  markers = names( which(m$classification == (which.min(m$parameters$mean)) ) )

  mutInfo = NULL
  if (mutualInfo == TRUE) {
    mutInfo = apply(gs$binary, 1, function(x) getMutInfo(x,gs$inputClass))
  }

  classEntropy = NULL
  if (classEnt == TRUE) {
    classEntropy = apply(scored, 2, getEntropy)
    names(classEntropy) = colnames(gs$specScore)
  }

  return(list(gene_shannon_index = ent, maxScaledSpecScore = maxi, markers = markers, mutInfo = mutInfo, classEntropy = classEntropy))

}
