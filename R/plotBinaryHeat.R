#' plotBinaryHeat
#'
#' Plot a heatmap binarized values for a select set of genes or peaks across all 
#' cells. Optionally also cluster those genes.
#'
#'\code{plotBinaryHeat} is similar to \code{plotMarkerHeat} but works on 
#' the binarized matrix, plotting only 0s and 1s (for scRNA-Seq, this means ignoring 
#' expression values and just plotting whether the gene is expressed or not). 
#' If \code{averageCells} is > 1, the fraction of the averaged cells where the gene 
#' or peak was detected will be plotted (a number between 0 and 1).
#'
#' \code{clusterGenesK} parameter controls the number of gene clusters when 
#' \code{clusterGenes} is \code{TRUE}. By default, it is equal to the number of
#' cell types/clusters. 
#'
#' By default, the heatmap plots every single cell in one column, this might
#' take forever if you have a lot of cells (would say >10k) or it can crash when 
#' you do not have enough RAM. If so, it might be good to set \code{averageCells} 
#' to \code{n} where \code{n} is the number of cells you 
#' want to average. For example, if \code{averageCells = 10}, every 10 cells will 
#' be averaged (without averaging across cell clusters) before plotting the heatmap.
#' if \code{averageCells =< 1}, no averaging happens. Hints: (1) If you want one column 
#' per cell cluster, set \code{averageCells} to a very high number (larger than the 
#' number of cells in the largest cluster). (2) Gene clustering occurs after cell 
#' averaging, so averaging will be useful in this case to prevent running k-means on 
#' binary values.
#'
#' @param exp A matrix of *binary* (0s and 1s) expression values. For example, the one 
#' found in  \code{sortGenes(...)$binary}.
#' @param classes A vector or factor of cell classes, whose length is equal to
#'   \code{ncol(exp)}.
#' @param markers A character vector of gene or peak names to plot in the heatmap.
#' @param colors Color palette used for the heatmap.
#' @param newOrder Reorder the clusters in the heatmap? See Examples.
#' @param clusterGenes Cluster genes before plotting?
#' @param clusterGenesK How many clusters should genes by clustered into? See
#'   Details.
#' @param averageCells Plot averages of cells instead of individual cells. You 
#' can use this when you have a large number of cells. See Details.
#' @param outs Should gene cluster output and pheatmap object be returned? FALSE
#'   by default.
#' @param plotheat Should the heatmap be drawn? TRUE by default.
#' @param gaps Should the heatmap have gaps between cell types and gene clusters? 
#' TRUE by default.
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
#' gs = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' mm = getMarkers(gs, quant = 0.999)
#'
#' #this plots a heatmap without reordering genes
#' plotBinaryHeat(gs$binary, gs$inputClass, mm$markers)
#'
#' @seealso
#' plotMarkerHeat
plotBinaryHeat = function(exp, classes, markers, colors = colorRampPalette(c("white","black"))(n=100), newOrder = 1:length(unique(classes)), clusterGenes = FALSE, clusterGenesK = length(unique(classes)), averageCells = 0, outs = FALSE, plotheat = TRUE, gaps = TRUE, seed = 10) {

	classes = as.integer(as.factor(classes))
	map = data.frame(oldO = 1:length(unique(classes)), newO = newOrder)
	classes = map$oldO[match(classes, map$newO)]
	
	
	temp = exp[which(rownames(exp) %in% markers),order(classes, decreasing = FALSE)]
	temp = temp[order(match(rownames(temp), markers)),] #catch error here? all clusters should be present after marker selection

	classesTemp = classes[order(classes, decreasing = FALSE)]
	if (averageCells > 1) {
		neh = lapply( 1:length(table(classesTemp)), function(i) rbind(apply(t(as.matrix(temp[,which(classesTemp == i)])), 2, function(x) binMean(x, averageCells))) )
		
		nah = unlist(sapply(1:length(neh), function(x) rep(x,nrow(neh[[x]]))))
		temp = t( do.call(rbind, neh) )
		colnames(temp) = nah
	} else {
		colnames(temp) = classesTemp
	}

	if (clusterGenes) {
		go = t(apply(temp, 1, function(x) x / (sqrt(sum(x^2)))))
		set.seed(seed)
		kk = kmeans(go, clusterGenesK, nstart = 10)$cluster
		temp = temp[order(kk, decreasing = FALSE),]
	}

	
	if (plotheat & (averageCells > 1)) {
		club = colnames(temp) #save column names
		cut = quantile(as.vector(temp), probs = 0.95, na.rm = TRUE)
		if (cut == 0) {
			cut = quantile(as.vector(temp), probs = 0.99, na.rm = TRUE)
		} 
		if (cut == 0) {
			cut = 1
			warning("Could not find an appropriate value to threshold color scale. Matrix may be too sparse. Please consider reporting this warning to mmibrahim@pm.me or at https://github.com/mahmoudibrahim/genesorteR/issues (preferred)")
		}
		temp[temp > cut] = cut
		colnames(temp) = club #reassign column names
	} else {cut = 1}

	if (clusterGenes) {
		if (plotheat) {
			if (gaps) {
				p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(0,cut,length.out = 101), gaps_col = cumsum(table(as.integer(colnames(temp)))), gaps_row =  cumsum(table(kk)), border_color = NA)
			} else {
				p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(0,cut,length.out = 101), border_color = NA)
			}
		} else {p = NULL}
		if (outs) {
			return(list(pheat = p, new_class_info = classes, gene_class_info = kk))
		}
    } else {
		if (plotheat) {
			if (gaps) {
				p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(0,cut,length.out = 101), gaps_col = cumsum(table(as.integer(colnames(temp)))), border_color = NA)
			} else {
				p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(0,cut,length.out = 101), border_color = NA)
			}
		} else {p = NULL}
		if (outs) {
			return(list(pheat = p, new_class_info = classes))
		}
	}

}
