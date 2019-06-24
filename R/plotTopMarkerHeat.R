#' plotTopMarkerHeat
#'
#' Plot a heatmap of expression values for the top genes in each cluster, 
#' defined by \code{sortGenes()}.
#'
#'
#' \code{plotTopMarkerHeat} is a convenience wrapper around \code{plotMarkerHeat}
#' that plots a heatmap of the top \code{top_n} (10 by default) genes in each cell 
#' cluster. Unlike \code{plotMarkerHeat}, \code{plotTopMarkerHeat} takes the output
#' of \code{sortGenes} as the only required input.
#'
#'
#' @param sg A list, typically the output of \code{sortGenes}.
#' @param top_n The number of top genes to plot for each cluster.
#' @param colors Color palette used for the heatmap.
#' @param newOrder Reorder the clusters in the heatmap? See Examples.
#' @param averageCells Plot averages of cells instead of individual cells. You 
#' can use this when you have a large number of cells. See Details.
#' @param gaps Should the heatmap have gaps between cell types and gene clusters? 
#' TRUE by default.
#' @param outs Should the top genes names be returned? FALSE by default.
#' @param plotheat Should the heatmap be drawn? TRUE by default.
#' @return If \code{outs} is TRUE, \code{plotMarkerHeat} returns a list containing 
#' the top \code{n} marker genes for each cluster.
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#' @export
#' @examples
#' data(kidneyTabulaMuris)
#' gs = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
#' plotTopMarkerHeat(gs) # plots the top 10 genes for each cluster
#'
#' #now plot the top 20 genes and average every 5 cells
#' plotTopMarkerHeat(gs, top_n= 20, averageCells=5)
#' 
#' #just identify the top 20 genes, do not make a plot
#' plotTopMarkerHeat(gs, top_n= 20, averageCells=5, outs = TRUE, plotheat = FALSE)
#'
#' @seealso
#' plotMarkerHeat
plotTopMarkerHeat = function(sg, top_n = 10, colors = colorRampPalette(rev(c("orangered4","orangered","gray90","dodgerblue","dodgerblue4")))(n=100), newOrder = 1:length(unique(sg$inputClass)), averageCells = 0, gaps = TRUE, outs = FALSE, plotheat = TRUE) {

	classes = as.integer(as.factor(sg$inputClass))
	map = data.frame(oldO = 1:length(unique(classes)), newO = newOrder)
	classes = map$oldO[match(classes, map$newO)]
	
	markers = list()
	for (i in 1:length(unique(classes))) {
		markers[[i]] = rownames(sg$specScore[order(sg$specScore[,map[i,2]], decreasing = T),])[1:top_n]
	}
	final_sel_odd = unique(unlist(markers))
	names(markers) = colnames(sg$specScore[,newOrder])

	
	temp = sg$inputMat[which(rownames(sg$inputMat) %in% final_sel_odd),order(classes, decreasing = FALSE)]
	temp = temp[order(match(rownames(temp), final_sel_odd)),] #catch error here? all clusters should be present after marker selection

	classesTemp = classes[order(classes, decreasing = FALSE)]
	if (averageCells > 1) {
		neh = lapply( 1:length(table(classesTemp)), function(i) rbind(apply(t(as.matrix(temp[,which(classesTemp == i)])), 2, function(x) binMean(x, averageCells))) )
		
		nah = unlist(sapply(1:length(neh), function(x) rep(x,nrow(neh[[x]]))))
		temp = t( do.call(rbind, neh) )
		colnames(temp) = nah
	} else {
		colnames(temp) = classesTemp
	}

		
	if (plotheat) {
		club = colnames(temp) #save column names
		temp = t(apply(temp, 1, scale))
		cut = quantile(as.vector(abs(temp)), probs = 0.95)
		temp[temp > cut] = cut
		temp[temp < -cut] = -cut
		colnames(temp) = club #reassign column names
	}
	

	if (plotheat) {
		if (gaps) {
			p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(-cut,cut,length.out = 101), gaps_col = cumsum(table(as.integer(colnames(temp)))), border_color = NA)
		} else {
			p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(-cut,cut,length.out = 101), border_color = NA)
		}
	}
	
	if (outs) {
		return(markers)
	}
	

}
