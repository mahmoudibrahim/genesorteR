#' plotTopBinaryHeat
#'
#' Plot a heatmap of binarized values for the top genes or peaks in each cluster, 
#' defined by \code{sortGenes()}.
#'
#'
#' \code{plotTopMarkerHeat} is similar to \code{plotTopMarkerHeat} but works on 
#' the binarized matrix, plotting only 0s and 1s (for scRNA-Seq, this means ignoring 
#' expression values and just plotting whether the gene is expressed or not). 
#' If \code{averageCells} is > 1, the fraction of the averaged cells where the gene 
#' or peak was detected will be plotted (a number between 0 and 1).
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
#' plotTopBinaryHeat(gs) # plots the top 10 genes for each cluster
#'
#' #now plot the top 20 genes and average every 5 cells
#' plotTopBinaryHeat(gs, top_n= 20, averageCells=5)
#' 
#' @seealso
#' plotTopMarkerHeat
plotTopBinaryHeat = function(sg, top_n = 10, colors = colorRampPalette(c("white","black"))(n=100), newOrder = 1:length(unique(sg$inputClass)), averageCells = 0, gaps = TRUE, outs = FALSE, plotheat = TRUE) {

	classes = as.integer(as.factor(sg$inputClass))
	map = data.frame(oldO = 1:length(unique(classes)), newO = newOrder)
	classes = map$oldO[match(classes, map$newO)]
	
	markers = list()
	for (i in 1:length(unique(classes))) {
		markers[[i]] = rownames(sg$specScore[order(sg$specScore[,map[i,2]], decreasing = T),])[1:top_n]
	}
	final_sel_odd = unique(unlist(markers))
	names(markers) = colnames(sg$specScore[,newOrder])

	
	temp = sg$binary[which(rownames(sg$binary) %in% final_sel_odd),order(classes, decreasing = FALSE)]
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

	
	
	if (plotheat) {
		if (gaps) {
			p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(0,cut,length.out = 101), gaps_col = cumsum(table(as.integer(colnames(temp)))), border_color = NA)
		} else {
			p=pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colors, display_numbers = F, fontsize = 5, show_colnames=FALSE, show_rownames=TRUE, breaks = seq(0,cut,length.out = 101), border_color = NA)
		}
	}
	
	if (outs) {
		return(markers)
	}
	

}
