#' write.files
#'
#' write.files saves gene-cluster sorting information obtained from
#' \code{sortGenes} to disk.
#'
#' Three files will be saved to disk, including the specificity score matrix,
#' the posterior cluster probability given the gene and the conditional
#' probability of observing the gene in a cluster. If \code{markers} is not
#' NULL, three additional files will be created and indicated with "_markers" in
#' their name. If \code{eachCluster} is \code{TRUE}, additional files (one per 
#' cluster will be written which includes all genes and their scaled specificity
#' score for each cluster, with genes sorted by their score.
#'
#' @param gs The output of \code{sortGenes()}.
#' @param prefix The prefix for saving the files.
#' @param markers Additionally, output files restricted to these genes. A
#'   character vector.
#' @param eachCluster Also write individual files for each cluster?
#' @author Mahmoud M Ibrahim <mmibrahim@pm.me>
#'
#' @examples
#' data(sim)
#' gs = sortGenes(sim$exp, sim$cellType)
#' \dontrun{write.files(gs)} #write all files for all genes.
#' \dontrun{write.files(gs, markers = c("g1","g2"))} #additionally write files that are restricted to genes g1 and g2.
write.files = function(gs, prefix = "genesorteROuts", markers = NULL, eachCluster = FALSE) {

	write.table(round(as.matrix(gs$condGeneProb),3), sep = "\t", file = paste0(prefix, "_condProbOfExp.csv"), quote = FALSE)
	write.table(round(as.matrix(gs$postClustProb),3), sep = "\t", file = paste0(prefix, "_postClustProb.csv"), quote = FALSE)
	write.table(round(as.matrix(gs$specScore),3), sep = "\t", file = paste0(prefix, "_specifictyScore.csv"),  quote = FALSE)

	if (eachCluster) {
		for (i in 1:ncol(gs$specScore)) {
			tt = as.matrix(gs$specScore[order(gs$specScore[,i], decreasing = TRUE),])
			write(paste(rownames(tt), round(score(tt[,i]),3), sep = ","), ncolumns = 1, file = paste0(prefix, "_genes_by_scaledSpecScore_", gsub("/", "_", colnames(gs$specScore)[i], fixed = TRUE), ".csv"))
		}
	}

	if (!(is.null(markers))) {
		ww = which(rownames(gs$condGeneProb) %in% markers)
		write.table(round(as.matrix(gs$condGeneProb[ww,]),3), sep = "\t", file = paste0(prefix, "_condProbOfExp_markers.csv"), quote = FALSE)
		write.table(round(as.matrix(gs$postClustProb[ww,]),3), sep = "\t", file = paste0(prefix, "_postClustProb_markers.csv"), quote = FALSE)
		write.table(round(as.matrix(gs$specScore[ww,]),3), sep = "\t", file = paste0(prefix, "_specificityScore_markers.csv"), quote = FALSE)
	}

}
