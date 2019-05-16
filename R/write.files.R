#' write.files
#'
#' write.files saves gene-cluster sorting information obtained from
#' \code{sortGenes} to disk.
#'
#' Three files will be saved to disk, including the specificity score matrix,
#' the posterior cluster probability given the gene and the conditional
#' probability of observing the gene in a cluster. If \code{markers} is not
#' NULL, three additional files will be created and indicated with "_markers" in
#' their name.
#'
#' @param gs A list containing \code{$specScore}, \code{$postClustProb} and
#'   \code{$specScore} sparse matrices. Typically the output of
#'   \code{sortGenes()}.
#' @param prefix The prefix for saving the files.
#' @param markers Additionally, output files restricted to these genes. A
#'   character vector.
#' @author Mahmoud M Ibrahim <m3i@@selfishscience.net>
#'
#' @examples
#' data(sim)
#' sg = sortGenes(sim$exp, sim$cellType)
#' \dontrun{write.files(sg)} #write all files for all genes.
#' \dontrun{write.files(sg, markers = c("g1","g2"))} #additionally write files that are restricted to g1 and g2.
write.files = function(gs, prefix = "genesorteROuts", markers = NULL) {

  write.table(round(as.matrix(gs$condGeneProb),3), sep = "\t", file = paste0(prefix, "_condProbOfExp.txt"), quote = FALSE)
  write.table(round(as.matrix(gs$postClustProb),3), sep = "\t", file = paste0(prefix, "_postClustProb.txt"), quote = FALSE)
  write.table(round(as.matrix(gs$specScore),3), sep = "\t", file = paste0(prefix, "_specifictyScore.txt"),  quote = FALSE)

  for (i in 1:ncol(gs$specScore)) {
    tt = as.matrix(gs$specScore[order(gs$specScore[,i], decreasing = T),])
    write(paste(rownames(tt), round(score(tt[,i]),3), sep = ","), ncolumns = 1, file = paste0(prefix, "_genes_by_scaledSpecScore_", colnames(gs$specScore)[i], ".txt"))
  }

  if (!(is.null(markers))) {
    ww = which(rownames(gs$condGeneProb) %in% markers)
    write.table(round(as.matrix(gs$condGeneProb[ww,]),3), sep = "\t", file = paste0(prefix, "_condProbOfExp_markers.txt"), quote = FALSE)
    write.table(round(as.matrix(gs$postClustProb[ww,]),3), sep = "\t", file = paste0(prefix, "_postClustProb_markers.txt"), quote = FALSE)
    write.table(round(as.matrix(gs$specScore[ww,]),3), sep = "\t", file = paste0(prefix, "_specifictyScore_markers.txt"), quote = FALSE)
  }

}
