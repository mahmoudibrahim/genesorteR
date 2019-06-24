#' A subset of the Tabula Muris SmartSeq2 data, restricted to three cell types from Kidney.
#'
#' A subset of the Tabula Muris data with 215 cells and 23341 genes. There are three
#' cell types, endothelial cells: 122 cells, kidney collecting duct epithelial cell: 78 cells and leukoocytes 15 cells.
#'
#' @format A list with two components:
#' \describe{\item{exp}{Gene expression sparse matrix of class dgCMatrix, with normalized log expression values. Rows are genes, columns are cells.} \item{cellType}{A character vector containing the cell types.}}
#' @references The Tabula Muris Consortium Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature, 562:7727, 2018
"kidneyTabulaMuris"
