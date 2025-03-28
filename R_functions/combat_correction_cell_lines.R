#' @import data.table
#' @importFrom sva ComBat
#' 
#' @export combat_correction_cell_lines

combat_correction_cell_lines <- \(matrix, method, 
                                  cell_line_covariates = covariates, 
                                  cell_line_column = "cell_line") {
  
  if (nrow(matrix) >= ncol(matrix)) {
    warning("The matrix has to have cell lines on rows and genes as columns!")}
  
  if (!method %in% names(cell_line_covariates)) {
    warning(str_c(method, " not found in the covariate file!"))}
  
  if (!all(rownames(matrix) %in% cell_line_covariates[, get(cell_line_column)])) {
    warning(str_c("Not all cell lines in covariate file!"))}
  
  .batch <- cell_line_covariates[rownames(matrix), on = cell_line_column, 
                                 get(method)]
  
  return(t(sva::ComBat(t(matrix), batch = .batch)))
}