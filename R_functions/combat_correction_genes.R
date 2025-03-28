#' @import data.table
#' @importFrom sva ComBat
#' @export combat_correction_genes


combat_correction_genes <- \(matrix, 
                             method, 
                             gene_covariates = gene_database, 
                             gene_column = "symbol") {
  
  if (nrow(matrix) >= ncol(matrix)) {
    warning("The matrix has to have cell lines on rows and genes as columns!")}
  
  if (method %in% c("chromosome", "chromosome_arm")) {
    .batch <- distinct(gene_covariates[, .SD, .SDcols = c(gene_column, method)])
    .batch <- .batch[colnames(matrix), on = gene_column]
    setnames(.batch, method, "batch")
    .batch[is.na(batch), batch := "Unknown"]
    
    if (.batch[, any(duplicated(get(gene_column)))]) {
      warning(str_c("Duplicated genes in batch information!"))}
    
    .batch <- .batch[, get("batch")]
    
    sva::ComBat(matrix, batch = .batch)}
}