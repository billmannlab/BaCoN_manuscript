#' @import data.table
#' @importFrom stringr str_c
#' @export collect_predictions

## ---- collect_predictions ----

collect_predictions <- \(matrix,
                         row.name = "expression_gene", col.name = "chronos_gene",
                         value.name = "score",
                         reference_matrix = NULL,
                         ref.name = "reference",
                         .which = "positive",
                         pairs,
                         verbose = F) {

  if (verbose) {
    message(stringr::str_c("Predicting", pairs, "pairs based on", .which, "scores...", sep = " "))
    if (is.null(reference_matrix)) {message("No reference matrix provided.")}
    }

  if (.which == "positive") {
    th <- sort(matrix[matrix > 0], decreasing = T)[pairs]
    i <- which(matrix >= th); i_arr_ind <- which(matrix >= th, arr.ind = T)
  }

  if (.which == "negative") {
    th <- sort(matrix[matrix < 0], decreasing = T)[pairs]
    i <- which(matrix <= th); i_arr_ind <- which(matrix <= th, arr.ind = T)
  }

  x <- data.table(i_arr_ind, score = matrix[i])
  if (!is.null(reference_matrix)) {x[, reference := reference_matrix[i]]}

  x[, `:=`(rowname = rownames(matrix)[row], colname = colnames(matrix)[col])]

  x[, `:=`(rank = .I, sorted_pair = mapply(sort_gene_pairs, rowname, colname))]

  if (is.null(reference_matrix)) {x <- x[order(score, decreasing = c("positive" = T, "negative" = F)[[.which]])]}
  if (!is.null(reference_matrix)) {x <- x[order(score, reference, decreasing = c("positive" = T, "negative" = F)[[.which]])]}

  data.table::setnames(x, old = c("rowname", "colname", "score"), new = c(row.name, col.name, value.name))
  if (!is.null(reference_matrix)) {data.table::setnames(x, old = "reference", new = ref.name)}

  return(x)
}
