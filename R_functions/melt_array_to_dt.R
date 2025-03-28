#' @import data.table
#' @export melt_array_to_dt

## ---- melt_array_to_dt function ----

melt_array_to_dt <- \(data, row.name = "row",
                      col.name = "col",
                      value.name = "value") {
  output <- list()

  if (!base::is.null(rownames(data))) {
    output[[row.name]] <- rep(rownames(data), ncol(data))
  } else {
    output[[row.name]] <- rep(as.character(1:nrow(data)), ncol(data))
  }

  if (!base::is.null(colnames(data))) {
    output[[col.name]] <- rep(colnames(data), each = nrow(data))
  } else {
    output[[col.name]] <- rep(as.character(1:ncol(data)), each = nrow(data))
  }

  output[[value.name]] <- base::as.vector(data)
  data.table(data.frame(output))
}
