#' @export empty_array

## ---- empty_array function ----

empty_array <- \(dnames, values = NA) {
  return(base::array(data = values,
                     dim = base::sapply(dnames, base::length),
                     dimnames = dnames))
}
