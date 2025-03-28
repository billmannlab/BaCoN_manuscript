#' @export mkdir

## ---- mkdir function ----

mkdir <- \(directory) {
  if (!base::dir.exists(directory)) {
    base::dir.create(directory, recursive = T, showWarnings = F)
  }
}
