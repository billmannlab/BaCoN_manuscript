#' @export cacheR

## ---- cacheR function ----

cacheR <- \(.,
            fname,
            dir = file.path(getwd(), "cache"),
            import = T,
            verbose = F) {

  if (!base::dir.exists(dir)) {
    base::dir.create(dir, recursive = T, showWarnings = F)
  }

  f <- base::file.path(dir, base::paste0(fname, ".rds"))

  messages <- list(a = base::paste0("Importing '", fname, "'."),
                   b = base::paste0("'", fname, "' exists. Skipping."),
                   c = base::paste0("Computing and caching '", fname, "'..."))

  if (file.exists(f)) {
    if (import) {
      if (verbose) {message(messages$a)}
      return(readRDS(f))}
    if (!import) {
      if (verbose) {message(messages$b)}
      return(NULL)
    }
  }

  if (!file.exists(f)) {
    if (verbose) {message(messages$c)}
    x <- base::eval(base::substitute(.), envir = parent.frame())
    saveRDS(x, file = base::file.path(dir, base::paste0(fname, ".rds")))
    return(x)
  }
}

