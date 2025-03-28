#' @importFrom readr read_lines
#' @export print_md

print_md <- \(file) {base::cat(readr::read_lines(file), sep = "\n")}
