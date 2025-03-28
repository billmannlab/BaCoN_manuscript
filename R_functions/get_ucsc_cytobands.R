#' @import data.table
#' @importFrom UCSC.utils UCSC_dbselect
#' @export get_ucsc_cytobands

get_ucsc_cytobands <- function(version = "hg38") {
  ucsc_cytobands <- UCSC.utils::UCSC_dbselect(dbname = version,
                                              from = "cytoBand")

  setDT(ucsc_cytobands)

  ucsc_cytobands <- ucsc_cytobands[!grepl("alt|random|fix|Un_|M", chrom),
                                   .(chr = gsub("chr", "", chrom),
                                     start = chromStart,
                                     end = chromEnd, name,
                                     arm = fcase(grepl("p", name), "p",
                                                 grepl("q", name), "q",
                                                 default = NA))]
  ucsc_cytobands <- ucsc_cytobands[, .(start = min(start, na.rm = T),
                                       end = max(end, na.rm = T)),
                                   by = .(chr, arm)]

  return(ucsc_cytobands)
}
