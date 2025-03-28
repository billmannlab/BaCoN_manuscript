#' @export get_ensembl_paralogs

get_ensembl_paralogs <- \(path, 
                          min_identity = NA, 
                          include_ensembl_families = F, 
                          only_pairs = T) {
  
  ensembl_paralogs <- sapply(list.files(path, full.names = T), 
                             \(.) {x <- readRDS(.); setDT(x)}, simplify = F) %>% 
    Reduce(\(x1, x2) {merge(x1, x2, by = c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene"))}, .)
  
  setnames(ensembl_paralogs, \(.) gsub("hsapiens_|paralog_", "", .))
  
  ensembl_paralogs <- ensembl_paralogs[
    !external_gene_name == "" & !associated_gene_name == "", 
    .(id1 = ensembl_gene_id, id2 = ensembl_gene, 
      name1 = external_gene_name, 
      name2 = associated_gene_name, 
      ident1 = perc_id_r1, 
      ident2 = perc_id, 
      subtype, orthology_type)] %>% 
    melt.data.table(measure.vars = c("ident1", "ident2"), value.name = "identity")
  
  
  ensembl_paralogs <- ensembl_paralogs[
    , .(ident_min = min(identity, na.rm = T), 
        ident_max = max(identity, na.rm = T)), 
    by = .(id1, id2, name1, name2, subtype, orthology_type)]
  
  
  ensembl_paralogs[, sorted_pair := sort_gene_pairs(name1, name2)]
  
  if (!is.na(min_identity)) {
    ensembl_paralogs <- ensembl_paralogs[ident_max > min_identity]
  }
  
  #if (include_ensembl_families) {
  #  families <- sapply(ensembl_paralogs[, 
  #                                      union(name1, name2)], 
  #                     \(g) {union(ensembl_paralogs[name1 == g, name2], 
  #                                 ensembl_paralogs[name2 == g, name1])}, 
  #                     simplify = F)
  #  ensembl_paralogs[, family_size := mapply(
  #    \(g1, g2) {length(union(families[[g1]], families[[g2]]))}, name1, name2)]
  
  #  d <- copy(ensembl_paralogs[, fam_id := "none"])
  #  
  #  families <- list()
  #  fam_counter <- 1
  #  pb <- progress::progress_bar$new(width = 50, total = ensembl_paralogs[, .N], force = T)
  #  while ("none" %in% d[, fam_id]) {
  #    n <- str_c("family_", fam_counter)
  #    families[[n]] <- d[fam_id == "none"][1, c(id1, id2)]
  #    new <- copy(d[fam_id == "none" & 
  #                    (id1 %in% families[[n]] | id2 %in% families[[n]])])
  #    while (new[, .N] != 0) {
  #      families[[n]] <- union(families[[n]], new[, c(id1, id2)])
  #      new <- copy(d[fam_id == "none" & 
  #                      (id1 %in% families[[n]] | id2 %in% families[[n]])])
  #      d[id1 %in% new[, c(id1, id2)] | 
  #          id2 %in% new[, c(id1, id2)], fam_id := n]}
  #    d <- d[fam_id == "none"]
  #    fam_counter <- fam_counter + 1
  #    pb$update(d[, .N] / ensembl_paralogs[, .N])}
  #
  #  for (n in names(families)) {
  #    ensembl_paralogs[id1 %in% families[[n]] | 
  #                       id2 %in% families[[n]], `:=`(fam_id = n)]}
  #}
  
  #if (!include_ensembl_families) {
  ensembl_paralogs <- copy(ensembl_paralogs)[, .(ident = max(ident_max, na.rm = T)), by = .(sorted_pair, orthology_type, subtype)]
  #  }
  
  #  if (include_ensembl_families) {
  #    ensembl_paralogs <- copy(ensembl_paralogs)[
  #      , .(ident = max(ident_max, na.rm = T)), 
  #      by = .(sorted_pair, orthology_type, subtype, family_size, fam_id)]
  #    }
  
  ensembl_paralogs[, c("gene1", "gene2") := tstrsplit(sorted_pair, "_")]  
  
  message(paste0(ensembl_paralogs[, .N], " ENSEMBL paralog pairs imported."))
  
    if (only_pairs) {return(ensembl_paralogs$sorted_pair)
  } else {
    return(ensembl_paralogs)}
}
