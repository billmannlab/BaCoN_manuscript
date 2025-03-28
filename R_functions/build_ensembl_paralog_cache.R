#' @export build_ensembl_paralog_cache

build_ensembl_paralog_cache <- \(
  path, 
  ensembl_infos = c("external_gene_name", 
                    str_c("hsapiens_paralog_", 
                          c("associated_gene_name", 
                            "perc_id_r1", 
                            "perc_id", 
                            "subtype", 
                            "orthology_type")))) {
  
  if (!dir.exists(path)) {dir.create(path)}
  
  all_files <- paste0("ensembl_info_", ensembl_infos)
  
  complete <- sapply(all_files, \(.) file.exists(file.path(path, str_c(., ".rds"))))
  
  print(paste0(sum(complete, na.rm = T), " / ", length(all_files)))
  
  # each information type is queued independently, but together with 
  # "ensembl_gene_id" and "hspapens_paralog_ensembl_gene". This allows 
  # to clearly connect all collected information types to the respective pairs
  
  job::job({
    ensembl <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                                   dataset = "hsapiens_gene_ensembl", 
                                   host = "https://www.ensembl.org")
    
    for (.info in sample(ensembl_infos)) {
      message(str_c(.info, "..."))
      d <- biomaRt::getBM(
        useCache = F, 
        attributes = c("ensembl_gene_id", 
                       "hsapiens_paralog_ensembl_gene", .info), 
        filters = c("with_hsapiens_paralog", 
                    "transcript_biotype"), 
        values = list(T, "protein_coding"), 
        mart = ensembl) %>% 
        cacheR(str_c("ensembl_info_", .info), path)}}, 
    packages = c("biomaRt", "stringr"), 
    import = c(cacheR, ensembl_infos, path), 
    title = str_c("Biomart call"))
}




#build_ensembl_paralog_cache_chromosome_wise <- \(
#  path, 
#  ensembl_infos = c("external_gene_name", 
#                    str_c("hsapiens_paralog_", 
#                          c("associated_gene_name", 
#                            "perc_id_r1", 
#                            "perc_id", 
#                            "subtype", 
#                            "orthology_type")))) {
#  
#  if (!dir.exists(path)) {dir.create(path)}
#  
#  all_files <- apply(expand.grid("ensembl_info", ensembl_infos, "chr", chromosomes), 
#                     1, paste, collapse = "_")
#  
#  complete <- sapply(all_files, \(.) file.exists(fp(path, str_c(., ".rds"))))
#  
#  print(paste0(sum(complete, na.rm = T), " / ", length(all_files)))
#  
#  ensembl <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
#                                 dataset = "hsapiens_gene_ensembl", 
#                                 host = "https://www.ensembl.org")
#  
#  # each information type is queued independently, but together with 
#  # "ensembl_gene_id" and "hspaoens_paralog_ensembl_gene". This allows 
#  # to clearly connect all collected information types to the respective pairs
#  job::job({
#    for (.info in sample(ensembl_infos)) {
#      
#      for (chr in sample(chromosomes)) {
#        message(str_c(.info, ", chromosome ", chr, "..."))
#        d <- biomaRt::getBM(
#          useCache = F, 
#          attributes = c("ensembl_gene_id", 
#                         "hsapiens_paralog_ensembl_gene", .info), 
#          filters = c("with_hsapiens_paralog", 
#                      "transcript_biotype", "chromosome_name"), 
#          values = list(T, "protein_coding", chr), 
#          mart = ensembl) %>% 
#          cacheR(str_c("ensembl_info_", .info, "_chr_", chr), path)}}}, 
#    packages = c("biomaRt", "stringr"), 
#    import = c(cacheR, ensembl_infos, chromosomes, ensembl, path), 
#    title = str_c("Biomart call"))
#}

