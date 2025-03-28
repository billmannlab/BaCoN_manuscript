
## ---- collect predictions ----

exp_gene_info <- copy(gene_database[, .(expression_gene = symbol, 
                                        exp_chr = chromosome, 
                                        exp_pos = position, 
                                        exp_chr_arm = chromosome_arm)])

eff_gene_info <- copy(gene_database[, .(effect_gene = symbol, 
                                        eff_chr = chromosome, 
                                        eff_pos = position, 
                                        eff_chr_arm = chromosome_arm)])

add_prediction_metadata <- \(
  predictions, 
  sorted_pair_colname = "sorted_pair", 
  sets = c("ensembl", "ohnologs", "anvar", "anvar_standard", 
           "dekegel", "dekegel_val", "gemini", "ito", "thompson")) {
  
  x <- copy(predictions)
  for (.ps in sets) {
    x[, c(.ps) := get(sorted_pair_colname) %in% paralog_pairs[[.ps]]$sorted_pair]
  }
  x}

pre_hoc_methods <- c(str_c("ComBat_", c("growthpattern", 
                                        "lineage", "chromosome", "chromosome_arm")), 
                     str_c("whitening_", methods_whiten))

for (.dmv in dm_versions_to_include) {
  
  message(.dmv)
  .fpaths <- sapply(c("predictions", "matrix_stats"), 
                    \(.) {file.path(cache, str_c(.dmv, "_results"), .)}, 
                    simplify = F)
  
  .ids <- gsub("_predictions.rds$", "", list.files(.fpaths$predictions))
  .info <- data.table(id = .ids)
  
  .checks <- data.table(
    id = .ids, 
    predictions = str_c(.ids, "_predictions.rds") %in% 
      list.files(.fpaths$predictions), 
    matrix_stats = str_c(.ids, "_matrix_stats.rds") %in% 
      list.files(.fpaths$matrix_stats))
  
  if (!.checks[, all(predictions & matrix_stats)]) {
    warning("No matrix stats files found for: ")
    warning(.checks[!(matrix_stats), str_c(id, collapse = ", ")])
  }
  
  # Read dataset characteristics out of ID string:
  for (. in c("22Q4", "23Q2", "24Q2")) {
    .info[grepl(., id), version := str_c("dm", .)]}
  
  for (. in c("default", dynamic_range_combinations[, id], 
              names(cl_subsets$random))) {
    .info[grepl(., id), setup := .]}
  
  for (. in c("CvE", "CuvE", "DEMvE", "CvCNV", "CvPE", 
              "DEMvPE_ms", "CvEnc", "CvMeth")) {
    .info[grepl(., id), datasets := .]}
  
  for (. in as.character(c(0, 0.001, 0.005, 0.01, 0.025, 0.05, 
                           0.075, 0.1, 0.125, 0.15))) {
    .info[grepl(str_c("BaCoN_", .), id), BaCoN_th := .]}
  
  for (. in pre_hoc_methods) {.info[grepl(., id), pre_hoc := .]}
  
  for (. in c("PCC", "SCC", "linreg")) {.info[grepl(., id), association := .]}
  
  .info[grepl("BaCoN", id), post_hoc := "BaCoN"]
  
  #.x$all_files <- lapply(.fpaths, list.files, full.names = T)
  
  .z_score_info <- list.files(.fpaths$matrix_stats) %>% 
    lapply(\(.) {data.table(data.frame(
      id = gsub("_matrix_stats.rds$", "", .), 
      readRDS(file.path(.fpaths$matrix_stats, .))))}) %>% 
    Reduce(bind_rows, .)
  
  .pb <- progbar(l(list.files(.fpaths$predictions)))
  
  .all_preds <- list.files(.fpaths$predictions) %>% 
    grep("dynamic", ., value = T, invert = T) %>% 
    lapply(\(.) {
      .pb$tick()
      x <- readRDS(file.path(.fpaths$predictions, .))
      if ("pcc_reference" %in% names(x)) {
        x <- x[order(score, pcc_reference, decreasing = T)]}
      else {x <- x[order(score, decreasing = T)]}
      x[, `:=`(id = gsub("_predictions.rds$", "", .), rank = 1:.N)]
      
      add_prediction_metadata(x)}) %>% Reduce(bind_rows, .)
  
  .everything <- .all_preds %>% 
    merge(.info, by = "id", all.x = T, sort = F) %>% 
    merge(.z_score_info, by = "id", all.x = T, sort = F) %>% 
    merge(exp_gene_info, by = "expression_gene", all.x = T, sort = F) %>% 
    merge(eff_gene_info, by = "effect_gene", all.x = T, sort = F)
  
  .everything[exp_chr_arm == eff_chr_arm, distance := abs(exp_pos - eff_pos)]
  .everything[, `:=`(
    self_association = expression_gene == effect_gene, 
    proximity = sorted_pair %in% important_pairs_complete$proximity, 
    proximity_class = fcase(
      expression_gene == effect_gene, "self-addiction", 
      exp_chr_arm == eff_chr_arm & distance <= 1e7, "same arm, <10mbp", 
      exp_chr_arm == eff_chr_arm, "same arm, >10mbp", 
      exp_chr == eff_chr & exp_chr_arm != eff_chr_arm, "different arm", 
      exp_chr != eff_chr, "different chromosomes", 
      default = "NA"))]
  
  .everything[, z_score := (pcc_reference - ref_mean) / ref_sd]
  .everything[, z_score_bin := factor(
    fcase(z_score < 0, "< 0", 
          z_score %between% c(0, 1), "0-1", 
          z_score %between% c(1, 2), "1-2", 
          z_score %between% c(2, 3), "2-3", 
          z_score > 3, "> 3"), 
    levels = rev(c("< 0", "0-1", "1-2", "2-3", "> 3")))]
  .everything
  results[[str_c("all_predictions_", .dmv)]] <- .everything
}
