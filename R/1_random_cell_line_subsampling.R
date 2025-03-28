.fpath <- file.path(large_cache, str_c(.dm, "_random_cell_line_subsets"))
mkdir(.fpath)
.pred <- T

for (.ss in names(cl_subsets$random)) {
  #Sys.sleep(20)
  message(.ss)
  .subset <- cl_subsets$random[[.ss]]
  .n <- str_c(.dm, .ss, "default_CvE", sep = "_")
  
  if (F) {
    .pcc_bacon <- BaCoN_pipeline(
      expression_matrix = gene_expression[.subset,], 
      effect_matrix = chronos[.subset,], 
      cache_path = file.path(.fpath, 
                             "PCC_and_PCC_BaCoN", 
                             str_c(.n, "PCC_BaCoN_0.05", sep = "_")))
    
    if (.pred) {
      .pcc_reference <- .pcc_bacon$correlation_matrix
      
      evaluation_routine(
        c(.pcc_bacon, list(ref = .pcc_reference)), 
        instructions = list(list(mat = "correlation_matrix"), 
                            list(mat = "bacon_matrix", ref = "ref")) %>% 
          setNames(str_c(.n, "_PCC", c("", "_BaCoN_0.05"))))}}
  
  # ComBat
  if (F) {
    for (.m in combat_methods) {
      .pcc_bacon <- BaCoN_pipeline(
        expression_matrix = combat_correction(gene_expression[.subset,], .m), 
        effect_matrix = combat_correction(chronos[.subset,], .m), 
        cache_path = file.path(
          .fpath, 
          str_c("ComBat", .m, "PCC_and_PCC_BaCoN", sep = "_"), 
          str_c(.n, "PCC_BaCoN_0.05_ComBat", "exp", .m, "chr", .m, sep = "_")))
      
      if (.pred) {
        evaluation_routine(
          c(.pcc_bacon, list(ref = .pcc_reference)), 
          instructions = list(list(mat = "correlation_matrix", ref = "ref"), 
                              list(mat = "bacon_matrix", ref = "ref")) %>% 
            setNames(str_c(.n, "_PCC_", c("", "BaCoN_0.05_"), "ComBat_exp_", .m, "_chr_", .m)))
      }}}
  
  # Cholesky whitening
  
  if (F) {
    .m <- "Cholesky"
    
    .pcc_bacon <- BaCoN_pipeline(
      expression_matrix = t(whitening::whiten(t(gene_expression[.subset,]), method = .m)), 
      effect_matrix = t(whitening::whiten(t(chronos[.subset,]), method = .m)), 
      cache_path = file.path(.fpath, 
                             str_c("whitening_", .m, "_PCC_and_PCC_BaCoN"), 
                             str_c(.n, "PCC_BaCoN_0.05_whitening", .m, sep = "_")))
    if (.pred) {
      
      evaluation_routine(
        c(.pcc_bacon, list(ref = .pcc_reference)), 
        instructions = list(list(mat = "correlation_matrix", ref = "ref"), 
                            list(mat = "bacon_matrix", ref = "ref")) %>% 
          setNames(str_c(.n, "_PCC_", c("", "BaCoN_0.05_"), "whitening_", .m)))
    }}
  
  ## ---- linear regression finalization random subsets ----
  
  if (F) { # requires the linreg cache files with random subsets to be computed. 
    linreg_array <- build_lm_matrix(
      rownames = expression_genes, 
      colnames = chronos_genes, 
      files_needed = chunk_definition[set == "default", name], 
      filepath = file.path(linreg_filepath, 
                           str_c("random_cl_subsets_", "eff1_exp2.lineage"), .ss)) %>% 
      cacheR(str_c(.dm, "default_CvE_linreg_matrix", "eff1_exp2.lineage", .ss, sep = "_"), 
             file.path(linreg_filepath, str_c("random_cl_subsets_", "eff1_exp2.lineage")))
    
    if (.pred) {
      list(linreg_matrix = linreg_array, ref = CvE_PCC) %>% 
        evaluation_routine(instructions = list(list(mat = "linreg_matrix", ref = "ref"#, add = "sample_scores"
        )) %>% 
          setNames(str_c(.dm, .ss, "default_CvE_linreg", "eff1_exp2.lineage", sep = "_")))
    }}}
