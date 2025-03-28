## ---- dynamic gene range ----

CvE_PCC_complete <- cor(
  gene_expression_complete[cl_subset_default,], 
  chronos_complete[cl_subset_default,], use = pco) %>% 
  cacheR(str_c(.dm, "_complete_CvE_PCC"), large_cache)

.pred <- F

ids_to_do <- dynamic_range_combinations[, .I]

for (.i in ids_to_do) {
  .id <- dynamic_range_combinations[.i, id]
  message(.id)
  exp_goi <- dynamic_range_combinations[.i, unlist(exp_genes)]
  chr_goi <- dynamic_range_combinations[.i, unlist(chr_genes)]
  
# dynamic gene range PCC + BaCoN
  
  if (T) {
    .pcc_bacon <- BaCoN_pipeline(
      expression_matrix = gene_expression_complete[cl_subset_default,exp_goi], 
      effect_matrix = chronos_complete[cl_subset_default,chr_goi], 
      cache_path = file.path(large_cache, 
                             str_c(.dm, "_dynamic_gene_space_PCC_BaCoN"), 
                             str_c(.dm, .id, "CvE_PCC_BaCoN_0.05", sep = "_")))
    
    if (.pred) {
      evaluation_routine(
        c(.pcc_bacon, list(ref = CvE_PCC_complete[exp_goi,chr_goi])), 
        instructions = list(list(mat = "correlation_matrix", pairs = 1500), 
                            list(mat = "bacon_matrix", ref = "ref", pairs = 1500)) %>% 
          setNames(str_c(.dm, .id, "CvE", c("PCC", "PCC_BaCoN_0.05"), sep = "_")))}}
  
  # dynamic gene space Cholesky PCC + BaCoN
  
  if (T) {
    
    .m <- "Cholesky"
    
    .PCC_m <- CvE_PCC_complete[exp_goi,chr_goi]
    
    if (dynamic_range_combinations[.i, chr_th != 1.5]) {
      # chronos threshold 12 does not work for whitening -> matrix too small!
      
      .pcc_bacon <- BaCoN_pipeline(
        expression_matrix = t(whiten(t(gene_expression_complete[cl_subset_default,exp_goi]), method = .m)), 
        effect_matrix = t(whiten(t(chronos_complete[cl_subset_default,chr_goi]), method = .m)), 
        cache_path = file.path(large_cache, 
                               str_c(.dm, "_dynamic_gene_space_whitening_", .m, "_PCC_BaCoN"), 
                               str_c(.dm, "_", .id, "_CvE_PCC_BaCoN_0.05_whitening_", .m)))
      
      if (.pred) {
        evaluation_routine(
          c(.pcc_bacon, list(ref = CvE_PCC_complete[exp_goi,chr_goi])), 
          instructions = list(list(mat = "correlation_matrix", ref = "ref", pairs = 1500), 
                              list(mat = "bacon_matrix", ref = "ref", pairs = 1500)) %>% 
            setNames(str_c(.dm, .id, "CvE", c("PCC", 
              "PCC_BaCoN_0.05"), "whitening", .m, sep = "_")))}}}
  
  # dynamic gene space, ComBat
  
  if (T) {
    for (.m in c("lineage", "chromosome", "chromosome_arm")) {
      .pcc_bacon <- BaCoN_pipeline(
        expression_matrix = combat_correction(gene_expression_complete[cl_subset_default,exp_goi], method = .m), 
        effect_matrix = combat_correction(chronos_complete[cl_subset_default,chr_goi], method = .m), 
        cache_path = file.path(
          large_cache, 
          str_c(.dm, "_dynamic_gene_space_ComBat.", .m, "_PCC_BaCoN"), 
          str_c(.dm, "_", .id, "_CvE_PCC_BaCoN_0.05_ComBat.", .m)))
      
      if (.pred) {
        
        evaluation_routine(
          c(.pcc_bacon, list(ref = CvE_PCC_complete[exp_goi,chr_goi])), 
          instructions = list(list(mat = "correlation_matrix", ref = "ref", pairs = 1500), 
                              list(mat = "bacon_matrix", ref = "ref", pairs = 1500)) %>% 
            setNames(str_c(.dm, "_", .id, "_CvE_PCC", c("", "_BaCoN_0.05"), "_ComBat_", .m)))}
    }}}
  
# dynamic gene range linreg

# The chunks for the entire ("complete") gene space are used to merge a large matrix, 
# which is then subsetted along the different levels of chronos and expression 
# gene set filtering. 

if (T) {
  lm_matrix_lineage_complete <- build_lm_matrix(
    rownames = colnames(gene_expression_complete), 
    colnames = colnames(chronos_complete), 
    files_needed = chunk_definition[, name], 
    filepath = file.path(large_cache, 
                         str_c(.dm, "_linear_regression_new_and_clean"), 
                         "eff1_exp2.lineage"))
  
  for (.i in 1:dynamic_range_combinations[, .N]
  ) {
    .id <- dynamic_range_combinations[.i, id]
    message(.id)
    exp_goi <- dynamic_range_combinations[.i, unlist(exp_genes)]
    chr_goi <- dynamic_range_combinations[.i, unlist(chr_genes)]
    
    if (.pred) {
      evaluation_routine(
        list(linreg_matrix = lm_matrix_lineage_complete[exp_goi,chr_goi], 
             ref = CvE_PCC_complete[exp_goi,chr_goi]), 
        instructions = list(list(mat = "linreg_matrix", ref = "ref", pairs = 1500)) %>% 
          setNames(str_c(.dm, .id, "CvE_linreg", "eff1_exp2.lineage", sep = "_")))
      }}}

