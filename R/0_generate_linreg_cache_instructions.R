
# generate cache instructions

for (.tag in linreg_models[, tag]) {
  .fpath <- file.path(lm_instructions_fpath, .tag); mkdir(.fpath)
  
  for (.cn in chunk_definition[set == "default", name]) {
    
    .fname <- str_c("lm_cache_instructions", .tag, .cn, sep = "_")
    
    .cl <- cl_subset_default
    .to_skip <- unexpressed_genes
    
    list(
      chunk_name = .cn, 
      genes_to_skip = .to_skip, 
      cell_lines = .cl, 
      to_compute = T, 
      lm_tag = .tag, 
      lm_formula = linreg_models[tag == .tag, formula], 
      expression_genes = chunk_definition[name == .cn, unlist(exp_genes)], 
      chronos_genes = chunk_definition[name == .cn, unlist(chr_genes)]) %>%
      cacheR(.fname, .fpath, F)
  }
  
  # for lineage, we need the entire gene space
  
  if (.tag == "eff1_exp2.lineage") {
    for (.cn in chunk_definition[, name]) {
      
      .fname <- str_c("lm_cache_instructions", .tag, .cn, sep = "_")
      list(
        chunk_name = .cn, 
        genes_to_skip = unexpressed_genes, 
        cell_lines = cl_subset_default, 
        to_compute = T, 
        lm_tag = .tag, 
        lm_formula = linreg_models[tag == .tag, formula], 
        expression_genes = chunk_definition[name == .cn, unlist(exp_genes)], 
        chronos_genes = chunk_definition[name == .cn, unlist(chr_genes)]) %>%
        cacheR(.fname, .fpath, F)
    }}}

# generate cache instructions for random cell line subsets:

.tag <- "eff1_exp2.lineage"

.pb <- progbar(l(cl_subsets$random))

for (.n in names(cl_subsets$random)) {
  .pb$tick()
  .fpath <- file.path(lm_instructions_fpath, 
                      str_c("random_cl_subsets", .tag, sep = "_"), .n); mkdir(.fpath)
  
  for (.cn in chunk_definition[set == "default", name]) {
    
    .fname <- str_c("lm_cache_instructions", .tag, .n, .cn, sep = "_")
    
    .exp_genes <- chunk_definition[name == .cn, unlist(exp_genes)]
    .chr_genes <- chunk_definition[name == .cn, unlist(chr_genes)]
    .subset <- cl_subsets$random[[.n]]
    
    if (!file.exists(file.path(.fpath, str_c(.fname, ".rds")))) {
      .pcc_m <- cor(gene_expression_complete[.subset,.exp_genes], 
                    chronos_complete[.subset,.chr_genes])
      
      list(
        chunk_name = .cn, 
        genes_to_skip = unexpressed_genes, 
        cell_lines = cl_subsets$random[[.n]], 
        to_compute = as.vector(.pcc_m >= mean(.pcc_m, na.rm = T) + 1 * sd(.pcc_m, na.rm = T)), 
        lm_tag = .tag, 
        lm_formula = linreg_models[tag == .tag, formula], 
        expression_genes = .exp_genes, 
        chronos_genes = .chr_genes) %>%
        cacheR(.fname, .fpath, F)
    }}}
