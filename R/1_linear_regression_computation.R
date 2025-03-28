## ---- linear regression computation ----

.pred <- F

if (F) {
  for (.tag in linreg_models[, tag]) {
    
    .fpath <- file.path(linreg_filepath, .tag)
    mkdir(.fpath)
    for (n_jobs in 1) {
      job({
        .shuffled_chunks <- chunk_definition[set == "default"][sample(1:.N)]
        
        
        if (.tag == "eff1_exp2.lineage") {
          # For Linear regression + lineage: Compute raw matrix for dynamic 
          # gene space. Required for dynamic gene space heatmaps
          .shuffled_chunks <- chunk_definition[sample(1:.N)]
        }
        
        for (.cn in .shuffled_chunks[, name]) {
          message(.cn)
          instruction_file <- file.path(
            lm_instructions_fpath, .tag, 
            str_c("lm_cache_instructions_",.tag, "_", .cn, ".rds"))
          
          target_file <- file.path(.fpath, str_c(.cn, ".rds"))
          
          if (!file.exists(target_file)) {
            work_on_cache_file(instruction_file = instruction_file, 
                               target_file = target_file)
          }}
        job::export("none")}, 
        packages = c("data.table", "progress", "job", "tidyverse"), 
        import = c(.fpath, lm_instructions_fpath, .n, 
                   expand_grid_dt, 
                   generate_model_data, 
                   get_lm_predictor_pval, 
                   get_lm_slope, 
                   work_on_cache_file, 
                   linreg_models, cl_subsets, 
                   covariates, 
                   chronos_complete, 
                   gene_expression_complete, 
                   loss_status_complete, 
                   copy_numbers_complete, 
                   methylation_complete, 
                   methylation_complete_cat, 
                   .tag, chunk_definition, unexpressed_genes), 
        title = str_c("linear models (", .tag, ")"))}
  }
}


## ---- linear regression finalization ----

## Build matrices using expression term p-values and slopes:

if (.pred) {
  for (.tag in linreg_models[, tag]) {
    message(.tag)
    linreg_array <- build_lm_matrix(
      rownames = expression_genes, 
      colnames = chronos_genes, 
      files_needed = chunk_definition[set == "default", name], 
      filepath = file.path(linreg_filepath, .tag)) %>% 
      cacheR(str_c(.dm, "_default_CvE_linreg_matrix_", .tag), linreg_filepath)
    
    list(linreg_matrix = linreg_array, ref = CvE_PCC) %>% 
      evaluation_routine(
        instructions = list(list(mat = "linreg_matrix", ref = "ref")) %>% 
          setNames(str_c(.dm, "_default_CvE_linreg_", .tag)))
  }}


## ---- linear regression computation random subsets ----

if (F) {
  for (.ss in names(cl_subsets$random)) {
    message(.ss)
    .fpath <- file.path(linreg_filepath, 
                        str_c("random_cl_subsets_", "eff1_exp2.lineage"), .n)
    mkdir(.fpath)
    
    for (n_jobs in 1) {
      job({
        .shuffled_chunks <- chunk_definition[set == "default"][sample(1:.N)]
        for (.cn in .shuffled_chunks[, name]) {
          message(.cn)
          instruction_file <- file.path(
            lm_instructions_fpath, 
            str_c("random_cl_subsets_", "eff1_exp2.lineage"), .ss, 
            str_c("lm_cache_instructions_", "eff1_exp2.lineage", "_", .ss, "_", 
                  .cn, ".rds"))
          
          target_file <- file.path(.fpath, str_c(.cn, ".rds"))
          
          if (!file.exists(target_file)) {
            work_on_cache_file(instruction_file = instruction_file, 
                               target_file = target_file)
          }
        }
        job::export("none")}, 
        packages = c("data.table", "progress", "job", "tidyverse"), 
        import = c(.fpath, lm_instructions_fpath, .ss, 
                   expand_grid_dt, 
                   generate_model_data, 
                   get_lm_predictor_pval, 
                   get_lm_slope, 
                   work_on_cache_file, 
                   linreg_models, cl_subsets, 
                   covariates, chronos_complete, gene_expression_complete, 
                   .tag, chunk_definition, unexpressed_genes), 
        title = str_c("random subset linear models (", "eff1_exp2.lineage", ")"))}
  }}

linreg_models[, .(tag, cache_files = sapply(
  tag, \(.tag) {l(list.files(file.path(linreg_filepath, .tag)))}))]
