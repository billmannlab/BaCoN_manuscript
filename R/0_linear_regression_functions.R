generate_model_data <- \(.chronos_gene, .expression_gene, .cl_set, .tag) {
  
  .x <- copy(covariates[.cl_set, on = "cell_line"])[
    , `:=`(g1_effect = chronos_complete[.cl_set,.chronos_gene], 
           g2_expression = gene_expression_complete[.cl_set,.expression_gene])]
  
.x}

get_lm_predictor_pval <- \(.model) {
  .model$coefficients[["g2_expression","Pr(>|t|)"]]}

get_lm_slope <- \(.model) {
  .model$coefficients[["g2_expression","Estimate"]]}



work_on_cache_file <- \(instruction_file, target_file) {
  
  instructions <- readRDS(instruction_file)
  
  .d <- expand_grid_dt(instructions$expression_genes, 
                       instructions$chronos_genes, "exp_gene", "chr_gene")
  .d[, to_compute := instructions$to_compute]
  .d[exp_gene %in% instructions$genes_to_skip, to_compute := F]
  
  .chronos_genes_oi <- .d[(to_compute), unique(get("chr_gene"))]
  pb <- progress_bar$new(
    format = "[:bar] :percent  (:current/:total, :tick_rate/s) ETA: :eta", 
    total = length(.chronos_genes_oi), 
    width = 75, force = T)
  
  for (.chr_goi in .chronos_genes_oi) {
    .i <- .d[, .I[(to_compute) & chr_gene == .chr_goi]]
    
    pb$tick()
    
    .models <- mapply(
      \(.cg, .eg) {
        lm(formula = instructions$lm_formula, 
           data = generate_model_data(.chronos_gene = .cg, 
                                      .expression_gene = .eg, 
                                      .cl_set = instructions$cell_lines, 
                                      .tag = instructions$lm_tag)) %>% summary()}, 
      .d[.i, get("chr_gene")], .d[.i, get("exp_gene")], SIMPLIFY = F)
    
    
    
    .d[.i, `:=`(predictor_pval = sapply(.models, get_lm_predictor_pval), 
                predictor_slope = sapply(.models, get_lm_slope))]
    
  }
  
  saveRDS(.d[, .(exp_gene, chr_gene, predictor_pval, predictor_slope)], target_file)}


build_lm_matrix <- \(rownames, colnames, 
                     files_needed, filepath, 
                     return_everything = F) {
  if (return_everything) {print("Returning everything.")}
  
  files_found <- gsub(".rds", "", basename(list.files(filepath, ".rds$", full.names = T)))
  files_available <- intersect(files_needed, files_found)
  
  if (all(files_needed %in% files_available)) {
    message("Cache complete!")
  } else {
    message(str_c(l(setdiff(files_needed, files_available)), 
                  " cache files are missing. Please compute them first."))
  }
  
  message("Generating list of arrays...")
  
  small_arrays <- list()
  pb <- progbar(l(files_available))
  
  for (.f in files_available) {
    .egenes <- chunk_definition[name == .f, unlist(exp_genes)]
    .cgenes <- chunk_definition[name == .f, unlist(chr_genes)]
    
    .data <- readRDS(file.path(filepath, str_c(.f, ".rds")))
    .check <- .data[, all(exp_gene == rep(.egenes, l(.cgenes))) & 
                      all(chr_gene == rep(.cgenes, each = l(.egenes)))]
    
    .pval_arr <- empty_array(list(.egenes, .cgenes))
    .slope_arr <- empty_array(list(.egenes, .cgenes))
    
    if (.check) {
      .pval_arr[] <- .data[, get("predictor_pval")]
      .slope_arr[] <- .data[, get("predictor_slope")]
    }
    if (!.check) {
      message(str_c("The file ", .f, " seems to be out of order."))
      
      for (.cg in .cgenes) {
        .pval_arr[,.cg] <- .data[chr_gene == .cg][
          .egenes, on = "exp_gene", get("predictor_pval")]
        .slope_arr[,.cg] <- .data[chr_gene == .cg][
          .egenes, on = "exp_gene", get("predictor_slope")]}}
    
    
    .output_arr <- -log10(.pval_arr)
    .i <- which(.slope_arr < 0)
    .output_arr[.i] <- -.output_arr[.i]
    
    small_arrays[[.f]] <- .output_arr
    pb$tick()
  }
  
  message("Fill complete array...")
  
  pb <- progbar(l(files_available))
  complete_array <- empty_array(list(rownames, colnames))
  
  for (.f in files_available) {
    complete_array[rownames(small_arrays[[.f]]),
                   colnames(small_arrays[[.f]])] <- small_arrays[[.f]]
    pb$tick()}
  
  
  
  if (return_everything) {
    return(list(small_arrays = small_arrays, 
                complete_array = complete_array))}
  if (!return_everything) {return(complete_array)}
}
