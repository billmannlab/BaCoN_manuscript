## ---- functions ----

## ---- combat_correction function ----


combat_correction <- \(matrix, method, 
                       cell_line_covariates = covariates, 
                       cell_line_column = "cell_line", 
                       gene_covariates = gene_database, 
                       gene_column = "symbol") {
  
  if (nrow(matrix) >= ncol(matrix)) {
    warning("The matrix has to have cell lines on rows and genes as columns!")}
  
  if (method %in% c("growthpattern", "lineage")) {
    if (!method %in% names(cell_line_covariates)) {
      warning(str_c(method, " not found in the covariate file!"))}
    
    if (!all(rownames(matrix) %in% cell_line_covariates[, get(cell_line_column)])) {
      warning(str_c("Not all cell lines in covariate file!"))}
    
    .batch <- cell_line_covariates[rownames(matrix), on = cell_line_column, 
                                   get(method)]
    
    return(t(sva::ComBat(t(matrix), batch = .batch)))
  }
  
  if (method %in% c("chromosome", "chromosome_arm")) {
    .batch <- distinct(gene_covariates[, .SD, .SDcols = c(gene_column, method)])
    .batch <- .batch[colnames(matrix), on = gene_column]
    setnames(.batch, method, "batch")
    .batch[is.na(batch), batch := "Unknown"]
    
    
    if (.batch[, any(duplicated(get(gene_column)))]) {
      warning(str_c("Duplicated genes in batch information!"))}
    
    .batch <- .batch[, get("batch")]
    
    sva::ComBat(matrix, batch = .batch)}
}

## ---- collect_predictions function ----

collect_predictions <- \(matrix, 
                         reference_matrix, 
                         posneg = "positive", 
                         set = "extreme", 
                         pairs = 5000) {
  
  message(paste("Predicting", pairs, "pairs based on", 
                set, posneg, "scores...", sep = " "))
  reference_given <- !missing(reference_matrix)
  
  if (!reference_given) {
    message(paste0("No reference matrix provided! ", 
                   "No reference PCC values will be collected."))}
  
  if (posneg == "positive") {
    th <- sort(matrix[matrix > 0], 
               decreasing = c("extreme" = T, "moderate" = F)[[set]])[pairs]
    i <- which(matrix >= th); i_arr_ind <- which(matrix >= th, arr.ind = T)
  }
  
  if (posneg == "negative") {
    th <- sort(matrix[matrix < 0], 
               decreasing = c("extreme" = F, "moderate" = T)[[set]])[pairs]
    i <- which(matrix <= th); i_arr_ind <- which(matrix <= th, arr.ind = T)
  }
  
  
  x <- data.table(i_arr_ind, score = matrix[i])
  if (reference_given) {x[, pcc_reference := reference_matrix[i]]}
  
  x[, `:=`(expression_gene = rownames(matrix)[row], 
           effect_gene = colnames(matrix)[col])]
  
  
  
  x[, `:=`(sorted_pair = sort_gene_pairs(expression_gene, effect_gene))]
  
  x <- x[, !c("row", "col")]
  
  .dec <- c("positive" = T, "negative" = F)[[posneg]]
  
  if (!reference_given) {
    x <- x[order(score, decreasing = .dec)]
  }
  if (reference_given) {
    x <- x[order(score, pcc_reference, decreasing = .dec)]
  }
  return(x)
}

## ---- evaluation_routine function ----


evaluation_routine <- \(matrices, instructions, out_dir = result_dir) {
  
  if (length(instructions) == 0) {
    warning("No instructions provided.")
  }
  
  
  for (.d in file.path(out_dir, c("predictions", 
                                  "matrix_stats", 
                                  #"spearman", 
                                  "sample_scores"))) {mkdir(.d)}
  for (id in names(instructions)) {
    message(id)
    
    instr <- instructions[[id]]
    
    n_pairs <- ifelse("pairs" %in% names(instr), instr$pairs, 5000)
    
    ref_given <- "ref" %in% names(instr)
    
    if (!ref_given) {
      collect_predictions(matrix = matrices[[instr$mat]], pairs = n_pairs) %>% 
        cacheR(str_c(id, "_predictions"), file.path(out_dir, "predictions"))
      
#      if (any(c("spearman") %in% instr$add)) {
#        warning("Instructions require reference matrix, but none is provided.")
#    }
    }
    
    #if ("matrix_stats" %in% instr$add) {
    message("Collecting matrix stats...")
    {.stats <- list(
      nrow = nrow(matrices[[instr$mat]]), 
      ncol = ncol(matrices[[instr$mat]]), 
      matrix_na_perc = na_perc(matrices[[instr$mat]]), 
      matrix_mean = mean(matrices[[instr$mat]], na.rm = T), 
      matrix_sd = sd(matrices[[instr$mat]], na.rm = T))
      
      if (ref_given) {
        .stats <- c(.stats, 
                    list(
                      ref_nrow = nrow(matrices[[instr$ref]]), 
                      ref_ncol = ncol(matrices[[instr$ref]]), 
                      ref_na_perc = na_perc(matrices[[instr$ref]]), 
                      ref_mean = mean(matrices[[instr$ref]], na.rm = T), 
                      ref_sd = sd(matrices[[instr$ref]], na.rm = T)))}
      .stats} %>% 
      cacheR(str_c(id, "_matrix_stats"), file.path(out_dir, "matrix_stats"))

    if (ref_given) {
      collect_predictions(matrix = matrices[[instr$mat]], 
                          reference_matrix = matrices[[instr$ref]], 
                          pairs = n_pairs) %>% 
        cacheR(str_c(id, "_predictions"), file.path(out_dir, "predictions"))
    }
    
    if ("sample_scores" %in% instr$add) {
      message("Collecting sample scores...")
      set.seed(692097)
      ss <- sample(1:length(matrices[[instr$mat]]), 10000000)
      .scores <- list(matrix = matrices[[instr$mat]][ss])
      if (ref_given) {
        .scores$reference <- matrices[[instr$ref]][ss]
      }
      .scores  %>% 
        cacheR(str_c(id, "_sample_scores"), file.path(out_dir, "sample_scores"))
    }}}


