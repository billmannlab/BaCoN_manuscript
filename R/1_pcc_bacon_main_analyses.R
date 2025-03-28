## ---- PCC and BaCoN ---

# toggle to switch on after the matrices are computed. 
# Initiates the collection of top predictions
.pred <- F

.pcc_bacon <- BaCoN_pipeline(
  expression_matrix = gene_expression, 
  effect_matrix = chronos, 
  cache_path = file.path(cache, 
                         str_c(.dm, "_default_CvE_PCC_BaCoN_0.05")))


if (.pred) {
  evaluation_routine(
    .pcc_bacon, 
    instructions = list(
      list(mat = "correlation_matrix", ref = "correlation_matrix", pairs = 10000), 
      list(mat = "bacon_matrix", ref = "correlation_matrix", 
           add = "sample_scores", pairs = 10000)) %>% 
      setNames(str_c(.dm, c("_default_CvE_PCC", 
                            "_default_CvE_PCC_BaCoN_0.05"))))
}


CvE_PCC <- .pcc_bacon$correlation_matrix
CvE_BaCoN <- .pcc_bacon$bacon_matrix

## ---- combat ----

if (F) {
  for (.m in combat_methods) {
    message(str_c("exp", .m, "chr", .m, sep = "_"))
    
    .pcc_bacon <- BaCoN_pipeline(
      expression_matrix = combat_correction(gene_expression, .m), 
      effect_matrix = combat_correction(chronos, .m), 
      cache_path = file.path(
        large_cache, 
        str_c(.dm, "_ComBat"), 
        str_c(.dm, "default_CvE_PCC_BaCoN_0.05_ComBat_exp", .m, "chr", .m, sep = "_")))
    
    if (.pred) {
      evaluation_routine(
        c(.pcc_bacon, list(ref = CvE_PCC)), 
        instructions = list(list(mat = "correlation_matrix", 
                                 ref = "ref", pairs = 10000), 
                            list(mat = "bacon_matrix", ref = "ref", 
                                 add = "sample_scores", pairs = 10000)) %>% 
          setNames(str_c(.dm, c("default_CvE_PCC", "default_CvE_PCC_BaCoN_0.05"), 
                         "ComBat_exp", .m, "chr", .m, sep = "_")))
    }}}

## ---- whitening ----

if (T) {
  for (.m in methods_whiten) {
    message(.m)
    .pcc_bacon <- BaCoN_pipeline(
      expression_matrix = t(whitening::whiten(t(gene_expression), method = .m)), 
      effect_matrix = t(whitening::whiten(t(chronos), method = .m)), 
      cache_path = file.path(large_cache, 
                             str_c(.dm, "_whitening"), 
                             str_c(.dm, "_default_CvE_PCC_BaCoN_0.05_whitening_", .m)))
    
    if (.pred) {
      
      evaluation_routine(
        c(.pcc_bacon, list(ref = CvE_PCC)), 
        instructions = list(
          list(mat = "correlation_matrix", ref = "ref", 
               add = list("sample_scores"), pairs = 10000), 
          list(mat = "bacon_matrix", ref = "ref", 
               add = "sample_scores", pairs = 10000)) %>% 
          setNames(str_c(.dm, c("default_CvE_PCC", "default_CvE_PCC_BaCoN_0.05"), 
                         "whitening", .m, sep = "_")))
    }}}

## ---- Spearman correlation ----

if (T) {
  .pcc_bacon <- BaCoN_pipeline(
    expression_matrix = gene_expression, 
    effect_matrix = chronos, 
    cor_method = "spearman", 
    cache_path = file.path(large_cache, str_c(.dm, "_default_CvE_SCC_BaCoN_0.05")))
  if (.pred) {
    
    evaluation_routine(
      c(.pcc_bacon, list(ref = CvE_PCC)), 
      instructions = list(list(mat = "correlation_matrix", ref = "ref", 
                               add = "sample_scores", pairs = 10000), 
                          list(mat = "bacon_matrix", ref = "ref", 
                               add = "sample_scores", pairs = 10000)) %>% 
        setNames(str_c(.dm, c("_default_CvE_SCC", "_default_CvE_SCC_BaCoN_0.05"))))
  }}

## ---- Other correction factors ----

if (F & .dm == "23Q2") {
  for (.i in c(0, 0.001, 0.005, 0.01, 0.025, 0.075, 0.1, 0.125, 0.15)) {
    .pcc_bacon <- BaCoN_pipeline(
      expression_matrix = gene_expression, 
      effect_matrix = chronos, 
      bacon_correction_factor = .i, 
      cache_path = file.path(reviewer_response_fpath, 
                             str_c(.dm, "_other_correction_factors"), 
                             str_c(.dm, "_default_CvE_PCC_BaCoN_", .i)))
    if (.pred) {
      evaluation_routine(
        .pcc_bacon, 
        instructions = list(
          list(mat = "bacon_matrix", ref = "correlation_matrix")) %>% 
          setNames(str_c(.dm, "default_CvE_PCC_BaCoN", .i, sep = "_")))
    }
  }
}
