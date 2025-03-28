## ---- reviewer response uncorrected scores computation ----

.fpath <- file.path(reviewer_response_fpath, str_c(.dm, "_uncorrected_chronos_scores"))
.pred <- T


# uncorrected Chronos scores:

.pcc_bacon <- BaCoN_pipeline(
  expression_matrix = gene_expression, 
  effect_matrix = chronos_uncorrected, 
  cache_path = file.path(.fpath, str_c(.dm, "_default_CuvE_PCC_BaCoN_0.05")))

if (.pred) {
  evaluation_routine(.pcc_bacon, 
                     instructions = list(list(mat = "correlation_matrix"), 
                                         list(mat = "bacon_matrix", ref = "correlation_matrix")) %>% 
                       setNames(str_c(.dm, c("_default_CuvE_PCC", 
                                             "_default_CuvE_PCC_BaCoN_0.05"
                       ))))}

CuvE_PCC <- .pcc_bacon$correlation_matrix

## with whitening

for (.m in methods_whiten) {
  message(.m)
  
  .pcc_bacon <- BaCoN_pipeline(
    expression_matrix = t(whitening::whiten(t(gene_expression), method = .m)), 
    effect_matrix = t(whitening::whiten(t(chronos_uncorrected), method = .m)), 
    cache_path = file.path(.fpath, str_c(.dm, "_default_CuvE_PCC_BaCoN_0.05_whitening_", .m)))
  if (.pred) {
    
    evaluation_routine(
      c(.pcc_bacon, list(ref = CuvE_PCC)), 
      instructions = list(list(mat = "correlation_matrix", ref = "ref"), 
                          list(mat = "bacon_matrix", ref = "ref")) %>% 
        setNames(str_c(.dm, c("default_CuvE_PCC", "default_CuvE_PCC_BaCoN_0.05"), "whitening", .m, sep = "_")))}
  }

## with ComBat

for (.m in combat_methods) {
  message(.m)

  .pcc_bacon <- BaCoN_pipeline(expression_matrix = combat_correction(gene_expression, .m), 
                               effect_matrix = combat_correction(chronos_uncorrected, .m), 
                               cache_path = file.path(.fpath, str_c(.dm, "_default_CuvE_PCC_BaCoN_0.05_ComBat_", .m)))
  
  if (.pred) {
    evaluation_routine(
      c(.pcc_bacon, list(ref = CuvE_PCC)), 
      instructions = list(list(mat = "correlation_matrix", ref = "ref"), 
                          list(mat = "bacon_matrix", ref = "ref")) %>% 
        setNames(str_c(.dm, c("default_CuvE_PCC", "default_CuvE_PCC_BaCoN_0.05"), "ComBat_exp", .m, "chr", .m, sep = "_")))}
  }
