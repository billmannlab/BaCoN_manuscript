
## ---- reviewer response other datatypes computation ----

.fpath <- file.path(reviewer_response_fpath, str_c(.dm, "_other_datatypes"))
.pred <- F

# Chronos ~ Methylation:

.pcc_bacon <- BaCoN_pipeline(
  expression_matrix = methylation_complete, 
  effect_matrix = chronos, 
  cache_path = file.path(.fpath, str_c(.dm, "_test_CvMeth_PCC_BaCoN_0.05")))

if (.pred) {
  evaluation_routine(
    .pcc_bacon, 
    instructions = list(list(mat = "correlation_matrix"), 
                        list(mat = "bacon_matrix", 
                             ref = "correlation_matrix")) %>% 
      setNames(str_c(.dm, c("_test_CvMeth_PCC", 
                            "_test_CvMeth_PCC_BaCoN_0.05"))))}

# Demeter ~ Protein expression:

.pcc_bacon <- BaCoN_pipeline(
  expression_matrix = gygi_pe, 
  effect_matrix = demeter, 
  cache_path = file.path(.fpath, str_c(.dm, "_test_DvPE_ms_PCC_BaCoN_0.05")))

if (.pred) {
  
  evaluation_routine(
    .pcc_bacon, 
    instructions = list(list(mat = "correlation_matrix"), 
                        list(mat = "bacon_matrix", ref = "correlation_matrix")) %>% 
      setNames(str_c(.dm, c("_test_DvPE_ms_PCC", 
                            "_test_DvPE_ms_PCC_BaCoN_0.05"))))}


# Demeter ~ Gene expression:

.pcc_bacon <- BaCoN_pipeline(
  expression_matrix = gene_expression, 
  effect_matrix = demeter, 
  cache_path = file.path(.fpath, str_c(.dm, "_test_DvE_PCC_BaCoN_0.05")))

if (.pred) {
  evaluation_routine(
    .pcc_bacon, 
    instructions = list(list(mat = "correlation_matrix"), 
                        list(mat = "bacon_matrix", 
                             ref = "correlation_matrix")) %>% 
      setNames(str_c(.dm, c("_test_DvE_PCC", 
                            "_test_DvE_PCC_BaCoN_0.05"))))}


# Chronos ~ Protein expression:

.pcc_bacon <- BaCoN_pipeline(
  expression_matrix = gygi_pe, 
  effect_matrix = chronos_complete, 
  cache_path = file.path(.fpath, str_c(.dm, "_test_CvPE_ms_PCC_BaCoN_0.05")))
if (.pred) {
  evaluation_routine(
    .pcc_bacon, 
    instructions = list(list(mat = "correlation_matrix"), 
                        list(mat = "bacon_matrix", ref = "correlation_matrix")) %>% 
      setNames(str_c(.dm, c("_test_CvPE_ms_PCC", 
                            "_test_CvPE_ms_PCC_BaCoN_0.05"))))}

# Demeter ~ Protein expression:

.pcc_bacon <- BaCoN_pipeline(
  expression_matrix = gygi_pe, 
  effect_matrix = demeter, 
  cache_path = file.path(.fpath, str_c(.dm, "_test_DvPE_ms_PCC_BaCoN_0.05")))
if (.pred) {
  
  evaluation_routine(
    .pcc_bacon, 
    instructions = list(list(mat = "correlation_matrix"), 
                        list(mat = "bacon_matrix", ref = "correlation_matrix")) %>% 
      setNames(str_c(.dm, c("_test_DvPE_ms_PCC", 
                            "_test_DvPE_ms_PCC_BaCoN_0.05"))))
  
}


# Chronos ~ CNV:

.pcc_bacon <- BaCoN_pipeline(
  expression_matrix = copy_numbers_complete, 
  effect_matrix = chronos_complete, 
  cache_path = file.path(.fpath, str_c(.dm, "_test_CvCNV_PCC_BaCoN_0.05")))
if (.pred) {
  
  evaluation_routine(
    .pcc_bacon, 
    instructions = list(list(mat = "correlation_matrix"), 
                        list(mat = "bacon_matrix", ref = "correlation_matrix")) %>% 
      setNames(str_c(.dm, c("_test_CvCNV_PCC", 
                            "_test_CvCNV_PCC_BaCoN_0.05"))))
  
}

