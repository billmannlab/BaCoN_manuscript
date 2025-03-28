.tag <- linreg_models[2, tag]

.fpath <- file.path(cache, 
                    str_c(.dm, "_random_subsets_linear_regression_", .tag))
mkdir(.fpath)

# Start jobs

for (i in 1) {
  job ({
    .reordered_chunks <- chunk_definition[set == "default"][sample(1:.N)]
    .genes_to_skip <- genes_excluded_from_lm
    .formula_to_use <- linreg_models[tag == .tag, formula]
    for (.n in .reordered_chunks[, name]) {
      .fname <- str_c("lm_cache_")}})}

##### ----

lm_random_subsets <- data.table(tag = "eff1_exp2.lineage", 
                                subset = names(cl_subsets$random))
lm_random_subsets[, complete := sapply(
  names(cl_subsets$random), 
  \(.ss) {length(list.files(file.path(
    cache2, 
    "linear_regression_random_subsets_eff1_exp2.lineage_z_score_th_1", 
    str_c(.ss, "_z_score_th_1")))) == 480})]

tag <- "eff1_exp2.lineage"

for (.ss in names(cl_subsets$random)) {
  message(.ss)
  
  # Using the PCC matrix to only compute pairs above z-score 1:
  
  .pcc_m <- readRDS(file.path(cache2, "random_subsets", "PCC", 
                              str_c("23Q2_CvE_PCC_", .ss, ".rds")))
  start_linreg_job(
    .jd = generate_jobdata(
      .chunk_definition = chunk_definition[set == "full"][sample(1:.N)], 
      .cl_set = cl_subsets$random[[.ss]], .tag = tag), 
    .tag = tag, 
    .coords_to_compute = .pcc_m >= (mean(.pcc_m, na.rm = T) + 
                                      sd(.pcc_m, na.rm = T)), 
    .fpath = file.path(
      cache2, 
      "linear_regression_random_subsets_eff1_exp2.lineage_z_score_th_1",  
      str_c(.ss, "_z_score_th_1")))
  
  Sys.sleep(1)
}
