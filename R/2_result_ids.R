## ---- result IDs ----

for (.dmv in c("22Q4", "23Q2", "24Q2")) {
  
  .n <- c(
    "default_CvE_PCC" = "PCC", 
    "default_CvE_PCC_BaCoN_0.05" = "PCC + BaCoN", 
    "default_CvE_SCC" = "SCC", 
    "default_CvE_SCC_BaCoN_0.05" = "SCC + BaCoN", 
    "default_CvE_PCC_ComBat_exp_lineage_chr_lineage" = "ComBat tissue", 
    "default_CvE_PCC_ComBat_exp_growthpattern_chr_growthpattern" = "ComBat growth", 
    "default_CvE_PCC_ComBat_exp_chromosome_chr_chromosome" = "ComBat chromosome", 
    "default_CvE_PCC_ComBat_exp_chromosome_arm_chr_chromosome_arm" = "ComBat chr arm", 
    "default_CvE_PCC_whitening_Cholesky" = "Cholesky white.", 
    "default_CvE_PCC_whitening_PCA" = "PCA white.",
    "default_CvE_linreg_eff1_exp2" = "Linreg", 
    "default_CvE_linreg_eff1_exp2.lineage" = "Linreg + tissue", 
    "default_CvE_linreg_eff1_exp2.growthpattern" = "Linreg + growth")
  
  results[[str_c("IDs_", .dmv)]]$benchmarking <- data.table(
    id = str_c(.dmv, "_", names(.n)), 
    name = factor(.n, levels = rev(.n)))
  
  .n <- c(
    "default_CvE_PCC" = "PCC", 
    "default_CvE_PCC_BaCoN_0.05" = "PCC + BaCoN", 
    "default_CvE_SCC" = "SCC", 
    "default_CvE_SCC_BaCoN_0.05" = "SCC + BaCoN", 
    "default_CvE_PCC_ComBat_exp_lineage_chr_lineage" = "ComBat tissue", 
    "default_CvE_PCC_BaCoN_0.05_ComBat_exp_lineage_chr_lineage" = 
      "ComBat tissue + BaCoN", 
    "default_CvE_PCC_ComBat_exp_growthpattern_chr_growthpattern" = "ComBat growth", 
    "default_CvE_PCC_BaCoN_0.05_ComBat_exp_growthpattern_chr_growthpattern" = 
      "ComBat growth + BaCoN", 
    "default_CvE_PCC_ComBat_exp_chromosome_chr_chromosome" = "ComBat chromosome", 
    "default_CvE_PCC_BaCoN_0.05_ComBat_exp_chromosome_chr_chromosome" = 
      "ComBat chromosome + BaCoN", 
    "default_CvE_PCC_ComBat_exp_chromosome_arm_chr_chromosome_arm" = "ComBat chr arm", 
    "default_CvE_PCC_BaCoN_0.05_ComBat_exp_chromosome_arm_chr_chromosome_arm" = 
      "ComBat chr arm + BaCoN", 
    "default_CvE_PCC_whitening_Cholesky" = "Cholesky white.", 
    "default_CvE_PCC_BaCoN_0.05_whitening_Cholesky" = "Cholesky white. + BaCoN")
  
  results[[str_c("IDs_", .dmv)]]$benchmarking_with_bacon <- data.table(
    id = str_c(.dmv, "_", names(.n)), 
    name = factor(.n, levels = rev(unique(.n))))
  
  results[[str_c("IDs_", .dmv)]]$hq_set <- data.table(
    id = str_c(.dmv, "_default_CvE_PCC_BaCoN_0.05_whitening_Cholesky"), 
    name = "Cholesky whitening + BaCoN")
  
  .n <- c(
    "default_CvE_PCC" = "PCC", 
    "default_CvE_PCC_BaCoN_0.05" = "PCC", 
    "default_CuvE_PCC" = "PCC", 
    "default_CuvE_PCC_BaCoN_0.05" = "PCC", 
    "default_CvE_PCC_ComBat_exp_lineage_chr_lineage" = "ComBat tissue", 
    "default_CvE_PCC_BaCoN_0.05_ComBat_exp_lineage_chr_lineage" = "ComBat tissue", 
    "default_CvE_PCC_ComBat_exp_chromosome_chr_chromosome" = "ComBat chromosome", 
    "default_CvE_PCC_BaCoN_0.05_ComBat_exp_chromosome_chr_chromosome" = 
      "ComBat chromosome", 
    "default_CvE_PCC_ComBat_exp_chromosome_arm_chr_chromosome_arm" = "ComBat chr arm", 
    "default_CvE_PCC_BaCoN_0.05_ComBat_exp_chromosome_arm_chr_chromosome_arm" = 
      "ComBat chr arm", 
    "default_CuvE_PCC_ComBat_exp_lineage_chr_lineage" = "ComBat tissue", 
    "default_CuvE_PCC_BaCoN_0.05_ComBat_exp_lineage_chr_lineage" = "ComBat tissue", 
    "default_CuvE_PCC_ComBat_exp_chromosome_chr_chromosome" = "ComBat chromosome", 
    "default_CuvE_PCC_BaCoN_0.05_ComBat_exp_chromosome_chr_chromosome" = 
      "ComBat chromosome", 
    "default_CuvE_PCC_ComBat_exp_chromosome_arm_chr_chromosome_arm" = "ComBat chr arm", 
    "default_CuvE_PCC_BaCoN_0.05_ComBat_exp_chromosome_arm_chr_chromosome_arm" = 
      "ComBat chr arm", 
    "default_CvE_PCC_whitening_Cholesky" = "Cholesky white.", 
    "default_CvE_PCC_BaCoN_0.05_whitening_Cholesky" = "Cholesky white.", 
    "default_CuvE_PCC_whitening_Cholesky" = "Cholesky white.", 
    "default_CuvE_PCC_BaCoN_0.05_whitening_Cholesky" = "Cholesky white.")
  
  results[[str_c("IDs_", .dmv)]]$all_methods <- data.table(
    id = str_c(.dmv, "_", names(.n)), 
    name = factor(.n, levels = rev(unique(.n))))
}

.n <- c(
  "PCC" = "PCC", 
  "PCC_BaCoN_0.05" = "PCC + BaCoN", 
  "PCC_whitening_Cholesky" = "Cholesky whitening + PCC", 
  "PCC_BaCoN_0.05_whitening_Cholesky" = "Cholesky whitening + PCC + BaCoN", 
  "PCC_ComBat_exp_lineage_chr_lineage" = "ComBat (tissue) + PCC", 
  "PCC_BaCoN_0.05_ComBat_exp_lineage_chr_lineage" = 
    "ComBat (tissue) + PCC + BaCoN", 
  "PCC_ComBat_exp_growthpattern_chr_growthpattern" = "ComBat (growth) + PCC", 
  "PCC_BaCoN_0.05_ComBat_exp_growthpattern_chr_growthpattern" = 
    "ComBat (growth) + PCC + BaCoN", 
  "PCC_ComBat_exp_chromosome_chr_chromosome"  = "ComBat (chromosome) + PCC", 
  "PCC_BaCoN_0.05_ComBat_exp_chromosome_chr_chromosome" = 
    "ComBat (chromosome) + PCC + BaCoN", 
  "PCC_ComBat_exp_chromosome_arm_chr_chromosome_arm" = 
    "ComBat (chromosome arm) + PCC", 
  "PCC_BaCoN_0.05_ComBat_exp_chromosome_arm_chr_chromosome_arm" = 
    "ComBat (chromosome arm) + PCC + BaCoN", 
  "linreg_eff1_exp2.lineage" = "Linear regression + tissue")

results$IDs_23Q2$random_cell_line_subsets <- c(
  sapply(names(.n), \(.) {
    data.table(id = str_c("23Q2", names(cl_subsets$random), 
                          "default_CvE", ., sep = "_"), 
               name = .n[[.]])}, simplify = F), 
  list(data.table(id = str_c("23Q2", "default_CvE", names(.n), sep = "_"), 
                  name = .n))) %>% rbindlist()


results$IDs_23Q2$random_cell_line_subsets[
  , name := factor(name, levels = rev(unique(name)))]

.n <- c("PCC" = "PCC", 
        "linreg_eff1_exp2.lineage" = "linreg", 
        "PCC_BaCoN_0.05" = "PCC + BaCoN", 
        "PCC_BaCoN_0.05_ComBat_lineage" = "ComBat (lineage) + PCC + BaCoN", 
        "PCC_BaCoN_0.05_ComBat_chromosome" = "ComBat (chromosome) + PCC + BaCoN", 
        "PCC_BaCoN_0.05_ComBat_chromosome_arm" = 
          "ComBat (chromosome arm) + PCC + BaCoN", 
        "PCC_whitening_Cholesky" = "Cholesky whitening + PCC", 
        "PCC_BaCoN_0.05_whitening_Cholesky" = "Cholesky whitening + PCC + BaCoN")

results$IDs_23Q2$dynamic_gene_space <- sapply(
  names(.n), 
  \(.) {data.table(
    id = str_c("23Q2", 
               dynamic_range_combinations[, unique(id)], "CvE", ., sep = "_"), 
    name = .n[[.]])}, simplify = F) %>% 
  rbindlist()

results$IDs_23Q2$other_correction_factors <- sapply(
  c(0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15), 
  \(.cf) {data.table(
    id = str_c("23Q2_default_CvE_PCC_BaCoN_", .cf), 
    name = str_c("PCC + BaCoN (correction factor ", .cf, ")"))}, simplify = F) %>% 
  rbindlist()

.n <- c("default_CvE_PCC" = "Chronos ~ Expression", 
        "default_CvE_PCC_BaCoN_0.05" = "Chronos ~ Expression + BaCoN", 
        "test_CvMeth_PCC" = "Chronos ~ Methylation", 
        "test_CvMeth_PCC_BaCoN_0.05" = "Chronos ~ Methylation + BaCoN", 
        "test_DvPE_ms_PCC" = "Demeter ~ Protein expression (MS)", 
        "test_DvPE_ms_PCC_BaCoN_0.05" = "Demeter ~ Protein expression (MS) + BaCoN",
        "test_DvE_PCC" = "Demeter ~ Expression", 
        "test_DvE_PCC_BaCoN_0.05" = "Demeter ~ Expression + BaCoN", 
        "test_CvPE_ms_PCC" = "Chronos ~ Protein Expression (MS)", 
        "test_CvPE_ms_PCC_BaCoN_0.05" = "Chronos ~ Protein Expression (MS) + BaCoN", 
        "test_CvCNV_PCC" = "Chronos ~ CNV", 
        "test_CvCNV_PCC_BaCoN_0.05" = "Chronos ~ CNV + BaCoN")

results$IDs_24Q2$other_datatypes <- data.table(
  id = str_c("24Q2", "_", names(.n)), 
  name = factor(.n, levels = rev(unique(.n))))


.n <- c("test_CvEnc_exp_1_PCC_BaCoN_0.05_whiten.Cholesky" = 
          "Cholesky whitening + PCC + BaCoN")

results$IDs_24Q2$ncRNAs <- data.table(
  id = str_c("24Q2", "_", names(.n)), 
  name = factor(.n, levels = rev(unique(.n))))

