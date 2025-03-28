## ---- generate EV1 set ----

.fpath <- file.path(large_cache, str_c(.dm, "_EV1")); mkdir(.fpath)

template <- {
  x <- data.table(expression_gene = rep(expression_genes, l(chronos_genes)), 
                  effect_gene = rep(chronos_genes, each = l(expression_genes)))
  
  x[, `:=`(ID = .I, sorted_pair = sort_gene_pairs(expression_gene, effect_gene))]
  x <- x %>% add_prediction_metadata()
  x[, `:=`(self_addiction = expression_gene == effect_gene, 
           neighbor = sorted_pair %in% important_pairs_default$proximity)]} %>% 
  cacheR(str_c(.dm, "_template"), file.path(large_cache, str_c(.dm, "_EV1")))


.i <- 5000
.fdr_cutoff <- 1000000

message(str_c("Collecting top ", .i, " scores for best performing methods..."))

EV1 <- list(PCC = file.path(
  cache, "23Q2_default_CvE_PCC_BaCoN_0.05", "correlation_matrix.rds"), 
  BaCoN = file.path(cache, "23Q2_default_CvE_PCC_BaCoN_0.05", 
                    "bacon_matrix.rds"), 
  #ComBat = file.path(
  #  large_cache, "23Q2_ComBat", 
  #  "23Q2_default_CvE_PCC_BaCoN_0.05_ComBat_exp_lineage_chr_lineage", 
  #  "correlation_matrix.rds"), 
  ComBat_BaCoN = file.path(
    large_cache, "23Q2_ComBat", 
    "23Q2_default_CvE_PCC_BaCoN_0.05_ComBat_exp_lineage_chr_lineage", 
    "bacon_matrix.rds"), 
  Cholesky = file.path(large_cache, "23Q2_whitening", 
                       "23Q2_default_CvE_PCC_BaCoN_0.05_whitening_Cholesky", 
                       "correlation_matrix.rds"), 
  Cholesky_BaCoN = file.path(large_cache, "23Q2_whitening", 
                             "23Q2_default_CvE_PCC_BaCoN_0.05_whitening_Cholesky", 
                             "bacon_matrix.rds"), 
  PCA = file.path(large_cache, "23Q2_whitening", 
                  "23Q2_default_CvE_PCC_BaCoN_0.05_whitening_PCA", 
                  "correlation_matrix.rds"), 
  PCA_BaCoN = file.path(large_cache, "23Q2_whitening", 
                        "23Q2_default_CvE_PCC_BaCoN_0.05_whitening_PCA", 
                        "bacon_matrix.rds"), 
  linreg_mlog10_pval = file.path(linreg_filepath, 
                                 "23Q2_default_CvE_linreg_matrix_eff1_exp2.lineage.rds")
)

EV1 <- cbind(template, 
             sapply(EV1, \(.n) {as.vector(readRDS(.n))}, simplify = F) %>% setDT)

message("PCC...")
EV1 <- EV1[order(PCC, decreasing = T)]
EV1[, PCC_rank := .I]

EV1[1:.fdr_cutoff, `:=`(PCC_FDR = prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  proximity_col = neighbor, 
  coex_z_score = 3))]


message("BaCoN...")
EV1 <- EV1[order(BaCoN, PCC, decreasing = T)]
EV1[, BaCoN_rank := .I]

EV1[1:.fdr_cutoff, `:=`(BaCoN_FDR = prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  proximity_col = neighbor, 
  coex_z_score = 3))]


message("ComBat + BaCoN...")
EV1 <- EV1[order(ComBat_BaCoN, PCC, decreasing = T)]
EV1[, ComBat_BaCoN_rank := .I]
EV1[1:.fdr_cutoff, `:=`(ComBat_BaCoN_FDR = prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  proximity_col = neighbor, 
  coex_z_score = 3))]



message("Cholesky...")
EV1 <- EV1[order(Cholesky, decreasing = T)]
EV1[, Cholesky_rank := .I]
EV1[1:.fdr_cutoff, `:=`(Cholesky_FDR = prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  proximity_col = neighbor, 
  coex_z_score = 3))]


message("Cholesky + BaCoN...")
EV1 <- EV1[order(Cholesky_BaCoN, PCC, decreasing = T)]
EV1[, Cholesky_BaCoN_rank := .I]
EV1[1:.fdr_cutoff, `:=`(Cholesky_BaCoN_FDR = prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  proximity_col = neighbor, 
  coex_z_score = 3))]



message("PCA...")
EV1 <- EV1[order(PCA, decreasing = T)]
EV1[, PCA_rank := .I]
EV1[1:.fdr_cutoff, `:=`(PCA_FDR = prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  proximity_col = neighbor, 
  coex_z_score = 3))]


message("PCA + BaCoN...")
EV1 <- EV1[order(PCA_BaCoN, PCC, decreasing = T)]
EV1[, PCA_BaCoN_rank := .I]
EV1[1:.fdr_cutoff, `:=`(PCA_BaCoN_FDR = prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  proximity_col = neighbor, 
  coex_z_score = 3))]

message("Linear regression...")
EV1 <- EV1[order(linreg_mlog10_pval, PCC, decreasing = T)]
EV1[, linreg_rank := .I]
EV1[, linreg_FDR_BH := p.adjust(10^(-linreg_mlog10_pval), method = "BH")]
EV1[1:.fdr_cutoff, `:=`(linreg_FDR = prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson | self_addiction), 
  proximity_col = neighbor, 
  coex_z_score = 3))]

saveRDS(EV1[(PCC_rank <= 50000 | 
               BaCoN_rank <= 50000 | 
               ComBat_BaCoN_rank <= 50000 | 
               Cholesky_rank <= 50000 | 
               Cholesky_BaCoN_rank <= 50000 | 
               PCA_rank <= 50000 | 
               PCA_BaCoN_rank <= 50000 | 
               linreg_rank <= 50000
)], file.path(cache, str_c(.dm, "_data_export"), 
              str_c("EV1_top_50000_predictions.rds")))
