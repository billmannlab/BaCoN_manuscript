
## ---- export datasets ----

library(openxlsx)
mkdir("data_export")

ev1_labels <- c(
  "buffered_gene" = "Buffered partner of the predicted gene pair", 
  "buffering_gene" = "Buffering partner of the predicted gene pair. The expression of the buffering gene is associated with a reduced cellular dependency on the buffered gene.", 
 # "sorted_pair" = "Alphabetically sorted concatenation of the buffering and the buffered gene", 
  
  "PCC" = "Pearson's correlation coefficient between the expression of the buffering gene and Chronos score of the buffered gene", 
  "BaCoN" = "BaCoN score of the prediction", 
  "ComBat_BaCoN" = "Score computed using ComBat (tissue) + PCC + BaCoN", 
  "Cholesky" = "Score computed using Cholesky whitening + PCC", 
  "Cholesky_BaCoN" = "Score computed using Cholesky whitening + PCC + BaCoN", 
  "PCA" = "Score computed using PCA whitening + PCC", 
  #"PCA_BaCoN" = "Score computed using PCA whitening + PCC + BaCoN", 
  "linreg_mlog10_pval" = "Negative logarithm of the p-value of the main predictor variable gene expression of the buffering gene, after the addiction between the Chronos scores of the buffered genes and the expression of the buffering gene has been modeled (buffered gene Chronos scores ~ buffering gene expression + Lineage)", 
  
  "PCC_z_score" = "Z-score of the Pearson's correlation coefficient between the expression of the buffering gene and Chronos score of the buffered gene", 
  
  "PCC_rank" = "Rank, defined by Pearson's Correlation", 
  "BaCoN_rank" = "Rank, defined by BaCoN score", 
  "ComBat_BaCoN_rank" = "Rank of the gene pair, defined by the score computed by ComBat (tissue) + PCC + BaCoN", 
  "Cholesky_rank" = "Rank of the gene pair, defined by the score computed by Cholesky whitening + PCC", 
  "Cholesky_BaCoN_rank" = "Rank of the gene pair, defined by the score computed by Cholesky whitening + PCC + BaCoN", 
  "PCA_rank" = "Rank of the gene pair, defined by the score computed by PCA whitening + PCC", 
  #"PCA_BaCoN_rank" = "Rank of the gene pair, defined by the score computed by PCA whitening + PCC + BaCoN", 
  "linreg_rank" = "Rank of the gene pair, defined by the negative logarithm of the linear model p-value", 
  
  "ensembl_paralog" = "The gene pair pair is annotated as Ensembl paralog with at least 20% sequence similarity", 
  "ohnolog" = "The pair is in the Ohnologs set at the standard threshold", 
  "anvar" = "The pair is in the Anvar et al. SL set", 
  "anvar_standard" = "The pair is in the Anvar et al. SL standard set (see methods for definition)", 
  "dekegel" = "The pair is in the DeKegel et al. SL set", 
  "dekegel_val" = "The pair is in the DeKegel et al. validated SL set", 
  "ito" = "The pair is in the Ito et al. SL set", 
  "thompson" = "The pair is in the Thompson et al. SL set", 
  
  "buffering_chr" = "Chromosome of the buffering gene", 
  "buffering_chr_arm" = "Chromosome arm of the buffering gene", 
  "buffering_pos" = "Location of the buffering gene on its chromosome arm", 
  "buffered_chr" = "Chromosome of the buffered gene", 
  "buffered_chr_arm" = "Chromosome arm of the buffered gene", 
  "buffered_pos" = "Location of the buffered gene on its chromosome arm", 
  "genomic_bp_distance" = "If both genes are on the same chromosome arm, genomic distance between the buffering and the buffered gene, in base pairs", 
  "proximity" = "Buffered and buffering genes are on the same arm and less than 10Mbp apart", 
  "self_addiction" = "The buffering and the buffered gene are identical", 
  
  "PCC_FDR" = "False discovery rate of the gene pair, computed on PCC scores", 
  "BaCoN_FDR" = "False discovery rate of the gene pair, computed on PCC + BaCoN", 
  "ComBat_BaCoN_FDR" = "False discovery rate of the gene pair, computed on ComBat (tissue) + PCC + BaCoN", 
  "Cholesky_FDR" = "False discovery rate of the gene pair, computed on Cholesky whitening + PCC", 
  "Cholesky_BaCoN_FDR" = "False discovery rate of the gene pair, computed on Cholesky whitening + PCC + BaCoN", 
  "PCA_FDR" = "False discovery rate of the gene pair, computed on PCA whitening + PCC", 
  "linreg_FDR" = "False discovery rate of the gene pair, computed on linear regression", 
  "high_confidence_prediction" = "The prediction is part of the 808 high-confidence predictions, based on Cholesky whitening + PCC + BaCoN, FDR < 0.5, after coexpressed and proximal pairs are removed."
 )

ev2_labels <- c(
  "buffered_gene" = "Buffered partner of the predicted gene pair", 
  "buffering_gene" = "Buffering partner of the predicted gene pair. The expression of the buffering gene is associated with a reduced cellular dependency on the buffered gene.", 
  "Cholesky_BaCoN" = "Score computed using Cholesky whitening + PCC + BaCoN", 
  "PCC" = "Pearson's correlation coefficient between the expression of the buffering gene and Chronos score of the buffered gene", 
  "PCC_z_score" = "Z-score of the Pearson's correlation coefficient between the expression of the buffering gene and Chronos score of the buffered gene", 
  "Cholesky_BaCoN_rank" = "Rank of the gene pair, defined by the score computed by Cholesky whitening + PCC + BaCoN", 
  "ensembl_paralog" = "The gene pair pair is annotated as Ensembl paralog with at least 20% sequence similarity", 
  "ohnolog" = "The pair is in the Ohnologs set at the standard threshold", 
  "anvar" = "The pair is in the Anvar et al. SL set", 
  "anvar_standard" = "The pair is in the Anvar et al. SL standard set (see methods for definition)", 
  "dekegel" = "The pair is in the DeKegel et al. SL set", 
  "dekegel_val" = "The pair is in the DeKegel et al. validated SL set", 
  "ito" = "The pair is in the Ito et al. SL set", 
  "thompson" = "The pair is in the Thompson et al. SL set", 
  "buffering_chr" = "Chromosome of the buffering gene", 
  "buffering_chr_arm" = "Chromosome arm of the buffering gene", 
  "buffered_chr" = "Chromosome of the buffered gene", 
  "buffered_chr_arm" = "Chromosome arm of the buffered gene", 
  "Cholesky_BaCoN_FDR" = "False discovery rate of the gene pair, computed on Cholesky whitening + PCC + BaCoN")

ev3_labels <- c(
  "buffered_gene" = "Buffered partner of the predicted gene pair", 
  "buffering_gene_ncRNA" = "Buffering partner of the predicted gene pair. The expression of the buffering gene or ncRNA is associated with a reduced cellular dependency on the buffered gene", 
  "prediction_class" = "Describes if the buffering and buffered partner are protein-coding genes or non-coding RNAs", 
  "Cholesky_BaCoN" = "Score computed using Cholesky whitening + PCC + BaCoN", 
  "PCC" = "Pearson's correlation coefficient between the expression of the buffering gene and Chronos score of the buffered gene",   
  "PCC_z_score" = "Z-score of the Pearson's correlation coefficient between the expression of the buffering gene and Chronos score of the buffered gene", 
  "Cholesky_BaCoN_rank" = "Rank of the gene pair, defined by the score computed by Cholesky whitening + PCC + BaCoN", 
  "buffering_chr" = "Chromosome of the buffering gene", 
  "buffering_chr_arm" = "Chromosome arm of the buffering gene", 
  "buffered_chr" = "Chromosome of the buffered gene", 
  "buffered_chr_arm" = "Chromosome arm of the buffered gene", 
  "genomic_bp_distance"  = "If both genes are on the same chromosome arm, genomic distance between the buffering and the buffered gene, in base pairs", 
  "proximity" = "Buffered and buffering genes are on the same arm and less than 10mbp apart")


data_export <- lapply(list(EV1_labels = list(column_label = names(ev1_labels), 
                                        description = ev1_labels), 
                    EV2_labels = list(column_label = names(ev2_labels), 
                                       description = ev2_labels), 
                    EV3_labels = list(column_label = names(ev3_labels), 
                                       description = ev3_labels)), setDT)


data_export$EV1 <- readRDS(file.path(cache, 
                                     str_c(.dm, "_data_export"), "EV1_top_50000_predictions.rds"))

data_export$EV1 <- data_export$EV1[PCC_rank <= 5000 | 
                                     BaCoN_rank <= 5000 |
                                     ComBat_BaCoN_rank <= 5000 |
                                     Cholesky_rank <= 5000 |
                                     Cholesky_BaCoN_rank <= 5000 |
                                     PCA_rank <= 5000 | 
                                     PCA_BaCoN_rank <= 5000 |
                                     linreg_rank <= 5000]

data_export$EV1 <- data_export$EV1 %>% 
  merge(exp_gene_info, by = "expression_gene", all.x = T, sort = F) %>% 
  merge(eff_gene_info, by = "effect_gene", all.x = T, sort = F)

data_export$EV1[exp_chr_arm == eff_chr_arm, distance := abs(exp_pos - eff_pos)]
data_export$EV1[, `:=`(
  self_addiction = expression_gene == effect_gene, 
  proximity = sorted_pair %in% important_pairs_complete$proximity, 
  proximity_class = fcase(
    expression_gene == effect_gene, "self-addiction", 
    exp_chr_arm == eff_chr_arm & distance <= 1e7, "same arm, <10mbp", 
    exp_chr_arm == eff_chr_arm, "same arm, >10mbp", 
    exp_chr == eff_chr & exp_chr_arm != eff_chr_arm, "different arm", 
    exp_chr != eff_chr, "different chromosomes", 
    default = "NA"))]

data_export$EV1[, PCC_z_score := (PCC - 0.0002780176) / 0.0553656]


data_export$EV1[, high_confidence_prediction := Cholesky_BaCoN_FDR <= 0.5 & 
                  Cholesky_BaCoN_rank %in% results$hq_predictions_23Q2_fdr_0.5$rank]

data_export$EV1 <- data_export$EV1[order(Cholesky_BaCoN_rank, PCC_rank), .(
  buffered_gene = effect_gene, 
  buffering_gene = expression_gene, 
  #sorted_pair, 
  PCC, 
  BaCoN, 
  ComBat_BaCoN, 
  Cholesky, 
  Cholesky_BaCoN, 
  PCA, 
  #PCA_BaCoN, 
  linreg_mlog10_pval, 
  
  PCC_z_score, 
  
  PCC_rank, 
  BaCoN_rank, 
  ComBat_BaCoN_rank, 
  Cholesky_rank, 
  Cholesky_BaCoN_rank, 
  PCA_rank, 
  #PCA_BaCoN_rank, 
  linreg_rank, 
  
  ensembl_paralog = as.numeric(ensembl), 
  ohnolog = as.integer(ohnologs), 
  anvar = as.integer(anvar), 
  anvar_standard = as.integer(anvar_standard),  
  dekegel = as.integer(dekegel), 
  dekegel_val = as.integer(dekegel_val), 
  ito = as.integer(ito), 
  thompson = as.integer(thompson), 
  
  buffering_chr = exp_chr, 
  buffering_chr_arm = exp_chr_arm, 
  buffering_pos = exp_pos, 
  buffered_chr = eff_chr, 
  buffered_chr_arm = eff_chr_arm, 
  buffered_pos = eff_pos, 
  genomic_bp_distance = distance, 
  proximity = as.integer(neighbor), 
  self_addiction = as.integer(self_addiction), 
  PCC_FDR, 
  BaCoN_FDR, 
  ComBat_BaCoN_FDR, 
  Cholesky_FDR, 
  Cholesky_BaCoN_FDR, 
  PCA_FDR, 
  linreg_FDR, 
  high_confidence_prediction = as.numeric(high_confidence_prediction)
  )]


# Derive EV2 from EV1:


data_export$EV2 <- data_export$EV1[high_confidence_prediction == 1][order(Cholesky_BaCoN_rank)]

data_export$EV2 <- data_export$EV2[, .(buffered_gene, 
                                       buffering_gene, 
                                       Cholesky_BaCoN, 
                                       PCC, 
                                       PCC_z_score, 
                                       Cholesky_BaCoN_rank, 
                                       ensembl_paralog, 
                                       ohnolog, 
                                       anvar, 
                                       anvar_standard, 
                                       dekegel, 
                                       dekegel_val, 
                                       ito, 
                                       thompson, 
                                       buffering_chr, 
                                       buffering_chr_arm, 
                                       buffered_chr, 
                                       buffered_chr_arm, 
                                       Cholesky_BaCoN_FDR)]

all(data_export$EV2[, sort_gene_pairs(buffered_gene, buffering_gene)] == 
      results$hq_predictions_23Q2_fdr_0.5$sorted_pair)


data_export$EV3 <- list(hq_set = results$hq_predictions_23Q2_fdr_0.5, 
                        ncRNAs = copy(readRDS(file.path(cache, "23Q2_data_export", "55_lncRNAs.rds"))))

data_export$EV3$hq_set[, prediction_class := "protein_coding_protein_coding"]
data_export$EV3$ncRNAs[, `:=`(proximity = eff_chr_arm == exp_chr_arm & distance <= proximity_threshold, 
                              prediction_class = str_c(effect_class, "_", expression_class))]

data_export$EV3 <- sapply(data_export$EV3, \(.) {.[
  , .SD, .SDcols = 
    intersect(names(data_export$EV3$hq_set), 
              names(data_export$EV3$ncRNAs))]}, simplify = F) %>% rbindlist()



data_export$EV3 <- data_export$EV3[order(score, pcc_reference, decreasing = T), 
                                   .(buffered_gene = effect_gene, 
                                     buffering_gene_ncRNA = expression_gene, 
                                     prediction_class, 
                                     Cholesky_BaCoN = score, 
                                     PCC = pcc_reference, 
                                     PCC_z_score = z_score, 
                                     Cholesky_BaCoN_rank = rank, 
                                     buffering_chr = exp_chr, 
                                     buffering_chr_arm = exp_chr_arm, 
                                     buffered_chr = eff_chr, 
                                     buffered_chr_arm = eff_chr_arm, 
                                     genomic_bp_distance = distance, 
                                     proximity = as.integer(proximity))]

for (.n in c("EV1", "EV2", "EV3")) {
  if (all(names(data_export[[.n]]) == data_export[[str_c(.n, "_labels")]]$column_label)) {
    message(str_c(.n, ": All columns match!"))
    
    fwrite(data_export[[.n]], file.path("data_export", str_c("Dataset_", .n, ".csv")))
    fwrite(data_export[[str_c(.n, "_labels")]], file.path("data_export", str_c("Dataset_", .n, "_labels.csv")), sep = ";")
    
    saveRDS(data_export[[.n]], file.path("data_export", str_c("Dataset_", .n, ".rds")))
    saveRDS(data_export[[str_c(.n, "_labels")]], file.path("data_export", str_c("Dataset_", .n, "_labels.rds")))
    
    .wb <- createWorkbook()
    addWorksheet(.wb, "data")
    writeData(.wb, "data", data_export[[.n]])
    addWorksheet(.wb, "legend")
    writeData(.wb, "legend", data_export[[str_c(.n, "_labels")]])
    saveWorkbook(.wb, file.path("data_export", str_c("Dataset_", .n, ".xlsx")), overwrite = T)
   
  } else {
    print(setdiff(names(data_export[[.n]]), data_export[[str_c(.n, "_labels")]]$column_label))
    print(setdiff(data_export[[str_c(.n, "_labels")]]$column_label, names(data_export[[.n]])))
  }}
