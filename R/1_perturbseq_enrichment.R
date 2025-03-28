## ---- angela preparation ----

library(reshape2)
library(FLEX)

# BaCoN's 1000 Top Prediction Pairs (TPP)
TPP_Bacon <- hq_predictions

## Import PerturbSeq datasets

Replogle_K562 <- readRDS(
  file.path(perturbseq_filepath, "PerturbSeq_QCE.rds"))
Replogle_RPE1 <- readRDS(
  file.path(perturbseq_filepath, "RPE1_singleCellData.RDS"))
Replogle_K562_essential <- readRDS(
  file.path(perturbseq_filepath, "K562essential_singleCellData.RDS"))


## Import DepMap Gene Effects version 23Q2

CRISPRGeneEffect23Q2 <- readRDS(
  file.path(perturbseq_filepath, "CRISPRGeneEffect23Q2.rds"))

DM_PCC <- readRDS(file.path(perturbseq_filepath, "DM_PCC.rds"))

## ---- angela functions ----

z_normalize <- \(x) {
  res <- (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
  return(res)
}

normalized_expression <- \(x) {
  common_genes <- intersect(rownames(x),colnames(x))
  
  for (gene in common_genes) {
    x[gene, gene] <- NA
  }
  
  result <- t(apply(x, MARGIN = 2, z_normalize))
  
  return(result)
}

TopPredEnrichment <- \(NE, TPP = TPP_BaCoN, th) {
  df <- reshape2::melt(NE, varnames = c("expression_gene","effect_gene"), 
                       value.name = "z_score", as.is = T) %>%
    left_join(TPP[,c("expression_gene", "effect_gene", "score")], 
              by = c("expression_gene","effect_gene"))
  
  df <- df[df$expression_gene != df$effect_gene,]
  
  q = sum(!is.na(df$z_score) & !is.na(df$score) & df$z_score > th)
  m = sum(!is.na(df$z_score) & df$z_score > th)
  n = sum(!is.na(df$z_score) & df$z_score < th)
  k = sum(!is.na(df$score))
  
  result <- phyper(q - 1, m, n, k, lower.tail = F)
  print(str_c(q, " of the top predicted BaCoN pairs also have PS z-score of > ", 
              th))
  print(str_c(m, " pairs are above the PS threshold of ", th))
  print(str_c(k, " of the top predicted BaCoN pairs can be found in PS"))
  print(str_c("The foldchange is ", round((q / k) / (m / (n+m)), 2)))
  print(str_c("The p-value is ", round(result, 6)))
  
  return(result)
}

TopPredEnrichment_DM <- \(x, th) {
  
  gene_overlap <- colnames(x)[colnames(x) %in% TPP_BaCoN$expression_gene | 
                                colnames(x) %in% TPP_BaCoN$effect_gene]
  
  diag(x) <- NA
  mat_norm <- z_normalize(x)
  
  mat_norm_overlap <- mat_norm[gene_overlap, gene_overlap]
  
  mat_norm_overlap[upper.tri(mat_norm_overlap)] <- NA
  
  mat_norm[upper.tri(mat_norm)] <- NA
  
  mat_melt <- reshape2::melt(mat_norm_overlap, 
                             varnames = c("expression_gene", "effect_gene"), 
                             value.name = "z_score", as.is = T, na.rm = T)
  
  sorted_gene_pairs <- apply(mat_melt[,c("expression_gene", "effect_gene")], 1, 
                             \(x) str_c(sort(x), collapse = "_"))
  
  df <- mat_melt %>%
    mutate("sorted_pair" = sorted_gene_pairs) %>%
    left_join(TPP_BaCoN[,c("sorted_pair", "score")], by = "sorted_pair")
  
  q = sum(!is.na(df$score) & df$z_score > th)
  m = sum(mat_norm > th, na.rm = T)
  n = sum(mat_norm < th, na.rm = T)
  k = sum(!is.na(df$score))
  
  result <- phyper(q - 1, m, n, k, lower.tail = F)
  
  print(str_c(
    q, " of the top predicted BaCoN pairs also have a DM z-score of > ", th))
  print(str_c(m, " pairs are above the DM threshold of ", th))
  print(str_c(k, " of the top predicted BaCoN pairs can be found in DM"))
  print(str_c("The foldchange is ", round((q / k)/(m / (n + m)), 2)))
  print(str_c("The p-value is ", result))
  
  return(result)
}

## ---- angela analysis perturbseq ----

# Preprocessing: Z-normalize PerturbSeq expression on gene-level

K562_gw_NE <- normalized_expression(x = t(Replogle_K562))

rpe1_NE <- normalized_expression(x = t(Replogle_RPE1))

K562_essential_NE <- normalized_expression(x = t(Replogle_K562_essential))

# Enrichment Analysis

## Perturb-Seq

### Z-threshold = 2

K562_gw_pvalue <- TopPredEnrichment(NE = K562_gw_NE, th = 2)
rpe1_pvalue <- TopPredEnrichment(NE = rpe1_NE, th = 2)
K562_ess_pvalue <- TopPredEnrichment(NE = K562_essential_NE, th = 2)

### Z-threshold = 3

K562_gw_pvalue <- TopPredEnrichment(NE = K562_gw_NE, th = 3)
rpe1_pvalue <- TopPredEnrichment(NE = rpe1_NE, th = 3)
K562_ess_pvalue <- TopPredEnrichment(NE = K562_essential_NE, th = 3)

## ---- angela analysis depmap ----

## DepMap

### Z-threshold = 4

DM_pvalue <- TopPredEnrichment_DM(x = DM_PCC, th = 4)

### Z-threshold = 5

DM_pvalue <- TopPredEnrichment_DM(x = DM_PCC, th = 5)


