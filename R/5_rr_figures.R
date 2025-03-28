



















## ---- supplementary figures cholesky against pca ----

buffering <- empty_array(list(expression_genes, chronos_genes, c("PCC", "Cholesky", "PCA", "ComBat", 
                                                                 "PCC_z", "Cholesky_z", "PCA_z", "ComBat_z"
)))

.d <- list(PCC = file.path(cache_filepath, "23Q2_default_CvE_PCC_BaCoN_0.05", "correlation_matrix.rds"), 
           Cholesky = file.path(large_cache_filepath, "23Q2_whitening", "23Q2_default_CvE_PCC_BaCoN_0.05_whitening_Cholesky", "correlation_matrix.rds"), 
           PCA = file.path(large_cache_filepath, "23Q2_whitening", "23Q2_default_CvE_PCC_BaCoN_0.05_whitening_PCA", "correlation_matrix.rds"), 
           ComBat = file.path(large_cache_filepath, "23Q2_ComBat", "23Q2_default_CvE_PCC_BaCoN_0.05_ComBat_exp_lineage_chr_lineage", "correlation_matrix.rds")
)

for (.n in names(.d)) {
  message(.n)
  buffering[,,.n] <- readRDS(.d[[.n]])
  buffering[,,str_c(.n, "_z")] <- (buffering[,,.n] - mean(buffering[,,.n], na.rm = T)) / sd(buffering[,,.n], na.rm = T)}


coessentiality <- empty_array(list(chronos_genes, chronos_genes, c("PCC", "Cholesky", "PCA", "ComBat", 
                                                                   "PCC_z", "Cholesky_z", "PCA_z", "ComBat_z"
)))

coessentiality[,,"PCC"] <- cor(chronos, use = pco) %>% cacheR(str_c(.dm, "_default_CoEss_PCC"), cache_filepath)
coessentiality[,,"Cholesky"] <- cor(t(whitening::whiten(t(chronos), method = .m)), use = pco) %>% cacheR(str_c(.dm, "_default_CoEss_PCC_whitening_Cholesky"), cache_filepath)
coessentiality[,,"PCA"] <- cor(t(whitening::whiten(t(chronos), method = .m)), use = pco) %>% cacheR(str_c(.dm, "_default_CoEss_PCC_whitening_PCA"), cache_filepath)
coessentiality[,,"ComBat"] <- cor(combat_correction(chronos, "lineage"), use = pco) %>% cacheR(str_c(.dm, "_default_CoEss_PCC_ComBat_lineage"), cache_filepath)

for (.n in c("PCC", "Cholesky", "PCA", "ComBat")) {
  message(.n)
  coessentiality[,,.n][lower.tri(coessentiality[,,.n], diag = T)] <- NA
  coessentiality[,,str_c(.n, "_z")] <- (coessentiality[,,.n] - mean(coessentiality[,,.n], na.rm = T)) / sd(coessentiality[,,.n], na.rm = T)
}

.d <- list(top_buffering_pairs = prepare_predictions_for_plotting(
  results$IDs_23Q2$benchmarking[1], 
  results$all_predictions_23Q2, 10)$predictions[, .(expression_gene, effect_gene, sorted_pair, score)], 
  top_coess_pairs = collect_predictions(coessentiality[,,"PCC"], pairs = 10)[, .(expression_gene, effect_gene, sorted_pair, score)])

.d$top_coess_pairs[, `:=`(sorted_pair = fct_reorder(sorted_pair, score, .desc = T))]
.d$top_buffering_pairs[, `:=`(sorted_pair = fct_reorder(sorted_pair, score, .desc = T))]

.d$top_buffering_pairs <- sapply(c("PCC", "Cholesky", "PCA", "ComBat"), \(.n) {
  .d$top_buffering_pairs[, .(expression_gene, 
                             effect_gene, 
                             sorted_pair, 
                             variable = factor(.n, levels = c("PCC", "ComBat", "Cholesky", "PCA")), 
                             value = mapply(\(.g1, .g2) {buffering[.g1,.g2,.n]}, expression_gene, effect_gene), 
                             z_score = mapply(\(.g1, .g2) {buffering[.g1,.g2,str_c(.n, "_z")]}, expression_gene, effect_gene))]
}, simplify = F) %>% rbindlist()


.d$top_coess_pairs <- sapply(c("PCC", "Cholesky", "PCA", "ComBat"), \(.n) {
  .d$top_coess_pairs[, .(expression_gene, 
                         effect_gene, 
                         sorted_pair, 
                         variable = factor(.n, levels = c("PCC", "ComBat", "Cholesky", "PCA")), 
                         value = mapply(\(.g1, .g2) {coessentiality[.g1,.g2,.n]}, expression_gene, effect_gene), 
                         z_score = mapply(\(.g1, .g2) {coessentiality[.g1,.g2,str_c(.n, "_z")]}, expression_gene, effect_gene))]
}, simplify = F) %>% rbindlist()

.d$buff_cor <- cor(data.frame(sapply(c("PCC", "ComBat", "Cholesky", "PCA"), \(.) {as.vector(buffering[,,.])}, simplify = F)), use = pco)
.d$coess_cor <- cor(data.frame(sapply(c("PCC", "ComBat", "Cholesky", "PCA"), \(.) {as.vector(coessentiality[,,.])}, simplify = F)), use = pco)

.fpath <- file.path(fig_directory, "supplementary_figures", "Cholesky against PCA")

.plt_elements <- list(
  geom_col(position = position_dodge()), 
  scale_fill_viridis_d(begin = 0.1, end = 0.9, option = "G"), 
  theme(legend.title = element_blank(), 
        legend.position = "bottom", 
        axis.text.x = element_text(angle = 20, hjust = 1))
)

.d$top_buffering_pairs %>% 
  ggplot(aes(factor(sorted_pair), value, fill = factor(variable, levels = c("PCC", "ComBat", "Cholesky", "PCA")))) + .plt_elements + 
  labs(x = "Top 10 buffering predictions", y = "Correlation")

plotsaver(.fpath, "Top_buffering_predictions_correlation")

.d$top_buffering_pairs %>% 
  ggplot(aes(sorted_pair, z_score, fill = variable)) + .plt_elements + 
  labs(x = "Top 10 buffering predictions", y = "Z-score")

plotsaver(.fpath, "Top_buffering_predictions_z_scores")

.d$top_coess_pairs %>% 
  ggplot(aes(sorted_pair, value, fill = variable)) + .plt_elements + 
  labs(x = "Top 10 coessentiality pairs", y = "Correlation")

plotsaver(.fpath, "Top_coessential_pares_correlation")


.d$top_coess_pairs %>% 
  ggplot(aes(sorted_pair, z_score, fill = variable)) + .plt_elements + 
  labs(x = "Top 10 coessentiality pairs", y = "Z-score")

plotsaver(.fpath, "Top_coessential_pairs_z_score")


mapply(\(.m1, .m2) {data.table(methods = str_c(.m1, " ~ ", .m2), 
                               cor = .d$buff_cor[.m1,.m2])}, "PCC", c("ComBat", "Cholesky", "PCA"), SIMPLIFY = F) %>% 
  rbindlist() %>% 
  ggplot(aes(cor, factor(methods, levels = str_c("PCC ~ ", c("ComBat", "Cholesky", "PCA"))))) + geom_col() + 
  labs(x = "Global Buffering Correlation") + 
  theme(axis.title.y = element_blank())

plotsaver(.fpath, "global_correlation_buffering", 3, 3)

mapply(\(.m1, .m2) {data.table(methods = str_c(.m1, " ~ ", .m2), 
                               cor = .d$coess_cor[.m1,.m2])}, "PCC", c("ComBat", "Cholesky", "PCA"), SIMPLIFY = F) %>% rbindlist() %>% 
  ggplot(aes(cor, factor(methods, levels = str_c("PCC ~ ", c("ComBat", "Cholesky", "PCA"))))) + geom_col() + 
  labs(x = "Global Coessentiality Correlation") + 
  theme(axis.title.y = element_blank())

plotsaver(.fpath, "global_correlation_coessentiality", 3, 3)
```