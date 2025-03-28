

results$hq_predictions_23Q2_fdr_0.5 %>% 
  ggplot() + 
  geom_histogram(
    mapping = aes(pcc_reference, 
                  fill = factor(z_score_bin, 
                                levels = c("< 0", "0-1", "1-2", "2-3", "> 3"))), 
    alpha = 0.9, bins = 100) + 
  scale_fill_manual(values = z_score_colors) + 
  labs(x = "PCC of top predictions", y = "Cholesky whitening + PCC + BaCoN", 
       fill = "z-score", 
       caption = results$hq_predictions_23Q2_fdr_0.5[
         sorted_pair %in% pairs_oi, 
         str_c(sorted_pair, " : ", round(pcc_reference, 5), collapse = ",\n")]) + 
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "bottom")

plotsaver(file.path(fig_directory, "figure_5"), "fig5a_hq_set_PCC_histogram", 5, 4)


list(pcc = data.table(
  name = "PCC, all", 
  pcc_reference = readRDS(file.path(
    cache_filepath, "23Q2_results", "sample_scores", 
    
    "23Q2_default_CvE_PCC_BaCoN_0.05_sample_scores.rds"))$reference[1:10000]),
  top_pcc = prepare_predictions_for_plotting(data.table(
    id = "23Q2_default_CvE_PCC", name = "Top 1k PCCs"), 
    results$all_predictions_23Q2, 1000)$predictions[, .(name, pcc_reference = score)], 
  high_conf = results$hq_predictions_23Q2_fdr_0.5[
    , .(pcc_reference, name = "Top predictions")]) %>% 
  Reduce(bind_rows, .) %>% 
  ggplot(aes(pcc_reference, color = name)) + 
  geom_density() + 
  labs(x = "PCC") + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        legend.title = element_blank(), legend.position = "top")

plotsaver(file.path(fig_directory, "figure_5"), "fig5a_hq_set_PCC_density", 5, 3)

## ---- plots hq set linreg density ----


fig5_data <- list()
fig5_data$linreg_top_preds <- prepare_predictions_for_plotting(
  results$IDs_23Q2$benchmarking[12], results$all_predictions_23Q2, 5000)$predictions

lm_matrix_lineage_complete <- build_lm_matrix(
  rownames = colnames(gene_expression_complete), 
  colnames = colnames(chronos_complete), 
  files_needed = chunk_definition[, name], 
  filepath = file.path(
    large_cache_filepath, 
    str_c(.dm, "_linear_regression_new_and_clean"), "eff1_exp2.lineage")) %>% 
  cacheR("23Q2_default_CvE_linreg_matrix_eff1_exp2.lineage", cache_filepath)

fig5_data$linreg_all_scores <- lm_matrix_lineage_complete %>% 
  melt_array_to_dt("expression_gene", "effect_gene", "mlog10_pval")

fig5_data$linreg_all_scores[, `:=`(pcc_reference = as.vector(CvE_PCC))]

fig5_data$linreg_all_scores <- fig5_data$linreg_all_scores[
  order(mlog10_pval, pcc_reference, decreasing = T)]

fig5_data$linreg_all_scores[, `:=`(
  rank = .I, 
  FDR_BH = p.adjust(10^(-mlog10_pval), method = "BH"), 
  FDR_bonf = p.adjust(10^(-mlog10_pval), method = "bonferroni"))]

fig5_data$linreg_FDR_comparison <- fig5_data$linreg_all_scores[1:100000]
fig5_data$linreg_FDR_comparison[, sorted_pair := sort_gene_pairs(expression_gene, effect_gene)]
fig5_data$linreg_FDR_comparison <- fig5_data$linreg_FDR_comparison %>% add_prediction_metadata()

fig5_data$linreg_FDR_comparison[, FDR := prediction_FDR(
  exp_gene_col = expression_gene, 
  eff_gene_col = effect_gene, 
  TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | 
              thompson | expression_gene == effect_gene), 
  save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | 
                thompson | expression_gene == effect_gene), 
  proximity_col = (sorted_pair %in% important_pairs_default$proximity), 
  coex_z_score = 3)]

fig5_data$FDR_cutoff_lines <- c(sapply(c(0.05, 0.1, 0.2), \(.fdr) {
  fig5_data$linreg_FDR_comparison[
    FDR <= .fdr, .(facet = "FDR", FDR = .fdr, 
                   mlog10_pval = min(mlog10_pval, na.rm = T))]}, simplify = F), 
  sapply(c(0.05, 0.1, 0.2), \(.fdr) {
    fig5_data$linreg_all_scores[
      FDR_BH <= .fdr, .(facet = "FDR_BH", 
                        FDR = .fdr, mlog10_pval = min(mlog10_pval, na.rm = T))]
  }, simplify = F), 
  sapply(c(0.05, 0.1, 0.2), \(.fdr) {
    fig5_data$linreg_all_scores[
      FDR_bonf <= .fdr, .(facet = "FDR_bonf", FDR = .fdr, 
                          mlog10_pval = min(mlog10_pval, na.rm = T))]}, simplify = F)
) %>% rbindlist()


.d <- list(top_linreg_preds = fig5_data$linreg_top_preds[
  , 
  .(effect_gene, expression_gene, mlog10_pval = score, pcc_reference, 
    facet = "Top linreg predictions")], 
  
  all_pairs = fig5_data$linreg_all_scores[sample(1:.N, 1000000), 
                                          .(effect_gene, 
                                            expression_gene, 
                                            mlog10_pval, 
                                            pcc_reference, 
                                            facet = "All pairs")], 
  
  hq_predictions = fig5_data$linreg_all_scores[
    which(expression_gene %in% 
            results$hq_predictions_23Q2_fdr_0.5[, unique(expression_gene)] & 
            effect_gene %in% 
            results$hq_predictions_23Q2_fdr_0.5[, unique(effect_gene)]), 
    .(effect_gene, expression_gene, sorted_pair = 
        sort_gene_pairs(expression_gene, effect_gene), 
      mlog10_pval, pcc_reference, facet = "Top predictions")])

.d$hq_predictions <- .d$hq_predictions[
  sorted_pair %in% results$hq_predictions_23Q2_fdr_0.5[, get("sorted_pair")]]

.d %>% Reduce(bind_rows, .) %>% 
  ggplot(aes(mlog10_pval, color = facet)) + 
  geom_density(linewidth = 1) + 
  scale_color_manual(values = c("All pairs" = "darkgrey", 
                                "Top linreg predictions" = "grey", 
                                "Top predictions" = "black")) + 
  geom_vline(aes(xintercept = mlog10_pval), 
             data = fig5_data$FDR_cutoff_lines[
               facet %in% c("FDR_BH", "FDR_bonf") & FDR == 0.05], 
             linetype = "dashed") + 
  scale_x_continuous(limits = c(-20, 50)) + 
  labs(x = "Multiple linear regression -log10 p-value of top predictions", 
       caption = fig5_data$FDR_cutoff_lines[
         FDR == 0.05, 
         str_c("-log10(pval) for FDR 5%:\n", 
               str_c(facet, ": ", round(mlog10_pval, 5), collapse = ",\n"))]) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank())

plotsaver(file.path(fig_directory, "figure_5"), 
          "fig5b_low_z_score_pairs_linreg_density", 6, 4)

# for reviewer response

for (.top in c(100, 1000)) {
  
  fig5_data$linreg_FDR_comparison[, coexpressed := F][
    effect_gene %in% unique(effect_gene[duplicated(effect_gene)]), 
    coexpressed := flag_coex(expression_gene), by = effect_gene]
  
  fig5_data$linreg_FDR_comparison[, `:=`(
    coexpressed = cumsum(coexpressed) / .I, 
    self_addiction = cumsum(effect_gene == expression_gene) / .I, 
    proximity = cumsum(sort_gene_pairs(expression_gene, effect_gene) %in% 
                         important_pairs_default$proximity) / .I)]
  
  
  fig5_data$linreg_FDR_comparison[1:.top] %>% 
    melt.data.table(measure.vars = c("FDR", "FDR_bonf", "coexpressed", "proximity")) %>% 
    ggplot(aes(rank, value)) + 
    geom_line(aes(color = variable), linewidth = 1) + 
    geom_vline(xintercept = 100) + 
    scale_y_continuous(sec.axis = sec_axis(~.*.top, name = "FPs")) + 
    labs(x = str_c("Top ", .top, " linreg predictions"), y = "FDR") + 
    theme(legend.position = "bottom", legend.title = element_blank())
  
  plotsaver(file.path(fig_directory), str_c("fig5_linreg_FDR_bonf_vs_FPa_top_", .top), 6, 4.5)
  
  fig5_data$linreg_FDR_comparison[1:.top] %>% 
    melt.data.table(measure.vars = c("FDR", "coexpressed", "proximity")) %>% 
    ggplot(aes(FDR_bonf, value)) + 
    geom_line(aes(color = variable), linewidth = 1) + 
    labs(x = str_c("Bonferrony FDR of top ", .top, " linreg predictions"), y = "") + 
    theme(legend.position = "bottom", legend.title = element_blank())
  
  plotsaver(file.path(fig_directory), str_c("fig5_linreg_FDR_bonf_vs_FPb_top_", .top), 6, 4.5)}

## ---- plots large prediction set identity ----


list(
  paralog_pairs$ensembl[, `:=`(facet = "All paralog pairs")], 
  paralog_pairs$ensembl[sorted_pair %in% 
                          results$hq_predictions_23Q2_fdr_0.5[(ensembl), sorted_pair]][
                            , `:=`(facet = "Buffering paralog pairs")]) %>% 
  Reduce(bind_rows, .) %>% 
  ggplot(aes(ident, color = facet)) + 
  geom_density(linewidth = 1.5) + 
  scale_color_manual(values = c("All paralog pairs" = "darkgrey", 
                                "Buffering paralog pairs" = "black")) + 
  labs(x = "Sequence Identity [%]", y = "Density") +
  theme(legend.position = c(0.7, 0.9), 
        strip.background = element_blank(), legend.title = element_blank())

plotsaver(file.path(fig_directory, "figure_5"), 
          "fig5e_paralog_prediction_identity_density", 4.5, 3)


## ---- plots large prediction set identity binned ----

.d <- copy(paralog_pairs$ensembl)[, ident_bin := fcase(
  ident %between% c(0, 10), "0 - 10", ident %between% c(10, 20), "10 - 20",
  ident %between% c(20, 30), "20 - 30", ident %between% c(30, 40), "30 - 40",
  ident %between% c(40, 50), "40 - 50", ident %between% c(50, 60), "50 - 60",
  ident %between% c(60, 70), "60 - 70", ident %between% c(70, 80), "70 - 80",
  ident %between% c(80, 90), "80 - 90",ident %between% c(90, 100), "90 - 100")]

.d[, ident_bin := factor(
  ident_bin, 
  levels = str_c(c("0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", 
                   "50 - 60", "60 - 70", "70 - 80", "80 - 90", "90 - 100")))]

.d2 <- merge(results$hq_predictions_23Q2_fdr_0.5, .d, by = "sorted_pair", all.x = T)


merge(.d[, .(total_pairs = .N), by = ident_bin], 
      .d2[, .(predicted = .N), by = ident_bin]) %>% 
  ggplot(aes(ident_bin, (predicted*100) / total_pairs)) + geom_col() + 
  labs(x = "Sequence Identity [%]", y = "% buffering paralogs") + 
  theme(axis.text.x = element_text(angle = 90))

plotsaver(file.path(fig_directory, "figure_5"), 
          "fig5f_paralog_identity_binned", 3, 2.5)

## ---- plots large prediction set chromosomal distribution ----

.d <- bind_rows(
  results$hq_predictions_23Q2_fdr_0.5[
    , .(facet = "buffering genes", n_pairs = .N), by = .(chr = exp_chr)], 
  results$hq_predictions_23Q2_fdr_0.5[
    , .(facet = "buffered genes", n_pairs = .N), by = .(chr = eff_chr)])

.d[!is.na(chr)] %>% 
  ggplot(aes(factor(chr, levels = c(1:22, "X", "Y")), n_pairs)) + 
  geom_col(aes(fill = facet), position = position_dodge2()) + 
  scale_fill_manual(values = c("buffered genes" = "#782E93", 
                               "buffering genes" = "#16A567")) + 
  labs(x = "chromosome", y = "# pairs") + 
  theme(legend.position = c(0.85, 0.85), legend.title = element_blank())

plotsaver(.fpath, "fig5i_chromosomal_distribution", 6, 4)


## ---- plots large prediction set chromosome pairs fc ----

chromosomes_oi <- c(1:22, "X", "Y")

.d <- list()

goi <- union(colnames(gene_expression), colnames(chronos))
gene_space <- ncol(gene_expression) * ncol(chronos)

for (.c1 in chromosomes_oi) {
  for (.c2 in chromosomes_oi) {
    .all_pairs <- prod(gene_database[symbol %in% goi & chromosome == .c1, .N], 
                       gene_database[symbol %in% goi & chromosome == .c2, .N])
    
    .all_bufferings <- prod(
      gene_database[symbol %in% colnames(gene_expression) & chromosome == .c1, .N], 
      gene_database[symbol %in% colnames(chronos) & chromosome == .c2, .N])
    .buffer_preds <- results$hq_predictions_23Q2_fdr_0.5[eff_chr == .c1 & exp_chr == .c2, .N]
    
    .d[[str_c(.c1, .c2)]] <- data.table(
      buffered_chr = .c1, 
      buffer_chr = .c2, 
      all_pairs = .all_pairs, 
      all_bufferings = .all_bufferings, 
      buffer_preds = .buffer_preds, 
      
      FC = (.buffer_preds/1000) / (.all_pairs/gene_space), 
      pval = phyper(q = .buffer_preds - 1, m = .all_bufferings, k = 1000, 
                    n = gene_space - .all_bufferings, lower.tail = F))}}

.d <- copy(Reduce(bind_rows, .d))
.d[, `:=`(fdr = p.adjust(pval, method = "BH"), LFC = log2(FC))]


.scale <- c(0.1, 1, 2.5, 5, 7.5, 10, 11)

for (.fdr in c(0.1, 0.2, 0.5)) {
  .colscale <- c("black", "grey")
  names(.colscale) <- c(str_c("FDR ", .fdr, "%"), "n.s.")
  
  .d %>% ggplot(aes(factor(buffer_chr, levels = chromosomes_oi), 
                    factor(buffered_chr, levels = rev(chromosomes_oi)))) + 
    scale_fill_manual(values = .colscale) + 
    geom_point(aes(size = LFC**2, 
                   fill = fcase(fdr < .fdr, str_c("FDR ", .fdr, "%"), default = "n.s.")), 
               shape = 21, color = "grey") + 
    scale_size_continuous(range = c(1, 10), 
                          limits = c(NA, NA), 
                          breaks = .scale**2, labels = .scale) + 
    scale_x_discrete(expand = expansion(mult = 0.04, add = 0)) + 
    scale_y_discrete(expand = expansion(mult = 0.04, add = 0)) + 
    labs(x = "buffering chromosome", y = "buffered chromosome") + 
    theme(legend.position = "right", #legend.box = "vertical", 
          legend.title = element_blank())
  
  plotsaver(.fpath, str_c("fig5j_chromosome_pairs_fc_fdr_", .fdr), 5.5, 4.5)
}


.d <- list(preds = results$hq_predictions_23Q2_fdr_0.5)

.d$enrichments$ensembl <- data.table(
  pair_set = "Paralogs", 
  draws = .d$preds[(ensembl), .N], # no selfs
  whites_drawn = .d$preds[eff_chr == exp_chr & ensembl, .N], 
  whites_in_urn = gene_pair_set_stats[
    gene_space == "default" & gene_pair_set == "syntenic_ensembl", get("hits")], 
  urn_size = gene_pair_set_stats[
    gene_space == "default" & gene_pair_set == "ensembl", get("hits")])

.d$enrichments$ohnologs <- data.table(
  pair_set = "Ohnologs", 
  draws = .d$preds[(ohnologs), .N], # no selfs
  whites_drawn = .d$preds[eff_chr == exp_chr & ohnologs, .N], 
  whites_in_urn = gene_pair_set_stats[
    gene_space == "default" & gene_pair_set == "syntenic_ohnologs", get("hits")], 
  urn_size = gene_pair_set_stats[
    gene_space == "default" & gene_pair_set == "ohnologs", get("hits")])

.d$enrichments$all_pairs <- data.table(
  pair_set = "All pairs", 
  draws = .d$preds[, .N], 
  whites_drawn = .d$preds[eff_chr == exp_chr, .N], 
  whites_in_urn = gene_pair_set_stats[
    gene_space == "default" & gene_pair_set == "same_chromosome", get("hits")], 
  urn_size = gene_pair_set_stats[
    gene_space == "default" & gene_pair_set == "same_chromosome", get("n_pairs")])


.d$enrichments$all_but_ensembl <- data.table(
  pair_set = "All non-Ensembl pairs", 
  draws = .d$preds[!(ensembl),.N], # no selfs
  whites_drawn = .d$preds[eff_chr == exp_chr & !ensembl, .N], 
  whites_in_urn = gene_pair_set_stats[
    gene_space == "default" & gene_pair_set == "syntenic_no_ensembl", get("hits")], 
  urn_size = gene_pair_set_stats[
    gene_space == "default" & gene_pair_set == "syntenic_no_ensembl", get("n_pairs")])

.d$enrichments$all_but_ensembl_and_ohnologs <- data.table(
  pair_set = "All pairs without Ensembl and Ohnologs", 
  draws = .d$preds[!(ensembl | ohnologs), .N], # no selfs
  whites_drawn = .d$preds[eff_chr == exp_chr & !(ensembl | ohnologs), .N], 
  whites_in_urn = gene_pair_set_stats[
    gene_space == "default" & 
      gene_pair_set == "syntenic_no_ensembl_no_ohnologs", get("hits")], 
  urn_size = gene_pair_set_stats[
    gene_space == "default" & 
      gene_pair_set == "syntenic_no_ensembl_no_ohnologs", get("n_pairs")])


.d$enrichments <- rbindlist(.d$enrichments)


.d$enrichments[, FC := (whites_drawn / draws) / (whites_in_urn / urn_size)]
.d$enrichments[whites_drawn != 0, `:=`(LFC = log2(FC), lower.tail = log2(FC) <= 0)]
.d$enrichments[, pval := phyper(
  # white balls are syntenic pairs among predictions, black balls are non-syntenic ones
  q = whites_drawn - 1, # white balls drawn
  m = whites_in_urn, # white balls in urn
  k = draws - whites_drawn, # black balls drawn
  n = (urn_size - whites_in_urn), # black balls in urn
  lower.tail = lower.tail, # lower.tail = T tests for depletion, lower.tail = F tests for enrichment
), by = .I][, sign := fcase(pval %between% c(0.01, 0.05), "*", 
                            pval %between% c(0.001, 0.01), "**", 
                            pval <= 0.001, "***", default = "")]

.pred_set_titles <- c("manuscript" = "Manuscript version", 
                      "old_hq_with_neighbors" = "Old, with neighbors", 
                      "old_hq_no_neighbors" = "Old, no neighbors", 
                      "new_hq_fdr_0.1" = "New, FDR 0.1", 
                      "new_hq_fdr_0.2" = "New, FDR 0.2", 
                      "new_hq_fdr_0.25" = "New, FDR 0.25", 
                      "new_hq_fdr_0.5" = "New, FDR 0.5")

#.d[, set := factor(.pred_set_titles[set], levels = .pred_set_titles)]

.d$enrichments %>% 
  ggplot(aes(LFC, pair_set, fill = pair_set)) + 
  geom_col() + 
  scale_fill_manual(values = c("All pairs" = "darkgrey", 
                               "All non-Ensembl pairs" = "darkgrey", 
                               "All pairs without Ensembl and Ohnologs" = "darkgrey", 
                               "Paralogs" = "#BC2C1A", 
                               "Ohnologs" = "#EC0868")) + 
  labs(x = "Syntenic buffering [LFC]", 
       caption = .d$enrichments[, str_c("P-values:\n", 
                                        str_c(pair_set, " : ", pval, collapse = ",\n"))]
  ) + 
  guides(fill = guide_legend(ncol = 2)) + 
  theme(legend.position = "bottom", axis.title.y = element_blank(), 
        legend.title = element_blank())

plotsaver(.fpath, "fig5h_synteny_enrichments", 6, 4)

## ---- plots circos plots circos plot ----

.d <- results$hq_predictions_23Q2_fdr_0.5[!is.na(eff_chr) & !is.na(exp_chr)]

pdf(file.path(.fpath, "fig5g_hq_set_circos_all.pdf"), 
    width = 3.5, height = 3.5)

draw_circos_plot(prepare_circos_data_from_predictions(.d), 
                 .d[, .(lwd = 0.75, 
                        color = fcase(ensembl & eff_chr == exp_chr, "#D73027", 
                                      ensembl & eff_chr != exp_chr, "#FDAE61", 
                                      !ensembl & eff_chr == exp_chr, "#4575B4", 
                                      !ensembl & eff_chr != exp_chr, "#ABD9E9"))])

dev.off()

pdf(file.path(.fpath, "fig5g_hq_set_circos_paralogs.pdf"), 
    width = 3.5, height = 3.5)


draw_circos_plot(
  prepare_circos_data_from_predictions(
    .d[(ensembl)]), 
  .d[(ensembl), .(lwd = 0.75, 
                  color = fcase(ensembl & eff_chr == exp_chr, "#D73027", 
                                ensembl & eff_chr != exp_chr, "#FDAE61"))])


dev.off()

pdf(file.path(.fpath, "fig5g_hq_set_circos_no_paralogs.pdf"), 
    width = 3.5, height = 3.5)

draw_circos_plot(prepare_circos_data_from_predictions(.d[!(ensembl)]), 
                 .d[!(ensembl), .(
                   lwd = 0.75, 
                   color = fcase(!ensembl & eff_chr == exp_chr, "#4575B4", 
                                 !ensembl & eff_chr != exp_chr, "#ABD9E9"))])

dev.off()


## ---- plots circos plots barplot ----

.d <- results$hq_predictions_23Q2_fdr_0.5[, .(gene1 = effect_gene, 
                                              gene2 = expression_gene, ensembl)] %>% 
  merge(gene_database[, .(gene1 = symbol, chr1 = chromosome, start1 = start, 
                          end1 = end)], by ="gene1", all.x = T, sort = F) %>% 
  merge(gene_database[, .(gene2 = symbol, chr2 = chromosome, start2 = start, 
                          end2 = end)], by = "gene2", all.x = T, sort = F)

.d <- .d[, .(type = fcase(
  ensembl & chr1 == chr2, "Syntenic & paralogs", 
  ensembl & chr1 != chr2, "Non-syntenic, paralogs", 
  !ensembl & chr1 == chr2, "Syntenic, no paralogs", 
  !ensembl & chr1 != chr2, "Non-syntenic, no paralogs"))]
.d <- .d[!is.na(type)]

.d[, .N, by = type] %>% ggplot(aes(N, type)) + 
  geom_col(aes(fill = type)) + 
  geom_text(aes(x = N + 20, label = N)) + 
  scale_fill_manual(values = c(
    "Syntenic & paralogs" = "#D73027", 
    "Non-syntenic, paralogs" = "#FDAE61", 
    "Syntenic, no paralogs" = "#4575B4", 
    "Non-syntenic, no paralogs" = "#ABD9E9")) + 
  guides(fill = guide_legend(ncol = 1)) + 
  theme(legend.position = "none", legend.title = element_blank(), 
        axis.title.y = element_blank()) + 
  labs(x = "# gene pairs")

plotsaver(.fpath, "fig5g_hq_set_bars", 4, 2)


