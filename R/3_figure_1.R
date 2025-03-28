## ---- figure 1 plots ----

mkdir(file.path(fig_directory, "figure_1"))

# Coexpression example CDS1 - CDS2

.d <- data.table(cds1_exp = gene_expression[,"CDS1"], 
                 cds2_eff = chronos[,"CDS2"])
.cor <- .d[, str_c("PCC = ", round(cor(cds1_exp, cds2_eff, use = pco), 2))]
.d %>% 
  ggplot(aes(cds1_exp, cds2_eff)) + 
  geom_point(shape = 21, color = "#121420", fill = "darkgrey", 
                            stroke = 0.2, size = 2) + 
  geom_smooth(method = "lm", color = "black", se = F) + 
  geom_text(aes(4.5, -0.65, label = .cor)) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, NA)) + 
  labs(x = "CDS1 expression\n[log2 TPM+1]", y = "CDS2 fitness\n[chronos]")
  

plotsaver(file.path(fig_directory, "figure_1"), "fig1b", 3, 3)


# buffering concept PCC Histogram

.d <- data.table(data.frame(
  readRDS(file.path(result_dir, "sample_scores", 
                    str_c(.dm, "_default_CvE_PCC_BaCoN_0.05_sample_scores.rds")))))
.d <- .d[sample(1:.N, 100000), .(PCC = reference)]

.top_100_pcc <- 0.5202841 # sort(CvE_PCC, decreasing = T)[100]
.n_pairs <- 147 # round(length(CvE_PCC)/1e6, 0)

pcc_histogram <- .d %>% 
  ggplot() + 
  geom_vline(xintercept = .top_100_pcc, linetype = "dashed") + 
  geom_histogram(aes(PCC), alpha = 0.75, bins = 100) + 
  scale_x_continuous(limits = c(-0.2, 0.6)) + 
  labs(x = str_c("PCC fitness ~ expression\n(", .n_pairs, 
                 " M gene pairs)")) + 
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank())

plotsaver(file.path(fig_directory, "figure_1"), "fig1c", 4, 2.5)

# co-expression bias heatmaps

heatmap_colors <- c("white", "#003958")

# tissue types represented with at least 10 cell lines
.lineages <- covariates[cell_line %in% rownames(gene_expression), .N, 
                        by = lineage][N > 10, lineage]

bias_coexpression <- list()

for (.b_gene in c("CDS2", "CHMP4B")) { 

    #find genes with interacting with biased genes in top 100 predictions
  .d <- results$all_predictions_23Q2[id == str_c("23Q2_default_CvE_PCC")][
    order(score, pcc_reference, decreasing = T)][1:100]
  .genes <- .d[effect_gene == .b_gene, unique(expression_gene)]
  
  .array <- empty_array(list(.lineages, .genes), 0)
  
  for (.l in .lineages) {
    for (.g in .genes) {
      cl_oi <- intersect(rownames(gene_expression), 
                         covariates[lineage == .l, cell_line])
      .array[.l,.g] <- mean(gene_expression[cl_oi,.g], na.rm = T)}}
  
  to_rename <- c("Ovary_Fallopian_Tube" = "Ovary", 
                 "Bladder_Urinary_Tract" = "Bladder urinary", 
                 "CNS_Brain" = "CNS Brain", 
                 "Soft_Tissue" = "Soft Tissue", 
                 "Biliary_Tract" = "Biliary Tract", 
                 "Esophagus_Stomach" = "Esoph. Stomach", 
                 "Head_and_Neck" = "Head & neck")
  
  for (.t in names(to_rename)) {
    rownames(.array)[rownames(.array) == .t] <- to_rename[.t]}
  
  .hm <- Heatmap(
    matrix = .array, 
    name = "mRNA expression [log2 TPM+1]", 
    col = circlize::colorRamp2(breaks = c(0, max(.array)), heatmap_colors), 
    row_title = "Cell line lineage", 
    row_title_side = "right", 
    column_title = str_c("Predicted ", .b_gene, 
                         " buffering partners\nin top 100 (uncorrected)"), 
    column_title_side = "bottom", 
    rect_gp = gpar(col = "white", lwd = 0.5), 
    heatmap_legend_param = list(direction = "horizontal"))
  
  bias_coexpression[[.b_gene]] <- list(
    lineages = .lineages, 
    genes = .genes, 
    array = .array, 
    heatmap = .hm)}


pwidth <- c("CDS2" = 7, "CHMP4B" = 5)

for (.b_gene in c("CDS2", "CHMP4B")) {# save heatmaps
  png(file.path(file.path(fig_directory, "figure_1"), str_c("fig1e_", .b_gene, ".png")), 
      pwidth[.b_gene], 5.5, units = "in", res = 300)
  draw(bias_coexpression[[.b_gene]]$heatmap, 
       heatmap_legend_side = "bottom")
  dev.off()
  
  pdf(file.path(file.path(fig_directory, "figure_1"), str_c("fig1e_", .b_gene, ".pdf")), 
      pwidth[.b_gene], 5.5)
  draw(bias_coexpression[[.b_gene]]$heatmap, 
       heatmap_legend_side = "bottom")
  dev.off()}

# co-expression bias histograms

.d <- prepare_predictions_for_plotting(results$IDs_23Q2$benchmarking[1], 
                                       results$all_predictions_23Q2, 100)$predictions

.d <- list(
  all_top_100 = .d[, expression_gene],
  CDS2 = .d[effect_gene == "CDS2", expression_gene], 
  CHMP4B = .d[effect_gene == "CHMP4B", expression_gene], 
  not_CDS2_CHMP4B = .d[!effect_gene %in% c("CDS2", "CHMP4B"), 
                       expression_gene])


.d <- sapply(names(.d), \(.n) {
  goi <- .d[[.n]]
  cormat <- cor(gene_expression[cl_subset_default,
                                intersect(goi, expression_genes)], use = pco)
  cormat[which(lower.tri(cormat, diag = T))] <- NA
  x <- data.table(which(!is.na(cormat), arr.ind = T), cormat[!is.na(cormat)])
  x <- x[, .(gene1 = rownames(cormat)[row], 
             gene2 = colnames(cormat)[col], 
             cor_val = V2, 
             set = .n)]}, simplify = F) %>% Reduce(bind_rows, .)

.d[, .(set = fcase(set == "all_top_100", "All top 100", 
                   set == "not_CDS2_CHMP4B", "no CDS2, CHMP4B", 
                   default = set), cor_val)] %>%
  ggplot(aes(cor_val)) + 
  geom_histogram(bins = 80) + 
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "dotted") + 
  labs(
    x = "Co-expression\namong buffering predictions") + 
  facet_wrap(. ~ set, ncol = 1, scales = "free_y") + 
  theme(legend.position = "none", 
        strip.background = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

plotsaver(file.path(fig_directory, "figure_1"), "fig1f", 2.5, 5)


.top <- 5000

.d <- prepare_predictions_for_plotting(
  results$IDs_23Q2$benchmarking[1], 
  results$all_predictions_23Q2, .top, c("predictions", "network"))

.d$predictions[, coexpressed := F][
  effect_gene %in% unique(effect_gene[duplicated(effect_gene)]), 
  coexpressed := flag_coex(expression_gene), by = effect_gene]
.d$network[, coexpressed := F][
  effect_gene %in% unique(effect_gene[duplicated(effect_gene)]), 
  coexpressed := flag_coex(expression_gene), by = effect_gene]


.plt_elements <- list(network_plt_theme(), 
                      geom_segment(aes(xend = xend, yend = yend)), 
                      geom_point(aes(x, y)), 
                      geom_point(aes(xend, yend)), 
                      coord_fixed(), 
                      labs(color = ""))

# Coloring by Ensembl

.d$network %>% 
  ggplot(aes(x, y, color = fcase(ensembl, "Ensembl paralog", !ensembl, "Non-paralog"))) + 
  .plt_elements + 
  scale_color_manual(values = c("Ensembl paralog" = paralog_colors[1], 
                                "Non-paralog" = paralog_colors[2]))

plotsaver(file.path(fig_directory, "figure_1"), "fig1d_PCC_ensembl_network", 3, 3)

# Coling by Co-expression

.d$predictions %>% 
  ggplot(aes(rank, cumsum(coexpressed))) + 
  geom_line(color = "#DD1C1A") + 
  geom_abline(slope = 1) + 
  labs(x = "Prediction rank", y = "# co-expressed", 
       caption = .d$predictions[
         1:5000, 
         str_c(sum(coexpressed), 
               " of ", .top, " predictions are coexpressed (", 
               sum(coexpressed)/.N*100, "%)")])

plotsaver(file.path(fig_directory, "figure_1"), "fig1_PCC_co_expressed_line", 3, 3)


.d$network %>% 
  ggplot(aes(x, y, color = fcase(coexpressed, "Co-expressed", 
                                 !coexpressed, "Not co-expressed"))) + 
  .plt_elements + 
  scale_color_manual(values = c("Co-expressed" = "#DD1C1A", 
                                "Not co-expressed" = "darkgrey"))

plotsaver(file.path(fig_directory, "figure_1"),
          "fig1_PCC_co_expressed_network", 3, 3)

