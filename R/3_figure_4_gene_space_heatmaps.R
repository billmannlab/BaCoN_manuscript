## ---- figure 4 plots heatmaps ----

# Heatmaps

.top <- 100
.fpath <- file.path(fig_directory, "fig4_dynamic_gene_space_heatmaps")
mkdir(.fpath)

.d <- list(stats = prepare_predictions_for_plotting(
  results$IDs_23Q2$dynamic_gene_space, 
  results$all_predictions_23Q2, .top, "stats", F)$stats %>% 
  merge(dynamic_range_combinations[, .(id, n_exp_genes, n_chr_genes)], 
        by.x = "setup", by.y = "id"))

.d$methods <- results$IDs_23Q2$dynamic_gene_space[, unique(name)]

.d$universe <- empty_array(
  list(rev(levels(dynamic_range_combinations$n_exp_genes)), 
       rev(levels(dynamic_range_combinations$n_chr_genes)), 
       c("universe_size", 
         str_c(rep(names(paralog_pairs), 2), "_", 
               rep(c("hits", "density"), each = length(paralog_pairs))))))


.d$universe[,,"universe_size"] <- distinct(gene_pair_set_stats[
  , .(gene_space, n_pairs)])[dynamic_range_combinations[, id], 
                             on = "gene_space", get("n_pairs")]


for (.n in names(paralog_pairs)) {
  .d$universe[,,str_c(.n, "_hits")] <- gene_pair_set_stats[gene_pair_set == .n][
    dynamic_range_combinations[, id], on = "gene_space", get("hits")]
  .d$universe[,,str_c(.n, "_density")] <- gene_pair_set_stats[gene_pair_set == .n][
    dynamic_range_combinations[, id], on = "gene_space", get("density")]
}

.d$hms <- empty_array(list(
  rev(levels(dynamic_range_combinations$n_exp_genes)), 
  rev(levels(dynamic_range_combinations$n_chr_genes)), 
  .d$methods, 
  c(#"n_genes", "diversity_index", 
    str_c("n_", names(paralog_pairs)), 
    str_c(names(paralog_pairs), "_fc"))
  
))

for (.m in .d$methods) {
  for (.n in names(paralog_pairs)) {
    .d$hms[,,.m,str_c("n_", .n)] <- .d$stats[name == .m][
      dynamic_range_combinations[, id], on = "setup", get(str_c("n_", .n))]
    .d$hms[,,.m,str_c(.n, "_fc")] <- .d$hms[,,.m,str_c("n_", .n)] / 
      (.d$universe[,,str_c(.n, "_density")]*.top)
  }
}

.d$hms %>% apply(3:4, na_perc)

# Plot heatmaps:

ha_row <- HeatmapAnnotation(
  "n exp genes" = anno_barplot(
    dynamic_range_combinations[, as.numeric(rev(levels(n_exp_genes)))], 
    add_numbers = F, 
    border = F, 
    axis = F, 
    axis_param = list(direction = "reverse")), 
  width = unit(1, "cm"), which = "row")

ha_col <- HeatmapAnnotation(
  "n chr genes" = anno_barplot(
    dynamic_range_combinations[, as.numeric(rev(levels(n_chr_genes)))], 
    add_numbers = F, 
    border = F, 
    axis = F), 
  height = unit(1, "cm"), which = "col")

shared_params <- list(cluster_rows = F, cluster_columns = F, 
                      rect_gp = gpar(col = "white", lwd = 1), 
                      column_names_gp = gpar(fontsize = 7), 
                      row_names_gp = gpar(fontsize = 7), 
                      column_title_gp = gpar(fontsize = 6))

hm_matrices <- list(
  method_comparison = sapply(
    dimnames(.d$hms)[[4]], 
    \(.v) {Reduce(cbind, lapply(.d$methods, \(i) .d$hms[,,i,.v]))}, simplify = F), 
  universe_stats = sapply(dimnames(.d$universe)[[3]], 
                          \(.v) {.d$universe[,,.v]}, simplify = F))

## Methods_comparison

pdf(file.path(fig_directory, "figure_4", "fig4fg_methods_comparison.pdf"), 
    width = 15, height = 3)

for (.n in names(hm_matrices$method_comparison)) {
  
  message(.n)
  
  .arr <- hm_matrices$method_comparison[[.n]]
  
  .col <- "#560A0A"
  if (grepl("_fc$", .n)) {.col <- "#CE4500"}
  colpalette <- circlize::colorRamp2(breaks = c(min(.arr, na.rm = T), 
                                                max(.arr, na.rm = T)), 
                                     c("white", .col))
  
  hm_plt <- do.call(Heatmap, c(list(.arr), 
                               list(row_title = .n, name = .n, 
                                    column_split = rep(.d$methods, each = 12), 
                                    col = colpalette), 
                               shared_params))
  draw(hm_plt)
}
dev.off()

## Universe stats:
                                                                                                                                                                                                                                                                                                                              

pdf(file.path(fig_directory, "figure_4", "fig4cde_universe_stats.pdf"), 
    width = 5, height = 3.5)

for (.n in names(hm_matrices$universe_stats)) {
  
  .arr <- hm_matrices$universe_stats[[.n]]
  colpalette <- circlize::colorRamp2(breaks = c(min(.arr, na.rm = T), 
                                                max(.arr, na.rm = T)), 
                                     c("white", "#003958"))
  
  hm_plt <- do.call(Heatmap, c(list(.arr), 
                               list(row_title = .n, name = .n, 
                                    left_annotation = ha_row, 
                                    top_annotation = ha_col, 
                                    col = colpalette#, 
                                    #left_annotation = ha_row, 
                                    #top_annotation = ha_col
                               ), shared_params))
  draw(hm_plt)
}

dev.off()

hm_matrices$universe_stats[c("universe_size", "ensembl_hits", "ohnologs_hits")]

