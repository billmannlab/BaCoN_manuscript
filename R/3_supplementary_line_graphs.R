
.top <- 5000

.d <- rbind(results$IDs_23Q2$benchmarking[12], 
            results$IDs_23Q2$benchmarking_with_bacon)

.d <- sapply(.d[, id], \(.id) {
  x <- prepare_predictions_for_plotting(
    .d[id == .id], results$all_predictions_23Q2, .top)$predictions
  
  x[, coexpressed := F][
    effect_gene %in% unique(effect_gene[duplicated(effect_gene)]), 
    coexpressed := flag_coex(expression_gene), by = .(id, name, effect_gene)]
}, simplify = F) %>% rbindlist()



.dx <- .d[, .(rank, ensembl, anvar, 
              z_l0 = z_score <= 0, 
              z_l1 = z_score <= 1, 
              z_l2 = z_score <= 2, 
              z_l3 = z_score <= 3, 
              self_association, 
              proximity, 
              coexpressed, 
              hq = z_score >= 3 & !ensembl & !coexpressed & !proximity), 
          by = .(id, name)] %>% 
  melt.data.table(
    measure.vars = c("ensembl", "anvar", str_c("z_l", 0:3), "self_association", 
                     "coexpressed", "proximity", "hq"), 
    variable.factor = F, value.factor = F)

.dx[, variable := c(
  "ensembl" = "Ensembl paralogs", 
  "anvar" = "Anvar et al. SL", 
  "self_association" = "Self-association", 
  "coexpressed" = "Coexpressed", 
  "proximity" = "Proximity", 
  "z_l0" = "z-score < 0", 
  "z_l1" = "z-score < 1", 
  "z_l2" = "z-score < 2", 
  "z_l3" = "z-score < 3", 
  "hq" = "z-score > 3 &\nnot Ensembl paralog &\nnot proximal &\nnot coexpressed")[variable]]

.colscale <- scale_color_manual(
  values = c(
    "Ensembl paralogs" = "#007EA7", 
    "Anvar et al. SL" = "#004E89", 
    "Self-association" = "#073B3A", 
    "Coexpressed" = "#F0386B", 
    "Proximity" = "#7F2982", 
    "z-score < 0" = z_score_colors[["< 0"]], 
    "z-score < 1" = z_score_colors[["0-1"]], 
    "z-score < 2" = z_score_colors[["1-2"]], 
    "z-score < 3" = z_score_colors[["2-3"]], 
    str_c("z-score > 3 &", 
    "not Ensembl paralog &", "not proximal &", "not coexpressed", 
    sep = "\n") = z_score_colors[["> 3"]]
  ))



.dx[, .(rank, value = cumsum(value)), by = .(id, name, variable)] %>% 
  ggplot(aes(rank, value, color = variable)) + 
  geom_line() + 
  .colscale + 
  geom_vline(xintercept = c(100, 1000), linetype = "dotted") + 
  geom_abline(slope = 1, linetype = "dotted") + 
  facet_wrap(. ~ fct_rev(name), scales = "free", ncol = 5) + 
  labs(x = "Top buffering predictions (gene pairs)", 
       y = "# Gene pairs in sub-category") + 
  #guides(color = guide_legend(ncol = 1)) + 
  theme(
        legend.position = "bottom", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 5), 
        axis.text.y = element_text(size = 5), 
        strip.text = element_text(size = 7), 
        strip.background = element_blank())

plotsaver(file.path(fig_directory, "supplementary_figures"), 
          str_c("lines_top_", .top), 9, 7)

