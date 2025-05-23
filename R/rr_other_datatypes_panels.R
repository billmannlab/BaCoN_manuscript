
.top <- 1000
.pheight <- 6

fig_data <- prepare_predictions_for_plotting(
  results$IDs_24Q2$other_datatypes, 
  results$all_predictions_24Q2, 
  .top, 
  list("predictions", "coexpression", "stats"))

fig_data$paralogs <- fig_data$stats %>% 
  melt.data.table(measure.vars = str_c("n_", names(paralog_set_titles)), 
                  variable.name = "paralog_set", value.name = "n_paralogs", 
                  variable.factor = F, value.factor = F)

fig_data$paralogs[, paralog_set := paralog_set_titles[gsub("n_", "", get("paralog_set"))]]

for (.n in names(paralog_set_titles)) {
  fig_data$paralogs[paralog_set == paralog_set_titles[[.n]], paralog_fc := n_paralogs / 
                      (gene_pair_set_stats[gene_space == "default" & gene_pair_set == .n, density*.top])]
}

fillscale <- scale_fill_manual(values = c("TRUE" = "navy", "FALSE" = "#434343"))

shared_elements <- list(
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(size = 10, angle = 90, hjust = 0))
)


pdf(file.path(fig_directory, "figure_3", str_c("all_paralog_sets_top_", .top, ".pdf")), 9, .pheight)

fig_data$paralogs %>%
  ggplot(aes(n_paralogs, name)) + 
  geom_col(aes(fill = grepl("BaCoN", name))) + 
  geom_vline(aes(xintercept = n_paralogs), 
             fig_data$paralogs[name == "PCC + BaCoN"], 
             color = "darkgrey", linetype = "dotted") + 
  shared_elements + 
  fillscale + 
  labs(x = "# paralog pairs among top100 buffering predictions", 
       y = "A priori correction + PCC") + 
  facet_wrap(. ~ paralog_set, nrow = 1, scales = "free_x")


plotsaver(file.path(fig_directory, "reviewer_response_other_datatypes"), 
          str_c("all_paralog_sets_paralog_hits_top_", .top), 9, .pheight)

fig_data$paralogs %>%
  ggplot(aes(paralog_fc/1000, name)) + 
  geom_col(aes(fill = grepl("BaCoN", name))) + 
  geom_vline(aes(xintercept = paralog_fc/1000), 
             fig_data$paralogs[name == "PCC + BaCoN"], 
             color = "darkgrey", linetype = "dotted") + 
  shared_elements + 
  fillscale + 
  scale_x_continuous(position = "top", 
                     #labels = \(x) format(x, scientific = T)#scales::comma
  ) + 
  labs(x = "Enrichment [foldchange x1000]", 
       y = "A priori correction + PCC") + 
  facet_wrap(. ~ paralog_set, nrow = 1, scales = "free_x")


plotsaver(file.path(fig_directory, "reviewer_response_other_datatypes"), 
          str_c("all_paralog_sets_paralog_fc_top_", .top), 9, .pheight)



dev.off()


fig_data$stats[, .(name, div_ind_0.5 = as.numeric(div_ind_0.5), 
                   n_neighbors = as.numeric(n_neighbors), 
                   n_self_addiction = as.numeric(n_self_addiction))] %>% 
  melt.data.table(measure.vars = c("div_ind_0.5", "n_neighbors", "n_self_addiction")) %>%
  ggplot(aes(value, name)) + 
  geom_col(aes(fill = grepl("BaCoN", name))) + 
  #shared_elements + 
  fillscale + 
  #labs() + 
  facet_grid(. ~ fcase(
    variable == "div_ind_0.5", "Div. index [...]", 
    variable == "n_neighbors", "# Neighbors", 
    variable == "n_self_addiction", "# Self-addictions"), scales = "free_x", switch = "x") + 
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(angle = 0, hjust = 0.5))

plotsaver(file.path(fig_directory, "reviewer_response_other_datatypes"), 
          str_c("div_ind_neighbors_selfs_top_", .top), 5, 4)




fig_data$coexpression %>% 
  ggplot(aes(coex, fill = grepl("BaCoN", name))) + 
  geom_histogram(bins = 80) + 
  fillscale + 
  geom_vline(xintercept = 0) + 
  facet_grid(fct_rev(name) ~ ., scales = "free_y") + 
  scale_x_continuous(limits = c(-1, 1)) + 
  labs(x = "Co-Expression") + coex_hist_theme()

plotsaver(file.path(fig_directory, "reviewer_response_other_datatypes"), 
          str_c("coex_histogram_top_", .top), 3, .pheight)


fig_data$z_scores <- fig_data$stats %>% 
  melt.data.table(measure.vars = grep("z_sc", names(fig_data$stats), value = T), 
                  variable.name = "z_score_bin")


fig_data$z_scores[, z_score_bin := factor(
  z_score_bin, 
  levels = rev(c("z_sc_l0", 
                 "z_sc_0_1", 
                 "z_sc_1_2", 
                 "z_sc_2_3", "z_sc_gr3")), 
  labels = rev(c("< 0", "0-1", "1-2", "2-3", "> 3")))]

fig_data$z_scores %>% 
  ggplot(aes(value, name, fill = z_score_bin)) + 
  geom_col() + 
  scale_fill_manual(values = z_score_colors) + 
  labs(x = str_c("% of top", .top, " predictions"), 
       fill = "PCC z-score") + 
  guides(fill = guide_legend(ncol = 1)) + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = "top", legend.direction = "vertical", 
        legend.key.size = unit(5, "mm"))

plotsaver(file.path(fig_directory, "reviewer_response_other_datatypes"), 
          str_c("z_scores_top_", .top), 2, .pheight)

