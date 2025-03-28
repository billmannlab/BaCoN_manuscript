## ---- figure 2 plots ----

.pheight <- 6

mkdir(file.path(fig_directory, "figure_2"))

.top <- 1000

.d <- prepare_predictions_for_plotting(
  results$IDs_23Q2$benchmarking, 
  results$all_predictions_23Q2, 
  .top, 
  list("predictions", "coexpression", "stats"))

.d$paralogs <- .d$stats %>% 
  melt.data.table(measure.vars = str_c("n_", names(paralog_set_titles)), 
                  variable.name = "paralog_set", 
                  value.name = "n_paralogs", variable.factor = F, 
                  value.factor = F)

.d$paralogs[, paralog_set := paralog_set_titles[gsub("n_", "", get("paralog_set"))]]

for (.n in names(paralog_set_titles)) {
  .d$paralogs[paralog_set == paralog_set_titles[[.n]], 
              paralog_fc := n_paralogs / 
                (gene_pair_set_stats[gene_space == "default" & 
                                       gene_pair_set == .n, density*.top])]
  }

colorscale <- scale_fill_manual(values = c("TRUE" = "navy", "FALSE" = "#434343"))


shared_elements <- list(
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(size = 10, angle = 90, hjust = 0)))

# Paralog pairs bar plot

.d$paralogs %>%
  ggplot(aes(n_paralogs, name)) + 
  geom_col(aes(fill = name == "PCC + BaCoN")) + 
  geom_vline(aes(xintercept = n_paralogs), 
             .d$paralogs[name == "PCC + BaCoN"], 
             color = "darkgrey", linetype = "dotted") + 
  shared_elements + colorscale +
  labs(x = "# paralog pairs among top100 buffering predictions", 
       y = "A priori correction\n+ PCC") + 
  facet_wrap(. ~ paralog_set, nrow = 1, scales = "free_x")

plotsaver(file.path(fig_directory, "figure_2"), 
          str_c("fig2a_paralogs_top_", .top), 6, .pheight)

# Paralog FC bar plot

.d$paralogs %>%
  ggplot(aes(paralog_fc/1000, name)) + 
  geom_col(aes(fill = name == "PCC + BaCoN")) + 
  geom_vline(aes(xintercept = paralog_fc/1000), 
             .d$paralogs[name == "PCC + BaCoN"], 
             color = "darkgrey", linetype = "dotted") + 
  shared_elements + colorscale + 
  scale_x_continuous(position = "top") + 
  labs(x = "Enrichment [foldchange x1000]", 
       y = "A priori correction\n+ PCC") + 
  facet_wrap(. ~ paralog_set, nrow = 1, scales = "free_x")

plotsaver(file.path(fig_directory, "figure_2"), 
          str_c("fig2a_paralog_fc_top_", .top), 
          6, .pheight)


# Other prediction statistics

.d$stats[, .(name, div_ind_0.5 = as.numeric(div_ind_0.5), 
             n_neighbors = as.numeric(n_neighbors), 
             n_self_addiction = as.numeric(n_self_addiction))] %>% 
  melt.data.table(measure.vars = c(
    "div_ind_0.5", "n_neighbors", "n_self_addiction")) %>%
  ggplot(aes(value, name)) + 
  geom_col(aes(fill = name == "PCC + BaCoN")) + 
  colorscale + 
  facet_grid(. ~ fcase(variable == "div_ind_0.5", "Div. index [...]", 
                       variable == "n_neighbors", "# Neighbors", 
                       variable == "n_self_addiction", "# Self-addictions"), 
             scales = "free_x", switch = "x") + 
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(angle = 0, hjust = 0.5))

plotsaver(file.path(fig_directory, "figure_2"), str_c("fig2bde_top_", .top), 5, 4)

# Coexpression histogram

.d$coexpression %>% 
  ggplot(aes(coex, fill = name == "PCC + BaCoN")) + 
  geom_histogram(bins = 80) + 
  colorscale + 
  geom_vline(xintercept = 0) + 
  facet_grid(fct_rev(name) ~ ., scales = "free_y") + 
  scale_x_continuous(limits = c(-1, 1)) + 
  labs(x = "Co-Expression") + coex_hist_theme() + 
  theme(legend.position = "none")

plotsaver(file.path(fig_directory, "figure_2"), 
          str_c("fig2c_", .top), 3, .pheight)

# Z-score histogram

.d$predictions %>% 
  ggplot(aes(pcc_reference, fill = z_score_bin)) + 
  geom_histogram(bins = 80) + 
  scale_fill_manual(values = z_score_colors) + 
  geom_vline(xintercept = pcc_z_scores, color = "darkgrey") + 
  geom_vline(xintercept = 0) + 
  facet_grid(fct_rev(name) ~ ., scales = "free_y") + 
  labs(x = "PCC z-score", fill = "PCC z-score") + 
  theme(legend.position = "top", 
        legend.direction = "vertical", 
        legend.key.size = unit(5, "mm")) + 
  coex_hist_theme()

plotsaver(file.path(fig_directory, "figure_2"), 
          str_c("fig2f_top_", .top), 2, .pheight)
