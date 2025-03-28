## ---- figure 4 plots cell line subsamples ----

.top <- 100

.d <- prepare_predictions_for_plotting(
    results$IDs_23Q2$random_cell_line_subsets, 
    results$all_predictions_23Q2, .top, c("stats"))


.d$stats[setup == "default", n_cell_lines := 1019]
.d$stats[setup %in% names(cl_subsets$random), 
         n_cell_lines := as.numeric(gsub("^r", "", str_split_i(setup, "\\.", 1)))]

.d$stats[, n_neighbors := as.numeric(n_neighbors)]

.d$stats2 <- .d$stats %>% melt.data.table(measure.vars = c("FDR", "n_neighbors"))

.d$stats2[, .(val_mean = mean(value, na.rm = T), 
             val_sd = sd(value, na.rm = T)), 
         by = .(name, n_cell_lines, variable)] %>% 
  ggplot(aes(val_mean, name)) + 
  geom_col(aes(fill = name)) + 
  geom_errorbar(aes(xmin = val_mean - val_sd, 
                    xmax = val_mean + val_sd), width = 0.3) +
  facet_grid(n_cell_lines ~ variable, scales = "free_x", switch = "y") + 
  guides(fill = guide_legend(ncol = 2)) + 
  theme(legend.position = "top", 
        strip.placement = "outside", 
        legend.title = element_blank(), 
        strip.background = element_blank())

plotsaver(file.path(fig_directory, "figure_4"), 
          str_c("fig4_FDR_proximity_top_", .top), 9, 11)


.d$stats <- .d$stats %>% 
  melt.data.table(measure.vars = str_c("n_", names(paralog_set_titles)), 
                  value.name = "n_paralogs", variable.name = "paralog_set")

.d$stats[, `:=`(paralog_set = paralog_set_titles[paralog_set])]

.d$stats[, .(mean_paralogs = mean(n_paralogs, na.rm = T), 
             sd_paralogs = sd(n_paralogs, na.rm = T), .N), 
         by = .(name, n_cell_lines, paralog_set)] %>% 
  ggplot(aes(mean_paralogs, name, fill = name, group = n_cell_lines)) + 
  geom_col() + 
  geom_errorbar(aes(xmin = mean_paralogs - sd_paralogs, 
                    xmax = mean_paralogs + sd_paralogs), width = 0.3) + 
  facet_grid(n_cell_lines ~ paralog_set, scales = "free_x", switch = "y") + 
  labs(x = "# paralog pairs among top100 buffering predictions", 
       y = "# cell lines") + 
    guides(fill = guide_legend(ncol = 2)) + 
  theme(legend.position = "top", 
        strip.placement = "outside", 
        legend.title = element_blank(), 
        strip.background = element_blank())

plotsaver(file.path(fig_directory, "figure_4"), 
          str_c("fig4ab_all_paralog_sets_top_", .top), 12, 11)


.fpath <- file.path(fig_directory, "supplementary_random_cell_line_samples")

for (.n_cls in c(100, 200, 500)) {
  .d <- prepare_predictions_for_plotting(
    results$IDs_23Q2$random_cell_line_subsets[grepl(str_c("r", .n_cls), id)], 
                                         results$all_predictions_23Q2, 100, 
                                         "network")
  
  shared_elements <- list(
    network_plt_theme(), 
    geom_segment(aes(xend = xend, yend = yend), linewidth = 0.1), 
    geom_point(size = 0.1),
    geom_point(aes(xend, yend), size = 0.1#, data = .d$nodes
    ), 
    coord_fixed(), 
    facet_grid(fct_rev(name) ~ setup, switch = "y"), 
    theme(legend.title = element_blank(), 
          plot.margin = unit(rep(0.001, 4), "cm"), 
          panel.spacing = unit(0.001, "cm"), 
          strip.text.y.left = element_text(angle = 0, hjust = 0), 
          strip.text.x = element_blank()))
  
  .d$network %>% 
    ggplot(aes(x, y, color = fcase(ensembl, "Ensembl paralog", 
                                   !ensembl, "Non-paralog"))) + 
    shared_elements + 
    scale_color_manual(values = c("Ensembl paralog" = paralog_colors[1], 
                                  "Non-paralog" = paralog_colors[2])) + 
  plotsaver(file.path(fig_directory, "figure_4"), 
            str_c("fig4c_networks_", .n_cls, "_cell_lines_paralogs"), 9, 10)
  
  
  .d$network %>% 
    ggplot(aes(x, y, color = z_score_bin)) + 
    shared_elements + 
    scale_color_manual(values = z_score_colors)
  
    plotsaver(file.path(fig_directory, "figure_4"), 
              str_c("fig4c_networks_", .n_cls, "_cell_lines_z_scores"), 9, 10)
  
  
}

# With all cell lines:

.d <- prepare_predictions_for_plotting(
  results$IDs_23Q2$random_cell_line_subsets[grepl("23Q2_default", id)], 
  results$all_predictions_23Q2, 100, 
  "network")

.d$network %>% 
  ggplot(aes(x, y, color = fcase(ensembl, "Ensembl paralog", 
                                 !ensembl, "Non-paralog"))) + 
  shared_elements + 
  scale_color_manual(values = c("Ensembl paralog" = paralog_colors[1], 
                                "Non-paralog" = paralog_colors[2])) + 
  plotsaver(file.path(fig_directory, "figure_4"), 
            str_c("fig4c_networks_", "all", "_cell_lines_paralogs"), 4, 10)


.d$network %>% 
  ggplot(aes(x, y, color = z_score_bin)) + 
  shared_elements + 
  scale_color_manual(values = z_score_colors)

  plotsaver(file.path(fig_directory, "figure_4"), 
            str_c("fig4c_networks_", "all", "_cell_lines_z_scores"), 4, 10)
