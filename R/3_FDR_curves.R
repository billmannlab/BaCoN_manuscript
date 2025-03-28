.top <- 1000

mkdir(file.path(fig_directory, "reviewer_response_FDR_curves"))

.d <- list(fig2 = prepare_predictions_for_plotting(
  results$IDs_23Q2$benchmarking, 
  results$all_predictions_23Q2, .top)$predictions[order(id, rank)], 
  fig3 = prepare_predictions_for_plotting(
    results$IDs_23Q2$benchmarking_with_bacon, 
    results$all_predictions_23Q2, .top)$predictions[order(id, rank)])



for (.n in names(.d)) {
  .d[[.n]][, `:=`(FDR = prediction_FDR(
    exp_gene_col = expression_gene, 
    eff_gene_col = effect_gene, 
    TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | 
                thompson | proximity == "self-addiction"), 
    save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | 
                  thompson | proximity == "self-addiction"), 
    proximity_col = (proximity == "same arm, <10mbp"), 
    coex_z_score = 3, verbose = F), 
    FDR_alternative = prediction_FDR(
      exp_gene_col = expression_gene, 
      eff_gene_col = effect_gene, 
      TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | thompson), 
      save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | 
                    thompson | proximity == "self-addiction"), 
      proximity_col = (proximity == "same arm, <10mbp"), 
      coex_z_score = 3, verbose = F), 
    FDR_all_non_FP = prediction_FDR(
      exp_gene_col = expression_gene, 
      eff_gene_col = effect_gene, 
      save_col = (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | 
                    thompson | proximity == "self-addiction"), 
      proximity_col = (proximity == "same arm, <10mbp"), 
      coex_z_score = 3, verbose = F)), by = id]}


# Fig 2

.linecol <- "#707070"

.captions <- str_c(
  "Top ", .top, " predictions\n", 
  c("Default FDR, computed using all but DeKegel paralogs as well as self-addictions as TP", 
    "FDR computed using all but DeKegel paralog pairs as TP", 
    "FDR computed using all non-FP as TP"))




.plt_elements <- list(
  geom_vline(xintercept = c(100, 1000), color = .linecol), 
  geom_line(linewidth = 1.25), 
  geom_hline(yintercept = c(0, 0.1, 0.2, 0.25, 0.5, 1), color = .linecol), 
  theme(legend.position = "bottom", strip.background = element_blank()))



.d$fig2 %>% 
  ggplot(aes(rank, FDR, color = name)) + 
  labs(color = "", caption = .captions[[1]]) + 
  .plt_elements

plotsaver(file.path(fig_directory, "reviewer_response_FDR_curves"), 
          str_c("fig2_sets_FPR_paralogs_and_selfs_top_", .top), 9, 6)

.d$fig2 %>% 
  ggplot(aes(rank, FDR_alternative, color = name)) + 
  labs(color = "", caption = .captions[[2]]) + 
  .plt_elements

plotsaver(file.path(fig_directory, "reviewer_response_FDR_curves"), 
          str_c("fig2_sets_FPR_paralogs_top_", .top), 9, 6)


.d$fig2 %>% 
  ggplot(aes(rank, FDR_all_non_FP, color = name)) + 
  labs(color = "", caption = .captions[[3]]) + 
  .plt_elements

plotsaver(file.path(fig_directory, "reviewer_response_FDR_curves"), 
          str_c("fig2_sets_FPR_all_non_FP_top_", .top), 9, 6)



.plt_elements <- list(
  geom_vline(xintercept = c(100, 1000), color = .linecol), 
  geom_line(linewidth = 1.25), 
  geom_hline(yintercept = c(0, 0.1, 0.2, 0.25, 0.5, 1), color = .linecol), 
  theme(legend.position = "bottom", strip.background = element_blank()))



.d$fig3 %>% 
  ggplot(aes(rank, FDR, color = name)) + 
  labs(color = "", caption = .captions[[1]]) + 
  .plt_elements

plotsaver(file.path(fig_directory, "reviewer_response_FDR_curves"), 
          str_c("fig3_sets_FPR_paralogs_and_selfs_top_", .top), 9, 6)

.d$fig3 %>% 
  ggplot(aes(rank, FDR_alternative, color = name)) + 
  labs(color = "", caption = .captions[[2]]) + 
  .plt_elements

plotsaver(file.path(fig_directory, "reviewer_response_FDR_curves"), 
          str_c("fig3_sets_FPR_paralogs_top_", .top), 9, 6)


.d$fig3 %>% 
  ggplot(aes(rank, FDR_all_non_FP, color = name)) + 
  labs(color = "", caption = .captions[[3]]) + 
  .plt_elements

plotsaver(file.path(fig_directory, "reviewer_response_FDR_curves"), 
          str_c("fig3_sets_FPR_all_non_FP_top_", .top), 9, 6)

