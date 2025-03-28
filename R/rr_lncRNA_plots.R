## ---- reviewer response ncRNA plots ----


.d <- list(hq_predictions = results$hq_predictions_23Q2_fdr_0.5, 
           ncRNAs = prepare_predictions_for_plotting(
             results$IDs_24Q2$ncRNAs, 
             results$all_predictions_24Q2, 5000, "predictions", 
             include_FDR = F)$predictions)

.d$ncRNAs <- .d$ncRNAs[score >= .d$hq_predictions[, tail(score, 1)]]

.d$ncRNAs <- .d$ncRNAs[, .SD, 
                       .SDcols = setdiff(
                         names(.d$ncRNAs), 
                         c("exp_chr", "exp_pos", "exp_chr_arm", "eff_chr", 
                           "eff_pos", "eff_chr_arm", "distance", 
                           "proximity", "proximity_class"))] %>% 
  merge(exp_gene_info_all, by = "expression_gene", all.x = T, sort = F) %>% 
  merge(eff_gene_info_all, by = "effect_gene", all.x = T, sort = F)

.d$ncRNAs[exp_chr_arm == eff_chr_arm, distance := abs(exp_pos - eff_pos)]
.d$ncRNAs[, `:=`(proximity = exp_chr_arm == eff_chr_arm & distance <= proximity_threshold)]


.d$ncRNAs <- .d$ncRNAs[expression_class == "lncRNA" & 
                         exp_chr %in% chromosomes & eff_chr %in% chromosomes]

saveRDS(.d$ncRNAs, file.path(cache, 
                             "23Q2_data_export", str_c(.d$ncRNAs[, .N], "_lncRNAs.rds")))

.d$chrom_dist <- bind_rows(
  .d$ncRNAs[, .(facet = "buffering ncRNAs", n_pairs = .N), by = .(chr = exp_chr)], 
  .d$ncRNAs[, .(facet = "buffered genes", n_pairs = .N), by = .(chr = eff_chr)])

.d$chrom_dist[!is.na(chr)] %>% 
  ggplot(aes(factor(chr, levels = c(1:22, "X", "Y")), n_pairs)) + 
  geom_col(aes(fill = facet), position = position_dodge2(preserve = "single")) + 
  scale_fill_manual(values = c("buffered genes" = "#782E93", 
                               "buffering ncRNAs" = "#16A567")) + 
  labs(x = "chromosome", y = "# pairs") + 
  theme(legend.position = c(0.25, 0.75), legend.title = element_blank())

plotsaver(file.path(fig_directory, "reviewer_response_ncRNAs"), 
          "ncRNA_chromosomal_distribution", 5, 3)


.d$ncRNAs[!is.na(eff_chr) & !is.na(exp_chr), .N, 
          by = .(type = fcase(exp_chr == eff_chr, "Syntenic", 
                              exp_chr != eff_chr, "Non-syntenic"))] %>% 
  ggplot(aes(N, type)) + 
  geom_col(aes(fill = type)) + 
  geom_text(aes(x = N + 5, label = N)) + 
  scale_fill_manual(values = c(
    "Syntenic" = "#4575B4", "Non-syntenic" = "#ABD9E9")) + 
  guides(fill = guide_legend(ncol = 1)) + 
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        axis.title.y = element_blank()) + 
  labs(x = "# gene pairs")

plotsaver(file.path(fig_directory, "reviewer_response_ncRNAs"), 
          "ncRNA_circos_bars", 4, 2)

pdf(file.path(fig_directory, "reviewer_response_ncRNAs", "ncRNA_circos_plots.pdf"), 
    width = 3.5, height = 3.5)

draw_circos_plot(prepare_circos_data_from_predictions(
  .d$ncRNAs[!is.na(eff_chr) & !is.na(exp_chr)]), 
  .d$ncRNAs[!is.na(eff_chr) & !is.na(exp_chr), 
            .(lwd = 0.75, 
              color = fcase(eff_chr == exp_chr, "#4575B4", 
                            eff_chr != exp_chr, "#ABD9E9"))], 
  add_labels = F)

dev.off()
