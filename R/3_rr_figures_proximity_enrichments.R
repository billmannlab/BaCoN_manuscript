
## ---- reviewer response proximity enrichments ----

.d <- list(preds = prepare_predictions_for_plotting(
  results$IDs_23Q2$benchmarking_with_bacon[c(1, 2, 5, 6, 13, 14)], 
  results$all_predictions_23Q2, 1000, "predictions")$predictions)

for (.id in .d$preds[, unique(id)]) {
  
  for (.i in 1:5) {
    .th1 <- c(0, 1, 5, 10, 15, 20)[.i]
    .th2 <- c(0, 1, 5, 10, 15, 20)[.i+1]
    .n <- str_c("same arm, ", .th1, "-", .th2, "mbp")
    
    .d$enrichments$with_coex[[str_c(.id, "_", .n)]] <- data.table(
      id = .id, 
      name = .d$preds[id == .id, unique(name)], 
      proximity = .n, 
      draws = .d$preds[id == .id & 
                         expression_gene != effect_gene, .N], # no selfs
      whites_drawn = .d$preds[id == .id & 
                                expression_gene != effect_gene & # no selfs
                                eff_chr_arm == exp_chr_arm & 
                                distance %between% c(.th1*1e6, .th2*1e6), .N], 
      
      whites_in_urn = gene_pair_set_stats[
        gene_space == "default" & 
          gene_pair_set == str_c("proximity", .th1, .th2, "mbp", sep = "_"), get("hits")], 
      urn_size = gene_pair_set_stats[
        gene_space == "default" & 
          gene_pair_set == str_c("proximity", .th1, .th2, "mbp", sep = "_"), get("n_pairs")])
  }
  
  .n <- str_c("same arm, >20mbp")
  
  .d$enrichments$with_coex[[str_c(.id, "_", .n)]] <- data.table(
    id = .id, 
    name = .d$preds[id == .id, unique(name)], 
    proximity = .n, 
    draws = .d$preds[id == .id & 
                       expression_gene != effect_gene, .N], # no selfs
    whites_drawn = .d$preds[id == .id & 
                              expression_gene != effect_gene & # no selfs
                              eff_chr_arm == exp_chr_arm & 
                              distance >= 20*1e6, .N], 
    whites_in_urn = gene_pair_set_stats[
      gene_space == "default" & gene_pair_set == "proximity_l20_mbp", get("hits")], 
    urn_size = gene_pair_set_stats[
      gene_space == "default" & gene_pair_set == "proximity_l20_mbp", get("n_pairs")])
  
  
  .n <- str_c("different arm")
  
  .d$enrichments$with_coex[[str_c(.id, "_", .n)]] <- data.table(
    id = .id, 
    name = .d$preds[id == .id, unique(name)], 
    proximity = .n, 
    draws = .d$preds[id == .id & 
                       expression_gene != effect_gene,.N], # no selfs
    whites_drawn = .d$preds[id == .id & 
                              expression_gene != effect_gene & # no selfs
                              eff_chr == exp_chr & eff_chr_arm != exp_chr_arm, # same chr, diff arm
                            .N], 
    whites_in_urn = gene_pair_set_stats[
      gene_space == "default" & gene_pair_set == "different_arm", get("hits")], 
    urn_size = gene_pair_set_stats[
      gene_space == "default" & gene_pair_set == "different_arm", get("n_pairs")])}



for (.id in .d$preds[, unique(id)]) {
  
  for (.i in 1:5) {
    .th1 <- c(0, 1, 5, 10, 15, 20)[.i]
    .th2 <- c(0, 1, 5, 10, 15, 20)[.i+1]
    .n <- str_c("same arm, ", .th1, "-", .th2, "mbp")
    
    .d$enrichments$no_coex[[str_c(.id, "_", .n)]] <- data.table(
      id = .id, 
      name = .d$preds[id == .id, unique(name)], 
      proximity = .n, 
      draws = .d$preds[
        id == .id & expression_gene != effect_gene & # no selfs
          !sorted_pair %in% important_pairs_default$coexpressed_z3, .N], # no coex predictions
      whites_drawn = .d$preds[
        id == .id & 
          expression_gene != effect_gene & # no selfs
          !sorted_pair %in% important_pairs_default$coexpressed_z3 & # no coex predictions
          eff_chr_arm == exp_chr_arm & 
          distance %between% c(.th1*1e6, .th2*1e6), .N], 
      
      whites_in_urn = gene_pair_set_stats[
        gene_space == "default" & 
          gene_pair_set == str_c("proximity", .th1, .th2, "mbp", 
                                 "no_coex_z3", sep = "_"), get("hits")], 
      
      urn_size = gene_pair_set_stats[
        gene_space == "default" & 
          gene_pair_set == str_c("proximity", .th1, .th2, "mbp", 
                                 "no_coex_z3", sep = "_"), get("n_pairs")])
    
  }
  
  .n <- str_c("same arm, >20mbp")
  
  .d$enrichments$no_coex[[str_c(.id, "_", .n)]] <- data.table(
    id = .id, 
    name = .d$preds[id == .id, unique(name)], 
    proximity = .n, 
    draws = .d$preds[id == .id & 
                       expression_gene != effect_gene & # no selfs
                       !sorted_pair %in% important_pairs_default$coexpressed_z3, .N], # no coex predictions
    whites_drawn = .d$preds[id == .id & 
                              expression_gene != effect_gene & # no selfs
                              !sorted_pair %in% important_pairs_default$coexpressed_z3 & # no coex predictions
                              eff_chr_arm == exp_chr_arm & 
                              distance >= 20*1e6, .N], 
    whites_in_urn = gene_pair_set_stats[
      gene_space == "default" & 
        gene_pair_set == "proximity_l20_mbp_no_coex_z3", get("hits")], 
    urn_size = gene_pair_set_stats[
      gene_space == "default" & gene_pair_set == "proximity_l20_mbp_no_coex_z3", get("n_pairs")])
  
  
  .n <- str_c("different arm")
  
  .d$enrichments$no_coex[[str_c(.id, "_", .n)]] <- data.table(
    id = .id, 
    name = .d$preds[id == .id, unique(name)], 
    proximity = .n, 
    draws = .d$preds[
      id == .id & 
        expression_gene != effect_gene & # no selfs
        !sorted_pair %in% important_pairs_default$coexpressed_z3, .N], # no coex predictions
    whites_drawn = .d$preds[
      id == .id & 
        expression_gene != effect_gene & # no selfs
        !sorted_pair %in% important_pairs_default$coexpressed_z3 & # no coex predictions
        eff_chr == exp_chr & eff_chr_arm != exp_chr_arm, # same chr, diff arm
      .N], 
    whites_in_urn = gene_pair_set_stats[
      gene_space == "default" & gene_pair_set == "different_arm", get("hits")], 
    urn_size = gene_pair_set_stats[
      gene_space == "default" & gene_pair_set == "different_arm", get("n_pairs")])}



.d$enrichments$with_coex <- rbindlist(.d$enrichments$with_coex)
.d$enrichments$no_coex <- rbindlist(.d$enrichments$no_coex)

.d$enrichments$with_coex[, facet := "with coexpressed pairs"]
.d$enrichments$no_coex[, facet := "no coexpressed pairs"]

.d$enrichments <- rbindlist(.d$enrichments)

.d$enrichments[, FC := (whites_drawn / draws) / (whites_in_urn / urn_size)]
.d$enrichments[whites_drawn != 0, `:=`(LFC = log2(FC), lower.tail = log2(FC) <= 0)]
# white balls are syntenic pairs among predictions, black balls are non-syntenic ones
.d$enrichments[, pval := phyper(
  q = whites_drawn - 1, # white balls drawn
  m = whites_in_urn, # white balls in urn
  k = draws - whites_drawn, # black balls drawn
  n = (urn_size - whites_in_urn), # black balls in urn
  # lower.tail = T tests for depletion, lower.tail = F tests for enrichment
  lower.tail = lower.tail, 
), by = .I][, sign := fcase(pval %between% c(0.01, 0.05), "*", 
                            pval %between% c(0.001, 0.01), "**", 
                            pval <= 0.001, "***", default = "")]

.d$enrichments[, proximity := factor(
  proximity, levels = c("same arm, 0-1mbp", "same arm, 1-5mbp", 
                        "same arm, 5-10mbp", "same arm, 10-15mbp", 
                        "same arm, 15-20mbp", "same arm, >20mbp", 
                        "different arm", "different chromosome"))]

for (.n in .d$enrichments[, unique(facet)]) {
  .d$enrichments[facet == .n] %>% 
    ggplot(aes(proximity, -log10(pval))) + 
    geom_col(position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") + 
    facet_grid(name ~ .) + 
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 40, hjust = 1), 
          strip.text.y = element_text(angle = 0), 
          strip.background = element_blank())
  
  plotsaver(file.path(fig_directory, "reviewer_response_proximity_enrichments"), 
            str_c("proximity_enrichments_", .n), 6, 8)
  
}
