## ---- generate hq set ----

.d <- prepare_predictions_for_plotting(
  results$IDs_23Q2$hq_set, 
  results$all_predictions_23Q2, 2000)$predictions[order(rank)]

.d[, coexpressed := F][
  effect_gene %in% unique(effect_gene[duplicated(effect_gene)]), 
  coexpressed := flag_coex(expression_gene), by = effect_gene]

.d[coexpressed & (ensembl | ohnologs | anvar | anvar_standard | ito | gemini | 
                    thompson | proximity == "self-addiction"), coexpressed := F]

for (.fdr in 0.5) {
  results[[str_c("hq_predictions_23Q2_fdr_", .fdr)]] <- .d[
    FDR <= .fdr & !coexpressed & expression_gene != effect_gene & !proximity]
}

