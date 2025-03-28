
## ---- prepare plotting ----

paralog_colors <- c("#D02A26", "#007EA7")

z_score_colors <- c("< 0" = "#5D0718", "0-1" = "#6E0D25", 
                    "1-2" = "#EA2B1F", "2-3" = "#FFAD05", 
                    "> 3" = "#25B06F")


paralog_set_titles <- factor(c("ensembl" = "Ensembl", 
                               "ohnologs" = "Ohnologs", 
                               "dekegel" = "DeKegel et al. SL", 
                               "dekegel_val" = "DeKegel et al. validated SL", 
                               "anvar" = "Anvar et al. SL", 
                               "anvar_standard" = "Anvar et al. SL standard", 
                               "ito" = "Ito et al. SL", 
                               "thompson" = "Thompson et al. SL"), 
                             levels = c("Ensembl", "Ohnologs", 
                                        "DeKegel et al. SL", 
                                        "DeKegel et al. validated SL", 
                                        "Anvar et al. SL", 
                                        "Anvar et al. SL standard", 
                                        "Ito et al. SL", 
                                        "Thompson et al. SL"))

pcc_z_scores <- sapply(
  0:3, \(.z) {mean(CvE_PCC, na.rm = T) + .z * sd(CvE_PCC, na.rm = T)}) %>% 
  cacheR(str_c(.dm, "_PCC_z_scores"), cache)



## ---- diversity_index function ----

diversity_index <- \(set1, set2, .i = 0.5) {
  if (base::length(set1) != base::length(set2)) {
    warning("Set size is not equal!")
  }
  x <- cumsum(base::sort(table(c(set1, set2)), decreasing = T))
  max_possible <- sum(base::length(set1), base::length(set2)) * .i
  return(base::length(x[which(x <= max_possible)]) / max_possible)}



## ---- FDR function ----

prediction_FDR <- \(exp_gene_col, 
                    eff_gene_col, 
                    TP_col, 
                    save_col, 
                    proximity_col, 
                    coex_z_score = 3, 
                    return_everything = F, 
                    verbose = T, 
                    show_progress = T) {
  
  if (verbose) {
    if (missing(TP_col)) {
      message("No values provided to use as TP. Using all that are not FP.")}}
  
  output <- data.table(exp_genes = exp_gene_col, 
                       eff_genes = eff_gene_col)
  
  if (missing(save_col)) {warning("No values to save provided!")} else {
    output[, save := save_col]}
  if (missing(proximity_col)) {warning("No proximity pairs provided.")} else {
    output[, FP_proximity := proximity_col]}
  
  output[, FP := F]
  duplicated_effect_genes <- output[, unique(eff_genes[duplicated(eff_genes)])]
  
  if (show_progress) {pb <- progbar(l(duplicated_effect_genes))}
  
  for (.g in duplicated_effect_genes) {
    if (show_progress) {pb$tick()}
    .i <- output[, .I[eff_genes == .g]]
    .exp_genes_in_cluster <- output[.i, exp_genes]
    
    .coex_m <- coexpression_default_z[.exp_genes_in_cluster,.exp_genes_in_cluster]
    .coex_m[lower.tri(.coex_m, diag = T)] <- NA
    
    .coex_m <- .coex_m >= coex_z_score
    
    output[.i, FP := colSums(.coex_m, na.rm = T) != 0]
  }
  
  output[, FP := FP | FP_proximity]
  
  output[(save & FP), FP := F]
  
  if (missing(TP_col)) {
    output[, TP := F]
    output[!(FP), TP := T]
  } else {
    output[, TP := TP_col]
  }
  
  #if (verbose) {message("Computing FDR...")}
  
  output[, `:=`(FP_cs = cumsum(FP), FP_TP_cs = cumsum(FP + TP))]
  output[, `:=`(FDR = FP_cs / FP_TP_cs, FPR = FP_cs / .I)]
  
  if (return_everything) {return(output)
  } else {return(output[, FDR])}}


## ---- prepare_predictions_for_plotting function ----

prepare_predictions_for_plotting <- \(ids, predictions, top, 
                                      type = list("predictions"), 
                                      include_FDR = T, 
                                      verbose = T) {
  
  if (verbose) {
    if (!all(ids[, id] %in% predictions[, id])) {
      message("No predictions found for: ")
      message(str_c(setdiff(ids[, id], predictions[, id]), collapse = ", \n"))
    }}
  
  output <- list()
  
  x <- ids %>% 
    merge(predictions[rank %in% 1:top], by = "id")
  
  if (include_FDR) {
    x[, FDR := prediction_FDR(
      exp_gene_col = expression_gene, 
      eff_gene_col = effect_gene, 
      TP_col = (ensembl | ohnologs | anvar | anvar_standard | ito | 
                  thompson | expression_gene == effect_gene), 
      save_col = (ensembl | ohnologs | anvar | anvar_standard | 
                    ito | thompson | expression_gene == effect_gene), 
      proximity_col = (sorted_pair %in% important_pairs_default$proximity), 
      coex_z_score = 3), by = .(id, name, setup)]
  }
  
  
  if ("predictions" %in% type) {output$predictions <- x}
  
  if ("stats" %in% type) {
    
    output$stats <- x[, .(
      n_predictions = .N, 
      n_genes = length(unique(c(expression_gene, effect_gene))), 
      n_ensembl = sum(ensembl), 
      n_ohnologs = sum(ohnologs), 
      n_dekegel = sum(dekegel), 
      n_dekegel_val = sum(dekegel_val), 
      n_anvar = sum(anvar), 
      n_anvar_standard = sum(anvar_standard),
      n_ito = sum(ito), 
      n_thompson = sum(thompson), 
      n_self_addiction = sum(expression_gene == effect_gene), 
      n_neighbors = sum(proximity, na.rm = T), 
      div_ind_0.5 = diversity_index(expression_gene, effect_gene, 0.5), 
      FDR = if (include_FDR) tail(FDR, 1) else NA_real_), 
      by = .(id, version, name, setup, datasets, pre_hoc, association, 
             post_hoc, BaCoN_th)]
  }
  
  if ("coexpression" %in% type) {
    
    output$coexpression <- sapply(
      x[, unique(id)], 
      \(.id) {
        goi <- x[id == .id, intersect(unique(expression_gene), expression_genes)]
        cormat <- cor(gene_expression[cl_subset_default,goi], use = pco)
        cormat[which(lower.tri(cormat, diag = T))] <- NA
        xx <- data.table(which(!is.na(cormat), arr.ind = T), cormat[!is.na(cormat)])
        xx <- xx[, .(gene1 = rownames(cormat)[row], 
                     gene2 = colnames(cormat)[col], coex = V2, 
                     id = .id)]}, simplify = F) %>% 
      rbindlist() %>% 
      merge(distinct(x[, .(id, name, version, setup, datasets, 
                           pre_hoc, association, post_hoc)]), by = "id", all.x = T)
  }
  
  if ("network" %in% type) {
    
    output$network <- sapply(x[, unique(id)], \(.id) {
      
      set.seed(12345)
      
      .edges <- copy(x[id == .id, .(expression_gene, effect_gene, score)])
      .g <- graph_from_data_frame(.edges, directed = F)
      .layout <- layout_with_fr(.g)
      .nodes <- data.table(node = V(.g)$name, x = .layout[,1], y = .layout[,2])
      .edges <- .edges %>% 
        merge(.nodes, by.x = c("effect_gene"), by.y = c("node"), all.x = T) %>% 
        merge(.nodes, by.x = c("expression_gene"), 
              by.y = c("node"), all.x = T, suffixes = c("end", "")) %>%
        merge(x[id == .id], by = c("expression_gene", "effect_gene", "score"))
      .edges}, simplify = F) %>% rbindlist()}
  return(output)
}

network_plt_theme <- \(.legend_pos = "bottom") {
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        strip.text = element_text(size = 6), 
        strip.background = element_blank(), 
        plot.caption = element_text(hjust = 0.5), 
        legend.position = .legend_pos)
}



.x_axis.text.size <- 5

coex_hist_theme <- \(.hide_y_axis = T) {
  .t <- c(list(theme(strip.background = element_blank(), 
                     strip.text.y.right = element_blank())), 
          ifelse(.hide_y_axis, 
                 list(theme(axis.text.x = element_text(size = .x_axis.text.size), 
                            axis.text.y = element_blank(), 
                            axis.title.y = element_blank(), 
                            axis.ticks.y = element_blank(), 
                            axis.line.y = element_blank())), list()))
  .t}


prepare_circos_data_from_predictions <- \(.data, 
                                          .show_all_chromosomes = F, 
                                          remove_unknown_chromosomes = T) {
  
  
  .chromosome_list <- c(1:22, "X", "Y")
  
  if (.data[, any(!exp_chr %in% .chromosome_list) | any(!eff_chr %in% .chromosome_list)]) {
    warning("Remove predictions with unknown chromosomes!")
    
    
    if (remove_unknown_chromosomes) {
      .data <- .data[exp_chr %in% .chromosome_list & eff_chr %in% .chromosome_list]
    }
  }
  .cols_needed <- c("effect_gene", "expression_gene", 
                    "exp_chr", "exp_pos", "eff_chr", "eff_pos")
  if (!all(.cols_needed) %in% names(.data)) {
    warning(str_c("Add column names : ", .cols_meeded[!.cols_needed %in% names(.data)], 
                  collapse = ", "))
  }
  #.data <- .data[exp_chr %in% .chromosome_list & eff_chr %in% .chromosome_list]
  
  .circos <- list(#predictions = .data
  )
  
  if(.show_all_chromosomes) {
    .circos$chr_oi <- str_c("chr", .chromosome_list)
  } else {
    .circos$chr_oi <- str_c(
      "chr", .chromosome_list[.chromosome_list %in% .data[, union(eff_chr, exp_chr)]])}
  
  hg38 <- BSgenome.Hsapiens.UCSC.hg38
  seqinf <- seqinfo(hg38)[.circos$chr_oi]
  gR_hg38 <- makeGRangesFromDataFrame(GRanges(seqinf), seqinfo = seqinf)
  
  .circos$links_1 <- .data[, .(chr = str_c("chr", eff_chr), 
                               start = eff_pos - 10000, 
                               end = eff_pos + 10000, 
                               gene = effect_gene)]
  .circos$links_2 <- .data[, .(chr = str_c("chr", exp_chr), 
                               start = exp_pos - 10000, 
                               end = exp_pos + 10000, 
                               gene = expression_gene)]
  
  .circos$labels <- distinct(bind_rows(.circos$links_1, .circos$links_2), 
                             .keep_all = T)
  return(.circos)
}

draw_circos_plot <- \(.data, lineinfo, .title = NA, 
                      add_genomicIdeogram = T, 
                      add_ideogram_coords = T, 
                      add_labels = T, add_links = T) {
  
  if (add_ideogram_coords) {
    .plttype <- c("axis", "labels")
  } else {
    .plttype <- c("labels")}
  
  circos.par("start.degree" = 90, 
             "gap.degree" = rep(2, length(.data$chr_oi)), 
             track.margin = c(0, 0))
  
  circos.initializeWithIdeogram(species = "hg38", # scale, axis
                                chromosome.index = .data$chr_oi, 
                                plotType = .plttype, 
                                track.height = 0.1, major.by = 1e8)
  
  if (add_genomicIdeogram) {circos.genomicIdeogram(track.height = 0.03)}
  
  if (add_labels) {
    
    circos.genomicLabels(.data$labels, 
                         labels.column = 4, connection_height = 0.03, 
                         side = "inside", 
                         cex = 0.5, 
                         labels_height = 0.25)} # gene pair labels
  
  if (add_links) {circos.genomicLink(region1 = .data$links_1, 
                                     region2 = .data$links_2, 
                                     lwd = lineinfo$lwd, 
                                     col = lineinfo$col)} # SL links
  
  if (!missing(.title)) {title(.ttl)}
  circos.clear()
}

flag_coex <- \(exp_genes) {
  coex_m <- coexpression_default_z[exp_genes,exp_genes]
  coex_m[lower.tri(coex_m, diag = T)] <- NA
  return(colSums(coex_m >= 3, na.rm = T) != 0)
}
