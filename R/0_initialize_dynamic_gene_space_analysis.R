
## ---- initialize dynamic gene space parameters ----

expression_ths <- c(10, 30, 60, 100, 300, 600, 900, 1000)
chronos_ths <- c(seq(0.2, 0.5, 0.05), c(0.6, 0.8, 1, 1.2, 1.5))

dynamic_range_combinations <- list(expression_gene_sets = data.table(
  exp_th = expression_ths, 
  exp_genes = sapply(expression_ths, \(.e_th) {
    names(which(apply(gene_expression_complete[cl_subset_default,] >= 3, 
                      2, sum, na.rm = T) >= .e_th))}, simplify = F)), 
  chronos_gene_sets = data.table(
    chr_th = chronos_ths, 
    chr_genes = sapply(
      chronos_ths, 
      \(.c_th) {names(
        which(apply(abs(chronos_complete) > .c_th, 2, sum) >= 30))}, 
      simplify = F)))

dynamic_range_combinations$all_combinations <- data.table(
  expand.grid(exp_th = expression_ths, chr_th = chronos_ths))
dynamic_range_combinations$all_combinations[, `:=`(
  id = str_c("dynamic", "exp_th", exp_th, "chr_th", chr_th, sep = "_"))]


dynamic_range_combinations <- dynamic_range_combinations$all_combinations %>%
  merge(dynamic_range_combinations$expression_gene_sets, by = "exp_th", all.x = T) %>%
  merge(dynamic_range_combinations$chronos_gene_sets, by = "chr_th", all.x = T)


dynamic_range_combinations[, `:=`(n_exp_genes = sapply(exp_genes, l), 
                                  n_chr_genes = sapply(chr_genes, l))]

dynamic_range_combinations[, `:=`(
  n_exp_genes = factor(n_exp_genes, levels = sort(unique(n_exp_genes))), 
  n_chr_genes = factor(n_chr_genes, levels = sort(unique(n_chr_genes))))]