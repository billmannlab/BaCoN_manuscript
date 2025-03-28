
## ---- important gene pairs ----

.fpath <- file.path(cache, str_c(.dm, "_important_gene_pairs"))
mkdir(.fpath)

important_pairs_default <- list()
important_pairs_complete <- list()

all_default_genes <- union(expression_genes, chronos_genes)
all_complete_genes <- union(colnames(gene_expression_complete), 
                            colnames(chronos_complete))

# Coexpressed genes

coexpression_default_z <- {
  x <- cor(gene_expression, use = pco)
  x <- x %>% `diag<-`(NA)
  x <- (x - mean(x, na.rm = T)) / sd(x, na.rm = T)} %>% 
  cacheR(str_c(.dm, "_coexpression_default_z_normalized"), 
         file.path(cache, str_c(.dm, "_coexpression")))


important_pairs_default$coexpressed_z3 <- {
  x <- coexpression_default_z %>% 
    melt_array_to_dt("gene1", "gene2", "coex_z")
  
  x <- x[coex_z >= 3]
  x[, sorted_pair := sort_gene_pairs(gene1, gene2)]
  x[gene1 != gene2, get("sorted_pair")]} %>% 
  cacheR("coexpressed_z3", .fpath)



important_pairs_default$ensembl_ident_20_50 <- paralog_pairs$ensembl[
  ident %between% c(20, 50), get("sorted_pair")]
important_pairs_default$ensembl_ident_50_100 <- paralog_pairs$ensembl[
  ident %between% c(20, 100), get("sorted_pair")]


# pairs in chromosomal proximity

armwise_gene_combinations <- list()

message("Same arm...")

armwise_gene_combinations$same_arm_default <- rbindlist(
  sapply(chromosome_arms, \(.arm) {
  .gdb <- gene_database[chromosome_arm == .arm]
  .all_genes <- gene_database[chromosome_arm == .arm, 
                              intersect(symbol, all_default_genes)]
  .all_combinations <- data.table(arm = .arm, 
                                  expand.grid(gene1 = .all_genes, 
                                              gene2 = .all_genes, 
                                              stringsAsFactors = F))  
  .all_combinations <- .all_combinations %>% 
    merge(.gdb[, .(gene1 = symbol, pos1 = position)], by = "gene1", all.x = T) %>%
    merge(.gdb[, .(gene2 = symbol, pos2 = position)], by = "gene2", all.x = T)
  
  .all_combinations[, `:=`(distance = abs(pos1 - pos2), 
                           sorted_pair = sort_gene_pairs(gene1, gene2))]
  
  .all_combinations <- distinct(
    .all_combinations[gene1 != gene2, .(sorted_pair, arm, distance)])
  .all_combinations}, simplify = F)) %>%
  cacheR("armwise_gene_combinations_same_arm_default", .fpath)


armwise_gene_combinations$same_arm_complete <- rbindlist(
  sapply(chromosome_arms, \(.arm) {
    .gdb <- gene_database[chromosome_arm == .arm]
    .all_genes <- gene_database[chromosome_arm == .arm, 
                                intersect(symbol, all_complete_genes)]
    .all_combinations <- data.table(
      arm = .arm, 
      expand.grid(gene1 = .all_genes, gene2 = .all_genes, stringsAsFactors = F))  
    .all_combinations <- .all_combinations %>% 
      merge(.gdb[, .(gene1 = symbol, pos1 = position)], by = "gene1", all.x = T) %>%
      merge(.gdb[, .(gene2 = symbol, pos2 = position)], by = "gene2", all.x = T)
    
  .all_combinations[, `:=`(distance = abs(pos1 - pos2), 
                           sorted_pair = sort_gene_pairs(gene1, gene2))]
  
  .all_combinations <- distinct(.all_combinations[gene1 != gene2, 
                                                  .(sorted_pair, arm, distance)])
  .all_combinations}, simplify = F)) %>%
  cacheR("armwise_gene_combinations_same_arm_complete", .fpath)


# cutoff at 10 Mbp

proximity_threshold <- 1e7

important_pairs_default$proximity <- armwise_gene_combinations$same_arm_default[
  distance <= proximity_threshold, get("sorted_pair")] %>%
  cacheR(str_c("proximity", proximity_threshold, "mbp", "default", sep = "_"), .fpath)

important_pairs_complete$proximity <- armwise_gene_combinations$same_arm_complete[
  distance <= proximity_threshold, get("sorted_pair")] %>%
  cacheR(str_c("proximity", proximity_threshold, "mbp", "complete", sep = "_"), .fpath)


message(ifelse(l(important_pairs_default$proximity) == 711272, 
               "Proximity standard correct.", "Proximity standard changed!"))


.ths <- c(0, 1, 5, 10, 15, 20)

for (.i in 1:5) {
  .th1 <- .ths[.i]
  .th2 <- .ths[.i+1]
  .x_default <- armwise_gene_combinations$same_arm_default[
    distance %between% c(.th1*1e6, .th2*1e6), get("sorted_pair")]
  .x_complete <- armwise_gene_combinations$same_arm_complete[
    distance %between% c(.th1*1e6, .th2*1e6), get("sorted_pair")]
  
  important_pairs_default[[str_c("proximity", .th1, .th2, "mbp", sep = "_")]] <- .x_default
  important_pairs_complete[[str_c("proximity", .th1, .th2, "mbp", sep = "_")]] <- .x_complete

  # to understand wether coexpression drives the proximity signal, 
  # we also want values for the proximity sets without coexpression
  
  #if (.n != "coexpressed_z3") {
    important_pairs_default[[
      str_c("proximity", .th1, .th2, "mbp", "no_coex_z3", sep = "_")
      ]] <- setdiff(.x_default, important_pairs_default$coexpressed_z3)
}
#}

.x_default <- armwise_gene_combinations$same_arm_default[distance >= 2e7, 
                                                         get("sorted_pair")]
important_pairs_default$proximity_l20_mbp <- .x_default

#if (.n != "coexpressed_z3") {
  important_pairs_default$proximity_l20_mbp_no_coex_z3 <- setdiff(
    .x_default, important_pairs_default$coexpressed_z3)
#}

important_pairs_complete$proximity_l20_mbp <- 
  armwise_gene_combinations$same_arm_complete[distance >= 2e7, get("sorted_pair")]
# largest distance: 147024160


message("Same chromosome...")

important_pairs_default$same_chromosome <- sapply(chromosomes[1:24], \(.chr) {
  .gdb <- gene_database[chromosome == .chr]
  .all_genes <- gene_database[chromosome == .chr, intersect(symbol, all_default_genes)]
  .all_combinations <- data.table(chromosome = .chr, 
                                  expand.grid(gene1 = .all_genes, gene2 = .all_genes, 
                                              stringsAsFactors = F))
  .all_combinations <- .all_combinations[gene1 != gene2] %>% 
    merge(.gdb[, .(gene1 = symbol, pos1 = position)], by = "gene1", all.x = T) %>%
    merge(.gdb[, .(gene2 = symbol, pos2 = position)], by = "gene2", all.x = T)
  
  .all_combinations[, sorted_pair := sort_gene_pairs(gene1, gene2)]
  
  .all_combinations[, get("sorted_pair")]}, simplify = F) %>% 
  Reduce(c, .) %>% 
  unique() %>% 
  cacheR("same_chromosome_default", .fpath)

important_pairs_default$syntenic_ensembl <- intersect(
  important_pairs_default$same_chromosome, 
  paralog_pairs$ensembl[, get("sorted_pair")])

important_pairs_default$syntenic_ohnologs <- intersect(
  important_pairs_default$same_chromosome, 
  paralog_pairs$ohnologs[, get("sorted_pair")])

important_pairs_default$syntenic_no_ensembl <- setdiff(
  important_pairs_default$same_chromosome, 
  paralog_pairs$ensembl[, get("sorted_pair")])

important_pairs_default$syntenic_no_ensembl_no_ohnologs <- setdiff(
  important_pairs_default$same_chromosome, 
  union(paralog_pairs$ensembl[, get("sorted_pair")], 
        paralog_pairs$ohnologs[, get("sorted_pair")]))


important_pairs_complete$same_chromosome <- sapply(chromosomes[1:24], \(.chr) {
  .gdb <- gene_database[chromosome == .chr]
  .all_genes <- gene_database[chromosome == .chr, 
                              intersect(symbol, all_complete_genes)]
  .all_combinations <- data.table(
    chromosome = .chr, 
    expand.grid(gene1 = .all_genes, gene2 = .all_genes, stringsAsFactors = F))
  .all_combinations <- .all_combinations[gene1 != gene2] %>% 
    merge(.gdb[, .(gene1 = symbol, pos1 = position)], by = "gene1", all.x = T) %>%
    merge(.gdb[, .(gene2 = symbol, pos2 = position)], by = "gene2", all.x = T)
  
  .all_combinations[, sorted_pair := sort_gene_pairs(gene1, gene2)]
  
  .all_combinations[, get("sorted_pair")]}, simplify = F) %>% 
  Reduce(c, .) %>% 
  unique() %>% 
  cacheR("same_chromosome_complete", .fpath)

message("Different arm...")


important_pairs_default$different_arm <- sapply(chromosomes, \(.chr) {
  .arm1 <- str_c("chr", .chr, "p")
  .arm2 <- str_c("chr", .chr, "q")
  .gdb <- gene_database[chromosome_arm == .arm1]
  .all_genes1 <- gene_database[chromosome_arm == .arm1, 
                               intersect(symbol, all_default_genes)]
  .all_genes2 <- gene_database[chromosome_arm == .arm2, 
                               intersect(symbol, all_default_genes)]
  
  .all_combinations1 <- data.table(same_arm = F, arm = .arm1, 
                                   expand.grid(gene1 = .all_genes1, 
                                               gene2 = .all_genes2, 
                                               stringsAsFactors = F))
  # in contrast to same-arm predictions, we include these predictions twice, 
  # as each gene on arm 1 can be the buffering or buffered partner to arm 2
  
  .arm1 <- str_c("chr", .chr, "q")
  .arm2 <- str_c("chr", .chr, "p")
  .gdb <- gene_database[chromosome_arm == .arm1]
  .all_genes1 <- gene_database[chromosome_arm == .arm1, 
                               intersect(symbol, all_default_genes)]
  .all_genes2 <- gene_database[chromosome_arm == .arm2, 
                               intersect(symbol, all_default_genes)]
  
  .all_combinations2 <- data.table(same_arm = F, arm = .arm1, 
                                   expand.grid(gene1 = .all_genes1, 
                                               gene2 = .all_genes2, 
                                               stringsAsFactors = F))
  
  .all_combinations <- rbindlist(list(.all_combinations1, .all_combinations2))
  
  .all_combinations <- .all_combinations[gene1 != gene2] %>% 
    merge(.gdb[, .(gene1 = symbol, pos1 = position)], by = "gene1", all.x = T) %>%
    merge(.gdb[, .(gene2 = symbol, pos2 = position)], by = "gene2", all.x = T)
  
  .all_combinations[, sorted_pair := sort_gene_pairs(gene1, gene2)]
  .all_combinations[, get("sorted_pair")]}, simplify = F) %>% 
  Reduce(c, .) %>% 
  unique() %>% 
  cacheR("different_arm_default", .fpath)




important_pairs_complete$different_arm <- sapply(chromosomes, \(.chr) {
  .arm1 <- str_c("chr", .chr, "p")
  .arm2 <- str_c("chr", .chr, "q")
  .gdb <- gene_database[chromosome_arm == .arm1]
  .all_genes1 <- gene_database[chromosome_arm == .arm1, 
                               intersect(symbol, all_complete_genes)]
  .all_genes2 <- gene_database[chromosome_arm == .arm2, 
                               intersect(symbol, all_complete_genes)]
  
  .all_combinations1 <- data.table(same_arm = F, arm = .arm1, 
                                   expand.grid(gene1 = .all_genes1, 
                                               gene2 = .all_genes2, 
                                               stringsAsFactors = F))
  
  # in contrast to same-arm predictions, we include these predictions twice, 
  # as each gene on arm 1 can be the buffering or buffered partner to arm 2
  # This does not affect the computed densities. 
  
  .arm1 <- str_c("chr", .chr, "q")
  .arm2 <- str_c("chr", .chr, "p")
  .gdb <- gene_database[chromosome_arm == .arm1]
  .all_genes1 <- gene_database[chromosome_arm == .arm1, 
                               intersect(symbol, all_complete_genes)]
  .all_genes2 <- gene_database[chromosome_arm == .arm2, 
                               intersect(symbol, all_complete_genes)]
  
  .all_combinations2 <- data.table(same_arm = F, arm = .arm1, 
                                   expand.grid(gene1 = .all_genes1, 
                                               gene2 = .all_genes2, 
                                               stringsAsFactors = F))
  
  
  .all_combinations <- rbindlist(list(.all_combinations1, .all_combinations2))
  
  .all_combinations <- .all_combinations[gene1 != gene2] %>% 
    merge(.gdb[, .(gene1 = symbol, pos1 = position)], by = "gene1", all.x = T) %>%
    merge(.gdb[, .(gene2 = symbol, pos2 = position)], by = "gene2", all.x = T)
  
  .all_combinations[, sorted_pair := sort_gene_pairs(gene1, gene2)]
  .all_combinations[, get("sorted_pair")]}, simplify = F) %>% 
  Reduce(c, .) %>% 
  unique() %>% 
  cacheR("different_arm_complete", .fpath)


if (F) {

message("Different chromosomes...")

important_pairs_default$different_chromosome <- sapply(seq_along(1:23), \(.i) {
  .chr <- chromosomes[.i]
  message(str_c("Chromosome ", .chr, "..."))
  .all_other_chromosomes <- chromosomes[(.i+1):24]

  .all_genes <- gene_database[chromosome == .chr, intersect(symbol, all_default_genes)]
  .all_on_other_chomosomes <- gene_database[
    chromosome %in% .all_other_chromosomes, intersect(symbol, all_default_genes)]
  
  .all_combinations <- data.table(same_arm = F, chr = str_c("chr", .chr), 
                                  expand.grid(gene1 = .all_genes, 
                                              gene2 = .all_on_other_chomosomes, 
                                              stringsAsFactors = F))
  # also here we count each pairing twice, as each chromosome can be buffering and buffered. 
  
  .all_combinations <- rbindlist(list(.all_combinations, .all_combinations))
  .all_combinations <- .all_combinations[gene1 != gene2]
  .all_combinations[, sorted_pair := sort_gene_pairs(gene1, gene2)]
  
  .all_combinations[, get("sorted_pair")]}, simplify = F) %>% 
  Reduce(c, .) %>% 
  unique() %>% 
  cacheR("different_chromosome_default", .fpath)




important_pairs_complete$different_chromosome <- sapply(seq_along(1:23), \(.i) {
  .chr <- chromosomes[.i]
  .all_other_chromosomes <- chromosomes[(.i+1):24]
  .all_genes <- gene_database[chromosome == .chr, intersect(symbol, all_complete_genes)]
  .all_on_other_chomosomes <- gene_database[chromosome %in% .all_other_chromosomes, 
                                            intersect(symbol, all_complete_genes)]
  
  .all_combinations <- data.table(same_arm = F, chr = str_c("chr", .chr), 
                                  expand.grid(gene1 = .all_genes, 
                                              gene2 = .all_on_other_chomosomes, 
                                              stringsAsFactors = F))
  
  # also here we count each pairing twice, 
  # as each chromosome can be buffering and buffered. 
  
  .all_combinations <- rbindlist(list(.all_combinations, .all_combinations))
  .all_combinations <- .all_combinations[gene1 != gene2]
  .all_combinations[, sorted_pair := sort_gene_pairs(gene1, gene2)]
  
  .all_combinations[, get("sorted_pair")]}, simplify = F) %>% 
  Reduce(c, .) %>% 
  unique() %>% 
  cacheR("different_chromosome_complete", .fpath)

}

