## ---- gene pair set stats ----

.fpath1 <- file.path(cache_filepath, str_c(.dm, "_gene_set_coordinates"))
.fpath2 <- file.path(cache_filepath, str_c(.dm, "_gene_pair_set_stats"))
mkdir(.fpath2)

gene_pair_set_stats <- list()

for (.n in names(paralog_pairs)) {
  .set_size <- paralog_pairs[[.n]][, .N]
  for (.s in c("default", "complete"
  )) {
    gene_pair_set_stats[[str_c(.s, "_", .n)]] <- {
      x <- readRDS(file.path(.fpath1, str_c(
        "gene_set_coordinates", .s, "paralogs", .n, .set_size, "pairs.rds", sep = "_")))
      
      data.table(gene_space = .s, 
                 nrow = nrow(x), 
                 ncol = ncol(x), 
                 gene_pair_set = .n, 
                 gene_pair_set_size = .set_size, 
                 hits = sum(x, na.rm = T))} %>% 
      cacheR(str_c("gene_pair_set_stats", .s, .n, sep = "_"), .fpath2)
    
    
    if (.s == "complete") {
      .set_size <- paralog_pairs[[.n]][, .N]
      
      gene_pair_set_stats[[str_c("dynamic_range_combinations", .n)]] <- 
        
        sapply(dynamic_range_combinations[, id], 
               \(.ss) {
                 x2 <- x[dynamic_range_combinations[id == .ss, unlist(exp_genes)],
                  dynamic_range_combinations[id == .ss, unlist(chr_genes)]]
                 
                 data.table(gene_space = .ss, 
                     nrow = nrow(x2), 
                     ncol = ncol(x2), 
                     gene_pair_set = .n, 
                     gene_pair_set_size = .set_size, 
                     hits = sum(x2, na.rm = T))}, simplify = F) %>% rbindlist() %>% 
          cacheR(str_c("gene_pair_set_stats", "dynamic_range_combinations", .n, sep = "_"), 
                 .fpath2)}
  }}

# Important pairs default

.s <- "default"

for (.n in names(important_pairs_default)) {
  message(.s, ", ", .n)
  .set_size <- l(important_pairs_default[[.n]])
  gene_pair_set_stats[[str_c(.s, "_", .n)]] <- {
    x <- readRDS(file.path(.fpath1, str_c("gene_set_coordinates", .s, .n, 
                                          .set_size, "pairs.rds", sep = "_")))
    
    data.table(gene_space = .s, 
               nrow = nrow(x), 
               ncol = ncol(x), 
               gene_pair_set = .n, 
               gene_pair_set_size = .set_size, 
               hits = sum(x, na.rm = T))} %>% 
    cacheR(str_c("gene_pair_set_stats", .s, .n, sep = "_"), .fpath2)
  }

# Important pairs complete

.s <- "complete"

for (.n in names(important_pairs_complete)) {
  message(.s, ", ", .n)  
  .set_size <- l(important_pairs_complete[[.n]])
  gene_pair_set_stats[[str_c(.s, "_", .n)]] <- {
    x <- readRDS(file.path(.fpath1, str_c("gene_set_coordinates", .s, .n, 
                                          .set_size, "pairs.rds", sep = "_")))
    
    data.table(gene_space = .s, 
               nrow = nrow(x), 
               ncol = ncol(x), 
               gene_pair_set = .n, 
               gene_pair_set_size = .set_size, 
               hits = sum(x, na.rm = T))} %>% 
    cacheR(str_c("gene_pair_set_stats", .s, .n, sep = "_"), .fpath2)}

gene_pair_set_stats <- gene_pair_set_stats %>% rbindlist()

gene_pair_set_stats[, `:=`(n_pairs = prod(nrow, ncol), 
                           density = hits / prod(nrow, ncol)), by = .I]

gene_pair_set_stats
