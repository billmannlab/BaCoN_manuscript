

## ---- linear regression models ----

linreg_models <- list(
  data.table("eff1_exp2", "g1_effect ~ g2_expression", ""), 
  data.table("eff1_exp2.lineage", "g1_effect ~ g2_expression + lineage", 
             " + Lineage"), 
  data.table("eff1_exp2.growthpattern", 
             "g1_effect ~ g2_expression + growthpattern", 
             " + Growthpattern")) %>% rbindlist()
linreg_models <- linreg_models[, .(tag = V1, formula = V2, description = V3)]

## ---- define linear regression chunks ----

split_gene_set <- \(vec, size) {split(vec, ceiling(seq_along(vec) / size))}
chunk_size <- 800

chunk_definition <- list()

expression_chunks <- split_gene_set(expression_genes, chunk_size)
chronos_chunks <- split_gene_set(chronos_genes, chunk_size)

expression_genes_only_complete <- symdiff(expression_genes, 
                                          colnames(gene_expression_complete))
chronos_genes_only_complete <- symdiff(chronos_genes, colnames(chronos_complete))

expression_chunks_complete <- split_gene_set(expression_genes_only_complete, chunk_size)
chronos_chunks_complete <- split_gene_set(chronos_genes_only_complete, chunk_size)

id <- 1

for (i in 1:l(expression_chunks)) {
  
  # Adding "default" gene space
  
  for (j in 1:l(chronos_chunks)) {
    n <- str_c("size", chunk_size, id, "default", sep = "_")
    chunk_definition[[n]] <- data.table(name = n, 
                                        size = chunk_size, 
                                        set = "default", 
                                        exp_genes = list(expression_chunks[[i]]), 
                                        chr_genes = list(chronos_chunks[[j]]))
    id <- id + 1}
  
  # Adding remaining area of the complete gene space
  
  for (j in 1:l(chronos_chunks_complete)) {
    n <- str_c("size", chunk_size, id, "complete", sep = "_")
    chunk_definition[[n]] <- data.table(
      name = n, 
      size = chunk_size, 
      set = "complete", 
      exp_genes = list(expression_chunks[[i]]), 
      chr_genes = list(chronos_chunks_complete[[j]]))
    id <- id + 1}
}

for (i in 1:l(expression_chunks_complete)) {
  for (j in 1:l(chronos_chunks)) {
    n <- str_c("size", chunk_size, id, "complete", sep = "_")
    chunk_definition[[n]] <- data.table(
      name = n, 
      size = chunk_size, 
      set = "complete", 
      exp_genes = list(expression_chunks_complete[[i]]), 
      chr_genes = list(chronos_chunks[[j]]))
    id <- id + 1}
  for (j in 1:l(chronos_chunks_complete)) {
    n <- str_c("size", chunk_size, id, "complete", sep = "_")
    chunk_definition[[n]] <- data.table(
      name = n, 
      size = chunk_size, 
      set = "complete", 
      exp_genes = list(expression_chunks[[i]]), 
      chr_genes = list(chronos_chunks_complete[[j]]))
    id <- id + 1}
}

chunk_definition <- rbindlist(chunk_definition)

chunks_default <- chunk_definition[set == "default", name]
