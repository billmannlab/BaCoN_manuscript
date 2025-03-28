
## ---- gene pair set coordinates ----

.fpath <- file.path(cache, str_c(.dm, "_gene_set_coordinates"))
mkdir(.fpath)



default_template <- empty_array(list(expression_genes, chronos_genes), F)
complete_template <- empty_array(list(colnames(gene_expression_complete), 
                                      colnames(chronos_complete)), F)


# Default sets

for (.n in names(paralog_pairs)) {
  message(.n)
  .d <- {
    x <- default_template
    x[coords_of_pairs(x, paralog_pairs[[.n]]$sorted_pair)] <- T
    x} %>% cacheR(str_c("gene_set_coordinates_default_paralogs", .n, 
                        paralog_pairs[[.n]][, .N], "pairs", sep = "_"), 
                  .fpath)}


for (.n in names(important_pairs_default)) {
  message(.n)
  .d <- {
    x <- default_template
    x[coords_of_pairs(x, important_pairs_default[[.n]])] <- T
    x} %>% cacheR(str_c("gene_set_coordinates_default", .n, 
                        l(important_pairs_default[[.n]]), "pairs", sep = "_"), 
                  .fpath)}


# Complete sets

for (.n in names(paralog_pairs)) {
  message(.n)
  .d <- {
    x <- complete_template
    x[coords_of_pairs(x, paralog_pairs[[.n]]$sorted_pair)] <- T
    x} %>% cacheR(str_c("gene_set_coordinates_complete_paralogs", .n, 
                        paralog_pairs[[.n]][, .N], "pairs", sep = "_"), 
                  .fpath)}


for (.n in names(important_pairs_complete)) {
  message(.n)
  .d <- {
    x <- complete_template
    x[coords_of_pairs(x, important_pairs_complete[[.n]])] <- T
    x} %>% cacheR(str_c("gene_set_coordinates_complete", .n, 
                        l(important_pairs_complete[[.n]]), "pairs", sep = "_"), 
                  .fpath)}