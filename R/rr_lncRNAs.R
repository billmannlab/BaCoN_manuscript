## ---- reviewer response ncRNAs ----

ge_with_RNAs_complete <- import_depmap(
  file.path("D:", "Promotion_databases", "depmap", "24Q2", 
            "24Q2_OmicsExpressionAllGenesTPMLogp1Profile.csv"))


gene_database_with_ncRNA <- feature_table[symbol %in% colnames(ge_with_RNAs_complete) & 
                                            x..feature %in% c("gene", "ncRNA")]


for (.chr in ucsc_cytobands[, unique(chr)]) {
  for (.arm in ucsc_cytobands[chr == .chr, arm]) {
    .x <- copy(ucsc_cytobands[chr == .chr & arm == .arm])
    gene_database_with_ncRNA[chromosome == .chr & start >= .x[, start] & end <= .x[, end], 
                             arm := .arm]}}

# for some genes / lcRNAs, the database contains multiple entries. 
# for these entries, the start and end position is approximated as median, 
# to be more outlier-resilient

gene_database_with_ncRNA <- gene_database_with_ncRNA[
  class %in% c("protein_coding", "lncRNA", "miRNA", "snRNA", "ncRNA"), 
  .(start = median(start, na.rm = T), 
    end = median(end, na.rm = T)), 
  by = .(class, chromosome, strand, arm, symbol, geneid)]

gene_database_with_ncRNA[, `:=`(chromosome_arm = str_c("chr", chromosome, arm), 
                                position = round((start + end)/2, 0))]


exp_gene_info_all <- copy(gene_database_with_ncRNA[, .(
  expression_gene = symbol, 
  expression_class = class, 
  exp_chr = chromosome, 
  exp_pos = position, 
  exp_start = start, 
  exp_end = end, 
  exp_chr_arm = chromosome_arm)])

eff_gene_info_all <- copy(gene_database_with_ncRNA[, .(
  effect_gene = symbol, 
  effect_class = class, 
  eff_chr = chromosome, 
  eff_pos = position, 
  eff_start = start, 
  eff_end = end, 
  eff_chr_arm = chromosome_arm)])

# translate cell line names:
omics_profiles <- fread(file.path("D:", "Promotion_databases", "depmap", "24Q2", 
                                  str_c("24Q2", "_OmicsProfiles.csv")))
rownames(ge_with_RNAs_complete) <- omics_profiles[rownames(ge_with_RNAs_complete), 
                                                  on = "ProfileID", get("ModelID")]

.exp_th <- 1

goi <- list(protein_coding = "protein_coding", 
            ncRNAs = c("lncRNA", "miRNA", "snRNA", "ncRNA"))
goi <- sapply(names(goi), \(.n) {gene_database_with_ncRNA[class %in% goi[[.n]], 
                                                          get("symbol")]}, simplify = F)

goi$protein_coding_expressed <- intersect(
  names(which(apply(ge_with_RNAs_complete, 2, mean) >= 2)), 
  goi$protein_coding)

goi$ncRNAs_expressed <- intersect(
  names(which(apply(ge_with_RNAs_complete > .exp_th, 2, sum) >= 100)), 
  goi$ncRNAs)
goi$all <- Reduce(c, goi[c("protein_coding_expressed", "ncRNAs_expressed")])
goi$chronos <- intersect(goi$all, colnames(chronos_complete))

expression_genes_with_ncRNA <- Reduce(c, goi[c("protein_coding_expressed", 
                                               "ncRNAs_expressed")])
expression_ncRNAs <- goi$ncRNAs_expressed
ge_with_RNAs <- ge_with_RNAs_complete[,expression_genes_with_ncRNA]

# Cell lines to the ones also found in the chronos score object
.cl_intersect <- intersect(rownames(ge_with_RNAs_complete), rownames(chronos_complete))
ge_with_RNAs <- ge_with_RNAs[.cl_intersect,]

if (T) {
  .pcc_bacon <- BaCoN_pipeline(
    expression_matrix = x1, 
    effect_matrix = x2, 
    cache_path = file.path(reviewer_response_fpath, str_c(.dm, "_ncRNAs"), 
                           str_c(.dm, "_test_CvEnc_exp_", .exp_th, "_PCC_BaCoN_0.05")))
  
  .pred <- T
  if (.pred) {
    evaluation_routine(
      .pcc_bacon, 
      instructions = list(list(mat = "correlation_matrix"), 
                          list(mat = "bacon_matrix", ref = "correlation_matrix")) %>% 
        setNames(str_c(.dm, "_test_CvEnc_exp_", .exp_th, "_PCC", c("", "_BaCoN_0.05"))))}
}

## with Cholesky whitening: 

if (.pred) {
  .pcc_bacon <- BaCoN_pipeline(
    expression_matrix = t(whiten(t(x1), method = "Cholesky")), 
    effect_matrix = t(whiten(t(x2), method = "Cholesky")), 
    cache_path = file.path(
      reviewer_response_fpath, str_c(.dm, "_ncRNAs"), 
      str_c(.dm, "_test_CvEnc_exp_", .exp_th, "_PCC_BaCoN_0.05_whiten.Cholesky")))
  
  
  if (T) {
    evaluation_routine(
      .pcc_bacon, 
      instructions = list(list(mat = "correlation_matrix"), 
                          list(mat = "bacon_matrix", ref = "correlation_matrix")) %>% 
        setNames(str_c(.dm, "_test_CvEnc_exp_", .exp_th, "_PCC_", 
                       c("whiten.Cholesky", "BaCoN_0.05_whiten.Cholesky"))))
  }}
