.fpath <- file.path(cache, "ensembl_paralogs")
mkdir(.fpath)

if (F) {
  #The Function executes a database query using Biomart. 
  #This might not always be successful. 
  #In these cases, this step has to be repeated, until a local cache of the 
  #paralog data types is generated. 
  build_ensembl_paralog_cache(.fpath)
}


# 20% minimum sequence identity
paralog_pairs$ensembl <- get_ensembl_paralogs(.fpath, 
                                              min_identity = 20, 
                                              only_pairs = F)


## ---- import paralogs curie ----

# download the relaxed ohnolog pairs file from the ohnologs.curie.fr database:
# "http://ohnologs.curie.fr/cgi-bin/DownloadBrowse.cgi?crit=[C]&org=hsapiens&opt=pairs&wgd=2R" -> 
# file.path(paralog_standards_filepath, "curie_ohnologs", "hsapiens.Pairs.Relaxed.2R.txt")

ohnolog_pairs_curie <- fread(
  file.path(paralog_standards_filepath, 
            "curie_ohnologs", 
            "hsapiens.Pairs.Relaxed.2R.txt"), 
  sep = "\t", check.names = T)
ohnolog_pairs_curie <- ohnolog_pairs_curie[
  Gene.type == "protein_coding", 
  .(Ohno1, Ohno2, 
    Symbol1, Symbol2, 
    Og.mean = Og.weighted.global.mean, 
    self.mean = Self.weighted.global.mean)]

paralog_pairs$ohnologs <- ohnolog_pairs_curie[
  Og.mean < 0.05 & self.mean < 0.3, 
  .(sorted_pair = sort_gene_pairs(Symbol1, Symbol2), 
    gene1 = Symbol1, 
    gene2 = Symbol2)]


## ---- import paralogs dekegel ----

# download "mmc2.csv" and mmc2.xlsx" from the supplementary files of 
# [DeKegel et al., 2021](https://doi.org/10.1016/j.cels.2021.08.006)
# to file.path(paralog_standards_filepath, "dekegel2021_supplementary").

paralog_pairs$dekegel_val <- data.table(
  read_excel(
    file.path(paralog_standards_filepath, 
              "dekegel2021_supplementary", 
              "mmc2.xlsx")))[
                !is.na(A1), 
                .(sorted_pair = sort_gene_pairs(A1, A2), 
                  gene1 = A1, gene2 = A2)]

paralog_pairs$dekegel <- fread(file.path(paralog_standards_filepath, 
                                         "dekegel2021_supplementary", 
                                         "mmc5.csv"), header = T)
paralog_pairs$dekegel[, sorted_pair := sort_gene_pairs(A1, A2)]


paralog_pairs$dekegel <- paralog_pairs$dekegel[
  A2_status_coef < 0 & A2_status_p_adj < 0.1, 
  .(sorted_pair, gene1 = A1, gene2 = A2)]

## ---- deprecated import paralogs dekegel ----

# run before paralog_pairs$dekegel is subsetted. The nonSL pairs are collected 
# based on other filter criteria. 
#nonSL_paralog_pairs$dekegel_nonSL <- paralog_pairs$dekegel[A2_status_p > 0.05, 
#.(sorted_pair, gene1 = A1, gene2 = A2)]

## ---- import paralogs ito ----


# download 41588_2021_967_MOESM2_ESM.xlsx from the supplementary files of 
# [Ito et al., 2021](https://doi.org/10.1038/s41588-021-00967-z)
# to file.path(paralog_standards_filepath, "Ito2021_supplementary").

files <- str_c("Supplementary table ", 8:10)
names(files) <- c("LFC", "pval", "FDR")
ito_cell_lines <- c("A549", "GI1", "HS936T", "HS944T", 
                    "HSC5", "IPC298", "MEL202", "MELJUSO", 
                    "MEWO", "PATU8988S", "PK1")

ito_predictions <- sapply(names(files), 
                          \(n) data.table(
                            read_excel(
                              file.path(paralog_standards_filepath, 
                                        "Ito2021_supplementary", 
                                        "41588_2021_967_MOESM2_ESM.xlsx"), 
                              sheet = files[[n]], skip = 2)) %>% 
                            melt.data.table(id.vars = "...1", 
                                            variable.name = "tissue", 
                                            value.name = n), simplify = F)

ito_predictions <- data.table(Reduce(
  \(x1, x2) {merge(x1, x2, by = c("...1", "tissue"))}, ito_predictions))
ito_predictions[, c("gene1", "gene2") := tstrsplit(...1, ";")]
for (cl in ito_cell_lines) {ito_predictions[grepl(cl, tissue), cell_line := cl]}
ito_predictions[, tissue := tolower(gsub(
  str_c(ito_cell_lines, "_", collapse = "|"), "", tissue))]
ito_predictions <- ito_predictions[, .(
  sorted_pair = sort_gene_pairs(gene1, gene2), 
  gene1, gene2, cell_line, tissue, LFC, pval, FDR)]
ito_predictions <- ito_predictions[order(cell_line, sorted_pair)]

ito_predictions[, n_sign_cell_lines := sum(LFC > 0 & FDR <= 0.05), 
                by = sorted_pair]

#nonSL_paralog_pairs$ito_nonSL <- ito_predictions[n_sign_cell_lines == 0, 
# .(sorted_pair, gene1, gene2)]
paralog_pairs$ito <- ito_predictions[n_sign_cell_lines >= 1, 
                                     .(sorted_pair, gene1, gene2)]




## ---- import paralogs anvar ----

# download "Supp_table7.xlsx" from the supplementary files of 
# [Anvar et al.](https://doi.org/10.6084/m9.figshare.24243832.v1)
# to file.path(paralog_standards_filepath, "Anvar2024_supplementary").

anvar_lfc <- fread(file.path(paralog_standards_filepath, 
                             "Anvar2024_supplementary", 
                             "Supp_table7.txt"))
anvar_lfc <- anvar_lfc[, .(Gene, 
                           A375 = A375CP1906_T21, 
                           A549 = A549_T21, 
                           Meljuso = MeljusoCP1905_T21, 
                           K562 = k562_T8)]

anvar_lfc[!grepl("_", Gene), gene1 := Gene]
anvar_lfc[grepl("_", Gene), str_c("gene", 1:4) := tstrsplit(Gene, "_")]

anvar_lfc <- anvar_lfc %>% 
  melt.data.table(id.vars = c("Gene", str_c("gene", 1:4)), 
                  value.name = "comb_LFC", 
                  variable.name = "cell_line", 
                  variable.factor = F)
anvar_lfc <- anvar_lfc %>% split(anvar_lfc$cell_line)

anvar_lfc <- sapply(names(anvar_lfc), \(.cl) {
  .d <- anvar_lfc[[.cl]][str_count(Gene, "_") == 1] %>% 
    merge(anvar_lfc[[.cl]][, .(Gene, gene1_LFC = comb_LFC)], 
          by.x = "gene1", by.y = "Gene", all.x = T) %>%
    merge(anvar_lfc[[.cl]][, .(Gene, gene2_LFC = comb_LFC)], 
          by.x = "gene2", by.y = "Gene", all.x = T) %>%
    merge(anvar_lfc[[.cl]][, .(Gene, gene3_LFC = comb_LFC)], 
          by.x = "gene3", by.y = "Gene", all.x = T) %>%
    merge(anvar_lfc[[.cl]][, .(Gene, gene4_LFC = comb_LFC)], 
          by.x = "gene4", by.y = "Gene", all.x = T)
  .d[, expected := apply(.SD, 1, sum, na.rm = T), 
     .SDcols = str_c("gene", 1:4, "_LFC")]
  .d[, c(.cl) := comb_LFC - expected]
  .d[, .SD, .SDcols = c("Gene", "gene1", "gene2", .cl)]
}, simplify = F) %>% 
  Reduce(\(.d1, .d2) {merge(.d1, .d2, by = c("Gene", "gene1", "gene2"))}, .)

anvar_lfc[, n_ess := apply(.SD <= -1, 1, sum, na.rm = T), 
          .SDcols = c("A375", "A549", "K562", "Meljuso")]
anvar_lfc[, sorted_pair := sort_gene_pairs(gene1, gene2)]

paralog_pairs$anvar <- copy(
  anvar_lfc[n_ess >= 1, 
            .(sorted_pair, gene1, gene2, n_ess, A375, A549, K562, Meljuso)])

## ---- import paralogs anvar standard ----

# download "Supp_table2.xlsx" from the supplementary files of 
# [Anvar et al.](https://doi.org/10.6084/m9.figshare.24243832.v1)
# to file.path(paralog_standards_filepath, "Anvar2024_supplementary").

paralog_pairs$anvar_standard <- data.table(read_xlsx(
  file.path(paralog_standards_filepath, 
            "Anvar2024_supplementary", 
            "Supp_table2.xlsx"), sheet = 1))
paralog_pairs$anvar_standard[, c("gene1", "gene2") := tstrsplit(...1, ":")]

paralog_pairs$anvar_standard <- paralog_pairs$anvar_standard[
  , .(sorted_pair = sort_gene_pairs(gene1, gene2), gene1, gene2, score)]

## ---- import paralogs thompson ----

# download "41467_2021_21478_MOESM8_ESM.xlsx" from the supplementary files of 
# [Thompson et al.](https://doi.org/10.1038/s41467-021-21478-9)
# to file.path(dirs$databases, "thompson2021_supplementary")


thompson_paralogs <- read_xlsx(
  file.path(paralog_standards_filepath, 
            "thompson2021_supplementary", 
            "41467_2021_21478_MOESM8_ESM.xlsx"), skip = 3) %>% 
  setDT
paralog_pairs$thompson <- thompson_paralogs[, unique(c(A375, Mewo, RPE))]
paralog_pairs$thompson <- data.table(
  sorted_pair = paralog_pairs$thompson[!is.na(paralog_pairs$thompson)])

paralog_pairs$thompson[, c("gene1", "gene2") := tstrsplit(sorted_pair, "_")]
paralog_pairs$thompson[, `:=`(A375 = sorted_pair %in% thompson_paralogs[, A375], 
                              Mewo = sorted_pair %in% thompson_paralogs[, Mewo], 
                              RPE = sorted_pair %in% thompson_paralogs[, RPE])]
paralog_pairs$thompson[, `:=`(sorted_pair = sort_gene_pairs(gene1, gene2), 
                              sign_cell_lines = apply(.SD, 1, sum, na.rm = T)), 
                       .SDcols = c("A375", "Mewo", "RPE")]
paralog_pairs$thompson <- paralog_pairs$thompson[sign_cell_lines > 1]

## ---- add chromosome info ----

# using the NCBI data, we add the chromosome of each gene to each paralog set
# of the list

paralog_pairs <- lapply(
  paralog_pairs, 
  \(x) {
    x %>% merge(feature_table_protein_coding[, .(gene = symbol, chr1 = chromosome)], 
                by.x = "gene1", by.y = "gene", all.x = T) %>% 
      merge(feature_table_protein_coding[, .(gene = symbol, chr2 = chromosome)], 
            by.x = "gene2", by.y = "gene", all.x = T)})


## ---- import paralogs schaeffer ----

#schaeffer_pairs <- read_xlsx(file.path(paralog_standards_filepath, 
#                                       "schaeffer2024_supplementary", 
#                                       "mmc2.xlsx")) %>% setDT
#schaeffer_pairs <- schaeffer_pairs[, .(sorted_pair = sort_gene_pairs(mut_gene, target_gene), gene1 = mut_gene, gene2 = target_gene)]


# export paralog sets

mkdir(file.path(cache, "important_gene_pairs"))
for (.n in names(paralog_pairs)) {
  saveRDS(paralog_pairs[[.n]], 
          file.path(cache, 
                    "important_gene_pairs", 
                    str_c("paralog_pairs_", .n, ".rds")))}
