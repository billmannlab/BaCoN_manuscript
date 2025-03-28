## ---- reviewer response import other datatypes ----

# Methylation

.d <- fread(file.path("D:", "Promotion_databases", "depmap",
                      "CCLE_RRBS_TSS_1kb_20180614.txt"))

.d <- .d %>% melt.data.table(measure.vars = names(.d)[8:850], 
                             variable.name = "cell_line", 
                             value.name = "methylation", variable.factor = F)

# subset the data to genes expressed in Depmap and cell lines in the depmap model info

.d <- .d[cell_line %in% covariates[, unique(ccle_code)] & 
           gene %in% colnames(gene_expression_complete)]
.d[, cell_line := covariates[cell_line, on = "ccle_code", get("cell_line")]]
.d <- .d[, .(methylation = mean(methylation, na.rm = T)), 
         by = .(cell_line, chr, gene)]
.d <- .d[order(cell_line, gene)]

.cl_oi <- .d[, unique(cell_line)]
.goi <- .d[, unique(gene)]
methylation_complete <- empty_array(list(.cl_oi, .goi))

for (cl in .cl_oi) {methylation_complete[cl,] <- .d[cell_line == cl][.goi, 
                                                                     on = "gene", 
                                                                     get("methylation")]}

# copy number data

copy_numbers_complete <- file.path(depmap_filepath, str_c(.dm, "_OmicsCNGene.csv")) %>% 
  import_depmap()
# A few genes with missing data are discarded. 
copy_numbers <- copy_numbers_complete[,apply(copy_numbers_complete, 2, na_perc) == 0]

expression_genes_with_cnv <- intersect(colnames(gene_expression), colnames(copy_numbers))

copy_numbers <- copy_numbers_complete[,expression_genes_with_cnv]

# RNAi based Demeter dependency scores

demeter_complete <- import_depmap(
  file.path("D:", "Promotion_databases", "CCLE_DEMETER2_v6", "D2_Achilles_gene_dep_scores.csv"))
rownames(demeter_complete) <- stringr::str_split_i(rownames(demeter_complete), " \\(", 1)

demeter_complete <- demeter_complete[,colnames(demeter_complete) %in% covariates[, ccle_code]]

colnames(demeter_complete) <- covariates[colnames(demeter_complete), on = "ccle_code", get("cell_line")]
demeter_complete <- t(demeter_complete)

demeter <- demeter_complete[,colnames(demeter_complete) %in% colnames(gene_expression_complete)]

# MS-based protein expression

## Import Gygi MS Protein quant data:

.d <- fread(
  file.path("D:", "Promotion_databases", "gygi_protein_expression", 
            "protein_quant_current_normalized.csv.gz"))
setnames(.d, c("Protein_Id", "Gene_Symbol", "Description", "Group_ID", "Uniprot", "Uniprot_Acc"), 
         c("protein_id", "gene_symbol", "description", "group_id", "uniprot", "uniprot_acc"))

.metadata_cols <- c("protein_id", "gene_symbol")
.cl_oi <- names(.d)[49:426]
.d <- .d[gene_symbol != "", .SD, .SDcols = c(.metadata_cols, .cl_oi)]

# some genes are measured multiple times. We mean-summarize these observations. 
.d <- .d[, lapply(.SD, mean, na.rm = T), .SDcols = .cl_oi, by = .(gene_symbol)]

# arrange values into an array:
genes_oi <- .d[, get("gene_symbol")]
.d <- .d[genes_oi, on = "gene_symbol"]
gygi_pe <- empty_array(list(.cl_oi, genes_oi))

for (.cl in .cl_oi) {gygi_pe[.cl,genes_oi] <- .d[, get(.cl)]}

# we are interested in genes represented in at least 200 cell lines:
gygi_pe <- gygi_pe[,names(which(apply(gygi_pe, 2, \(.x) {sum(!is.na(.x))}) > 200))]

rownames(gygi_pe) <- str_split_i(rownames(gygi_pe), pattern = "_TenPx", 1)

# Translate Gygi cell lines:

gygi_to_dm <- data.table(read_xlsx(
  file.path("D:", "Promotion_databases", 
            "gygi_protein_expression", "Table_S1_Sample_information.xlsx"), 
  sheet = 2), check.names = T)
gygi_to_dm[Cell.Line == "CNI-H1568", Cell.Line := "NCI-H1568"]
gygi_to_dm <- gygi_to_dm[, .(ccle_code = CCLE.Code)] %>% distinct()

gygi_to_dm <- merge(gygi_to_dm, covariates, by = "ccle_code")

rownames(gygi_pe) <- gygi_to_dm[rownames(gygi_pe), on = "ccle_code", get("cell_line")]



# 23Q2 DepMap mutation status

message(str_c("Importing mutation data for version ", .dm, "..."))

if (.dm == "23Q2") {
  x <- fread(file.path(depmap_filepath, str_c(.dm, "_OmicsSomaticMutations.csv")))
  
  loss_mutation_types <- toupper(c("frame_shift_del", "frame_shift_ins", 
                                   "intron", "missense", "nonsense", "nonstop", 
                                   "start_codon_ins"))
  
  x <- x[VariantInfo %in% loss_mutation_types, 
         .(ModelID,
           HugoSymbol, 
           VariantInfo, 
           LikelyLoF)]
  
  .cl_oi <- x[, unique(ModelID)]
  .goi <- x[, unique(HugoSymbol)]
  
  loss_status_complete <- empty_array(list(.cl_oi, .goi), 0)
  
  .pb <- progbar(l(.cl_oi))
  
  for (cl in .cl_oi) {
    .pb$tick()
    loss_status_complete[cl,x[ModelID == cl, get("HugoSymbol")]] <- 1}
}


if (.dm == "24Q2") {
  loss_status_complete <- import_depmap(
    file.path(depmap_filepath, 
              str_c(.dm, "_OmicsSomaticMutationsMatrixDamaging.csv")))
}
