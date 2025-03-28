## ---- generate gene database ----

ncbi <- list(version = "GCF_000001405.40_GRCh38.p14")

ncbi$url <- str_c("https://ftp.ncbi.nih.gov", 
                  "genomes/refseq/vertebrate_mammalian/Homo_sapiens", 
                  "all_assembly_versions", ncbi$version, "/", sep = "/")

ncbi$fpath <- file.path(cache, 
                        "NCBI", str_c("assembly_", 
                                      ncbi$version)); mkdir(ncbi$fpath)

.fname <- str_c(ncbi$version, "_feature_table.txt.gz")
.file <- file.path(ncbi$fpath, .fname)
if (!file.exists(.file)) {
  utils::download.file(str_c(ncbi$url, .fname), .file, method = "curl", quiet = T)}

feature_table <- fread(.file, check.names = T)
setnames(feature_table, tolower)
feature_table <- feature_table[assembly_unit == "Primary Assembly" & 
                                 chromosome != ""]

.dup_xy_genes <- feature_table[class == "protein_coding" & 
                                 chromosome %in% c("X", "Y") & 
                                 duplicated(symbol), get("symbol")]

feature_table[x..feature == "gene" & class == "protein_coding" & 
                symbol %in% .dup_xy_genes, `:=`(chromosome = "X/Y")]


feature_table_protein_coding <- feature_table[x..feature == "gene" & 
                                                class == "protein_coding"]

feature_table_protein_coding <- feature_table_protein_coding[
  , .(geneid, symbol, chromosome, 
      start, end, strand, 
      description = name, class)] %>% distinct()


saveRDS(feature_table, file.path(cache, "ncbi_feature_table.rds"))
saveRDS(feature_table_protein_coding, 
        file.path(cache, "ncbi_feature_table_protein_coding.rds"))




# Retrieve UCSC cytoband and chromosome information
ucsc_cytobands <- get_ucsc_cytobands() %>% cacheR("UCSC_cytobands_hg38", cache)



chr_sizes <- ucsc_cytobands[, .(chr_size = max(end, na.rm = T)), by = chr]

# The gathered infos are merged into a gene database, 
# containing the chromosome arm location of each gene. 

gene_database <- copy(feature_table_protein_coding)

for (.chr in ucsc_cytobands[, unique(chr)]) {
  for (.arm in ucsc_cytobands[chr == .chr, arm]) {
    .x <- ucsc_cytobands[chr == .chr & arm == .arm]
    gene_database[chromosome == .chr & start >= .x[, start] & end <= .x[, end], 
                  arm := .arm]}}
gene_database[, `:=`(chromosome_arm = str_c("chr", chromosome, arm), 
                     position = round((start + end)/2, 0))]
gene_database



