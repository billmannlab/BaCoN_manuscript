# The following files are required for the script to run. 
# If they are not located in the depmap database directory, 
# please manually download them from the depmap portal 
# (https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q2)
# and rename them with the "23Q2_" prefix: 

# Model data: "Model.csv" -> "23Q2_Model.csv"
# Chronos scores: "CRISPRGeneEffect.csv" -> "23Q2_CRISPRGeneEffect.csv"
# Gene expression: "OmicsExpressionProteinCodingGenesTPMLogp1.csv" -> 
#  "23Q2_OmicsExpressionProteinCodingGenesTPMLogp1.csv"
# Copy Number Variation: "OmicsCNGene.csv" -> "23Q2_OmicsCNGene.csv"
# Assay Conditions: "ModelCondition.csv" -> "23Q2_ModelCondition.csv"

## ---- import DepMap ----

chronos_complete <- file.path(
  depmap_filepath, str_c(.dm, "_CRISPRGeneEffect.csv")) %>% 
  import_depmap(impute = T) # impute missing values with colmeans

gene_expression_complete <- file.path(
  depmap_filepath, 
  str_c(.dm, "_OmicsExpressionProteinCodingGenesTPMLogp1.csv")) %>% 
  import_depmap()


## ---- filter DepMap ----

# To correlate fitness effect with expression data, 
# it is necessary to equalize cell lines
cl_subset_default <- intersect(rownames(gene_expression_complete), 
                            rownames(chronos_complete))

#expression_genes_complete <- colnames(gene_expression_complete)
#chronos_genes_complete <- colnames(chronos_complete)


# To be able to buffer another gene, a gene needs to be expressed. 
# A mean expression threshold is applied on the raw expression gene set. 
expression_genes <- names(which(apply(gene_expression_complete, 2, mean) > 1))
gene_expression <- gene_expression_complete[cl_subset_default,expression_genes]

# We only keep effect genes that we have expression data for.
chronos_genes <- intersect(expression_genes, colnames(chronos_complete))
chronos <- chronos_complete[cl_subset_default,chronos_genes]


# Some genes (USP17L5, USP17L10, USP17L19) are not expressed in any of the 
# intersecting cell lines. 
# They are removed from some analyses like linear regression.  

unexpressed_genes <- names(which(
  apply(gene_expression_complete[cl_subset_default,] == 0, 2, all)))


if (.dm == "24Q2") {
  chronos_uncorrected_complete <- file.path(
    depmap_filepath, 
    str_c(.dm, "_CRISPRGeneEffectUncorrected.csv")) %>% import_depmap(impute = T)
  
  chronos_uncorrected <- chronos_uncorrected_complete[cl_subset_default,
                                                      chronos_genes]
  }

