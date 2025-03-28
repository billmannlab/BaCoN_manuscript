## ---- sample information ----

covariates <- fread(file.path(depmap_filepath, str_c(.dm, "_Model.csv")))

covariates <- covariates[, .(ccle_code = CCLEName, 
                             cell_line = ModelID, 
                             #media = OnboardedMedia, 
                             GrowthPattern, 
                             OncotreePrimaryDisease, 
                             OncotreeLineage)]
setnames(covariates, tolower)

covariates[, lineage := fcase(
  # renaming Lung cancers to NSCLC and SCLC
  oncotreeprimarydisease == "Non-Small Cell Lung Cancer", "NSCLC", 
  oncotreeprimarydisease == "Lung Neuroendocrine Tumor", "SCLC", 
  oncotreelineage %in% c("Lung", "Skin", "Peripheral Nervous System"), 
  oncotreeprimarydisease)]
covariates[is.na(lineage), lineage := oncotreelineage]

# Organoid and Neurosphere growthpatterns are set to unknown
# covariates[growthpattern %in% c("Organoid", "Neurosphere"), 
# growthpattern := "Unknown"]

# tidy the covariates dataset by replacing problematic strings with "_"
strs_to_replace <- c("-", "/", " ")
cols <- c("growthpattern", "oncotreeprimarydisease", 
          "oncotreelineage", #"media", 
          "lineage")
for (s in strs_to_replace) {
  covariates[, (cols) := lapply(.SD, \(.) gsub(s, "_", .)), .SDcols = cols]}

saveRDS(covariates, file.path(cache, str_c("depmap_", .dm, "_covariates.rds")))