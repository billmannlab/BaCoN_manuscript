## ---- initialize ----

# analysis parameters

chromosomes <- c(1:22, "X", "Y", "MT")
chromosome_arms <- str_c("chr", rep(chromosomes[1:24], each = 2), 
                         rep(c("p", "q"), 24))

cl_subsets <- list()
paralog_pairs <- list()

methods_whiten <- c("PCA", "Cholesky")
combat_methods <- c("lineage", "growthpattern", "chromosome", "chromosome_arm")

# filepaths

depmap_filepath <- file.path("D:", "Promotion_databases", "depmap", .dm)
cache <- file.path("D:", "Promotion_cache", "BaCoN")
large_cache <- file.path("G:", "Promotion_cache_large", "BaCoN")

reviewer_response_fpath <- file.path(large_cache, "reviewer_response")
mkdir(reviewer_response_fpath)

result_dir <- file.path(cache, paste0(.dm, "_results")); mkdir(result_dir)
fig_directory <- "figures_new"; mkdir(fig_directory)

perturbseq_filepath <- file.path("D:", "Promotion_databases", "perturbseq_datasets")
paralog_standards_filepath <- file.path("D:", "Promotion_databases", "paralog_standards")

lm_instructions_fpath <- file.path(large_cache, 
                                   str_c(.dm, "_linear_regression_cache_instructions"))
mkdir(lm_instructions_fpath)

linreg_filepath <- file.path(large_cache, 
                             str_c(.dm, "_linear_regression"))
mkdir(linreg_filepath)

