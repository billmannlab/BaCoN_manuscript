# load essential packages

library(data.table)
library(readxl)

#library(tidyverse); 

library(magrittr)
library(stringr)
library(dplyr)
library(ggplot2); theme_set(theme_classic())

library(sva) # ComBat
library(whitening) # whitening
library(job) # run multiple R commands in separate sessions
library(progress) # display progress bars
library(biomaRt) # request data from NCBI database
#library(UCSC.utils) # Cytoband information

# load packages required for plotting

library(BSgenome.Hsapiens.UCSC.hg38) # circos plots
library(circlize) # circos plots
library(ComplexHeatmap) # heatmaps
library(scattermore) # large scatter plots
#library(pheatmap)
library(viridis)
library(igraph) # networks

# load utility package
#if (F) {devtools::install_github("thorohde/TR.R.utils")}
#library(TR.R.utils)
#devtools::load_all("D:/Promotion_projects_github/TR.R.utils/")

# load BaCoN package
if (F) {devtools::install_github("billmannlab/BaCoN")}

# load BaCoN library. The package is loaded locally here. 
devtools::load_all("D:/Promotion_projects_github/BaCoN_package/")
#library(BaCoN)
