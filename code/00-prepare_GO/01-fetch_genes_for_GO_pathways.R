library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
setwd("/Users/melis/Documents/GitHub/MRIxST/code")
dirs <- list.files(path = "../processed_data/00-prepare_GO/00-fetch_GO_offspring/.", pattern = "\\.csv$")

for (i in 1:6) {
  data = read.csv(paste0("../processed_data/00-prepare_GO/00-fetch_GO_offspring/", dirs[i]))
  go_terms <- data$x
  genes <- AnnotationDbi::select(org.Hs.eg.db, keys = go_terms, columns = "SYMBOL", keytype = "GO")
  write.csv(genes[!is.na(genes$SYMBOL),], file = paste0('../processed_data/00-prepare_GO/01-fetch_genes_for_GO_pathways/genes_', dirs[i]))
}
