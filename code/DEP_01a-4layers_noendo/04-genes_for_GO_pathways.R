library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
setwd("/Users/melis/Documents/MRI_cortical_layers/MRI_layers/code")
dirs <- list.files(path = "../data/GO pathways/GO_pathways_with_offspring/.", pattern = "\\.csv$")

for (i in 1:6) {
  data = read.csv(paste0("../data/GO pathways/GO_pathways_with_offspring/", dirs[i]))
  go_terms <- data$x
  genes <- AnnotationDbi::select(org.Hs.eg.db, keys = go_terms, columns = "SYMBOL", keytype = "GO")
  write.csv(genes[!is.na(genes$SYMBOL),], file = paste0('../data/GO pathways/GO_pathways_with_offspring/genes_', dirs[i]))
}
