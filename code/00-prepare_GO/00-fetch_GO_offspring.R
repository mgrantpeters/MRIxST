library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(biomaRt)
library(GO.db)
setwd("/Users/melis/Documents/GitHub/MRIxST/code")
dirs <- list.files(path = "../raw_data/keyword_search_amiGO/.", pattern = "\\.csv$")

for (i in 1:6) {
  # Read in amiGO pathways based on searched keywords
  data = read.csv(paste0("../raw_data/keyword_search_amiGO/", dirs[i]))
  print(length(data$GO_term))
  # Fetch offspring pathways
  for (j in 1:length(data$GO_term)) {
    if (is.na(Ontology(data$GO_term[j]))==TRUE){
      print("Is NA") # some pathway codes appear to not fit any ontology
                     # skip fetching ontology for these
    } else if (Ontology(data$GO_term[j])[[1]] == "BP") {
      offspring = GOBPOFFSPRING[[data$GO_term[j]]]
    } else if (Ontology(data$GO_term[j])[[1]] == "CC") {
      offspring = GOCCOFFSPRING[[data$GO_term[j]]]
    } else if (Ontology(data$GO_term[j])[[1]] == "MF") {
      offspring = GOMFOFFSPRING[[data$GO_term[j]]]
    }   
    if (j==1){
      all_pathways = c(data$GO_term[j], offspring)
    } else {
      all_pathways = c(all_pathways, data$GO_term[j], offspring)
    }
    
  }
  # Filter to avoid pathway duplication
  print(paste0("length pre unique filtration:", length(all_pathways)))
  final = unique(all_pathways)
  print(paste0("length post unique filtration:", length(final)))
  # Save expanded GO pathway list
  write.csv(final, file = paste0('../processed_data/00-prepare_GO/00-fetch_GO_offspring/', dirs[i]))
}
