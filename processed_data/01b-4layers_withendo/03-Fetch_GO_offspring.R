library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
setwd("/Users/melis/Documents/MRI_cortical_layers/MRI_layers/code")
dirs <- list.files(path = "../data/GO pathways/GO_pathways_with_offspring/.", pattern = "\\.csv$")

for (i in 1:6) {
  data = read.csv(paste0("../data/GO pathways/GO_pathways_with_offspring/", dirs[i]))
  #data = data %>% separate(GO_term, into = c("GO", "pathway"), sep = "\t")
  #interest = data[data$Exclude!='Exclude',]
  print(length(data$GO_term))
  for (j in 1:length(data$GO_term)) {
    if (is.na(Ontology(data$GO_term[j]))==TRUE){
      print("Is NA")
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
  print(paste0("length pre unique filtration:", length(all_pathways)))
  final = unique(all_pathways)
  print(paste0("length post unique filtration:", length(final)))
  write.csv(final, file = paste0('../data/GO pathways/GO_pathways_with_offspring_', dirs[i]))
}

  #interest$GO
  # Specify the GO term
  go_term <- data$GO_term[1]

  # Retrieve the genes associated with the GO term
  genes <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term, columns = "SYMBOL", keytype = "GO")
  
  write.csv(genes[!is.na(genes$SYMBOL),], file = paste0('../data/GO_pathways/pathways_with_offspring_', dirs[i]))
  }