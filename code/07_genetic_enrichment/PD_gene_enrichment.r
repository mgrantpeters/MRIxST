library(tidyverse)
library(vroom)

# Define enrichment function using fisher's test
gene_enrichment<-function(sig_genes, background, test_gene_list){
  
  sig_in<-length(sig_genes[sig_genes %in% test_gene_list])
  sig_out<-length(sig_genes[!sig_genes %in% test_gene_list])
  back_in<-length(background[background %in% test_gene_list])
  back_out<- length(background[!background %in% test_gene_list])
  
  dat <- data.frame("In Gene Set" = c(sig_in, sig_out),
                    "Not In Gene Set" = c(back_in, back_out),
                    row.names = c("Significantly Spliced Genes", "Background Genes"),
                    stringsAsFactors = FALSE
  )
  
  test <- fisher.test(dat, alternative = "greater")
  
  
  res<-tibble(p_value=test$p.value, odds_ratio=test$estimate, sig_ratio=sig_in/sig_out, bg_ratio=back_in/back_out)
  
  return(res)
  
}

disease_genes = vroom::vroom("~/MinaRyten/ASAP/downstream_analyses/check_for_enrichment/PD_commonGWAS_nalls2019_kim2023_GP2_2024_union_with_ens.csv", delim = ",")
de_res = vroom::vroom("../processed_data/07_genetic_enrichment/background_genes_DLPFC.csv", delim = ",")
category_list = c('myelin','iron','lipid','inflammation')

for (category in category_list){
    path_gene = vroom::vroom(paste0("~/MinaRyten/Melissa/MRIxST/processed_data/07_genetic_enrichment/filtered_genes_", category,"_GO_term.csv"), delim = ",")
    pathway_list = unique(path_gene$GO)
    for (pathway in pathway_list){
        subset = path_gene[path_gene$GO == pathway,]
        gene_list = subset$SYMBOL
        if (length(intersect(gene_list, disease_genes$gene_sym))>0) {
            intersect_lists=intersect(gene_list, disease_genes$gene_sym)
            n_intersect = length(intersect_lists)
            perc_overlap = (100*length(intersect_lists))/length(gene_list)
        } else {
            intersect_lists=NA   
            n_intersect=0
            perc_overlap=0
        }
        results = gene_enrichment(gene_list, unique(de_res$gene_name), disease_genes$gene_sym) %>% 
            as.data.frame() %>% 
            dplyr::mutate(pathway=pathway, category=category, n_intersect=n_intersect, perc_overlap=perc_overlap, list_of_genes = paste(unlist(intersect_lists),collapse=", ")) 
        if ((category == category_list[1]) & (pathway == pathway_list[1])) {
            all_results = results
        } else {
            all_results =  dplyr::bind_rows(all_results, results)
        }
        
    }
}

all_results = all_results %>% 
        mutate(p_adjust = p.adjust(p_value, method = "fdr")) %>% 
        dplyr::arrange(p_adjust) %>% 
        vroom::vroom_write(x=., file=paste0("../processed_data/07_genetic_enrichment/PD_gene_enrichment.csv"), delim = ",")