library(tidyverse)

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

# read in DEGs and create extra column (up or down reg)
de_res = readRDS("~/MinaRyten/Aine/wood_full/data/dreamlet/de_X_sum_all_DE.rds") %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(direction = case_when(
    log_fc<0 ~ "Down-regulated", 
    log_fc>0 ~ "Up-regulated"
  )) %>% 
  tidyr::separate(coef, into=c("tissue", "comparison"), "_") #%>% 
  #dplyr::mutate(entrezId = ensembl_to_entrez(id))

for (test in c('myelin_structural_proteins', 'oligo_maturation', 'PI3K')){
    if (test=='myelin_structural_proteins'){
        test_genes = vroom::vroom("~/MinaRyten/Melissa/MRIxST/raw_data/myelin_gene_lists/myelin_structural_proteins.csv", delim = ",")
    } else if(test=='oligo_maturation'){
        test_genes = vroom::vroom("~/MinaRyten/Melissa/MRIxST/raw_data/myelin_gene_lists/oligo_maturation.csv", delim = ",")
    }  else if(test=='PI3K'){
        test_genes = vroom::vroom("~/MinaRyten/Melissa/MRIxST/raw_data/PI3K_gene_lists/PI3K.csv", delim = ",")
    }
        
    results = de_res %>% 
      dplyr::filter(annot_level=="1") %>% 
      dplyr::select(annot_level, tissue, assay, comparison, direction) %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(annot_level, tissue, assay, comparison, direction) %>% 
      dplyr::group_split() %>% 
      parallel::mclapply(X=., mc.cores=10, FUN=function(x){

        message(x)

        l=x$annot_level
        t=x$tissue
        a=x$assay
        c=x$comparison
        d=x$direction

        all = de_res %>% 
          dplyr::filter(annot_level==l, tissue==t, assay==a, comparison==c, direction==d) # %>% 
          #tidyr::drop_na(entrezId)

        sig = all %>% 
          dplyr::filter(adj_p_val<0.05)

        message(length(intersect(unique(sig$gene_sym), test_genes$gene_sym)))

        if( (intersect(unique(sig$gene_sym), test_genes$gene_sym) %>% length()) != 0 ){

            gene_enrichment(unique(sig$gene_sym), unique(de_res$gene_sym), test_genes$gene_sym) %>% 
            as.data.frame() %>% 
            dplyr::mutate(annot_level=l, tissue=t, assay=a, comparison=c, direction=d) %>% 
            return()

        }else{
          return(NULL)
          message("none sig.")
        }
      }) %>% 
      purrr::discard(is.null)

    # collate
    if (length(results) != 0){
    results = results %>% 
      dplyr::bind_rows() %>% 
      janitor::clean_names() %>% 
      dplyr::select(annot_level, tissue, assay, comparison, direction, p_value) %>% 
      dplyr::as_tibble() %>% 
      dplyr::group_by(annot_level, tissue, assay, comparison, direction) %>% 
      dplyr::mutate(p_adjust = p.adjust(p_value)) %>% 
      as.data.frame()
    results$test_enrichment = test
    
    results = results %>%
        group_by(test_enrichment) %>% 
        mutate(p_adjust = p.adjust(p_value, method = "fdr")) %>% 
        #ungroup() %>%
        dplyr::arrange(p_adjust)
    
        if (test=='myelin_structural_proteins'){
            final_results = results
        } else{
            final_results = rbind(final_results, results)
        }
    }

}    

#final_results <- final_results %>% 
#        group_by(assay) %>% 
#        mutate(p_adjust = p.adjust(p_value, method = "fdr")) %>% 
        #ungroup() %>%
#        dplyr::arrange(p_adjust)

final_results %>% 
  vroom::vroom_write(x=., file="../../processed_data/11_transcriptomics/annot_l1_enrich_test_lists_fdr_all.csv", delim = ",")

