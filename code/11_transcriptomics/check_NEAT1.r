data = readRDS("~/MinaRyten/Aine/wood_full/data/dreamlet/de_X_sum_all_DE.rds")
sig = data[data$adj.P.Val<0.05,]
#level1 = sig[sig$annot_level=='2',]
sig[sig$gene_sym=='NEAT1',]

level1 = sig[sig$annot_level=='2',]
level1[level1$gene_sym=='NEAT1',]