data = readRDS("~/MinaRyten/Aine/wood_full/data/dreamlet/de_X_sum_all_DE.rds")
sig = data[data$adj.P.Val<0.05,]
#level1 = sig[sig$annot_level=='2',]
sig[sig$gene_sym=='NEAT1',]

level1 = data[data$annot_level=='2',]
level1[level1$gene_sym=='NEAT1',]

level1 = sig[sig$annot_level=='3',]
level1[level1$gene_sym=='UGCG',]


c('CDIPT', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIK3CA', 'PIK3CB', 'PIK3CD')


data = readRDS("~/MinaRyten/Aine/wood_full/data/dreamlet/de_X_sum_all_DE.rds")
sig = data[data$adj.P.Val<0.05,]
#level1 = sig[sig$annot_level=='2',]
sig[sig$gene_sym=='CDIPT',]

#level1 = sig[sig$annot_level=='2',]
#level1[level1$gene_sym=='CDIPT',]

sig = data[data$adj.P.Val<0.05,]

sig[sig$gene_sym=='PIK3R1',]

sig[sig$gene_sym=='PIK3R2',]

sig[sig$gene_sym=='PIK3R3',]

sig[sig$gene_sym=='PIK3CA',]

sig[sig$gene_sym=='PIK3CB',]

sig[sig$gene_sym=='PIK3CD',]