## Nov 2025
## Louise Huuki-Myers
## Preform spatial registration on 4 layers from k=7 spatialDLPFC vs. human pilot layers

#### prep & load data ####

library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

## prep dirs
data_dir <- here("processed_data", "07_spatialDLPFC_registration", "01_spatialDLPFC_4layer_spatial_reg")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "07_spatialDLPFC_registration", "01_spatialDLPFC_4layer_spatial_reg")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data
spe <- fetch_data("spatialDLPFC_Visium")

#### Add 4 layer classifications to k=7 BayesSpace clusters ####

colData(spe)
table(spe$BayesSpace_harmony_07)

vis_clus_k7 <- vis_clus(spe,
         clustervar = "BayesSpace_harmony_07",
         sampleid = "Br2743_mid")

ggsave(vis_clus_k7, filename = here(plot_dir, "vis_clus_k7-Br2743_mid.png"))

layer_def <- tibble(harmony_07 = 1:7,
                    SpD = c("SpD1", "SpD1", "SpD2", "SpD3", "WM", "SpD4", "WM"))

layer_def |> arrange(SpD)
# harmony_07 SpD
# <int> <chr>
# 1          1 SpD1
# 2          2 SpD1
# 3          3 SpD2
# 4          4 SpD3
# 5          6 SpD4
# 6          5 WM
# 7          7 WM

spe$SpD <- layer_def$SpD[match(spe$BayesSpace_harmony_07, layer_def$harmony_07)]

table(spe$SpD)

SpD_colors = c('SpD1' = '#ED5E70',
               'SpD2' = '#CCBB44',
               'SpD3' = '#238934',
               'SpD4' = '#66CCEE',
               'WM' = "grey90")

vis_clus_SpD_l4 <- vis_clus(spe,
                            clustervar = "SpD",
                            sampleid = "Br2743_mid",
                            colors = SpD_colors,
                            guide_point_size = 4)

ggsave(vis_clus_SpD_l4, filename = here(plot_dir, "vis_clus_SpD_l4-Br2743_mid.png"))


#### run registration wrapper ####

colData(spe)[,c("sample_id", "SpD", "position", "age", "sex")]
rowData(spe)

modeling_results <- registration_wrapper(
    spe,
    var_registration = "SpD",
    var_sample_id = "sample_id",
    covars = c("position", "age", "sex"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    pseudobulk_rds_file = here(data_dir, "spatialDLPFC_spe_pb_l4.rds")
)

saveRDS(modeling_results, file = here(data_dir, "modeling_results_spatialDLPFC_l4.Rds"))

#### register to pilot data ####

human_pilot_model <- fetch_data(type = "modeling_results")

cor_layer <- layer_stat_cor(
    stats = modeling_results$enrichment,
    modeling_results = fetch_data(type = "modeling_results"),
    model_type = "enrichment",
    top_n = 100
)

## annotate

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

#   cluster layer_confidence   layer_label layer_label_simple
# 1    SpD1             good        Layer2                 L2
# 2    SpD2             good Layer4/Layer3               L4/3
# 3    SpD3             good        Layer5                 L5
# 4    SpD4             good        Layer6                 L6
# 5      WM             good            WM                 WM

write.csv(cor_layer, file = here(data_dir, "cor_layer_spatialDLPFC_l4_v_human_pilot.csv"))
write.csv(anno, file = here(data_dir, "cor_layer_anno_spatialDLPFC_l4_v_human_pilot.csv"))

saveRDS(list(cor_layer = cor_layer, anno = anno), file = here(data_dir, "cor_layer_anno_spatialDLPFC_l4_v_human_pilot.Rds"))

#### plot ####

cor_layer <- cor_layer[,sort(colnames(cor_layer))]

pdf(here(plot_dir, "spatialDLPFC_l4_v_human_pilot_cor_plot.pdf"), height = 6)
layer_stat_cor_plot(cor_layer,
                    reference_colors = libd_layer_colors,
                    query_colors = SpD_colors,
                    annotation = anno,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE
)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

