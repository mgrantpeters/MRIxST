## Nov 2025
## Louise Huuki-Myers
## Preform spatial registration on k=7

#### prep & load data ####

library("spatialLIBD")
library("tidyverse")
library("here")

data_dir <- here()

## load data
spe <- fetch_data("spatialDLPFC_Visium")

#### Explore k=7  - add 4 layer classifications ####

colData(spe)
table(spe$BayesSpace_harmony_07)

vis_clus(spe,
         clustervar = "BayesSpace_harmony_07",
         sampleid = "Br2743_mid")

layer_def <- tibble(harmony_07 = 1:7,
                    layer_def = c("SpD1", "SpD1", "SpD2", "SpD3", "WM", "SpD4", "WM"))

spe$MRI_layer <- layer_def$layer_def[match(spe$BayesSpace_harmony_07, layer_def$harmony_07)]

table(spe$MRI_layer)

MRI_colors = c('SpD1' = '#ED5E70',
               'SpD2' = '#CCBB44',
               'SpD3' = '#238934',
               'SpD4' = '#66CCEE',
               'WM' = "grey90")

vis_clus(spe,
         clustervar = "MRI_layer",
         sampleid = "Br2743_mid",
         colors = MRI_colors)

#### run registration wrapper ####

colData(spe)[,c("sample_id", "MRI_layer", "position", "age", "sex")]
rowData(spe)

modeling_results <- registration_wrapper(
    spe,
    var_registration = "MRI_layer",
    var_sample_id = "sample_id",
    covars = c("position", "age", "sex"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

#### register to pilot data ####

human_pilot_model <- fetch_data(type = "modeling_results")

cor_stats_layer <- layer_stat_cor(
    stats = modeling_results$enrichment,
    modeling_results = fetch_data(type = "modeling_results"),
    model_type = "enrichment",
    top_n = 100
)


#### plot ####

layer_stat_cor_plot(cor_stats_layer,
                    reference_colors = libd_layer_colors)

