library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:49.941401+00:00 shinyapps[5096402]: Listening on http://127.0.0.1:32863
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’
# 2021-11-11T05:30:50.222437+00:00 shinyapps[5096402]:   111: <Anonymous>
library("Polychrome")

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the spe object
load("spatialDLPFC_equivolume_domains.Rdata", verbose = TRUE)
# load the pseudobulked object spe_pseudo
spe_pseudo <- readRDS("spatialDLPFC_equivolume_domains_pseudobulk.Rds")
# load modeling results for pd domains
modeling_results <- readRDS("spatialDLPFC_equivolume_domains_modeling_results.Rds")
sig_genes <- readRDS("spatialDLPFC_equivolume_domains_sig_genes.Rds")

vars <- colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = spe_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "spatialDLPFC, Visium, Sp09",
    spe_discrete_vars = c( # this is the variables for the spe object not the spe_pseudo object
        "BayesSpace",
        "ManualAnnotation",
        vars[grep("^SpaceRanger_|^scran_", vars)],
        vars[grep("^BayesSpace_harmony", vars)],
        vars[grep("^BayesSpace_pca", vars)],
        "graph_based_PCA_within",
        "PCA_SNN_k10_k7",
        "Harmony_SNN_k10_k7",
        "manual_layer_label",
        "wrinkle_type",
        "BayesSpace_colors"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        vars[grep("^VistoSeg_", vars)],
        vars[grep("^layer_", vars)],
        vars[grep("^broad_", vars)]
    ),
    default_cluster = "BayesSpace",
    docs_path = "www"
)
