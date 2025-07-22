## With codde from:
# https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/01_build_spe/03_add_deconvolution.R#L250
# https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/07_layer_differential_expression/01_create_pseudobulk_data.R
# https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/07_layer_differential_expression/03_model_BayesSpace.R

library("spatialLIBD")
library("here")
library("lobstr")
library("jaffelab")
library("scater")
library("sessioninfo")

dir_rdata <- here::here("processed_data", "05_spatialDLPFC_equivolume_domains")
dir.create(dir_rdata, showWarnings = FALSE)

## Load spatialDLPFC data
spe <- fetch_data("spatialDLPFC_Visium")
lobstr::obj_size(spe)
# 6.97 GB

## Create the 4 equivolume spatial domains + WM using spatialDLPFC Sp07 clusters
spe$BayesSpace <- NA
spe$BayesSpace[spe$BayesSpace_harmony_07 == 1 | spe$BayesSpace_harmony_07 == 2] <- "SpD1"
spe$BayesSpace[spe$BayesSpace_harmony_07 == 3] <- "SpD2"
spe$BayesSpace[spe$BayesSpace_harmony_07 == 4] <- "SpD3"
spe$BayesSpace[spe$BayesSpace_harmony_07 == 6] <- "SpD4"
spe$BayesSpace[spe$BayesSpace_harmony_07 == 5 | spe$BayesSpace_harmony_07 == 7] <- "WM"
table(spe$BayesSpace, useNA = "ifany")
#  SpD1  SpD2  SpD3  SpD4    WM
# 29766 27320 21793 20758 14290

## Pseudo-bulk data
spe_pseudo <- registration_pseudobulk(spe,
    var_registration = "BayesSpace",
    var_sample_id = "sample_id",
    min_ncells = 10
)
# 2025-01-27 11:59:27.811734 make pseudobulk object
# 2025-01-27 11:59:38.935869 dropping 0 pseudo-bulked samples that are below 'min_ncells'.
# 2025-01-27 11:59:38.966455 drop lowly expressed genes
# 2025-01-27 11:59:39.098746 normalize expression

summary(spe_pseudo$ncells)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 72.0   474.8   732.0   759.5   936.8  1930.0

tapply(spe_pseudo$ncells, spe_pseudo$BayesSpace, summary)
# $SpD1
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   209.0   671.5   984.0   992.2  1349.0  1930.0
#
# $SpD2
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   247.0   728.5   897.5   910.7  1083.8  1585.0
#
# $SpD3
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   225.0   586.8   724.5   726.4   793.5  1723.0
#
# $SpD4
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    72.0   471.0   732.0   691.9   831.8  1153.0
#
# $WM
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   165.0   269.8   391.0   476.3   502.2  1295.0

## Simplify the colData()  for the pseudo-bulked data
colData(spe_pseudo) <- colData(spe_pseudo)[, sort(c(
    "age",
    "sample_id",
    "BayesSpace",
    "subject",
    "sex",
    "position",
    "diagnosis",
    "ncells"
))]

## Compute PCs
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/f47daafa19b02e6208c7e0a9bc068367f806206c/code/analysis/09_region_differential_expression/preliminary_analysis.R#L60-L68
set.seed(20250127)
pca <- prcomp(t(assays(spe_pseudo)$logcounts))
message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
# $PCA_var_explained
#  [1] 16.30  9.24  4.78  4.22  3.57  3.36  2.84  2.45  2.36  2.16
# [11]  2.05  2.01  1.87  1.70  1.67  1.58  1.54  1.39  1.30  1.16
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

## Compute some reduced dims
set.seed(20250127)
spe_pseudo <- scater::runMDS(spe_pseudo, ncomponents = 20)
spe_pseudo <- scater::runPCA(spe_pseudo, name = "runPCA")

## Drop stuff we don't need for the spatialLIBD app
imgData(spe) <- imgData(spe)[imgData(spe)$image_id == "lowres", ]
## Check the size in GB
lobstr::obj_size(spe)
# 4.59 GB

## Drop the counts which take quite a bit of space
counts(spe) <- NULL
## Check the size in GB
lobstr::obj_size(spe)
# 2.45 GB

## Other app related changes
colors_BayesSpace <- Polychrome::palette36.colors(28)
names(colors_BayesSpace) <- unique(spe$BayesSpace)
m <- match(spe$BayesSpace, names(colors_BayesSpace))
stopifnot(all(!is.na(m)))
spe$BayesSpace_colors <- colors_BayesSpace[m]

m2 <- match(spe_pseudo$BayesSpace, names(colors_BayesSpace))
stopifnot(all(!is.na(m2)))
spe_pseudo$BayesSpace_colors <- colors_BayesSpace[m2]

## Save data
save(spe, file = file.path(dir_rdata, "spatialDLPFC_equivolume_domains.Rdata"))
saveRDS(spe_pseudo, file = file.path(dir_rdata, "spatialDLPFC_equivolume_domains_pseudobulk.Rds"))



## To avoid having to change parameters later on
spe_pseudo$registration_variable <- spe_pseudo$BayesSpace
spe_pseudo$registration_sample_id <- spe_pseudo$sample_id

## Set arguments used in spatialLIBD::registration_wrapper()
covars <- c("position", "age", "sex")
gene_ensembl <- "gene_id"
gene_name <- "gene_name"
suffix <- "all"

## Taken from spatialLIBD::registration_wrapper()
## https://github.com/LieberInstitute/spatialLIBD/blob/master/R/registration_wrapper.R
registration_mod <-
    registration_model(spe_pseudo, covars = covars)

block_cor <-
    registration_block_cor(spe_pseudo, registration_model = registration_mod)

results_enrichment <-
    registration_stats_enrichment(
        spe_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
results_pairwise <-
    registration_stats_pairwise(
        spe_pseudo,
        registration_model = registration_mod,
        block_cor = block_cor,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
results_anova <-
    registration_stats_anova(
        spe_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name,
        suffix = suffix
    )

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_enrichment,
    "pairwise" = results_pairwise
)


## Save the modeling results
saveRDS(
    modeling_results,
    file = file.path(
        dir_rdata,
        "spatialDLPFC_equivolume_domains_modeling_results.Rds"
    )
)


## Check that we have the right number of tests
k <- length(unique(spe_pseudo$BayesSpace))
tests <- lapply(modeling_results, function(x) {
    colnames(x)[grep("stat", colnames(x))]
})
stopifnot(length(tests$anova) == 1) ## assuming only "all"
stopifnot(length(tests$enrichment) == k)
stopifnot(length(tests$pairwise) == choose(k, 2))

## For sig_genes_extract_all() to work https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L37-L44
spe_pseudo$spatialLIBD <- spe_pseudo$BayesSpace

sig_genes <- sig_genes_extract_all(
    n = nrow(spe_pseudo),
    modeling_results = modeling_results,
    sce_layer = spe_pseudo
)

## Check that we have the right number of tests.
## the + 1 at the end assumes only "all"
stopifnot(length(unique(sig_genes$test)) == choose(k, 2) * 2 + k + 1)

lobstr::obj_size(sig_genes)
# 72.26 MB

dim(sig_genes)
# [1] 345592     13

## Drop parts we don't need to reduce the memory
sig_genes$in_rows <- NULL
sig_genes$in_rows_top20 <- NULL
lobstr::obj_size(sig_genes)
# 30.78 MB

## Extract FDR < 5%
## From
## https://github.com/LieberInstitute/brainseq_phase2/blob/be2b7f972bb2a0ede320633bf06abe1d4ef2c067/supp_tabs/create_supp_tables.R#L173-L181
fix_csv <- function(df) {
    for (i in seq_len(ncol(df))) {
        if (any(grepl(",", df[, i]))) {
            message(paste(Sys.time(), "fixing column", colnames(df)[i]))
            df[, i] <- gsub(",", ";", df[, i])
        }
    }
    return(df)
}
z <- fix_csv(as.data.frame(subset(sig_genes, fdr < 0.05)))
dim(z)
# [1] 136960     11
dim(subset(z, top <= 25))
# [1] 650  11
write.csv(
    subset(z, top <= 25),
    file = file.path(
        dir_rdata,
        "spatialDLPFC_equivolume_domains_model_results_FDR5perc_top25.csv"
    )
)

saveRDS(sig_genes, file = file.path(dir_rdata, "spatialDLPFC_equivolume_domains_sig_genes.Rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.4.2 Patched (2024-12-07 r87428)
#  os       macOS Sonoma 14.5
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2025-01-27
#  rstudio  2024.12.0+467 Kousa Dogwood (desktop)
#  pandoc   3.2 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version    date (UTC) lib source
#  abind                  1.4-8      2024-09-12 [1] CRAN (R 4.4.1)
#  AnnotationDbi          1.68.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  AnnotationHub          3.14.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  attempt                0.3.1      2020-05-03 [1] CRAN (R 4.4.0)
#  beachmat               2.22.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  beeswarm               0.4.0      2021-06-01 [1] CRAN (R 4.4.0)
#  benchmarkme            1.0.8      2022-06-12 [1] CRAN (R 4.4.0)
#  benchmarkmeData        1.0.4      2020-04-23 [1] CRAN (R 4.4.0)
#  Biobase              * 2.66.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocFileCache          2.14.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocGenerics         * 0.52.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocIO                 1.16.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocManager            1.30.25    2024-08-28 [1] CRAN (R 4.4.1)
#  BiocNeighbors          2.0.1      2024-11-28 [1] Bioconductor 3.20 (R 4.4.2)
#  BiocParallel           1.40.0     2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocSingular           1.22.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocVersion            3.20.0     2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
#  Biostrings             2.74.1     2024-12-16 [1] Bioconductor 3.20 (R 4.4.2)
#  bit                    4.5.0.1    2024-12-03 [1] CRAN (R 4.4.1)
#  bit64                  4.6.0-1    2025-01-16 [1] CRAN (R 4.4.1)
#  bitops                 1.0-9      2024-10-03 [1] CRAN (R 4.4.1)
#  blob                   1.2.4      2023-03-17 [1] CRAN (R 4.4.0)
#  brio                   1.1.5      2024-04-24 [1] CRAN (R 4.4.0)
#  bslib                  0.8.0      2024-07-29 [1] CRAN (R 4.4.0)
#  cachem                 1.1.0      2024-05-16 [1] CRAN (R 4.4.0)
#  circlize               0.4.16     2024-02-20 [1] CRAN (R 4.4.0)
#  cli                    3.6.3      2024-06-21 [1] CRAN (R 4.4.0)
#  clue                   0.3-66     2024-11-13 [1] CRAN (R 4.4.1)
#  cluster                2.1.8      2024-12-11 [1] CRAN (R 4.4.1)
#  codetools              0.2-20     2024-03-31 [1] CRAN (R 4.4.2)
#  colorout             * 1.3-2      2024-11-21 [1] Github (jalvesaq/colorout@2a5f214)
#  colorspace             2.1-1      2024-07-26 [1] CRAN (R 4.4.0)
#  ComplexHeatmap         2.22.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  config                 0.3.2      2023-08-30 [1] CRAN (R 4.4.0)
#  cowplot                1.1.3      2024-01-22 [1] CRAN (R 4.4.0)
#  crayon                 1.5.3      2024-06-20 [1] CRAN (R 4.4.0)
#  curl                   6.2.0      2025-01-23 [1] CRAN (R 4.4.2)
#  data.table             1.16.4     2024-12-06 [1] CRAN (R 4.4.1)
#  DBI                    1.2.3      2024-06-02 [1] CRAN (R 4.4.0)
#  dbplyr                 2.5.0      2024-03-19 [1] CRAN (R 4.4.0)
#  DelayedArray           0.32.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  devtools             * 2.4.5      2022-10-11 [1] CRAN (R 4.4.0)
#  digest                 0.6.37     2024-08-19 [1] CRAN (R 4.4.1)
#  doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.4.0)
#  dotCall64              1.2        2024-10-04 [1] CRAN (R 4.4.1)
#  dplyr                  1.1.4      2023-11-17 [1] CRAN (R 4.4.0)
#  DT                     0.33       2024-04-04 [1] CRAN (R 4.4.0)
#  edgeR                  4.4.1      2024-12-02 [1] Bioconductor 3.20 (R 4.4.2)
#  ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.4.0)
#  ExperimentHub          2.14.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  fastmap                1.2.0      2024-05-15 [1] CRAN (R 4.4.0)
#  fields                 16.3       2024-09-30 [1] CRAN (R 4.4.1)
#  filelock               1.0.3      2023-12-11 [1] CRAN (R 4.4.0)
#  foreach                1.5.2      2022-02-02 [1] CRAN (R 4.4.0)
#  fs                     1.6.5      2024-10-30 [1] CRAN (R 4.4.1)
#  gargle                 1.5.2      2023-07-20 [1] CRAN (R 4.4.0)
#  generics               0.1.3      2022-07-05 [1] CRAN (R 4.4.0)
#  GenomeInfoDb         * 1.42.1     2024-11-28 [1] Bioconductor 3.20 (R 4.4.2)
#  GenomeInfoDbData       1.2.13     2024-11-21 [1] Bioconductor
#  GenomicAlignments      1.42.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  GenomicRanges        * 1.58.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  GetoptLong             1.0.5      2020-12-15 [1] CRAN (R 4.4.0)
#  ggbeeswarm             0.7.2      2023-04-29 [1] CRAN (R 4.4.0)
#  ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.4.0)
#  ggrepel                0.9.6      2024-09-07 [1] CRAN (R 4.4.1)
#  GlobalOptions          0.1.2      2020-06-10 [1] CRAN (R 4.4.0)
#  glue                   1.8.0      2024-09-30 [1] CRAN (R 4.4.1)
#  golem                  0.5.1      2024-08-27 [1] CRAN (R 4.4.1)
#  googledrive            2.1.1      2023-06-11 [1] CRAN (R 4.4.0)
#  gridExtra              2.3        2017-09-09 [1] CRAN (R 4.4.0)
#  gtable                 0.3.6      2024-10-25 [1] CRAN (R 4.4.1)
#  here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.4.0)
#  hms                    1.1.3      2023-03-21 [1] CRAN (R 4.4.0)
#  htmltools              0.5.8.1    2024-04-04 [1] CRAN (R 4.4.0)
#  htmlwidgets            1.6.4      2023-12-06 [1] CRAN (R 4.4.0)
#  httpuv                 1.6.15     2024-03-26 [1] CRAN (R 4.4.0)
#  httr                   1.4.7      2023-08-15 [1] CRAN (R 4.4.0)
#  IRanges              * 2.40.1     2024-12-05 [1] Bioconductor 3.20 (R 4.4.2)
#  irlba                  2.3.5.1    2022-10-03 [1] CRAN (R 4.4.0)
#  iterators              1.0.14     2022-02-05 [1] CRAN (R 4.4.0)
#  jaffelab             * 0.99.34    2024-05-03 [1] Github (LieberInstitute/jaffelab@2710795)
#  jquerylib              0.1.4      2021-04-26 [1] CRAN (R 4.4.0)
#  jsonlite               1.8.9      2024-09-20 [1] CRAN (R 4.4.1)
#  KEGGREST               1.46.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  later                  1.4.1      2024-11-27 [1] CRAN (R 4.4.1)
#  lattice                0.22-6     2024-03-20 [1] CRAN (R 4.4.2)
#  lazyeval               0.2.2      2019-03-15 [1] CRAN (R 4.4.0)
#  lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.4.0)
#  limma                  3.62.2     2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
#  lobstr                 1.1.2      2022-06-22 [1] CRAN (R 4.4.0)
#  locfit                 1.5-9.10   2024-06-24 [1] CRAN (R 4.4.0)
#  lubridate              1.9.4      2024-12-08 [1] CRAN (R 4.4.1)
#  magick                 2.8.5      2024-09-20 [1] CRAN (R 4.4.1)
#  magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.4.0)
#  maps                   3.4.2.1    2024-11-10 [1] CRAN (R 4.4.1)
#  MASS                   7.3-64     2025-01-04 [1] CRAN (R 4.4.1)
#  Matrix                 1.7-1      2024-10-18 [1] CRAN (R 4.4.1)
#  MatrixGenerics       * 1.18.1     2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
#  matrixStats          * 1.5.0      2025-01-07 [1] CRAN (R 4.4.1)
#  memoise                2.0.1      2021-11-26 [1] CRAN (R 4.4.0)
#  mime                   0.12       2021-09-28 [1] CRAN (R 4.4.0)
#  miniUI                 0.1.1.1    2018-05-18 [1] CRAN (R 4.4.0)
#  munsell                0.5.1      2024-04-01 [1] CRAN (R 4.4.0)
#  nlme                   3.1-166    2024-08-14 [1] CRAN (R 4.4.2)
#  paletteer              1.6.0      2024-01-21 [1] CRAN (R 4.4.0)
#  pillar                 1.10.1     2025-01-07 [1] CRAN (R 4.4.1)
#  pkgbuild               1.4.6      2025-01-16 [1] CRAN (R 4.4.1)
#  pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.4.0)
#  pkgload                1.4.0      2024-06-28 [1] CRAN (R 4.4.0)
#  plotly                 4.10.4     2024-01-13 [1] CRAN (R 4.4.0)
#  png                    0.1-8      2022-11-29 [1] CRAN (R 4.4.0)
#  Polychrome             1.5.1      2022-05-03 [1] CRAN (R 4.4.0)
#  prettyunits            1.2.0      2023-09-24 [1] CRAN (R 4.4.0)
#  profvis                0.4.0      2024-09-20 [1] CRAN (R 4.4.1)
#  promises               1.3.2      2024-11-28 [1] CRAN (R 4.4.1)
#  prompt                 1.0.2.9000 2024-05-03 [1] Github (gaborcsardi/prompt@17bd0e1)
#  purrr                  1.0.2      2023-08-10 [1] CRAN (R 4.4.0)
#  R6                     2.5.1      2021-08-19 [1] CRAN (R 4.4.0)
#  rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.4.0)
#  rappdirs               0.3.3      2021-01-31 [1] CRAN (R 4.4.0)
#  RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.4.0)
#  Rcpp                   1.0.14     2025-01-12 [1] CRAN (R 4.4.1)
#  RCurl                  1.98-1.16  2024-07-11 [1] CRAN (R 4.4.0)
#  rematch2               2.1.2      2020-05-01 [1] CRAN (R 4.4.0)
#  remotes                2.5.0      2024-03-17 [1] CRAN (R 4.4.0)
#  restfulr               0.0.15     2022-06-16 [1] CRAN (R 4.4.0)
#  rjson                  0.2.23     2024-09-16 [1] CRAN (R 4.4.1)
#  rlang                  1.1.4      2024-06-04 [1] CRAN (R 4.4.0)
#  rprojroot              2.0.4      2023-11-05 [1] CRAN (R 4.4.0)
#  Rsamtools              2.22.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  RSQLite                2.3.9      2024-12-03 [1] CRAN (R 4.4.1)
#  rsthemes               0.5.0      2024-11-21 [1] Github (gadenbuie/rsthemes@96fb91a)
#  rstudioapi             0.17.1     2024-10-22 [1] CRAN (R 4.4.1)
#  rsvd                   1.0.5      2021-04-16 [1] CRAN (R 4.4.0)
#  rtracklayer            1.66.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  S4Arrays               1.6.0      2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  S4Vectors            * 0.44.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  sass                   0.4.9.9000 2024-05-03 [1] Github (rstudio/sass@9228fcf)
#  ScaledMatrix           1.14.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  scales                 1.3.0      2023-11-28 [1] CRAN (R 4.4.0)
#  scater               * 1.34.0     2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
#  scatterplot3d          0.3-44     2023-05-05 [1] CRAN (R 4.4.0)
#  scuttle              * 1.16.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  segmented              2.1-3      2024-10-25 [1] CRAN (R 4.4.1)
#  sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.4.0)
#  shape                  1.4.6.1    2024-02-23 [1] CRAN (R 4.4.0)
#  shiny                  1.10.0     2024-12-14 [1] CRAN (R 4.4.1)
#  shinyWidgets           0.8.7      2024-09-23 [1] CRAN (R 4.4.1)
#  SingleCellExperiment * 1.28.1     2024-11-10 [1] Bioconductor 3.20 (R 4.4.1)
#  spam                   2.11-0     2024-10-03 [1] CRAN (R 4.4.1)
#  SparseArray            1.6.1      2025-01-20 [1] Bioconductor 3.20 (R 4.4.2)
#  SpatialExperiment    * 1.16.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  spatialLIBD          * 1.19.0     2024-12-12 [1] Bioconductor
#  statmod                1.5.0      2023-01-06 [1] CRAN (R 4.4.0)
#  stringi                1.8.4      2024-05-06 [1] CRAN (R 4.4.0)
#  stringr                1.5.1      2023-11-14 [1] CRAN (R 4.4.0)
#  SummarizedExperiment * 1.36.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  suncalc                0.5.1      2022-09-29 [1] CRAN (R 4.4.0)
#  testthat             * 3.2.3      2025-01-13 [1] CRAN (R 4.4.1)
#  tibble                 3.2.1      2023-03-20 [1] CRAN (R 4.4.0)
#  tidyr                  1.3.1      2024-01-24 [1] CRAN (R 4.4.0)
#  tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.4.0)
#  timechange             0.3.0      2024-01-18 [1] CRAN (R 4.4.0)
#  UCSC.utils             1.2.0      2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  urlchecker             1.0.1      2021-11-30 [1] CRAN (R 4.4.0)
#  usethis              * 3.1.0      2024-11-26 [1] CRAN (R 4.4.1)
#  vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.4.0)
#  vipor                  0.4.7      2023-12-18 [1] CRAN (R 4.4.0)
#  viridis                0.6.5      2024-01-29 [1] CRAN (R 4.4.0)
#  viridisLite            0.4.2      2023-05-02 [1] CRAN (R 4.4.0)
#  withr                  3.0.2      2024-10-28 [1] CRAN (R 4.4.1)
#  XML                    3.99-0.18  2025-01-01 [1] CRAN (R 4.4.1)
#  xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.4.0)
#  XVector                0.46.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  yaml                   2.3.10     2024-07-26 [1] CRAN (R 4.4.0)
#  zlibbioc               1.52.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#
#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
