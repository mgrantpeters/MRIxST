library("rsconnect")
library("here")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
# load(here("code", "deploy_app_k09", ".deploy_info.Rdata"), verbose = TRUE)
# rsconnect::setAccountInfo(
#     name = deploy_info$name,
#     token = deploy_info$token,
#     secret = deploy_info$secret
# )

## Run only once to create the symbolic links
withr::with_dir(
    here("code", "06_spatialDLPFC_equivolume_domains_spatiaLIBD_app"),
    sapply(
        c(
            "ln -s ../../processed_data/05_spatialDLPFC_equivolume_domains/spatialDLPFC_equivolume_domains.Rdata spatialDLPFC_equivolume_domains.Rdata",
            "ln -s ../../processed_data/05_spatialDLPFC_equivolume_domains/spatialDLPFC_equivolume_domains_pseudobulk.Rds spatialDLPFC_equivolume_domains_pseudobulk.Rds",
            "ln -s ../../processed_data/05_spatialDLPFC_equivolume_domains/spatialDLPFC_equivolume_domains_modeling_results.Rds spatialDLPFC_equivolume_domains_modeling_results.Rds",
            "ln -s ../../processed_data/05_spatialDLPFC_equivolume_domains/spatialDLPFC_equivolume_domains_sig_genes.Rds spatialDLPFC_equivolume_domains_sig_genes.Rds"
        ), system
    )
)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = here("code", "06_spatialDLPFC_equivolume_domains_spatiaLIBD_app"),
    appFiles = c(
        "app.R",
        "spatialDLPFC_equivolume_domains.Rdata",
        "spatialDLPFC_equivolume_domains_pseudobulk.Rds",
        "spatialDLPFC_equivolume_domains_modeling_results.Rds",
        "spatialDLPFC_equivolume_domains_sig_genes.Rds",
        withr::with_dir(here("code", "06_spatialDLPFC_equivolume_domains_spatiaLIBD_app"), dir("www", full.names = TRUE))
    ),
    appName = "spatialDLPFC_equivolume_domains",
    account = "libd",
    server = "shinyapps.io"
)
