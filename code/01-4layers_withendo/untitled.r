library("spatialLIBD")
#library("escheR")
library("dplyr")
library("tidyr")
library("stringr")
library("tibble")
library("ggplot2")


# Load spatial data
spe <- fetch_data(type = "spatialDLPFC_Visium")

# Recode cluster labels
tmp <- as.data.frame(colData(spe)) %>%
  mutate(BayesSpace_harmony_07 = recode(BayesSpace_harmony_07,
                                        `1` = "SpD1",
                                        `2` = "SpD1",
                                        `3` = "SpD2",
                                        `4` = "SpD3",
                                        `5` = "WM",
                                        `6` = "SpD4",
                                        `7` = "WM"))

# Assign back to colData
colData(spe) <- DataFrame(tmp)

# Optional: Set factor levels for consistent legend order
colData(spe)$BayesSpace_harmony_07 <- factor(
  colData(spe)$BayesSpace_harmony_07,
  levels = c("SpD1", "SpD2", "SpD3", "SpD4", "WM")
)

spD_colors <- c(
  "SpD1" = "#ED5E70",  # blue
  "SpD2" = "#CCBB44",  # orange
  "SpD3" = "#238934",  # green
  "SpD4" = "#66CCEE",  # red
  "WM"   = "#C3B3A1"   # purple
)

vis_clus(
  spe = spe,
  clustervar = "BayesSpace_harmony_07",
  sampleid = "Br8667_mid",
  colors = spD_colors
)

ggsave(
  filename=paste0("~/MinaRyten/Melissa/MRIxST/plots/01-4layers_withendo/00-spatial_domain_definition/Br8667_mid_4layers.png"),
  device ="png",
  scale = 1,
  width = 7,
  height = 7,
  units = c("in"),
  dpi = 300
)

vis_clus(
  spe = spe,
  clustervar = "BayesSpace_harmony_07",
  sampleid = "Br6432_ant",
  colors = spD_colors
)

ggsave(
  filename=paste0("~/MinaRyten/Melissa/MRIxST/plots/01-4layers_withendo/00-spatial_domain_definition/Br6432_ant_4layers.png"),
  device ="png",
  scale = 1,
  width = 7,
  height = 7,
  units = c("in"),
  dpi = 300
)

vis_clus(
  spe = spe,
  clustervar = "BayesSpace_harmony_07",
  sampleid = "Br6522_ant",
  colors = spD_colors
)

ggsave(
  filename=paste0("~/MinaRyten/Melissa/MRIxST/plots/01-4layers_withendo/00-spatial_domain_definition/Br6522_ant_4layers.png"),
  device ="png",
  scale = 1,
  width = 7,
  height = 7,
  units = c("in"),
  dpi = 300
)
