# script to do expanded pathway analysis

# load in some libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(hydroGOF)
library(matrixStats)
library(dplyr)
library(scales)
library(viridis)
library(gridExtra)

# define directories
d0 <- getwd();
qMRI.data.dir <- paste(d0, "/../../processed_data/02-4layers_qMRI/", sep = "")
cell.data.dir <- paste(d0, "/../../processed_data/01a-4layers_noendo/02-top3_cells_counts/", sep = "")
GO.data.dir <- paste(d0, "/../../processed_data/01a-4layers_noendo/05-pathway_scoring/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir <- paste(d0, "/../../plots/02a-4layers_qMRIvsST_noendo", sep = "")

# make a function to rename layers like so:
# L4 --> SpD1
# L3 --> SpD2
# L2 --> SpD3
# L1 --> SpD4
rename_layers <- function(x) {
  x <- gsub("_L1", "_X", x)  
  x <- gsub("_L2", "_Y", x)  
  x <- gsub("_L3", "_Z", x)  
  x <- gsub("_L4", "_SpD1", x) 
  x <- gsub("_Z", "_SpD2", x) 
  x <- gsub("_Y", "_SpD3", x)  
  x <- gsub("_X", "_SpD4", x)  
  return(x)
}

#### Load in & reorganise qMRI data ####

# load in the qMRI data and demographics for the DLPFC sub-regions
qMRI.dems <- read.csv(paste(dems.dir,"anon_DLPFC_qMRI_demographics.csv",sep = ""))
qMRI.data <-  read.csv(paste(qMRI.data.dir,"anon_DLPFC_AMP_qMRI_medians_4_equivol_layers.csv",sep = ""))

# set age and diagnosis as factors
qMRI.dems$sex <- factor(qMRI.dems$sex)
qMRI.dems$diagnosis <- factor(qMRI.dems$diagnosis,levels = c("Control","PD"))

# re-format these so that the sub-regions are in different columns (for regression analysis)
# anterior
qMRI.data.ant <- qMRI.data[grepl("ant", qMRI.data$subject_region),2:ncol(qMRI.data)]
colnames(qMRI.data.ant) <- paste(colnames(qMRI.data.ant),"_ant", sep = "")
# middle
qMRI.data.mid <- qMRI.data[grepl("mid", qMRI.data$subject_region),2:ncol(qMRI.data)]
colnames(qMRI.data.mid) <- paste(colnames(qMRI.data.mid),"_mid", sep = "")
# posterior
qMRI.data.pos <- qMRI.data[grepl("pos", qMRI.data$subject_region),2:ncol(qMRI.data)]
colnames(qMRI.data.pos) <- paste(colnames(qMRI.data.pos),"_pos", sep = "")

# merge them together and add the dems
qMRI.data.dems <- cbind(qMRI.dems,qMRI.data.ant, qMRI.data.mid, qMRI.data.pos)
rownames(qMRI.data.dems) <- 1:nrow(qMRI.data.dems)

# load in the QC info sheet 
QC.info <- read.csv(paste(dems.dir,"anon_7T_QC_info.csv",sep = ""))

# specify the metrics we are interested in
qMRI.metrics <- c("QSM","R1","R2s","PD","MTsat")
MOI <- c("QSM","R1","R2s","MTsat")

# match up QC and qMRI data frames
QC.info <- subset(QC.info, QC.info$subject %in% qMRI.data.dems$subject)
rownames(QC.info) <- 1:nrow(QC.info)

# remove data from demographics and qMRI data using QC info
# if they failed MPRAGE QC or FS QC, then remove from everywhere
qMRI.data.dems <- subset(qMRI.data.dems, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
QC.info <- subset(QC.info, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)

# rename columns
colnames(qMRI.data.dems) <- lapply(colnames(qMRI.data.dems), rename_layers)
# get all ROI names
qMRI.ROI.names <- colnames(qMRI.data.dems)[5:ncol(qMRI.data.dems)]

# need to create some new columns in the QC.info dataframe
QC.info$QSM_QC <- QC.info$QSM_QC == 1 & QC.info$has_QSM == 1
QC.info$R1_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$R2s_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$PD_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$MTsat_QC <- QC.info$MPM_QC == 1 & QC.info$has_MT == 1

#### adjust qMRI values ####

# initialize a data frame for the adjusted values
adj.qMRI.data.dems <- qMRI.data.dems

# get mean and diff from mean age and sex for all subjects
mn.qMRI.age <- mean(qMRI.data.dems$age)
diff.qMRI.age <- mn.qMRI.age - qMRI.data.dems$age
mn.qMRI.sex <- mean(as.numeric(qMRI.data.dems$sex))
diff.qMRI.sex <- mn.qMRI.sex - as.numeric(qMRI.data.dems$sex)

# adjust all ROIs for age and sex using the covariance method
for (ROI in qMRI.ROI.names) {
  mdl.qMRI.adj <- lm(qMRI.data.dems[,ROI] ~ age + sex, data = qMRI.data.dems)
  adj.qMRI.vals <- qMRI.data.dems[,ROI] +
    diff.qMRI.age * mdl.qMRI.adj$coefficients[2] +
    diff.qMRI.sex * mdl.qMRI.adj$coefficients[3]
  adj.qMRI.data.dems[,ROI] <- adj.qMRI.vals
}


#### load in & reorganise GO data ####

# load in ST patient dems
ST.dems <- read.csv(paste(dems.dir,"ST_patient_info.csv",sep = ""))
ST.dems <- ST.dems[order(ST.dems$sample_id),]
ST.dems$sex <- factor(ST.dems$sex)

# get names of 4 layer pathway files
GO.fnames <- dir(GO.data.dir,"scores_4layers_genes_(.+)")

# use this to define some GO type names
GO.codes <- substr(GO.fnames,nchar(GO.fnames)-12,nchar(GO.fnames)-4)
GO.descs <- tolower(substr(GO.fnames, 22, nchar(GO.fnames) - 23))
GO.names <- paste(GO.descs,GO.codes,sep = "_")

# put these all into one data frame with the dems
GO.data.dems <- ST.dems

# do some reorganising and put them into wide format for adjustment
for (i in 1:length(GO.fnames)){
  GO.fname <- GO.fnames[i]
  GO.data <- read.csv(paste(GO.data.dir, GO.fname, sep = ""))
  GO.data <- GO.data[,2:ncol(GO.data)]
  colnames(GO.data) <- c("sample_id", "region", GO.names[i])
  GO.data <- GO.data[order(
    GO.data$sample_id,
    factor(GO.data$region, levels = c("WM", "SpD4","SpD3","SpD2","SpD1"))),]
  GO.data <- pivot_wider(GO.data,names_from = region, values_from = GO.names[i],
                         names_prefix = paste(GO.names[i],"_", sep = ""))
  GO.data.dems <- cbind(GO.data.dems,GO.data[,2:ncol(GO.data)])
}
colnames(GO.data.dems) <- gsub("_WM", "_SWM", colnames(GO.data.dems))
#colnames(GO.data.dems) <- lapply(colnames(GO.data.dems), rename_SpDs)

# reorganise in the same way as the qMRI data
GO.data.ant <- GO.data.dems[ST.dems$region == "anterior", 7:ncol(GO.data.dems)]
colnames(GO.data.ant) <- paste(colnames(GO.data.ant),"_ant", sep = "")
GO.data.mid <- GO.data.dems[ST.dems$region == "middle", 7:ncol(GO.data.dems)]
colnames(GO.data.mid) <- paste(colnames(GO.data.mid),"_mid", sep = "")
GO.data.pos <- GO.data.dems[ST.dems$region == "posterior", 7:ncol(GO.data.dems)]
colnames(GO.data.pos) <- paste(colnames(GO.data.pos),"_pos", sep = "")

# merge together with the demographics
GO.data.dems <- cbind(ST.dems[seq(1,nrow(ST.dems),3),], GO.data.ant, GO.data.mid, GO.data.pos)
rownames(GO.data.dems) <- 1:nrow(GO.data.dems)
GO.data.dems <- subset(GO.data.dems, select = -c(sample_id, region))

# remove columns with all zero or all nan values
GO.data.dems <- 
  GO.data.dems[, colSums(GO.data.dems == 0) != nrow(GO.data.dems) & 
                 colSums(is.na(GO.data.dems)) != nrow(GO.data.dems)]
# from this, work out which GO.names should be excluded
GO.names <- unique((regmatches(
  colnames(GO.data.dems), gregexpr("\\w*GO\\w{7}", colnames(GO.data.dems))
)))
GO.names <- GO.names[2:length(GO.names)]

#### adjust ST data ####

# initialize data frame for the adjusted values
adj.GO.data.dems <- GO.data.dems

# get mean and diff from mean age and sex for all subjects
mn.ST.age <- mean(GO.data.dems$age)
diff.ST.age <- mn.ST.age - GO.data.dems$age
mn.ST.sex <- mean(as.numeric(GO.data.dems$sex))
diff.ST.sex <- mn.ST.sex - as.numeric(GO.data.dems$sex)

# adjust pathway ST data
GO.ROI.names <- colnames(GO.data.dems)[5:ncol(GO.data.dems)]
for (ROI in GO.ROI.names) {
  mdl.GO.adj <- lm(GO.data.dems[,ROI] ~ age + sex, data = GO.data.dems)
  adj.GO.vals <- GO.data.dems[,ROI] +
    diff.ST.age * mdl.GO.adj$coefficients[2] +
    diff.ST.sex * mdl.GO.adj$coefficients[3]
  adj.GO.data.dems[,ROI] <- adj.GO.vals
}

#### calculate average across subjects ####

# making columns
dom <- c("SWM", "SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

adj.qMRI.GO.medians <- expand.grid(SpD = dom, region = reg)

# go through each qMRI metric and remove QC-failed subjects / nan values before averaging
qMRI.metrics = c("QSM","R1","R2s","MTsat")

# will still want to know pathways of interest
POI.names = c("inflammation_GO0006954", "iron_GO0060586", 
              "lipid_GO0006629", "myelin_GO0043209")

# do this for controls and PD separately
# (although only looking at Ctrl vs Ctrl for now)
adj.qMRI.data.dems.Ctrl <- subset(adj.qMRI.data.dems, diagnosis == "Control")
QC.info.Ctrl <- subset(QC.info, adj.qMRI.data.dems$diagnosis == "Control" )
adj.qMRI.data.dems.PD <- subset(adj.qMRI.data.dems, diagnosis == "PD")
QC.info.PD <- subset(QC.info, adj.qMRI.data.dems$diagnosis == "PD" )

for (m in qMRI.metrics){
  # remove QC failed subs before calculating mean
  metric.data.Ctrl <- adj.qMRI.data.dems.Ctrl[QC.info.Ctrl[,paste(m,"_QC",sep = "")], 
                                              grepl(m, colnames(adj.qMRI.data.dems))]
  metric.medians.Ctrl <- colMedians(as.matrix(metric.data.Ctrl),na.rm = TRUE)
  metric.data.PD <- adj.qMRI.data.dems.PD[QC.info.PD[,paste(m,"_QC",sep = "")], 
                                          grepl(m, colnames(adj.qMRI.data.dems))]
  metric.medians.PD <- colMedians(as.matrix(metric.data.PD),na.rm = TRUE)
  adj.qMRI.GO.medians[,paste(m, "_Control", sep = "")] <- metric.medians.Ctrl
  adj.qMRI.GO.medians[,paste(m, "_PD", sep = "")] <- metric.medians.PD
}

# removing neuronal_death and mineral (for now) as these were not in the main analysis
adj.GO.data.dems <- adj.GO.data.dems[,
                                     !(grepl("neuron_death", colnames(adj.GO.data.dems)) | 
                                         grepl("mineral",colnames(adj.GO.data.dems)))]
GO.names <- GO.names[!(grepl("neuron_death", GO.names) | 
                         grepl("mineral",GO.names))]

for (c in GO.names){
  GO.data <- adj.GO.data.dems[,grepl(paste(c,"_",sep = ""), colnames(adj.GO.data.dems))]
  GO.medians <- colMedians(as.matrix(GO.data), na.rm = TRUE)
  adj.qMRI.GO.medians[,paste(c, "_Control", sep = "")] <- GO.medians
}

#### associations with SWM ####

# calculate associations between all GO pathways and all qMRI metrics
all.rho.vals <- data.frame(matrix(nrow = length(GO.names), ncol = length(qMRI.metrics)))
rownames(all.rho.vals) <- GO.names
colnames(all.rho.vals) <- qMRI.metrics
all.R2.vals <- all.rho.vals
all.pvals <- all.rho.vals

for (metric in qMRI.metrics) {
  for (GO in GO.names) {
    # Compute correlation
    namemetric <- paste(metric, "_Control", sep = "")
    nameGO <- paste(GO, "_Control", sep = "")
    res <- cor.test(adj.qMRI.GO.medians[[namemetric]], adj.qMRI.GO.medians[[nameGO]], method = "pearson")
    all.rho.vals[GO, metric] <- res$estimate
    all.R2.vals[GO,metric] <- res$estimate^2
    all.pvals[GO,metric] <- res$p.value
  }
}

# add in a column indicating GO term category
all.rho.vals$GO.category <- substr(GO.names, 1, nchar(GO.names) - 10)
all.R2.vals$GO.category <- all.rho.vals$GO.category

# put these into long form
all.rho.vals.melt <- melt(all.rho.vals)
all.R2.vals.melt <- melt(all.R2.vals)
colnames(all.rho.vals.melt) <- c("GO.category", "qMRI.metric", "rho")
colnames(all.R2.vals.melt) <- c("GO.category", "qMRI.metric", "R2")

GO.categories <- unique(all.R2.vals.melt$GO.category)

# save some stats
write.csv(all.R2.vals, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_R2.csv", sep = ""), 
  row.names = TRUE)
write.csv(all.rho.vals, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_rho.csv", sep = ""), 
  row.names = TRUE)
write.csv(all.pvals, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_pval.csv", sep = ""), 
  row.names = TRUE)

# NOW DO A FANCY PLOT

# create 2x2 grid of plots for each qMRI metric
plot.list.by.qMRI <- list()

# find the specific R2 values to be plotted (for POI)

# Define the function
plot.density.by.qMRI <- function(qMRI.factor, all.R2.melt, all.R2, ST) {
  
  # Filter data for the current qMRI metric
  df.qMRI <- subset(all.R2.melt, qMRI.metric == qMRI.factor)
  
  # Calculate density for each GO category
  qMRI.density.data <- df.qMRI %>%
    group_by(GO.category) %>%
    do({
      dens <- density(.$R2, from = 0, to = 1)
      data.frame(
        R2 = dens$x,  # R2 values from density calculation
        density = rescale(dens$y, to = c(0, 1)),  # density values (normalized)
        ST = rep(unique(.$GO.category), length(dens$x))  # add ST variable info 
      )
    }) %>%
    ungroup()
  
  # get 
  POI.data <- data.frame(
    x = 1:4, 
    y = sapply(1:4, function(i) all.R2[ST[i], qMRI.factor])
  )
  
  # Create the plot
  ggplot(qMRI.density.data, aes(x = GO.category, y = R2, fill = density)) +
    geom_tile(width = 0.8, height = 0.02) +  # bar with constant size
    scale_fill_gradient(low = "blue", high = "red", limits = c(0, 1)) +  # color by (normalized) density value
    geom_point(data = POI.data, aes(x = x, y = y), col = "white", inherit.aes = FALSE) +
    labs(title = qMRI.factor, x = element_blank(), y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "none")
}


plot.list.by.qMRI <- lapply(1:length(qMRI.metrics), function(i) 
  plot.density.by.qMRI(qMRI.factor = qMRI.metrics[i],
                       all.R2.melt = all.R2.vals.melt,
                       all.R2 = all.R2.vals,
                       ST = POI.names))

png(filename = paste(results.dir, "/Plot_expanded_pathway_analysis_by_qMRI.png", sep = ""),
    width = 800, height = 800, units = "px", res = 120)
multi.plot.by.qMRI <- ggarrange(plotlist = plot.list.by.qMRI, 
                   nrow = 2, 
                   ncol = 2, 
                   common.legend = TRUE, 
                   legend = "bottom")
plot(multi.plot.by.qMRI)
annotate_figure(multi.plot.by.qMRI, 
                top = text_grob(expression("Density of qMRI vs ST R"^2*" values")),
                left = text_grob(expression("R"^2), rot = 90))
dev.off()


# create 2x2 grid of plots for each GO category
plot.list.by.GO <- list()

# find the specific R2 values to be plotted (for POI)

# Define the function
plot.density.by.GO <- function(GO.factor, all.R2.melt, all.R2, qMRI, POI) {
  
  # Filter data for the current qMRI metric
  df.GO <- subset(all.R2.melt, GO.category == GO.factor)
  
  # Calculate density for each GO category
  GO.density.data <- df.GO %>%
    group_by(qMRI.metric) %>%
    do({
      dens <- density(.$R2, from = 0, to = 1)
      data.frame(
        R2 = dens$x,  # R2 values from density calculation
        density = rescale(dens$y, to = c(0, 1)),  # density values (normalized)
        qMRI = rep(unique(.$qMRI.metric), length(dens$x))  # add ST variable info 
      )
    }) %>%
    ungroup()
  
  # get 
  POI.data <- data.frame(
    x = 1:4, 
    y = sapply(1:4, function(i) all.R2[POI, qMRI[i]])
  )
  
  # Create the plot
  ggplot(GO.density.data, aes(x = qMRI.metric, y = R2, fill = density)) +
    geom_tile(width = 0.8, height = 0.02) +  # bar with constant size
    scale_fill_gradient(low = "blue", high = "red", limits = c(0, 1)) +  # color by (normalized) density value
    geom_point(data = POI.data, aes(x = x, y = y), col = "white", inherit.aes = FALSE) +
    labs(title = GO.factor, x = element_blank(), y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "none")
}


plot.list.by.GO <- lapply(1:length(GO.categories), function(i) 
  plot.density.by.GO(GO.factor = GO.categories[i],
                     all.R2.melt = all.R2.vals.melt,
                     all.R2 = all.R2.vals,
                     qMRI = qMRI.metrics, 
                     POI = POI.names[i]))

png(filename = paste(results.dir, "/Plot_expanded_pathway_analysis_by_GO.png", sep = ""),
    width = 800, height = 800, units = "px", res = 120)
multi.plot.by.GO <- ggarrange(plotlist = plot.list.by.GO, 
                                nrow = 2, 
                                ncol = 2, 
                                common.legend = TRUE, 
                                legend = "bottom")
plot(multi.plot.by.GO)
annotate_figure(multi.plot.by.GO, 
                top = text_grob(expression("Density of qMRI vs ST R"^2*" values")),
                left = text_grob(expression("R"^2), rot = 90))
dev.off()

#### associations without SWM ####

# remove SWM
adj.qMRI.GO.medians.no.SWM <- subset(adj.qMRI.GO.medians, SpD != "SWM")

# calculate associations between all GO pathways and all qMRI metrics
all.rho.vals.no.SWM <- data.frame(matrix(nrow = length(GO.names), ncol = length(qMRI.metrics)))
rownames(all.rho.vals.no.SWM) <- GO.names
colnames(all.rho.vals.no.SWM) <- qMRI.metrics
all.R2.vals.no.SWM <- all.rho.vals.no.SWM
all.pvals.no.SWM <- all.rho.vals.no.SWM

for (metric in qMRI.metrics) {
  for (GO in GO.names) {
    # Compute correlation
    namemetric <- paste(metric, "_Control", sep = "")
    nameGO <- paste(GO, "_Control", sep = "")
    res <- cor.test(adj.qMRI.GO.medians.no.SWM[[namemetric]], adj.qMRI.GO.medians.no.SWM[[nameGO]], method = "pearson")
    all.rho.vals.no.SWM[GO, metric] <- res$estimate
    all.R2.vals.no.SWM[GO,metric] <- res$estimate^2
    all.pvals.no.SWM[GO,metric] <- res$p.value
  }
}

# add in a column indicating GO term category
all.rho.vals.no.SWM$GO.category <- substr(GO.names, 1, nchar(GO.names) - 10)
all.R2.vals.no.SWM$GO.category <- all.rho.vals.no.SWM$GO.category

# put these into long form
all.rho.vals.no.SWM.melt <- melt(all.rho.vals.no.SWM)
all.R2.vals.no.SWM.melt <- melt(all.R2.vals.no.SWM)
colnames(all.rho.vals.no.SWM.melt) <- c("GO.category", "qMRI.metric", "rho")
colnames(all.R2.vals.no.SWM.melt) <- c("GO.category", "qMRI.metric", "R2")

# save some stats
write.csv(all.R2.vals.no.SWM, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_R2_no_SWM.csv", sep = ""), 
  row.names = TRUE)
write.csv(all.rho.vals.no.SWM, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_rho_no_SWM.csv", sep = ""), 
  row.names = TRUE)
write.csv(all.pvals.no.SWM, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_pval_no_SWM.csv", sep = ""), 
  row.names = TRUE)


# NOW DO A FANCY PLOT

# create 2x2 grid of plots for each qMRI metric
plot.list.by.qMRI.no.SWM <- list()

plot.list.by.qMRI.no.SWM <- lapply(1:length(qMRI.metrics), function(i) 
  plot.density.by.qMRI(qMRI.factor = qMRI.metrics[i],
                       all.R2.melt = all.R2.vals.no.SWM.melt,
                       all.R2 = all.R2.vals.no.SWM,
                       ST = POI.names))

png(filename = paste(results.dir, "/Plot_expanded_pathway_analysis_by_qMRI_no_SWM.png", sep = ""),
    width = 800, height = 800, units = "px", res = 120)
multi.plot.by.qMRI.no.SWM <- ggarrange(plotlist = plot.list.by.qMRI.no.SWM, 
                   nrow = 2, 
                   ncol = 2, 
                   common.legend = TRUE, 
                   legend = "bottom")
plot(multi.plot.by.qMRI.no.SWM)
annotate_figure(multi.plot.by.qMRI.no.SWM, 
                top = text_grob(expression("Density of qMRI vs ST R"^2*" values")),
                left = text_grob(expression("R"^2), rot = 90))
dev.off()


# create 2x2 grid of plots for each GO category
plot.list.by.GO.no.SWM <- list()

plot.list.by.GO.no.SWM <- lapply(1:length(GO.categories), function(i) 
  plot.density.by.GO(GO.factor = GO.categories[i],
                     all.R2.melt = all.R2.vals.no.SWM.melt,
                     all.R2 = all.R2.vals.no.SWM,
                     qMRI = qMRI.metrics, 
                     POI = POI.names[i]))

png(filename = paste(results.dir, "/Plot_expanded_pathway_analysis_by_GO_no_SWM.png", sep = ""),
    width = 800, height = 800, units = "px", res = 120)
multi.plot.by.GO.no.SWM <- ggarrange(plotlist = plot.list.by.GO.no.SWM, 
                                nrow = 2, 
                                ncol = 2, 
                                common.legend = TRUE, 
                                legend = "bottom")
plot(multi.plot.by.GO.no.SWM)
annotate_figure(multi.plot.by.GO.no.SWM, 
                top = text_grob(expression("Density of qMRI vs ST R"^2*" values")),
                left = text_grob(expression("R"^2), rot = 90))
dev.off()


#### associations without SWM or SpD1 ####

# remove SpD1
adj.qMRI.GO.medians.no.SWMSpD1 <- subset(adj.qMRI.GO.medians.no.SWM, SpD != "SpD1")

# calculate associations between all GO pathways and all qMRI metrics
all.rho.vals.no.SWMSpD1 <- data.frame(matrix(nrow = length(GO.names), ncol = length(qMRI.metrics)))
rownames(all.rho.vals.no.SWMSpD1) <- GO.names
colnames(all.rho.vals.no.SWMSpD1) <- qMRI.metrics
all.R2.vals.no.SWMSpD1 <- all.rho.vals.no.SWMSpD1
all.pvals.no.SWMSpD1 <- all.rho.vals.no.SWMSpD1

for (metric in qMRI.metrics) {
  for (GO in GO.names) {
    # Compute correlation
    namemetric <- paste(metric, "_Control", sep = "")
    nameGO <- paste(GO, "_Control", sep = "")
    res <- cor.test(adj.qMRI.GO.medians.no.SWMSpD1[[namemetric]], adj.qMRI.GO.medians.no.SWMSpD1[[nameGO]], method = "pearson")
    all.rho.vals.no.SWMSpD1[GO, metric] <- res$estimate
    all.R2.vals.no.SWMSpD1[GO,metric] <- res$estimate^2
    all.pvals.no.SWMSpD1[GO,metric] <- res$p.value
  }
}

# add in a column indicating GO term category
all.rho.vals.no.SWMSpD1$GO.category <- substr(GO.names, 1, nchar(GO.names) - 10)
all.R2.vals.no.SWMSpD1$GO.category <- all.rho.vals.no.SWMSpD1$GO.category

# put these into long form
all.rho.vals.no.SWMSpD1.melt <- melt(all.rho.vals.no.SWMSpD1)
all.R2.vals.no.SWMSpD1.melt <- melt(all.R2.vals.no.SWMSpD1)
colnames(all.rho.vals.no.SWMSpD1.melt) <- c("GO.category", "qMRI.metric", "rho")
colnames(all.R2.vals.no.SWMSpD1.melt) <- c("GO.category", "qMRI.metric", "R2")

# save some stats
write.csv(all.R2.vals.no.SWMSpD1, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_R2_no_SWM_SpD1.csv", sep = ""), 
  row.names = TRUE)
write.csv(all.rho.vals.no.SWMSpD1, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_rho_no_SWM_SpD1.csv", sep = ""), 
  row.names = TRUE)
write.csv(all.pvals.no.SWMSpD1, file = paste(
  results.dir, "/Stats_expanded_pathway_analysis_pval_no_SWM_SpD1.csv", sep = ""), 
  row.names = TRUE)

# NOW DO A FANCY PLOT

# create 2x2 grid of plots for each qMRI metric
plot.list.by.qMRI.no.SWMSpD1 <- list()

plot.list.by.qMRI.no.SWMSpD1 <- lapply(1:length(qMRI.metrics), function(i) 
  plot.density.by.qMRI(qMRI.factor = qMRI.metrics[i],
                              all.R2.melt = all.R2.vals.no.SWMSpD1.melt,
                              all.R2 = all.R2.vals.no.SWMSpD1,
                              ST = POI.names))

png(filename = paste(results.dir, "/Plot_expanded_pathway_analysis_by_qMRI_no_SWM_SpD1.png", sep = ""),
    width = 800, height = 800, units = "px", res = 120)
multi.plot.by.qMRI.no.SWMSpD1 <- ggarrange(plotlist = plot.list.by.qMRI.no.SWMSpD1, 
                                       nrow = 2, 
                                       ncol = 2, 
                                       common.legend = TRUE, 
                                       legend = "bottom")
plot(multi.plot.by.qMRI.no.SWMSpD1)
annotate_figure(multi.plot.by.qMRI.no.SWMSpD1, 
                top = text_grob(expression("Density of qMRI vs ST R"^2*" values")),
                left = text_grob(expression("R"^2), rot = 90))
dev.off()

# create 2x2 grid of plots for each GO category
plot.list.by.GO.no.SWMSpD1 <- list()

plot.list.by.GO.no.SWMSpD1 <- lapply(1:length(GO.categories), function(i) 
  plot.density.by.GO(GO.factor = GO.categories[i],
                            all.R2.melt = all.R2.vals.no.SWMSpD1.melt,
                            all.R2 = all.R2.vals.no.SWMSpD1,
                            qMRI = qMRI.metrics, 
                            POI = POI.names[i]))

png(filename = paste(results.dir, "/Plot_expanded_pathway_analysis_by_GO_no_SWM_SpD1.png", sep = ""),
    width = 800, height = 800, units = "px", res = 120)
multi.plot.by.GO.no.SWMSpD1 <- ggarrange(plotlist = plot.list.by.GO.no.SWMSpD1, 
                                     nrow = 2, 
                                     ncol = 2, 
                                     common.legend = TRUE, 
                                     legend = "bottom")
plot(multi.plot.by.GO.no.SWMSpD1)
annotate_figure(multi.plot.by.GO.no.SWMSpD1, 
                top = text_grob(expression("Density of qMRI vs ST R"^2*" values")),
                left = text_grob(expression("R"^2), rot = 90))
dev.off()