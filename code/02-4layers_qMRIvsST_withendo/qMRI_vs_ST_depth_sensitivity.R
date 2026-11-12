# Master script to look at qMRI vs ST across layers

# load in some libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(hydroGOF)
library(matrixStats)

# define directories
d0 <- getwd();
qMRI.data.dir <- paste(d0, "/../../raw_data/4layers_qMRI/", sep = "")
cell.data.dir <- paste(d0, "/../../processed_data/01-4layers_withendo/02-top3_cells_counts/", sep = "")
GO.data.dir <- paste(d0, "/../../processed_data/01-4layers_withendo/05-pathway_scoring/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir <- paste(d0, "/../../plots/02-4layers_qMRIvsST_withendo/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/02-4layers_qMRIvsST_withendo/", sep = "")
results.dir.qMRIonly <- paste(d0, "/../../plots/00-4layers_qMRI/", sep = "")
results.dir.cellonly <- paste(d0, "/../../plots/01-4layers_withendo/02-top3_cells_counts/", sep = "")
results.dir.GO.only <- paste(d0, "/../../plots/01-4layers_withendo/05-pathway_scoring/", sep = "")

# separate output folders for the additional depth-focused review analyses
results.dir.depth <- paste(d0, "/../../plots/02c-4layers_qMRIvsST_withendo_depth_sensitivity/", sep = "")
results.dir.stats.depth <- paste(d0, "/../../stats/02c-4layers_qMRIvsST_withendo_depth_sensitivity/", sep = "")
dir.create(results.dir.depth, recursive = TRUE, showWarnings = FALSE)
dir.create(results.dir.stats.depth, recursive = TRUE, showWarnings = FALSE)

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

qMRI.data <- qMRI.data[,!grepl("_SWM",colnames(qMRI.data))]

# load the additional continuous depth measure for the review analyses
depth.file <- paste(qMRI.data.dir, "anon_DLPFC_AMP_mean_depth_4_equivol_layers.csv", sep = "")
depth.data <- read.csv(depth.file)
depth.data <- depth.data[,!grepl("_SWM",colnames(depth.data))]
stopifnot(identical(qMRI.data$subject_region, depth.data$subject_region))

# reorganise depth in exactly the same anterior/middle/posterior order as qMRI
depth.data.ant <- depth.data[grepl("ant", depth.data$subject_region),2:ncol(depth.data)]
colnames(depth.data.ant) <- paste(colnames(depth.data.ant),"_ant",sep = "")
depth.data.mid <- depth.data[grepl("mid", depth.data$subject_region),2:ncol(depth.data)]
colnames(depth.data.mid) <- paste(colnames(depth.data.mid),"_mid",sep = "")
depth.data.pos <- depth.data[grepl("pos", depth.data$subject_region),2:ncol(depth.data)]
colnames(depth.data.pos) <- paste(colnames(depth.data.pos),"_pos",sep = "")
depth.data.dems <- cbind(subject = qMRI.dems$subject, depth.data.ant, depth.data.mid, depth.data.pos)

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
qMRI.metrics <- c("MTsat","PrD", "QSM", "R1","R2s")
MOI <- qMRI.metrics

# match up QC and qMRI data frames
QC.info <- subset(QC.info, QC.info$subject %in% qMRI.data.dems$subject)
rownames(QC.info) <- 1:nrow(QC.info)

# remove data from demographics and qMRI data using QC info
# if they failed MPRAGE QC or FS QC, then remove from everywhere
keep.structural.QC <- QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1
qMRI.data.dems <- qMRI.data.dems[keep.structural.QC,]
depth.data.dems <- depth.data.dems[keep.structural.QC,]
QC.info <- QC.info[keep.structural.QC,]

# rename columns
colnames(qMRI.data.dems) <- lapply(colnames(qMRI.data.dems), rename_layers)
colnames(depth.data.dems) <- lapply(colnames(depth.data.dems), rename_layers)
# get all ROI names
qMRI.ROI.names <- colnames(qMRI.data.dems)[7:ncol(qMRI.data.dems)]

# need to create some new columns in the QC.info dataframe
QC.info$QSM_QC <- QC.info$QSM_QC == 1 & QC.info$has_QSM == 1
QC.info$R1_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$R2s_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$PrD_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$MTsat_QC <- QC.info$MPM_QC == 1 & QC.info$has_MT == 1

#### adjust qMRI values ####

# initialize a data frame for the adjusted values
adj.qMRI.data.dems <- qMRI.data.dems

# get mean and diff from mean age and sex for all subjects
mn.qMRI.age <- mean(qMRI.data.dems$age)
diff.qMRI.age <- mn.qMRI.age - qMRI.data.dems$age
mn.qMRI.sex <- mean(as.numeric(qMRI.data.dems$sex))
diff.qMRI.sex <- mn.qMRI.sex - as.numeric(qMRI.data.dems$sex)
mn.qMRI.CT <- mean(as.numeric(qMRI.data.dems$DLPFC_CT))
diff.qMRI.CT <- mn.qMRI.CT - as.numeric(qMRI.data.dems$DLPFC_CT)

# adjust all ROIs for age and sex using the covariance method
for (ROI in qMRI.ROI.names) {
  mdl.qMRI.adj <- lm(qMRI.data.dems[,ROI] ~ age + sex + DLPFC_CT, data = qMRI.data.dems)
  adj.qMRI.vals <- qMRI.data.dems[,ROI] +
    diff.qMRI.age * mdl.qMRI.adj$coefficients[2] +
    diff.qMRI.sex * mdl.qMRI.adj$coefficients[3] +
    diff.qMRI.CT  * mdl.qMRI.adj$coefficients[4]
  adj.qMRI.data.dems[,ROI] <- adj.qMRI.vals
}

#### load in & reorganise cell data ####

# load in ST patient dems
ST.dems <- read.csv(paste(dems.dir,"ST_patient_info.csv",sep = ""))
ST.dems <- ST.dems[order(ST.dems$sample_id),]
ST.dems$sex <- factor(ST.dems$sex)

# specify the cell types we're interested in
COI <- c("Astro","Micro", "OPC", "Oligo")

# get names of 4 layer cell count files
all.cell.fnames <- dir(cell.data.dir,"count_(.+)_4SpD_withEndo")

# select only the interesting ones (preserve COI order)
COI.fnames <- unlist(lapply(COI, function(coi) {
  grep(coi, all.cell.fnames, value = TRUE)
}))

# put these all into one data frame with the dems
cell.data.dems <- ST.dems

# go through the fnames, load in and merge together
for (i in 1:length(COI.fnames)){
  COI.fname <- COI.fnames[i]
  cell.data <- read.csv(paste(cell.data.dir, COI.fname, sep = ""))
  cell.data <- cell.data[,2:ncol(cell.data)]
  colnames(cell.data) <- paste(COI[i],"_",colnames(cell.data),sep="")
  cell.data <- cell.data[,c(paste(COI[i], "_SpD4",sep = ""),
                            paste(COI[i], "_SpD3",sep = ""),
                            paste(COI[i], "_SpD2",sep = ""),
                            paste(COI[i], "_SpD1",sep = ""))]
  cell.data.dems <- cbind(cell.data.dems,cell.data)
}

# reorganise in the same way as the qMRI data
cell.data.ant <- cell.data.dems[ST.dems$region == "anterior", 7:ncol(cell.data.dems)]
colnames(cell.data.ant) <- paste(colnames(cell.data.ant),"_ant", sep = "")
cell.data.mid <- cell.data.dems[ST.dems$region == "middle", 7:ncol(cell.data.dems)]
colnames(cell.data.mid) <- paste(colnames(cell.data.mid),"_mid", sep = "")
cell.data.pos <- cell.data.dems[ST.dems$region == "posterior", 7:ncol(cell.data.dems)]
colnames(cell.data.pos) <- paste(colnames(cell.data.pos),"_pos", sep = "")

# merge together with the demographics
cell.data.dems <- cbind(ST.dems[seq(1,nrow(ST.dems),3),], cell.data.ant, cell.data.mid, cell.data.pos)
rownames(cell.data.dems) <- 1:nrow(cell.data.dems)
cell.data.dems <- subset(cell.data.dems, select = -c(sample_id, region))

#### load in & reorganise pathway data ####

POI.noSWM <- as.character(read.delim("pathways_of_interest_no_SWM.txt", sep = ",", header = FALSE, quote = "'"))
POI.names <- POI.noSWM
POI <- gsub(".*_","",POI.noSWM)

# get names of 4 layer pathway files
all.GO.fnames <- dir(GO.data.dir,"scores_4layers_genes_(.+)")

# select only the a priori ones
POI.fnames <- all.GO.fnames[grep(paste(POI, collapse = "|"), all.GO.fnames)]

# put these all into one data frame with the dems
GO.data.dems <- ST.dems

# do some reorganising and put them into wide format for adjustment
for (i in 1:length(POI.fnames)){
  POI.fname <- POI.fnames[i]
  GO.data <- read.csv(paste(GO.data.dir, POI.fname, sep = ""))
  GO.data <- GO.data[,2:ncol(GO.data)]
  colnames(GO.data) <- c("sample_id", "region", POI.names[i])
  GO.data <- GO.data[order(
    GO.data$sample_id,
    factor(GO.data$region, levels = c("SpD4","SpD3","SpD2","SpD1"))),]
  GO.data <- pivot_wider(GO.data,names_from = region, values_from = POI.names[i],
                         names_prefix = paste(POI.names[i],"_", sep = ""))
  GO.data.dems <- cbind(GO.data.dems,GO.data[,2:ncol(GO.data)])
}

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

GO.data.dems <- GO.data.dems[,!grepl("_WM_",colnames(GO.data.dems))]

#### adjust ST data ####

# initialize data frames for the adjusted values
adj.cell.data.dems <- cell.data.dems
adj.GO.data.dems <- GO.data.dems

# get mean and diff from mean age and sex for all subjects
mn.ST.age <- mean(cell.data.dems$age)
diff.ST.age <- mn.ST.age - cell.data.dems$age
mn.ST.sex <- mean(as.numeric(cell.data.dems$sex))
diff.ST.sex <- mn.ST.sex - as.numeric(cell.data.dems$sex)

# adjust cell ST data
cell.ROI.names <- colnames(cell.data.dems)[5:ncol(cell.data.dems)]
for (ROI in cell.ROI.names) {
  mdl.cell.adj <- lm(cell.data.dems[,ROI] ~ age + sex, data = cell.data.dems)
  adj.cell.vals <- cell.data.dems[,ROI] +
    diff.ST.age * mdl.cell.adj$coefficients[2] +
    diff.ST.sex * mdl.cell.adj$coefficients[3]
  adj.cell.data.dems[,ROI] <- adj.cell.vals
}

# adjust pathway ST data
GO.ROI.names <- colnames(GO.data.dems)[5:ncol(GO.data.dems)]
for (ROI in GO.ROI.names) {
  mdl.GO.adj <- lm(GO.data.dems[,ROI] ~ age + sex, data = GO.data.dems)
  adj.GO.vals <- GO.data.dems[,ROI] +
    diff.ST.age * mdl.GO.adj$coefficients[2] +
    diff.ST.sex * mdl.GO.adj$coefficients[3]
  adj.GO.data.dems[,ROI] <- adj.GO.vals
}

#### qMRI boxplots ####

# subjects
subj <- adj.qMRI.data.dems$subject

# define spatial domains / subregions
dom <- c("SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

# include in boxplots?
adj.qMRI.data.dems$inc_QSM <- QC.info$QSM_QC == 1 & QC.info$has_QSM == 1
adj.qMRI.data.dems$inc_R1 <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
adj.qMRI.data.dems$inc_R2s <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
adj.qMRI.data.dems$inc_PrD <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
adj.qMRI.data.dems$inc_MTsat <- QC.info$MPM_QC == 1 & QC.info$has_MT == 1

# do this for qMRI
adj.qMRI.data.long <- expand.grid(subject = subj, SpD = dom, region = reg)
adj.qMRI.data.dems.melt <- melt(adj.qMRI.data.dems,id.vars=c("subject","age","sex","diagnosis",
                                                             "inc_QSM","inc_R1","inc_R2s","inc_PrD","inc_MTsat"))

for (m in qMRI.metrics){
  metric.data <- adj.qMRI.data.dems.melt[grepl(m, adj.qMRI.data.dems.melt$variable),11]
  QC.data <- adj.qMRI.data.dems.melt[grepl(m, adj.qMRI.data.dems.melt$variable),paste("inc_",m,sep = "")]
  adj.qMRI.data.long[,paste(m, sep = "")] <- metric.data
  adj.qMRI.data.long[,paste("inc_", m, sep = "")] <- QC.data
}

adj.qMRI.data.long$diagnosis <- adj.qMRI.data.dems.melt$diagnosis[1:nrow(adj.qMRI.data.long)]

# plot for controls & PD
adj.qMRI.data.long.Ctrl <- subset(adj.qMRI.data.long, diagnosis == "Control")
adj.qMRI.data.long.PD <- subset(adj.qMRI.data.long, diagnosis == "PD")


qMRI.boxplot <- function(metric, dframe) {
  ggplot(
    dframe[dframe[,paste("inc_",metric,sep = "")],],
    aes(SpD, .data[[metric]], fill = SpD, shape = region)
  ) +
    geom_point(position = position_jitterdodge(), aes(colour = SpD), alpha = 0.9) +
    geom_boxplot(alpha = 0) +
    theme_minimal() +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    ggtitle(metric,) +
    scale_fill_manual(values = SpD.cols) +
    scale_color_manual(values = SpD.cols)
}
# run without SWM

adj.qMRI.data.long.no.SWM <- subset(adj.qMRI.data.long, SpD != "SWM")

# plot for controls & PD
adj.qMRI.data.long.Ctrl.no.SWM <- subset(adj.qMRI.data.long.no.SWM, diagnosis == "Control")
adj.qMRI.data.long.PD.no.SWM <- subset(adj.qMRI.data.long.no.SWM, diagnosis == "PD")

dom <- c("SpD4", "SpD3","SpD2","SpD1")
# colours for spatial domains
SpD.cols <- c("#66CCEE","#228833","#CCBB44","#EE6677")

# qMRI boxplots for Ctrl
qMRI.boxplotlist.Ctrl.no.SWM <- lapply(1:length(MOI), function(i) 
  qMRI.boxplot(metric = MOI[i], adj.qMRI.data.long.Ctrl.no.SWM))

png(filename = paste(results.dir.qMRIonly, "/Control_qMRI_boxplots_no_SWM.png", sep = ""),
    width = 1400, height = 600, units = "px", res = 120)
qMRI.multiboxplot.Ctrl.no.SWM <- ggarrange(plotlist = qMRI.boxplotlist.Ctrl.no.SWM, 
                                    nrow = 1, ncol = 5, 
                                    common.legend = TRUE, 
                                    legend = "right")
plot(qMRI.multiboxplot.Ctrl.no.SWM)
dev.off()

# qMRI boxplots for PD
qMRI.boxplotlist.PD.no.SWM <- lapply(1:length(MOI), function(i) 
  qMRI.boxplot(metric = MOI[i], adj.qMRI.data.long.PD.no.SWM))

png(filename = paste(results.dir.qMRIonly, "/PD_qMRI_boxplots_no_SWM.png", sep = ""),
    width = 1400, height = 600, units = "px", res = 120)
qMRI.multiboxplot.PD.no.SWM <- ggarrange(plotlist = qMRI.boxplotlist.PD.no.SWM, 
                                  nrow = 1, ncol = 5, 
                                  common.legend = TRUE, 
                                  legend = "right")
plot(qMRI.multiboxplot.PD.no.SWM)
dev.off()


# remove info about boxplot inclusion
adj.qMRI.data.dems<- subset(adj.qMRI.data.dems, 
                            select = -c(inc_QSM, 
                                        inc_R1,
                                        inc_R2s,
                                        inc_PrD,
                                        inc_MTsat))

#### cell type boxplots ####

# subjects
subj <- adj.cell.data.dems$subject

# define spatial domains / subregions
dom <- c("SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

# do this for cell types
adj.cell.data.long <- expand.grid(subject = subj, SpD = dom, region = reg)
adj.cell.data.dems.melt <- melt(adj.cell.data.dems,id.vars=c("subject","age","sex","diagnosis"))

for (c in COI){
  cell.data <- adj.cell.data.dems.melt[grepl(c, adj.cell.data.dems.melt$variable),6]
  adj.cell.data.long[,paste(c, sep = "")] <- cell.data
}

cell.boxplot <- function(ST, dframe) {
  ggplot(dframe,aes(SpD, .data[[ST]], fill = SpD, shape = region)) +
    geom_point(position = position_jitterdodge(), aes(colour = SpD), alpha = 0.9) +
    geom_boxplot(alpha = 0) +
    theme_minimal() +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    ggtitle(ST,) +
    scale_fill_manual(values = SpD.cols) +
    scale_color_manual(values = SpD.cols)
}

# run without SWM
dom <- c("SpD4", "SpD3","SpD2","SpD1")
adj.cell.data.long.no.SWM <- subset(adj.cell.data.long, SpD != "SWM")

# cell boxplots for Ctrl
cell.boxplotlist.no.SWM <- lapply(1:length(COI), function(i) 
  cell.boxplot(ST = COI[i], adj.cell.data.long.no.SWM))

png(filename = paste(results.dir.cellonly, "/Control_cell_boxplots_no_SWM.png", sep = ""),
    width = 1400, height = 600, units = "px", res = 120)
cell.multiboxplot.no.SWM <- ggarrange(plotlist = cell.boxplotlist.no.SWM, 
                                    nrow = 1, ncol = 4, 
                                    common.legend = TRUE, 
                                    legend = "right")
plot(cell.multiboxplot.no.SWM)
dev.off()

#### POI boxplots ####

# define spatial domains / subregions
dom <- c("SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

# long format
adj.GO.data.long <- expand.grid(subject = subj, SpD = dom, region = reg)
adj.GO.data.dems.melt <- melt(adj.GO.data.dems,id.vars=c("subject","age","sex","diagnosis"))

for (p in POI.names){
  GO.data <- adj.GO.data.dems.melt[grepl(p, adj.GO.data.dems.melt$variable),6]
  adj.GO.data.long[,paste(p, sep = "")] <- GO.data
}

# remove SWM
adj.GO.data.long.no.SWM <- subset(adj.GO.data.long, SpD != "SWM")

GO.boxplot <- function(ST, dframe) {
  ggplot(dframe,aes(SpD, .data[[ST]], fill = SpD, shape = region)) +
    geom_point(position = position_jitterdodge(), aes(colour = SpD), alpha = 0.9) +
    geom_boxplot(alpha = 0) +
    theme_minimal() +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    ggtitle(ST,) +
    scale_fill_manual(values = SpD.cols) +
    scale_color_manual(values = SpD.cols)
}

# GO boxplots for Ctrl
GO.boxplotlist <- lapply(1:length(POI.names), function(i) 
  GO.boxplot(ST = POI.names[i], adj.GO.data.long.no.SWM))

png(filename = paste(results.dir.GO.only, "/Control_POI_boxplots.png", sep = ""),
    width = 1600, height = 1200, units = "px", res = 120)
GO.multiboxplot <- ggarrange(plotlist = GO.boxplotlist, 
                             nrow = 3, ncol = 5, 
                             common.legend = TRUE, 
                             legend = "right")
plot(GO.multiboxplot)
dev.off()

#### qMRI vs ST correlation analysis set-up ####

# define spatial domains / subregions
dom <- c("SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

# set up table containing qMRI / ST median values

adj.qMRI.ST.medians <- expand.grid(SpD = dom,region = reg)
#adj.qMRI.GO.medians <- adj.qMRI.cell.medians

# go through each qMRI metric and remove QC-failed subjects / nan values before averaging

# do this for controls and PD separately
# (although only looking at Ctrl vs Ctrl for now)
adj.qMRI.data.dems.Ctrl <- subset(adj.qMRI.data.dems, diagnosis == "Control")
QC.info.Ctrl <- subset(QC.info, adj.qMRI.data.dems$diagnosis == "Control" )
adj.qMRI.data.dems.PD <- subset(adj.qMRI.data.dems, diagnosis == "PD")
QC.info.PD <- subset(QC.info, adj.qMRI.data.dems$diagnosis == "PD" )

for (m in MOI){
  # remove QC failed subs before calculating mean
  metric.data.Ctrl <- adj.qMRI.data.dems.Ctrl[QC.info.Ctrl[,paste(m,"_QC",sep = "")], 
                                              grepl(m, colnames(adj.qMRI.data.dems))]
  metric.medians.Ctrl <- colMedians(as.matrix(metric.data.Ctrl),na.rm = TRUE)
  metric.data.PD <- adj.qMRI.data.dems.PD[QC.info.PD[,paste(m,"_QC",sep = "")], 
                                          grepl(m, colnames(adj.qMRI.data.dems))]
  metric.medians.PD <- colMedians(as.matrix(metric.data.PD),na.rm = TRUE)
  adj.qMRI.ST.medians[,paste(m, "_Control", sep = "")] <- metric.medians.Ctrl
  adj.qMRI.ST.medians[,paste(m, "_PD", sep = "")] <- metric.medians.PD
}

# add in median cell values
for (c in COI){
  cell.data <- adj.cell.data.dems[,grepl(paste(c,"_",sep = ""), colnames(adj.cell.data.dems))]
  cell.medians <- colMedians(as.matrix(cell.data), na.rm = TRUE)
  adj.qMRI.ST.medians[,paste(c, "_Control", sep = "")] <- cell.medians
}

# add in median GO values
for (p in POI.names){
  GO.data <- adj.GO.data.dems[,grepl(paste(p,"_",sep = ""), colnames(adj.GO.data.dems))]
  GO.medians <- colMedians(as.matrix(GO.data), na.rm = TRUE)
  adj.qMRI.ST.medians[,paste(p, "_Control", sep = "")] <- GO.medians
}

# add the median continuous depth measure using the same control subjects and spatial order
depth.data.dems.Ctrl <- depth.data.dems[depth.data.dems$subject %in% adj.qMRI.data.dems.Ctrl$subject,]
depth.medians.Ctrl <- colMedians(as.matrix(depth.data.dems.Ctrl[,2:ncol(depth.data.dems.Ctrl)]),na.rm = TRUE)
stopifnot(length(depth.medians.Ctrl) == nrow(adj.qMRI.ST.medians))
adj.qMRI.ST.medians$depth <- depth.medians.Ctrl

write.csv(adj.qMRI.ST.medians, file = paste(
  results.dir.stats, "/adj_qMRI_ST_medians.csv", sep = ""), 
  row.names = FALSE)

# define a function to plot the data
plot.ST.qMRI <- function(dframe ,MRIST, group.string, sts) {
  MRI <- MRIST[1]
  ST <- MRIST[2]
  
  x.data <- dframe[[paste(MRI, "_", group.string, sep = "")]]
  y.data <- dframe[[paste(ST, "_Control", sep = "")]]
  
  R2 <- sts[sts[,1] == MRI[[1]] & sts[,2] == ST[[1]], ]$R2
  adj.p <- sts[sts[,1] == MRI[[1]] & sts[,2] == ST[[1]], ]$corr_p
  
  ggplot(data = dframe, aes(x = x.data, y = y.data)) +
    theme_minimal() +
    {
      if (ST == COI[length(COI)] || 
          ST == POI.names[length(POI.names)])
        xlab(MRI)
      else
        theme(axis.title.x = element_blank())
    } +
    {
      if (MRI == MOI[1])
        ylab(ST)
      else
        theme(axis.title.y = element_blank())
    } +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      linetype = "dotted",
      colour = "black",
      fill = "lightgrey"
    ) +
    geom_point(aes(colour = SpD, shape = region)) +
    scale_color_manual(values = SpD.cols) +
    scale_x_continuous(labels = function(x) {
      ifelse(
        x == 0, "0", # zero shows as 0
        ifelse(
          abs(x) < 1e-3, 
          format(x, scientific = TRUE, digits = 2), # Scientific format for small values
          signif(x, 2)                    # otherwise 2 sig figs
        )
      )
    }) +
    scale_y_continuous(labels = function(y) {
      ifelse(
        y == 0, 
        "0", 
        ifelse(
          abs(y) < 1e-3, 
          format(y, scientific = TRUE, digits = 2), 
          signif(y, 2)                              
        )
      )
    }) +
    annotate("text",
             x = -Inf,
             y = Inf,
             hjust = 0,
             vjust = 1,
             label = as.expression(bquote(
               R^2 == .(signif(R2, 2)) ~ "," ~ p[FDR] == .(if (adj.p < 1e-3) {
                 format(adj.p, digits = 2, scientific = TRUE)
               } else {
                 signif(adj.p, 2)
               })
             )),
             size = 3.25,
             parse = TRUE
    )
}

plot.ST.qMRI.cooksd <- function(dframe ,MRIST, group.string) {
  # sort by SpD
  dframe <- dframe[order(dframe$SpD),]
  
  MRI <- MRIST[1]
  ST <- MRIST[2]
  
  fit <- lm(as.formula(paste0(paste(MRI,group.string, sep = "_"),"~",
                              paste(ST,"Control", sep = "_"))),
            data = dframe)
  res <- summary(fit)
  dframe[,paste(MRI,ST,"cooksd", sep = "_")] <- cooks.distance(fit)
  
  # define cutoff as N-k-1
  cutoff <- 4/(nrow(dframe) - 1 - 1)
  
  x.data <- 1:nrow(dframe)
  y.data <- dframe[[paste(MRI,ST,"cooksd", sep = "_")]]
  
  ggplot(data = dframe, aes(x = x.data, y = y.data)) +
    theme_minimal() +
    {
      if (ST == COI[length(COI)] || 
          ST == POI.names[length(POI.names)])
        xlab(MRI)
      else
        theme(axis.title.x = element_blank())
    } +
    {
      if (MRI == MOI[1])
        ylab(ST)
      else
        theme(axis.title.y = element_blank())
    } +
    theme(axis.text.x = element_blank()
    ) +
    geom_point(aes(colour = SpD, shape = region)) +
    geom_hline(yintercept = cutoff, linetype = "dotted", size = 1) +
    scale_color_manual(values = SpD.cols)
}

# and define a function to generate some stats tables
stats.ST.qMRI <- function(dframe, MRIST, group.string) {
  MRI <- MRIST[1]
  ST <- MRIST[2]
  
  x.data <- dframe[[paste(MRI, "_", group.string, sep = "")]]
  y.data <- dframe[[paste(ST, "_Control", sep = "")]]
  
  res <- cor.test(x.data,y.data)
  
  Tstat <- res$statistic
  rho <- res$estimate
  R2 <- res$estimate^2
  p <- res$p.value
  
  return(c(Tstat,rho,R2,p))
  
}

# expand out the qMRI vs ST combos
MOICOI <- expand.grid(MOI,COI,stringsAsFactors = FALSE)
MOIPOI <- expand.grid(MOI,POI.names,stringsAsFactors = FALSE)

#### qMRI vs ST LOO set-up ####

# define a function to run the LOO
run.qMRI.ST.LOO <- function(dframe, ST, LOO, stats, group.string) {
  
  # loop though ST metrics
  for (ii in 1:length(ST)) {
    
    st <- ST[ii]
    
    # loop through each rregion-depth row
    for (jj in 1:nrow(dframe)) {
      
      # define train and test data
      train <- dframe[-jj, ]
      test <- dframe[jj, ]
      
      # actual data
      loo.actual <- test[, paste(st, "Control", sep = "_")]
      
      loo.fm <- as.formula(paste0(paste(st, "Control", sep = "_"), "~",
                                  paste("QSM", group.string, sep = "_"), "+",
                                  paste("R1", group.string, sep = "_"), "+",
                                  paste("R2s", group.string, sep = "_"), "+",
                                  paste("PrD", group.string, sep = "_"), "+",
                                  paste("MTsat", group.string, sep = "_")))
      # run LOO
      loo.fit <- lm(loo.fm, data = train)
      loo.pred <- predict(loo.fit, newdata = test)
      loo.SE <- (loo.pred - loo.actual) ^ 2
      
      # store results
      LOO[jj, paste(st, "actual", sep = "_")] <- loo.actual
      LOO[jj, paste(st, "predicted", sep = "_")] <- loo.pred
      LOO[jj, paste(st, "SE", sep = "_")] <- loo.SE
      
    }
    
    fit <- lm(as.formula(paste0(paste(st,"actual", sep = "_"),"~",
                                paste(st,"predicted", sep = "_"))),
              data = LOO)
    res <- summary(fit)
    LOO[,paste(st, "cooksd", sep = "_")]<-cooks.distance(fit)
    
    # Save stats
    stats[ii,1] <- st
    stats[ii,2] <- nrmse(LOO[, paste(st, "predicted", sep = "_")],
                         LOO[, paste(st, "actual", sep = "_")])
    stats[ii,3] <- res$coefficients[2,3]
    stats[ii,4] <- cor.test(LOO[, paste(st, "predicted", sep = "_")],
                            LOO[, paste(st, "actual", sep = "_")])$estimate
    stats[ii,5] <- res$r.squared
    stats[ii,6] <- res$coefficients[2, 4]
    
  }
  
  # return updated NRMSE and LOO
  return(list(LOO = LOO, stats = stats))
  
}

run.qMRI.ST.LOO.SpD <- function(dframe, ST, LOO, stats, group.string) {
  
  # loop though ST metrics
  for (ii in 1:length(ST)) {
    
    st <- ST[ii]
    
    # loop through each SpD and leave-one-SpD out
    for (spd in unique(dframe$SpD)) {
      
      # get SpD rows
      test.idx <- which(dframe$SpD == spd)
      
      # define train and test data
      train <- dframe[-test.idx, ]
      test <- dframe[test.idx, ]
      
      # actual data
      loo.actual <- test[, paste(st, "Control", sep = "_")]
      
      loo.fm <- as.formula(paste0(paste(st, "Control", sep = "_"), "~",
                                  paste("QSM", group.string, sep = "_"), "+",
                                  paste("R1", group.string, sep = "_"), "+",
                                  paste("R2s", group.string, sep = "_"), "+",
                                  paste("PrD", group.string, sep = "_"), "+",
                                  paste("MTsat", group.string, sep = "_")))
      # run LOO
      loo.fit <- lm(loo.fm, data = train)
      loo.pred <- predict(loo.fit, newdata = test)
      loo.SE <- (loo.pred - loo.actual) ^ 2
      
      # store results
      LOO[test.idx, paste(st, "actual", sep = "_")] <- loo.actual
      LOO[test.idx, paste(st, "predicted", sep = "_")] <- loo.pred
      LOO[test.idx, paste(st, "SE", sep = "_")] <- loo.SE
      
    }
    
    fit <- lm(as.formula(paste0(paste(st,"actual", sep = "_"),"~",
                                paste(st,"predicted", sep = "_"))),
              data = LOO)
    res <- summary(fit)
    LOO[,paste(st, "cooksd", sep = "_")]<-cooks.distance(fit)
    
    full.coef <- coef(fit)
    full.mse <- sum(residuals(fit)^2) / df.residual(fit)
    p <- length(full.coef)
    
    # SpD-level cooks distance
    
    full.coef <- coef(fit)
    full.mse <- sum(residuals(fit)^2) / df.residual(fit)
    
    for (spd in unique(LOO$SpD)) {
      
      fit.minus <- lm(
        as.formula(
          paste0(paste(st, "actual", sep = "_"),
                 " ~ ", paste(st, "predicted", sep = "_"))
        ),
        data = subset(LOO, SpD != spd)
      )
      
      beta.diff <- coef(fit.minus) - full.coef
      
      spd.cook <- sum(beta.diff^2) / (length(full.coef) * full.mse)
      
      # assign same value to all rows from this SpD
      LOO[LOO$SpD == spd,paste(st, "SpD_cooksd", sep = "_")] <- spd.cook
    }
    
    # Save stats
    stats[ii,1] <- st
    stats[ii,2] <- nrmse(LOO[, paste(st, "predicted", sep = "_")],
                         LOO[, paste(st, "actual", sep = "_")])
    stats[ii,3] <- res$coefficients[2,3]
    stats[ii,4] <- cor.test(LOO[, paste(st, "predicted", sep = "_")],
                            LOO[, paste(st, "actual", sep = "_")])$estimate
    stats[ii,5] <- res$r.squared
    stats[ii,6] <- res$coefficients[2, 4]
    
  }
  
  # return updated NRMSE and LOO
  return(list(LOO = LOO, stats = stats))
  
}

# define a function to plot the data
plot.LOO <- function(ST, LOO,sts) {
  x.data <- LOO[, paste(ST, "actual", sep = "_")]
  y.data <- LOO[[paste(ST, "predicted", sep = "_")]]
  
  R2 <- sts[sts[,1] == ST[[1]], ]$R2
  adj.p <- sts[sts[,1] == ST[[1]], ]$corr_p
  
  ggplot(data = LOO, aes(x = x.data, y = y.data)) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    ggtitle(ST) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      linetype = "dotted",
      colour = "black",
      fill = "lightgrey"
    ) +
    geom_point(aes(colour = SpD, shape = region)) +
    scale_color_manual(values = SpD.cols) +
    scale_x_continuous(labels = function(x) {
      ifelse(
        x == 0, "0", # zero shows as 0
        ifelse(
          abs(x) < 1e-3, 
          format(x, scientific = TRUE, digits = 2), # Scientific format for small values
          signif(x, 2)                    # otherwise 2 sig figs
        )
      )
    }) +
    scale_y_continuous(labels = function(y) {
      ifelse(
        y == 0, 
        "0", 
        ifelse(
          abs(y) < 1e-3, 
          format(y, scientific = TRUE, digits = 2), 
          signif(y, 2)                              
        )
      )
    }) +
    annotate("text",
             x = -Inf,
             y = Inf,
             hjust = 0,
             vjust = 1,
             label = as.expression(bquote(
               R^2 == .(signif(R2, 2)) ~ "," ~ p[FDR] == .(if (adj.p < 1e-3) {
                 format(adj.p, digits = 2, scientific = TRUE)
               } else {
                 signif(adj.p, 2)
               })
             )),
             size = 3.25,
             parse = TRUE
    )
}

# function to plot the LOO Cook's distance
plot.LOO.cooksd <- function(ST, LOO) {
  # sort by SpD
  LOO <- LOO[order(LOO$SpD),]
  x.data <- 1:nrow(LOO)
  y.data <- LOO[[paste(ST, "cooksd", sep = "_")]]
  
  # define cutoff as N-k-1
  cutoff <- 4/(nrow(LOO) - 1 - 1)
  
  ggplot(data = LOO, aes(x = x.data, y = y.data)) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    ggtitle(ST) +
    geom_point(aes(colour = SpD, shape = region)) +
    geom_hline(yintercept = cutoff, linetype = "dotted", size = 1) +
    scale_color_manual(values = SpD.cols)
}

# function to plot the LOO Cook's distance
plot.LOO.SpD.cooksd <- function(ST, LOO) {
  # sort by SpD
  LOO <- LOO[order(LOO$SpD),]
  x.data <- 1:nrow(LOO)
  y.data <- LOO[[paste(ST, "SpD_cooksd", sep = "_")]]
  
  ggplot(data = LOO, aes(x = x.data, y = y.data)) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    ggtitle(ST) +
    geom_point(aes(colour = SpD, shape = region)) +
    scale_color_manual(values = SpD.cols)
}


#### qMRI vs ST correlation no SWM ####

# remove the SWM ROIs
adj.qMRI.ST.medians.no.SWM <- subset(adj.qMRI.ST.medians, SpD != "SWM")

# CELL COUNTS

# generate a stats table
cell.stats.Ctrl.no.SWM <- lapply(1:nrow(MOICOI), function(i) stats.ST.qMRI(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOICOI[i,], group.string = "Control"))
cell.stats.Ctrl.no.SWM <- as.data.frame(do.call(rbind, cell.stats.Ctrl.no.SWM))
cell.stats.Ctrl.no.SWM <- cbind(MOICOI,cell.stats.Ctrl.no.SWM)
colnames(cell.stats.Ctrl.no.SWM) <- c("qMRI", "cell_type", "T", "rho", "R2","p")
cell.stats.Ctrl.no.SWM$corr_p <- p.adjust(cell.stats.Ctrl.no.SWM$p, method = "fdr")
write.csv(cell.stats.Ctrl.no.SWM, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_cell_count_no_SWM.csv", sep = ""), 
  row.names = FALSE)

# lapply the plot function to create a multi plot
cell.plot.list.Ctrl.no.SWM <- lapply(1:nrow(MOICOI), function(i) plot.ST.qMRI(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOICOI[i, ], group.string = "Control", sts = cell.stats.Ctrl.no.SWM))

# Ctrl
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_cell_count_no_SWM.png", sep = ""),
    width = 1200, height = 800, units = "px", res = 120)
cell.multi.plot.Ctrl.no.SWM <- ggarrange(
  plotlist = cell.plot.list.Ctrl.no.SWM,
  nrow = 4,
  ncol = 5,
  common.legend = TRUE,
  legend = "bottom"
)
plot(cell.multi.plot.Ctrl.no.SWM)
annotate_figure(cell.multi.plot.Ctrl.no.SWM, top = text_grob("Control qMRI versus control cell count"))
dev.off()

# lapply for cook's distances
cell.plot.list.cooksd.Ctrl.no.SWM <- lapply(1:nrow(MOICOI), function(i) plot.ST.qMRI.cooksd(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOICOI[i, ], group.string = "Control"))

# Ctrl
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_cell_count_cooksd_no_SWM.png", sep = ""),
    width = 1200, height = 800, units = "px", res = 120)
cell.multi.plot.cooksd.Ctrl.no.SWM <- ggarrange(
  plotlist = cell.plot.list.cooksd.Ctrl.no.SWM,
  nrow = 4,
  ncol = 5,
  common.legend = TRUE,
  legend = "right"
)
plot(cell.multi.plot.cooksd.Ctrl.no.SWM)
annotate_figure(cell.multi.plot.cooksd.Ctrl.no.SWM, 
                top = text_grob("Effect of regions on Control qMRI vs. Control cell count correlation"),
                left = text_grob("Cook's distance", rot = 90),
                bottom = text_grob("ROI"))
dev.off()

# PATHWAYS OF INTEREST

# generate a stats table
GO.stats.Ctrl.no.SWM <- lapply(1:nrow(MOIPOI), function(i) stats.ST.qMRI(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOIPOI[i,], group.string = "Control"))
GO.stats.Ctrl.no.SWM <- as.data.frame(do.call(rbind, GO.stats.Ctrl.no.SWM))
GO.stats.Ctrl.no.SWM <- cbind(MOIPOI,GO.stats.Ctrl.no.SWM)
colnames(GO.stats.Ctrl.no.SWM) <- c("qMRI", "pathway_score", "T", "rho", "R2","p")
GO.stats.Ctrl.no.SWM$corr_p <- p.adjust(GO.stats.Ctrl.no.SWM$p, method = "fdr")
write.csv(GO.stats.Ctrl.no.SWM, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_GO_POI_no_SWM.csv", sep = ""), 
  row.names = FALSE)

# lapply the plot function to create a multi plot
GO.plot.list.Ctrl.no.SWM <- lapply(1:nrow(MOIPOI), function(i) plot.ST.qMRI(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOIPOI[i, ], group.string = "Control", sts = GO.stats.Ctrl.no.SWM))

# Ctrl
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_GO_POI_no_SWM.png", sep = ""),
    width = 1200, height = 2500, units = "px", res = 120)
GO.multi.plot.Ctrl.no.SWM <- ggarrange(
  plotlist = GO.plot.list.Ctrl.no.SWM,
  nrow = 14,
  ncol = 5,
  common.legend = TRUE,
  legend = "bottom"
)
plot(GO.multi.plot.Ctrl.no.SWM)
annotate_figure(GO.multi.plot.Ctrl.no.SWM, top = text_grob("Control qMRI versus control pathway score"))
dev.off()

# lapply for cook's distances
GO.plot.list.cooksd.Ctrl.no.SWM <- lapply(1:nrow(MOIPOI), function(i) plot.ST.qMRI.cooksd(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOIPOI[i, ], group.string = "Control"))

# Ctrl
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_GO_POI_cooksd_no_SWM.png", sep = ""),
    width = 1200, height = 2500, units = "px", res = 120)
GO.multi.plot.cooksd.Ctrl.no.SWM <- ggarrange(
  plotlist = GO.plot.list.cooksd.Ctrl.no.SWM,
  nrow = 14,
  ncol = 5,
  common.legend = TRUE,
  legend = "right"
)
plot(GO.multi.plot.cooksd.Ctrl.no.SWM)
annotate_figure(GO.multi.plot.cooksd.Ctrl.no.SWM, 
                top = text_grob("Effect of regions on Control qMRI vs. Control GO pathway correlation"),
                left = text_grob("Cook's distance", rot = 90),
                bottom = text_grob("ROI"))
dev.off()


#### qMRI vs ST LOO no SWM ####

# CELL COUNTS

# define empty data frames for LOO results
cell.LOO.Ctrl.colnames <- paste(rep(COI,each = 3),c("actual","predicted", "SE"),sep = "_")
cell.LOO.Ctrl.no.SWM <- data.frame(matrix(nrow = nrow(adj.qMRI.ST.medians.no.SWM), ncol = length(cell.LOO.Ctrl.colnames)))
colnames(cell.LOO.Ctrl.no.SWM) <- cell.LOO.Ctrl.colnames
cell.LOO.Ctrl.no.SWM <-cbind(adj.qMRI.ST.medians.no.SWM[,1:2],cell.LOO.Ctrl.no.SWM)

cell.LOO.stats.Ctrl.no.SWM <- data.frame(matrix(nrow = length(COI), ncol = 6))
colnames(cell.LOO.stats.Ctrl.no.SWM) <- c("cell_type", "NRMSE", "T","rho", "R2", "p")

# run the LOO
cell.LOO.Ctrl.output.no.SWM <- run.qMRI.ST.LOO(dframe = adj.qMRI.ST.medians.no.SWM,
                                               ST = COI, 
                                               LOO = cell.LOO.Ctrl.no.SWM, 
                                               stats = cell.LOO.stats.Ctrl.no.SWM, 
                                               group.string = "Control")
cell.LOO.Ctrl.no.SWM <- cell.LOO.Ctrl.output.no.SWM$LOO
cell.LOO.stats.Ctrl.no.SWM <- cell.LOO.Ctrl.output.no.SWM$stats

# save some stats
cell.LOO.stats.Ctrl.no.SWM$corr_p <- p.adjust(cell.LOO.stats.Ctrl.no.SWM$p, method = "fdr")
write.csv(cell.LOO.stats.Ctrl.no.SWM, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_cell_count_LOO_no_SWM.csv", sep = ""), 
  row.names = FALSE)

# save the LOO predictions
write.csv(cell.LOO.Ctrl.no.SWM, file = paste(
  results.dir.stats, "/LOO_predictions_control_qMRI_vs_cell_count_no_SWM.csv", sep = ""),
  row.names = FALSE)

# plot the LOO
plot.list.cell.LOO.Ctrl.no.SWM <- lapply(1:length(COI), function(i) 
  plot.LOO(ST = COI[i], LOO = cell.LOO.Ctrl.no.SWM, sts = cell.LOO.stats.Ctrl.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_cell_count_LOO_no_SWM.png", sep = ""),
    width = 600, height = 500, units = "px", res = 120)
multi.plot.cell.LOO.Ctrl.no.SWM <- ggarrange(
  plotlist = plot.list.cell.LOO.Ctrl.no.SWM,
  nrow = 2,
  ncol = 2,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.cell.LOO.Ctrl.no.SWM)
annotate_figure(multi.plot.cell.LOO.Ctrl.no.SWM, bottom = text_grob("Control SRT count"),
                left = text_grob("qMRI-predicted cell count", rot = 90),
                top = text_grob("Control qMRI medians predicting control cell count"))
dev.off()

# plot the LOO cook's distance
plot.list.cell.LOO.cooksd.Ctrl.no.SWM <- lapply(1:length(COI), function(i) plot.LOO.cooksd(ST = COI[i], cell.LOO.Ctrl.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_cell_count_LOO_cooksd_no_SWM.png", sep = ""),
    width = 600, height = 500, units = "px", res = 120)
multi.plot.cell.LOO.cooksd.Ctrl.no.SWM <- ggarrange(
  plotlist = plot.list.cell.LOO.cooksd.Ctrl.no.SWM,
  nrow = 2,
  ncol = 2,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.cell.LOO.cooksd.Ctrl.no.SWM)
annotate_figure(multi.plot.cell.LOO.cooksd.Ctrl.no.SWM,  
                left = text_grob("Cook's distance", rot = 90),
                top = text_grob("Influence of regions on cell LOO correlations"),
                bottom = text_grob("ROI"))
dev.off()

# now run a LOO leaving out whole SpDs instead
cell.LOO.Ctrl.SpD.colnames <- paste(rep(COI,each = 3),c("actual","predicted", "SE"),sep = "_")
cell.LOO.Ctrl.SpD.no.SWM <- data.frame(matrix(nrow = nrow(adj.qMRI.ST.medians.no.SWM), ncol = length(cell.LOO.Ctrl.SpD.colnames)))
colnames(cell.LOO.Ctrl.SpD.no.SWM) <- cell.LOO.Ctrl.SpD.colnames
cell.LOO.Ctrl.SpD.no.SWM <-cbind(adj.qMRI.ST.medians.no.SWM[,1:2],cell.LOO.Ctrl.SpD.no.SWM)

cell.LOO.stats.Ctrl.SpD.no.SWM <- data.frame(matrix(nrow = length(COI), ncol = 6))
colnames(cell.LOO.stats.Ctrl.SpD.no.SWM) <- c("cell_type", "NRMSE", "T","rho", "R2", "p")

# run the LOO
cell.LOO.Ctrl.SpD.output.no.SWM <- run.qMRI.ST.LOO.SpD(dframe = adj.qMRI.ST.medians.no.SWM,
                                               ST = COI, 
                                               LOO = cell.LOO.Ctrl.SpD.no.SWM, 
                                               stats = cell.LOO.stats.Ctrl.SpD.no.SWM, 
                                               group.string = "Control")
cell.LOO.Ctrl.SpD.no.SWM <- cell.LOO.Ctrl.SpD.output.no.SWM$LOO
cell.LOO.stats.Ctrl.SpD.no.SWM <- cell.LOO.Ctrl.SpD.output.no.SWM$stats

# save some stats
cell.LOO.stats.Ctrl.SpD.no.SWM$corr_p <- p.adjust(cell.LOO.stats.Ctrl.SpD.no.SWM$p, method = "fdr")
write.csv(cell.LOO.stats.Ctrl.SpD.no.SWM, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_cell_count_LOO_SpD_no_SWM.csv", sep = ""), 
  row.names = FALSE)

# save the LOO predictions
write.csv(cell.LOO.Ctrl.SpD.no.SWM, file = paste(
  results.dir.stats, "/LOO_SpD_predictions_control_qMRI_vs_cell_count_no_SWM.csv", sep = ""),
  row.names = FALSE)

# plot the LOO
plot.list.cell.LOO.Ctrl.SpD.no.SWM <- lapply(1:length(COI), function(i) 
  plot.LOO(ST = COI[i], LOO = cell.LOO.Ctrl.SpD.no.SWM, sts = cell.LOO.stats.Ctrl.SpD.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_cell_count_LOO_SpD_no_SWM.png", sep = ""),
    width = 600, height = 500, units = "px", res = 120)
multi.plot.cell.LOO.Ctrl.SpD.no.SWM <- ggarrange(
  plotlist = plot.list.cell.LOO.Ctrl.SpD.no.SWM,
  nrow = 2,
  ncol = 2,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.cell.LOO.Ctrl.SpD.no.SWM)
annotate_figure(multi.plot.cell.LOO.Ctrl.SpD.no.SWM, bottom = text_grob("Control SRT count"),
                left = text_grob("qMRI-predicted cell count", rot = 90),
                top = text_grob("Control qMRI medians predicting control cell count"))
dev.off()

# plot the LOO cook's distance
plot.list.cell.LOO.cooksd.Ctrl.SpD.no.SWM <- lapply(1:length(COI), function(i) plot.LOO.SpD.cooksd(ST = COI[i], cell.LOO.Ctrl.SpD.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_cell_count_LOO_SpD_cooksd_no_SWM.png", sep = ""),
    width = 600, height = 500, units = "px", res = 120)
multi.plot.cell.LOO.cooksd.Ctrl.SpD.no.SWM <- ggarrange(
  plotlist = plot.list.cell.LOO.cooksd.Ctrl.SpD.no.SWM,
  nrow = 2,
  ncol = 2,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.cell.LOO.cooksd.Ctrl.SpD.no.SWM)
annotate_figure(multi.plot.cell.LOO.cooksd.Ctrl.SpD.no.SWM,  
                left = text_grob("Cook's distance", rot = 90),
                top = text_grob("Influence of regions on cell LOO correlations"),
                bottom = text_grob("ROI"))
dev.off()

# PATHWAYS OF INTEREST

# define empty data frames for LOO results
GO.LOO.Ctrl.colnames <- paste(rep(POI.names,each = 3),c("actual","predicted", "SE"),sep = "_")
GO.LOO.Ctrl.no.SWM <- data.frame(matrix(nrow = nrow(adj.qMRI.ST.medians.no.SWM), ncol = length(GO.LOO.Ctrl.colnames)))
colnames(GO.LOO.Ctrl.no.SWM) <- GO.LOO.Ctrl.colnames
GO.LOO.Ctrl.no.SWM <-cbind(adj.qMRI.ST.medians.no.SWM[,1:2],GO.LOO.Ctrl.no.SWM)

GO.LOO.stats.Ctrl.no.SWM <- data.frame(matrix(nrow = length(POI.names), ncol = 6))
colnames(GO.LOO.stats.Ctrl.no.SWM) <- c("GO_type", "NRMSE", "T","rho", "R2", "p")

# run the LOO
GO.LOO.Ctrl.no.SWM.output <- run.qMRI.ST.LOO(dframe = adj.qMRI.ST.medians.no.SWM,
                                             ST = POI.names,
                                             LOO = GO.LOO.Ctrl.no.SWM,
                                             stats = GO.LOO.stats.Ctrl.no.SWM,
                                             group.string = "Control")
GO.LOO.Ctrl.no.SWM <- GO.LOO.Ctrl.no.SWM.output$LOO
GO.LOO.stats.Ctrl.no.SWM <- GO.LOO.Ctrl.no.SWM.output$stats

# save some stats
GO.LOO.stats.Ctrl.no.SWM$corr_p <- p.adjust(GO.LOO.stats.Ctrl.no.SWM$p, method = "fdr")
write.csv(GO.LOO.stats.Ctrl.no.SWM, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_GO_POI_LOO_no_SWM.csv", sep = ""),
  row.names = FALSE)

# save the LOO predictions
write.csv(GO.LOO.Ctrl.no.SWM, file = paste(
  results.dir.stats, "/LOO_predictions_control_qMRI_vs_GO_POI_no_SWM.csv", sep = ""),
  row.names = FALSE)

# plot the LOO
plot.list.GO.LOO.Ctrl.no.SWM <- lapply(1:length(POI.names), function(i) 
  plot.LOO(ST = POI.names[i], LOO = GO.LOO.Ctrl.no.SWM, sts = GO.LOO.stats.Ctrl.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_GO_POI_LOO_no_SWM.png", sep = ""),
    width = 1400, height = 700, units = "px", res = 120)
multi.plot.GO.LOO.Ctrl.no.SWM <- ggarrange(
  plotlist = plot.list.GO.LOO.Ctrl.no.SWM,
  nrow = 3,
  ncol = 5,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.GO.LOO.Ctrl.no.SWM)
annotate_figure(multi.plot.GO.LOO.Ctrl.no.SWM, bottom = text_grob("Control SRT GO pathway score"),
                left = text_grob("qMRI-predicted GO pathway score", rot = 90),
                top = text_grob("Control qMRI medians predicting control GO pathway score"))
dev.off()

# plot the LOO cook's distance
plot.list.GO.LOO.cooksd.Ctrl.no.SWM <- lapply(1:length(POI.names), function(i) plot.LOO.cooksd(ST = POI.names[i], GO.LOO.Ctrl.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_GO_POI_LOO_cooksd_no_SWM.png", sep = ""),
    width = 1200, height = 700, units = "px", res = 120)
multi.plot.GO.LOO.cooksd.Ctrl.no.SWM <- ggarrange(
  plotlist = plot.list.GO.LOO.cooksd.Ctrl.no.SWM,
  nrow = 3,
  ncol = 5,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.GO.LOO.cooksd.Ctrl.no.SWM)
annotate_figure(multi.plot.GO.LOO.cooksd.Ctrl.no.SWM,  
                left = text_grob("Cook's distance", rot = 90),
                top = text_grob("Influence of regions on GO LOO correlations"),
                bottom = text_grob("ROI"))
dev.off()


# now run a LOO leaving out whole SpDs instead
GO.LOO.Ctrl.SpD.colnames <- paste(rep(POI.names,each = 3),c("actual","predicted", "SE"),sep = "_")
GO.LOO.Ctrl.SpD.no.SWM <- data.frame(matrix(nrow = nrow(adj.qMRI.ST.medians.no.SWM), ncol = length(GO.LOO.Ctrl.SpD.colnames)))
colnames(GO.LOO.Ctrl.SpD.no.SWM) <- GO.LOO.Ctrl.SpD.colnames
GO.LOO.Ctrl.SpD.no.SWM <-cbind(adj.qMRI.ST.medians.no.SWM[,1:2],GO.LOO.Ctrl.SpD.no.SWM)

GO.LOO.stats.Ctrl.SpD.no.SWM <- data.frame(matrix(nrow = length(POI.names), ncol = 6))
colnames(GO.LOO.stats.Ctrl.SpD.no.SWM) <- c("GO_type", "NRMSE", "T","rho", "R2", "p")

# run the LOO
GO.LOO.Ctrl.SpD.output.no.SWM <- run.qMRI.ST.LOO.SpD(dframe = adj.qMRI.ST.medians.no.SWM,
                                                       ST = POI.names, 
                                                       LOO = GO.LOO.Ctrl.SpD.no.SWM, 
                                                       stats = GO.LOO.stats.Ctrl.SpD.no.SWM, 
                                                       group.string = "Control")
GO.LOO.Ctrl.SpD.no.SWM <- GO.LOO.Ctrl.SpD.output.no.SWM$LOO
GO.LOO.stats.Ctrl.SpD.no.SWM <- GO.LOO.Ctrl.SpD.output.no.SWM$stats

# save some stats
GO.LOO.stats.Ctrl.SpD.no.SWM$corr_p <- p.adjust(GO.LOO.stats.Ctrl.SpD.no.SWM$p, method = "fdr")
write.csv(GO.LOO.stats.Ctrl.SpD.no.SWM, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_POI_LOO_SpD_no_SWM.csv", sep = ""), 
  row.names = FALSE)

# save the LOO predictions
write.csv(GO.LOO.Ctrl.SpD.no.SWM, file = paste(
  results.dir.stats, "/LOO_SpD_predictions_control_qMRI_vs_POI_no_SWM.csv", sep = ""),
  row.names = FALSE)

# plot the LOO
plot.list.GO.LOO.Ctrl.SpD.no.SWM <- lapply(1:length(POI.names), function(i) 
  plot.LOO(ST = POI.names[i], LOO = GO.LOO.Ctrl.SpD.no.SWM, sts = GO.LOO.stats.Ctrl.SpD.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_POI_LOO_SpD_no_SWM.png", sep = ""),
    width = 1400, height = 700, units = "px", res = 120)
multi.plot.GO.LOO.Ctrl.SpD.no.SWM <- ggarrange(
  plotlist = plot.list.GO.LOO.Ctrl.SpD.no.SWM,
  nrow = 3,
  ncol = 5,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.GO.LOO.Ctrl.SpD.no.SWM)
annotate_figure(multi.plot.GO.LOO.Ctrl.SpD.no.SWM, bottom = text_grob("Control SRT count"),
                left = text_grob("qMRI-predicted cell count", rot = 90),
                top = text_grob("Control qMRI medians predicting control GO pathway score"))
dev.off()

# plot the LOO cook's distance
plot.list.GO.LOO.cooksd.Ctrl.SpD.no.SWM <- lapply(1:length(POI.names), function(i) plot.LOO.SpD.cooksd(ST = POI.names[i], GO.LOO.Ctrl.SpD.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_POI_LOO_SpD_cooksd_no_SWM.png", sep = ""),
    width = 1200, height = 700, units = "px", res = 120)
multi.plot.GO.LOO.cooksd.Ctrl.SpD.no.SWM <- ggarrange(
  plotlist = plot.list.GO.LOO.cooksd.Ctrl.SpD.no.SWM,
  nrow = 3,
  ncol = 5,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.GO.LOO.cooksd.Ctrl.SpD.no.SWM)
annotate_figure(multi.plot.GO.LOO.cooksd.Ctrl.SpD.no.SWM,  
                left = text_grob("Cook's distance", rot = 90),
                top = text_grob("Influence of regions on GO LOO correlations"),
                bottom = text_grob("ROI"))
dev.off()


#### Additional depth-focused analyses ####

# These analyses are added after the complete original analysis. They do not
# alter or replace any of the original correlations, LOO models or plots.

depth.analysis <- subset(adj.qMRI.ST.medians, SpD != "SWM")
qMRI.predictors.depth <- paste(MOI, "Control", sep = "_")
cell.outcomes.depth <- paste(COI, "Control", sep = "_")
GO.outcomes.depth <- paste(POI.names, "Control", sep = "_")
all.ST.outcomes.depth <- c(cell.outcomes.depth, GO.outcomes.depth)

# correlation matrix of depth and qMRI metrics
qMRI.depth.cor <- cor(depth.analysis[,c("depth",qMRI.predictors.depth)],
                      use = "pairwise.complete.obs", method = "pearson")
write.csv(qMRI.depth.cor,
          paste(results.dir.stats.depth,"/qMRI_depth_correlation_matrix.csv",sep = ""))

qMRI.depth.cor.long <- melt(qMRI.depth.cor)
p.qMRI.depth.cor <- ggplot(qMRI.depth.cor.long,aes(Var1,Var2,fill = value)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1)) +
  geom_tile() +
  geom_text(aes(label = round(value,2)),size = 3) +
  scale_fill_gradient2(low = "#2166AC",mid = "white",high = "#B2182B",
                       midpoint = 0,limits = c(-1,1),name = "Pearson r")
ggsave(paste(results.dir.depth,"/qMRI_depth_correlation_matrix.png",sep = ""),
       p.qMRI.depth.cor,width = 7,height = 6,dpi = 300)

# variance inflation factors for the original qMRI predictor set, with and without depth
calculate.VIF.depth <- function(dframe,predictors,model.name) {
  do.call(rbind,lapply(predictors,function(pred) {
    fit <- lm(reformulate(setdiff(predictors,pred),response = pred),data = dframe)
    r2 <- summary(fit)$r.squared
    data.frame(model = model.name,predictor = pred,VIF = 1/(1-r2),
               tolerance = 1-r2)
  }))
}
VIF.depth.results <- rbind(
  calculate.VIF.depth(depth.analysis,qMRI.predictors.depth,"qMRI_only"),
  calculate.VIF.depth(depth.analysis,c("depth",qMRI.predictors.depth),"depth_plus_qMRI")
)
write.csv(VIF.depth.results,
          paste(results.dir.stats.depth,"/variance_inflation_factors.csv",sep = ""),
          row.names = FALSE)

# variance in each qMRI/ST measure explained by depth alone
variance.by.depth <- function(dframe,variables) {
  do.call(rbind,lapply(variables,function(v) {
    fit <- lm(reformulate("depth",response = v),data = dframe)
    sm <- summary(fit)
    data.frame(variable = v,depth_beta = coef(fit)["depth"],
               depth_R2 = sm$r.squared,depth_adjusted_R2 = sm$adj.r.squared,
               depth_p = coef(sm)["depth","Pr(>|t|)"])
  }))
}
qMRI.variance.by.depth <- variance.by.depth(depth.analysis,qMRI.predictors.depth)
qMRI.variance.by.depth$depth_p_FDR <- p.adjust(qMRI.variance.by.depth$depth_p,"fdr")
ST.variance.by.depth <- variance.by.depth(depth.analysis,all.ST.outcomes.depth)
ST.variance.by.depth$family <- ifelse(ST.variance.by.depth$variable %in% cell.outcomes.depth,
                                      "cell_count","pathway")
ST.variance.by.depth$depth_p_FDR <- ave(ST.variance.by.depth$depth_p,
                                        ST.variance.by.depth$family,
                                        FUN = function(x) p.adjust(x,"fdr"))
write.csv(qMRI.variance.by.depth,
          paste(results.dir.stats.depth,"/variance_in_qMRI_explained_by_depth.csv",sep = ""),
          row.names = FALSE)
write.csv(ST.variance.by.depth,
          paste(results.dir.stats.depth,"/variance_in_ST_explained_by_depth.csv",sep = ""),
          row.names = FALSE)

# depth-adjusted residual correlations, retaining the original Pearson approach
run.depth.adjusted.correlations <- function(dframe,ST.variables) {
  out <- expand.grid(qMRI = qMRI.predictors.depth,ST = ST.variables,
                     stringsAsFactors = FALSE)
  out$r <- out$R2 <- out$p <- NA_real_
  for (i in 1:nrow(out)) {
    qMRI.res <- residuals(lm(reformulate("depth",response = out$qMRI[i]),data = dframe))
    ST.res <- residuals(lm(reformulate("depth",response = out$ST[i]),data = dframe))
    test <- cor.test(qMRI.res,ST.res)
    out$r[i] <- test$estimate
    out$R2[i] <- test$estimate^2
    out$p[i] <- test$p.value
  }
  out$corr_p <- p.adjust(out$p,"fdr")
  out
}
cell.depth.adjusted.stats <- run.depth.adjusted.correlations(depth.analysis,cell.outcomes.depth)
GO.depth.adjusted.stats <- run.depth.adjusted.correlations(depth.analysis,GO.outcomes.depth)
write.csv(cell.depth.adjusted.stats,
          paste(results.dir.stats.depth,"/Stats_depth_adjusted_qMRI_vs_cell_count.csv",sep = ""),
          row.names = FALSE)
write.csv(GO.depth.adjusted.stats,
          paste(results.dir.stats.depth,"/Stats_depth_adjusted_qMRI_vs_POI.csv",sep = ""),
          row.names = FALSE)

#### Raw depth and depth-adjusted correlation plots ####

# Match the statistics displayed and saved for the original qMRI correlations:
# T, Pearson rho, R2, p and FDR-adjusted p.
run.depth.raw.correlations <- function(dframe,ST.variables,outcome.name) {
  out <- data.frame(
    depth = rep("depth",length(ST.variables)),
    outcome = ST.variables,
    T = NA_real_,rho = NA_real_,R2 = NA_real_,p = NA_real_
  )
  for (i in 1:length(ST.variables)) {
    test <- cor.test(dframe$depth,dframe[[ST.variables[i]]])
    out$T[i] <- test$statistic
    out$rho[i] <- test$estimate
    out$R2[i] <- test$estimate^2
    out$p[i] <- test$p.value
  }
  out$corr_p <- p.adjust(out$p,method = "fdr")
  colnames(out)[colnames(out) == "outcome"] <- outcome.name
  out
}

cell.depth.raw.stats <- run.depth.raw.correlations(
  depth.analysis,cell.outcomes.depth,"cell_type"
)
GO.depth.raw.stats <- run.depth.raw.correlations(
  depth.analysis,GO.outcomes.depth,"pathway_score"
)
write.csv(
  cell.depth.raw.stats,
  paste(results.dir.stats.depth,"/Stats_depth_vs_cell_count_no_SWM.csv",sep = ""),
  row.names = FALSE
)
write.csv(
  GO.depth.raw.stats,
  paste(results.dir.stats.depth,"/Stats_depth_vs_GO_POI_no_SWM.csv",sep = ""),
  row.names = FALSE
)

# Use the same axis-number formatting as the original correlation plots.
format.correlation.axis <- function(x) {
  ifelse(
    x == 0,"0",
    ifelse(abs(x) < 1e-3,
           format(x,scientific = TRUE,digits = 2),
           signif(x,2))
  )
}

# Extract one legend and remove legends from all component plots. This guarantees
# exactly one legend in every arranged depth-sensitivity figure.
save.single.legend.figure <- function(plot.list,nrow,ncol,filename,
                                      width,height,top = NULL,
                                      bottom = NULL,left = NULL) {
  legend <- get_legend(plot.list[[1]] + theme(legend.position = "right"))
  plot.list.no.legend <- lapply(plot.list,function(p) {
    p + theme(legend.position = "none")
  })
  plot.grid <- ggarrange(
    plotlist = plot.list.no.legend,
    nrow = nrow,ncol = ncol,
    common.legend = FALSE
  )
  plot.with.legend <- ggarrange(
    plot.grid,legend,
    ncol = 2,widths = c(1,0.10)
  )
  labelled <- annotate_figure(
    plot.with.legend,
    top = if (!is.null(top)) text_grob(top) else NULL,
    bottom = if (!is.null(bottom)) text_grob(bottom) else NULL,
    left = if (!is.null(left)) text_grob(left,rot = 90) else NULL
  )
  png(filename = filename,width = width,height = height,units = "px",res = 120)
  print(labelled)
  dev.off()
}

# Raw depth-versus-ST plots, styled like plot.ST.qMRI().
plot.depth.raw <- function(dframe,ST,sts,outcome.column,last.ST) {
  x.data <- dframe$depth
  y.data <- dframe[[ST]]
  R2 <- sts[sts[[outcome.column]] == ST,]$R2
  adj.p <- sts[sts[[outcome.column]] == ST,]$corr_p

  ggplot(data = dframe,aes(x = x.data,y = y.data)) +
    theme_minimal() +
    {
      if (ST == last.ST) xlab("Depth")
      else theme(axis.title.x = element_blank())
    } +
    ylab(gsub("_Control$","",ST)) +
    geom_smooth(
      method = "lm",formula = y ~ x,linetype = "dotted",
      colour = "black",fill = "lightgrey"
    ) +
    geom_point(aes(colour = SpD,shape = region)) +
    scale_color_manual(values = SpD.cols) +
    scale_x_continuous(labels = format.correlation.axis) +
    scale_y_continuous(labels = format.correlation.axis) +
    annotate(
      "text",x = -Inf,y = Inf,hjust = 0,vjust = 1,
      label = as.expression(bquote(
        R^2 == .(signif(R2,2)) ~ "," ~ p[FDR] == .(
          if (adj.p < 1e-3) format(adj.p,digits = 2,scientific = TRUE)
          else signif(adj.p,2)
        )
      )),
      size = 3.25,parse = TRUE
    )
}

cell.depth.raw.plots <- lapply(cell.outcomes.depth,function(st) {
  plot.depth.raw(
    depth.analysis,st,cell.depth.raw.stats,"cell_type",
    cell.outcomes.depth[length(cell.outcomes.depth)]
  )
})
save.single.legend.figure(
  cell.depth.raw.plots,4,1,
  paste(results.dir.depth,"/Plot_depth_vs_cell_count_no_SWM.png",sep = ""),
  500,900,
  top = "Depth versus control cell count"
)

GO.depth.raw.plots <- lapply(GO.outcomes.depth,function(st) {
  plot.depth.raw(
    depth.analysis,st,GO.depth.raw.stats,"pathway_score",
    GO.outcomes.depth[length(GO.outcomes.depth)]
  )
})
save.single.legend.figure(
  GO.depth.raw.plots,length(GO.outcomes.depth),1,
  paste(results.dir.depth,"/Plot_depth_vs_GO_POI_no_SWM.png",sep = ""),
  600,max(1200,170*length(GO.outcomes.depth)),
  top = "Depth versus control pathway score"
)

# Depth-adjusted qMRI-versus-ST plots. Both variables are residualised against
# depth before plotting
plot.depth.adjusted.qMRI <- function(dframe,qMRI.variable,ST.variable,sts,
                                     first.qMRI,last.ST) {
  qMRI.res <- residuals(lm(reformulate("depth",response = qMRI.variable),data = dframe))
  ST.res <- residuals(lm(reformulate("depth",response = ST.variable),data = dframe))
  plot.data <- data.frame(
    SpD = dframe$SpD,region = dframe$region,
    qMRI.residual = qMRI.res,ST.residual = ST.res
  )
  R2 <- sts[sts$qMRI == qMRI.variable & sts$ST == ST.variable,]$R2
  adj.p <- sts[sts$qMRI == qMRI.variable & sts$ST == ST.variable,]$corr_p

  ggplot(plot.data,aes(x = qMRI.residual,y = ST.residual)) +
    theme_minimal() +
    {
      if (ST.variable == last.ST) xlab(gsub("_Control$","",qMRI.variable))
      else theme(axis.title.x = element_blank())
    } +
    {
      if (qMRI.variable == first.qMRI) ylab(gsub("_Control$","",ST.variable))
      else theme(axis.title.y = element_blank())
    } +
    geom_smooth(
      method = "lm",formula = y ~ x,linetype = "dotted",
      colour = "black",fill = "lightgrey"
    ) +
    geom_point(aes(colour = SpD,shape = region)) +
    scale_color_manual(values = SpD.cols) +
    scale_x_continuous(labels = format.correlation.axis) +
    scale_y_continuous(labels = format.correlation.axis) +
    annotate(
      "text",x = -Inf,y = Inf,hjust = 0,vjust = 1,
      label = as.expression(bquote(
        R^2 == .(signif(R2,2)) ~ "," ~ p[FDR] == .(
          if (adj.p < 1e-3) format(adj.p,digits = 2,scientific = TRUE)
          else signif(adj.p,2)
        )
      )),
      size = 3.25,parse = TRUE
    )
}

cell.depth.adjusted.plots <- lapply(cell.outcomes.depth,function(st) {
  lapply(qMRI.predictors.depth,function(qmri) {
    plot.depth.adjusted.qMRI(
      depth.analysis,qmri,st,cell.depth.adjusted.stats,
      qMRI.predictors.depth[1],cell.outcomes.depth[length(cell.outcomes.depth)]
    )
  })
})
cell.depth.adjusted.plots <- unlist(cell.depth.adjusted.plots,recursive = FALSE)
save.single.legend.figure(
  cell.depth.adjusted.plots,length(cell.outcomes.depth),length(qMRI.predictors.depth),
  paste(results.dir.depth,"/Plot_depth_adjusted_qMRI_vs_cell_count_no_SWM.png",sep = ""),
  1200,800,
  top = "Depth-adjusted control qMRI versus control cell count"
)

GO.depth.adjusted.plots <- lapply(GO.outcomes.depth,function(st) {
  lapply(qMRI.predictors.depth,function(qmri) {
    plot.depth.adjusted.qMRI(
      depth.analysis,qmri,st,GO.depth.adjusted.stats,
      qMRI.predictors.depth[1],GO.outcomes.depth[length(GO.outcomes.depth)]
    )
  })
})
GO.depth.adjusted.plots <- unlist(GO.depth.adjusted.plots,recursive = FALSE)
save.single.legend.figure(
  GO.depth.adjusted.plots,length(GO.outcomes.depth),length(qMRI.predictors.depth),
  paste(results.dir.depth,"/Plot_depth_adjusted_qMRI_vs_GO_POI_no_SWM.png",sep = ""),
  1200,max(2500,170*length(GO.outcomes.depth)),
  top = "Depth-adjusted control qMRI versus control pathway score"
)

# Run the depth comparison with the same two validation schemes used in the
# original qMRI analysis:
#   region-level LOO: one region-depth row held out at a time
#   SpD-level LOO: all three regions from one SpD held out together
#
# The three models retain the original qMRI predictor specification:
#   Depth:     depth only
#   qMRI:      the original five qMRI metrics
#   Combined:  depth plus the original five qMRI metrics

run.depth.LOO <- function(dframe,ST.variables,LOO.level = c("region","SpD")) {
  LOO.level <- match.arg(LOO.level)
  model.predictors <- list(
    Depth = "depth",
    qMRI = qMRI.predictors.depth,
    Combined = c("depth",qMRI.predictors.depth)
  )
  predictions <- data.frame()
  stats <- data.frame()

  for (st in ST.variables) {
    st.predictions <- data.frame()
    held.out.groups <- if (LOO.level == "region") seq_len(nrow(dframe)) else unique(dframe$SpD)

    for (held.out in held.out.groups) {
      test.index <- if (LOO.level == "region") {
        seq_len(nrow(dframe)) == held.out
      } else {
        dframe$SpD == held.out
      }
      train <- dframe[!test.index,]
      test <- dframe[test.index,]

      for (model.name in names(model.predictors)) {
        loo.fm <- reformulate(model.predictors[[model.name]],response = st)
        loo.fit <- lm(loo.fm,data = train)
        loo.pred <- predict(loo.fit,newdata = test)
        st.predictions <- rbind(st.predictions,data.frame(
          ST = st,
          model = model.name,
          SpD = test$SpD,
          region = test$region,
          actual = test[[st]],
          predicted = as.numeric(loo.pred),
          SE = (as.numeric(loo.pred)-test[[st]])^2
        ))
      }
    }

    st.predictions$cooksd <- NA_real_
    for (model.name in names(model.predictors)) {
      model.predictions <- st.predictions[st.predictions$model == model.name,]
      fit <- lm(actual ~ predicted,data = model.predictions)
      res <- summary(fit)
      model.predictions$cooksd <- cooks.distance(fit)
      st.predictions$cooksd[st.predictions$model == model.name] <- model.predictions$cooksd

      stats <- rbind(stats,data.frame(
        ST = st,
        model = model.name,
        NRMSE = nrmse(model.predictions$predicted,model.predictions$actual),
        T = res$coefficients[2,3],
        rho = cor.test(model.predictions$predicted,model.predictions$actual)$estimate,
        R2 = res$r.squared,
        p = res$coefficients[2,4]
      ))
    }
    predictions <- rbind(predictions,st.predictions)
  }

  # FDR across outcomes within each model.
  stats$corr_p <- ave(stats$p,stats$model,FUN = function(x) p.adjust(x,method = "fdr"))

  # Direct comparison with the depth-only baseline, kept as extra columns after
  # the original NRMSE/T/rho/R2/p/corr_p statistics.
  stats$delta_NRMSE_vs_depth <- NA_real_
  stats$delta_R2_vs_depth <- NA_real_
  for (st in ST.variables) {
    depth.NRMSE <- stats$NRMSE[stats$ST == st & stats$model == "Depth"]
    depth.R2 <- stats$R2[stats$ST == st & stats$model == "Depth"]
    stats$delta_NRMSE_vs_depth[stats$ST == st] <-
      depth.NRMSE-stats$NRMSE[stats$ST == st]
    stats$delta_R2_vs_depth[stats$ST == st] <-
      stats$R2[stats$ST == st]-depth.R2
  }

  return(list(LOO = predictions,stats = stats))
}

# Plot 
plot.depth.model.comparison <- function(ST,LOO,sts) {
  plot.data <- LOO[LOO$ST == ST,]
  plot.data$model <- factor(plot.data$model,
                            levels = c("Depth","qMRI","Combined"))
  annotation.data <- sts[sts$ST == ST,]
  annotation.data$model <- factor(annotation.data$model,
                                  levels = c("Depth","qMRI","Combined"))
  annotation.data$label <- mapply(function(R2,adj.p) {
    p.text <- if (adj.p < 1e-3) {
      format(adj.p,digits = 2,scientific = TRUE)
    } else {
      signif(adj.p,2)
    }
    paste0("R^2 == ",signif(R2,2)," ~ ',' ~ p[FDR] == ",p.text)
  },annotation.data$R2,annotation.data$corr_p)

  ggplot(plot.data,aes(x = actual,y = predicted)) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank()) +
    geom_smooth(method = "lm",formula = y ~ x,linetype = "dotted",
                colour = "black",fill = "lightgrey") +
    geom_point(aes(colour = SpD,shape = region)) +
    scale_color_manual(values = SpD.cols) +
    scale_x_continuous(labels = function(x) {
      ifelse(x == 0,"0",ifelse(abs(x) < 1e-3,
        format(x,scientific = TRUE,digits = 2),signif(x,2)))
    }) +
    scale_y_continuous(labels = function(y) {
      ifelse(y == 0,"0",ifelse(abs(y) < 1e-3,
        format(y,scientific = TRUE,digits = 2),signif(y,2)))
    }) +
    facet_wrap(~model,nrow = 1,scales = "free") +
    geom_text(data = annotation.data,aes(x = -Inf,y = Inf,label = label),
              inherit.aes = FALSE,hjust = 0,vjust = 1,size = 3.25,parse = TRUE)
}

# Arrange outcomes like the original figures. The outcome name is a separate row
# title, while each row contains the short model headings Depth/qMRI/Combined.
save.depth.comparison.plot <- function(result,ST.variables,prefix,LOO.label,
                                       width,height) {
  base.plots <- lapply(ST.variables,function(st) {
    plot.depth.model.comparison(st,result$LOO,result$stats)
  })
  legend <- get_legend(base.plots[[1]] + theme(legend.position = "right"))
  outcome.rows <- lapply(seq_along(ST.variables),function(i) {
    annotate_figure(
      base.plots[[i]] + theme(legend.position = "none"),
      top = text_grob(gsub("_Control$","",ST.variables[i]),face = "bold",size = 10)
    )
  })
  combined.rows <- ggarrange(
    plotlist = outcome.rows,nrow = length(ST.variables),ncol = 1,
    common.legend = FALSE
  )
  combined <- ggarrange(
    combined.rows,legend,ncol = 2,widths = c(1,0.10)
  )
  labelled <- annotate_figure(
    combined,
    bottom = text_grob("Control SRT count"),
    left = text_grob("Predicted cell count / pathway score",rot = 90),
    top = text_grob(paste("Depth and qMRI model comparison:",LOO.label),
                    face = "bold")
  )
  png(paste(results.dir.depth,"/",prefix,"_model_comparison.png",sep = ""),
      width = width,height = height,units = "px",res = 120)
  print(labelled)
  dev.off()
}
# Format comparison statistics to preserve the original CSV order and names.
format.depth.stats <- function(sts,outcome.name) {
  colnames(sts)[colnames(sts) == "ST"] <- outcome.name
  sts[,c(outcome.name,"NRMSE","T","rho","R2","p","corr_p","model",
         "delta_NRMSE_vs_depth","delta_R2_vs_depth")]
}

# CELL COUNTS: region-level and SpD-level LOO
cell.depth.LOO.region <- run.depth.LOO(depth.analysis,cell.outcomes.depth,"region")
cell.depth.LOO.SpD <- run.depth.LOO(depth.analysis,cell.outcomes.depth,"SpD")
write.csv(format.depth.stats(cell.depth.LOO.region$stats,"cell_type"),
          paste(results.dir.stats.depth,"/Stats_depth_qMRI_vs_cell_count_LOO_no_SWM.csv",sep = ""),
          row.names = FALSE)
write.csv(cell.depth.LOO.region$LOO,
          paste(results.dir.stats.depth,"/LOO_predictions_depth_qMRI_vs_cell_count_no_SWM.csv",sep = ""),
          row.names = FALSE)
write.csv(format.depth.stats(cell.depth.LOO.SpD$stats,"cell_type"),
          paste(results.dir.stats.depth,"/Stats_depth_qMRI_vs_cell_count_LOO_SpD_no_SWM.csv",sep = ""),
          row.names = FALSE)
write.csv(cell.depth.LOO.SpD$LOO,
          paste(results.dir.stats.depth,"/LOO_SpD_predictions_depth_qMRI_vs_cell_count_no_SWM.csv",sep = ""),
          row.names = FALSE)

save.depth.comparison.plot(
  cell.depth.LOO.region,cell.outcomes.depth,
  "Plot_depth_qMRI_vs_cell_count_LOO_no_SWM",
  "region-level LOO",1200,1100
)
save.depth.comparison.plot(
  cell.depth.LOO.SpD,cell.outcomes.depth,
  "Plot_depth_qMRI_vs_cell_count_LOO_SpD_no_SWM",
  "SpD-level LOO",1200,1100
)

# PATHWAYS: region-level and SpD-level LOO
GO.depth.LOO.region <- run.depth.LOO(depth.analysis,GO.outcomes.depth,"region")
GO.depth.LOO.SpD <- run.depth.LOO(depth.analysis,GO.outcomes.depth,"SpD")
write.csv(format.depth.stats(GO.depth.LOO.region$stats,"GO_type"),
          paste(results.dir.stats.depth,"/Stats_depth_qMRI_vs_POI_LOO_no_SWM.csv",sep = ""),
          row.names = FALSE)
write.csv(GO.depth.LOO.region$LOO,
          paste(results.dir.stats.depth,"/LOO_predictions_depth_qMRI_vs_POI_no_SWM.csv",sep = ""),
          row.names = FALSE)
write.csv(format.depth.stats(GO.depth.LOO.SpD$stats,"GO_type"),
          paste(results.dir.stats.depth,"/Stats_depth_qMRI_vs_POI_LOO_SpD_no_SWM.csv",sep = ""),
          row.names = FALSE)
write.csv(GO.depth.LOO.SpD$LOO,
          paste(results.dir.stats.depth,"/LOO_SpD_predictions_depth_qMRI_vs_POI_no_SWM.csv",sep = ""),
          row.names = FALSE)

save.depth.comparison.plot(
  GO.depth.LOO.region,GO.outcomes.depth,
  "Plot_depth_qMRI_vs_POI_LOO_no_SWM",
  "region-level LOO",1200,3000
)
save.depth.comparison.plot(
  GO.depth.LOO.SpD,GO.outcomes.depth,
  "Plot_depth_qMRI_vs_POI_LOO_SpD_no_SWM",
  "SpD-level LOO",1200,3000
)

#### Descriptive comparison of cross-validated model performance ####

# The existing LOO statistics tables report the performance of the individual
# Depth, qMRI and Combined models. The analyses below add direct descriptive
# comparisons of held-out prediction error for all three model pairs.
#
# No significance tests are performed.
#
# For region-level LOO, the comparison unit is the held-out region-depth row.
#
# For SpD-level LOO, prediction errors are first averaged within each held-out
# SpD. Therefore, the comparison unit is the held-out depth domain rather than
# an individual regional observation.
#
# Comparisons are labelled using:
#
#   first-named model over comparison model
#
# For example:
#
#   qMRI over Depth
#
# Delta error values are calculated as:
#
#   error_first_named_model - error_comparison_model
#
# Therefore:
#
#   negative value = first-named model has lower error
#   positive value = comparison model has lower error
#   zero           = models have equal error


compare.LOO.models <- function(
    LOO,
    LOO.level,
    family) {
  
  # model.1 is the comparison model
  # model.2 is the first-named model in the displayed comparison
  model.pairs <- list(
    c("Depth", "qMRI"),
    c("Depth", "Combined"),
    c("qMRI", "Combined")
  )
  
  results <- data.frame()
  
  for (st in unique(LOO$ST)) {
    
    st.data <- LOO[
      LOO$ST == st,
    ]
    
    for (pair in model.pairs) {
      
      model.1 <- pair[1]
      model.2 <- pair[2]
      
      data.1 <- st.data[
        st.data$model == model.1,
      ]
      
      data.2 <- st.data[
        st.data$model == model.2,
      ]
      
      # Ensure that predictions from both models are placed in the same
      # SpD-region order before their prediction errors are compared
      order.1 <- order(
        data.1$SpD,
        data.1$region
      )
      
      order.2 <- order(
        data.2$SpD,
        data.2$region
      )
      
      data.1 <- data.1[
        order.1,
      ]
      
      data.2 <- data.2[
        order.2,
      ]
      
      stopifnot(
        identical(
          as.character(data.1$SpD),
          as.character(data.2$SpD)
        ),
        identical(
          as.character(data.1$region),
          as.character(data.2$region)
        ),
        isTRUE(
          all.equal(
            data.1$actual,
            data.2$actual
          )
        )
      )
      
      # For SpD-level LOO, calculate one mean squared error and one
      # mean absolute error for each held-out SpD
      if (LOO.level == "SpD") {
        
        SpD.levels <- unique(
          as.character(data.1$SpD)
        )
        
        paired.loss <- data.frame(
          unit = SpD.levels,
          
          MSE.model.1 = vapply(
            SpD.levels,
            function(x) {
              
              mean(
                data.1$SE[
                  as.character(data.1$SpD) == x
                ]
              )
            },
            numeric(1)
          ),
          
          MSE.model.2 = vapply(
            SpD.levels,
            function(x) {
              
              mean(
                data.2$SE[
                  as.character(data.2$SpD) == x
                ]
              )
            },
            numeric(1)
          ),
          
          MAE.model.1 = vapply(
            SpD.levels,
            function(x) {
              
              mean(
                abs(
                  data.1$predicted[
                    as.character(data.1$SpD) == x
                  ] -
                    data.1$actual[
                      as.character(data.1$SpD) == x
                    ]
                )
              )
            },
            numeric(1)
          ),
          
          MAE.model.2 = vapply(
            SpD.levels,
            function(x) {
              
              mean(
                abs(
                  data.2$predicted[
                    as.character(data.2$SpD) == x
                  ] -
                    data.2$actual[
                      as.character(data.2$SpD) == x
                    ]
                )
              )
            },
            numeric(1)
          )
        )
        
      } else {
        
        # For region-level LOO, each of the 12 held-out region-depth
        # observations is one comparison unit
        paired.loss <- data.frame(
          unit = seq_len(
            nrow(data.1)
          ),
          
          MSE.model.1 =
            data.1$SE,
          
          MSE.model.2 =
            data.2$SE,
          
          MAE.model.1 =
            abs(
              data.1$predicted -
                data.1$actual
            ),
          
          MAE.model.2 =
            abs(
              data.2$predicted -
                data.2$actual
            )
        )
      }
      
      RMSE.1 <- sqrt(
        mean(data.1$SE)
      )
      
      RMSE.2 <- sqrt(
        mean(data.2$SE)
      )
      
      MAE.1 <- mean(
        abs(
          data.1$predicted -
            data.1$actual
        )
      )
      
      MAE.2 <- mean(
        abs(
          data.2$predicted -
            data.2$actual
        )
      )
      
      NRMSE.1 <- nrmse(
        data.1$predicted,
        data.1$actual
      )
      
      NRMSE.2 <- nrmse(
        data.2$predicted,
        data.2$actual
      )
      
      # First-named model is model.2.
      #
      # Therefore, negative differences favour the first-named model:
      #
      #   error_model.2 - error_model.1
      #
      # Example:
      #
      #   qMRI over Depth
      #
      #   delta_NRMSE = NRMSE_qMRI - NRMSE_Depth
      #
      # Negative = qMRI has lower NRMSE
      # Positive = Depth has lower NRMSE
      
      MSE.difference <-
        paired.loss$MSE.model.2 -
        paired.loss$MSE.model.1
      
      MAE.difference <-
        paired.loss$MAE.model.2 -
        paired.loss$MAE.model.1
      
      results <- rbind(
        results,
        data.frame(
          family = family,
          
          LOO_level =
            LOO.level,
          
          ST = st,
          
          comparison = paste(
            model.2,
            "over",
            model.1
          ),
          
          first_named_model =
            model.2,
          
          comparison_model =
            model.1,
          
          n_comparison_units =
            nrow(paired.loss),
          
          RMSE_first_named =
            RMSE.2,
          
          RMSE_comparison =
            RMSE.1,
          
          # Negative values favour the first-named model
          delta_RMSE =
            RMSE.2 -
            RMSE.1,
          
          MAE_first_named =
            MAE.2,
          
          MAE_comparison =
            MAE.1,
          
          # Negative values favour the first-named model
          delta_MAE =
            MAE.2 -
            MAE.1,
          
          NRMSE_first_named =
            NRMSE.2,
          
          NRMSE_comparison =
            NRMSE.1,
          
          # Negative values favour the first-named model
          delta_NRMSE =
            NRMSE.2 -
            NRMSE.1,
          
          # Negative values indicate lower squared prediction error
          # for the first-named model
          mean_MSE_difference =
            mean(MSE.difference),
          
          median_MSE_difference =
            median(MSE.difference),
          
          # Negative values indicate lower absolute prediction error
          # for the first-named model
          mean_MAE_difference =
            mean(MAE.difference),
          
          median_MAE_difference =
            median(MAE.difference),
          
          first_named_model_better_units =
            sum(MSE.difference < 0),
          
          comparison_model_better_units =
            sum(MSE.difference > 0),
          
          tied_units =
            sum(MSE.difference == 0)
        )
      )
    }
  }
  
  return(results)
}


#### Run pairwise model comparisons ####

cell.LOO.model.comparisons <- rbind(
  
  compare.LOO.models(
    cell.depth.LOO.region$LOO,
    "region",
    "cell_count"
  ),
  
  compare.LOO.models(
    cell.depth.LOO.SpD$LOO,
    "SpD",
    "cell_count"
  )
)


GO.LOO.model.comparisons <- rbind(
  
  compare.LOO.models(
    GO.depth.LOO.region$LOO,
    "region",
    "pathway"
  ),
  
  compare.LOO.models(
    GO.depth.LOO.SpD$LOO,
    "SpD",
    "pathway"
  )
)


LOO.model.comparisons <- rbind(
  cell.LOO.model.comparisons,
  GO.LOO.model.comparisons
)


#### Save descriptive model-comparison statistics ####

write.csv(
  cell.LOO.model.comparisons,
  paste(
    results.dir.stats.depth,
    "/Stats_pairwise_LOO_model_comparisons_cell_count.csv",
    sep = ""
  ),
  row.names = FALSE
)


write.csv(
  GO.LOO.model.comparisons,
  paste(
    results.dir.stats.depth,
    "/Stats_pairwise_LOO_model_comparisons_POI.csv",
    sep = ""
  ),
  row.names = FALSE
)


write.csv(
  LOO.model.comparisons,
  paste(
    results.dir.stats.depth,
    "/Stats_pairwise_LOO_model_comparisons_all.csv",
    sep = ""
  ),
  row.names = FALSE
)


#### Plot pairwise LOO model comparisons ####

# Comparison labels name the model being evaluated first:
#
#   qMRI over Depth
#   Combined over Depth
#   Combined over qMRI
#
# Delta NRMSE is calculated as:
#
#   NRMSE_first_named_model - NRMSE_comparison_model
#
# Therefore:
#
#   negative values favour the first-named model
#   positive values favour the comparison model
#   zero indicates equal NRMSE
plot.LOO.model.comparisons <- function(
    comparison.data,
    family.name,
    filename,
    width,
    height) {
  
  plot.data <- comparison.data[
    comparison.data$family ==
      family.name,
  ]
  
  plot.data$ST_label <- gsub(
    "_Control$",
    "",
    plot.data$ST
  )
  
  plot.data$comparison <- factor(
    plot.data$comparison,
    levels = c(
      "qMRI over Depth",
      "Combined over Depth",
      "Combined over qMRI"
    )
  )
  
  plot.data$LOO_level <- factor(
    plot.data$LOO_level,
    levels = c(
      "region",
      "SpD"
    ),
    labels = c(
      "Region-level LOO",
      "SpD-level LOO"
    )
  )
  
  # If delta_NRMSE is negative, the first-named model has lower NRMSE.
  #
  # If delta_NRMSE is positive, the comparison model has lower NRMSE.
  #
  # The label therefore reports the number of held-out units won by
  # whichever model has the lower overall NRMSE.
  plot.data$better.units <- ifelse(
    plot.data$delta_NRMSE < 0,
    
    plot.data$first_named_model_better_units,
    
    ifelse(
      plot.data$delta_NRMSE > 0,
      
      plot.data$comparison_model_better_units,
      
      plot.data$tied_units
    )
  )
  
  plot.data$units.label <- ifelse(
    plot.data$delta_NRMSE == 0,
    
    "Equal NRMSE",
    
    paste0(
      plot.data$better.units,
      "/",
      plot.data$n_comparison_units
      )
  )
  
  # Place labels just beyond the end of each bar.
  #
  # Extra expansion is added below so that labels near the upper or
  # lower edge of a subplot are not cut off.
  plot.data$label.vjust <- ifelse(
    plot.data$delta_NRMSE < 0,
    1.35,
    -0.35
  )
  
  # Centre the label on the zero line if both models have equal NRMSE.
  plot.data$label.vjust[
    plot.data$delta_NRMSE == 0
  ] <- -0.35
  
  p <- ggplot(
    plot.data,
    aes(
      x = comparison,
      y = delta_NRMSE,
      fill = comparison
    )
  ) +
    
    theme_minimal() +
    
    theme(
      axis.title.x =
        element_blank(),
      
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        ),
      
      legend.position =
        "none",
      
      strip.background =
        element_blank(),
      
      # Additional outer space ensures that unclipped labels remain
      # inside the complete saved figure.
      plot.margin =
        margin(
          t = 12,
          r = 10,
          b = 12,
          l = 10
        )
    ) +
    
    # Zero indicates equal NRMSE.
    geom_hline(
      yintercept = 0,
      linetype = "dotted"
    ) +
    
    geom_col(
      width = 0.7
    ) +
    
    # Report the number of held-out units won by the model with the
    # lower overall NRMSE.
    geom_text(
      aes(
        label = units.label
      ),
      vjust = plot.data$label.vjust,
      size = 3
    ) +
    
    facet_grid(
      ST_label ~ LOO_level,
      scales = "free_y"
    ) +
    
    ylab(
      paste(
        expression(Delta * "NRMSE"),
        "Negative values favour the first-named model",
        sep = "\n"
      )
    ) +
    
    # Add space beyond the most extreme positive and negative bars
    # so that labels remain visible.
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.18,
          0.18
        )
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "qMRI over Depth" =
          "#66CCEE",
        
        "Combined over Depth" =
          "#228833",
        
        "Combined over qMRI" =
          "#EE6677"
      )
    ) +
    
    # Prevent labels positioned beyond the bar ends from being
    # clipped at the panel boundary.
    coord_cartesian(
      clip = "off"
    )
  
  png(
    filename = filename,
    width = width,
    height = height,
    units = "px",
    res = 120
  )
  
  print(p)
  
  dev.off()
}


#### Save cell-count model-comparison plot ####

plot.LOO.model.comparisons(
  LOO.model.comparisons,
  "cell_count",
  paste(
    results.dir.depth,
    "/Plot_pairwise_LOO_model_comparisons_cell_count.png",
    sep = ""
  ),
  900,
  900
)


#### Save pathway model-comparison plot ####

plot.LOO.model.comparisons(
  LOO.model.comparisons,
  "pathway",
  paste(
    results.dir.depth,
    "/Plot_pairwise_LOO_model_comparisons_POI.png",
    sep = ""
  ),
  900,
  max(
    1800,
    170 *
      length(GO.outcomes.depth)
  )
)