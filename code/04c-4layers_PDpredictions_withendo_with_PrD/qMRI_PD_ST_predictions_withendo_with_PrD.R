# Using control qMRI data to fit LMs and project PD ST scores

# load in some libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(hydroGOF)
library(matrixStats)

# define directories
d0 <- getwd();
qMRI.data.dir <- paste(d0, "/../../processed_data/02-4layers_qMRI/", sep = "")
cell.data.dir <- paste(d0, "/../../processed_data/01b-4layers_withendo/02-top3_cells_counts/", sep = "")
GO.data.dir <- paste(d0, "/../../processed_data/01b-4layers_withendo/05-pathway_scoring/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir.qMRIonly <- paste(d0, "/../../plots/00-4layers_qMRI/", sep = "")
results.dir <- paste(d0, "/../../plots/04c-4layers_PDpredictions_withendo_with_PrD/", sep = "")
results.dir2 <- paste(d0, "/../../plots/02c-4layers_qMRIvsST_withendo_with_PrD/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/04c-4layers_PDpredictions_withendo_with_PrD/", sep = "")
results.dir2.stats <- paste(d0, "/../../stats/02c-4layers_qMRIvsST_withendo_with_PrD/", sep = "")

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

# may not need this - spatial domain naming is different from "layer" naming
# # make a function to rename spatial domains like so:
# # SpD4 --> SpD1
# # SpD3 --> SpD2
# # SpD2 --> SpD3
# # SpD1 --> SpD4
# # WM --> SWM
# rename_SpDs <- function(x) {
#   x <- gsub("_WM", "_W", x)
#   x <- gsub("_SpD1", "_X", x)
#   x <- gsub("_SpD2", "_Y", x)
#   x <- gsub("_SpD3", "_Z", x)
#   x <- gsub("_SpD4", "_SpD1", x)
#   x <- gsub("_Z", "_SpD2", x)
#   x <- gsub("_Y", "_SpD3", x)
#   x <- gsub("_X", "_SpD4", x)
#   x <- gsub("_W", "_SWM", x)
#   return(x)
# }

#### Load in & reorganise qMRI data ####

# load in the qMRI data and demographics for the DLPFC sub-regions
qMRI.dems <- read.csv(paste(dems.dir,"anon_DLPFC_qMRI_demographics.csv",sep = ""))
qMRI.data <-  read.csv(paste(qMRI.data.dir,"anon_DLPFC_AMP_qMRI_medians_4_equivol_layers.csv",sep = ""))
qMRI.ACG.data <- read.csv(paste(qMRI.data.dir,"anon_ACG_AMP_qMRI_medians_4_equivol_layers.csv",sep = ""))

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

# do this also for the ACG data
qMRI.ACG.data.ant <- qMRI.ACG.data[grepl("ant", qMRI.ACG.data$subject_region),2:ncol(qMRI.ACG.data)]
colnames(qMRI.ACG.data.ant) <- paste(colnames(qMRI.ACG.data.ant),"_ant", sep = "")
# middle
qMRI.ACG.data.mid <- qMRI.ACG.data[grepl("mid", qMRI.ACG.data$subject_region),2:ncol(qMRI.ACG.data)]
colnames(qMRI.ACG.data.mid) <- paste(colnames(qMRI.ACG.data.mid),"_mid", sep = "")
# posterior
qMRI.ACG.data.pos <- qMRI.ACG.data[grepl("pos", qMRI.ACG.data$subject_region),2:ncol(qMRI.ACG.data)]
colnames(qMRI.ACG.data.pos) <- paste(colnames(qMRI.ACG.data.pos),"_pos", sep = "")

# merge them together and add the dems
qMRI.data.dems <- cbind(qMRI.dems,qMRI.data.ant, qMRI.data.mid, qMRI.data.pos)
rownames(qMRI.data.dems) <- 1:nrow(qMRI.data.dems)
qMRI.ACG.data.dems <- cbind(qMRI.dems,qMRI.ACG.data.ant, qMRI.ACG.data.mid, qMRI.ACG.data.pos)
rownames(qMRI.ACG.data.dems) <- 1:nrow(qMRI.ACG.data.dems)

# load in the QC info sheet 
QC.info <- read.csv(paste(dems.dir,"anon_7T_QC_info.csv",sep = ""))

# specify the metrics we are interested in
qMRI.metrics <- c("QSM","R1","R2s","PrD","MTsat")
MOI <- c("QSM","R1","R2s","PrD","MTsat")

# match up QC and qMRI data frames
QC.info <- subset(QC.info, QC.info$subject %in% qMRI.data.dems$subject)
rownames(QC.info) <- 1:nrow(QC.info)

# remove data from demographics and qMRI data using QC info
# if they failed MPRAGE QC or FS QC, then remove from everywhere
qMRI.data.dems <- subset(qMRI.data.dems, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
qMRI.ACG.data.dems <- subset(qMRI.ACG.data.dems, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
QC.info <- subset(QC.info, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)

# rename columns
colnames(qMRI.data.dems) <- lapply(colnames(qMRI.data.dems), rename_layers)
colnames(qMRI.ACG.data.dems) <- lapply(colnames(qMRI.ACG.data.dems), rename_layers)
# get all ROI names
qMRI.ROI.names <- colnames(qMRI.data.dems)[5:ncol(qMRI.data.dems)]

# need to create some new columns in the QC.info dataframe
QC.info$QSM_QC <- QC.info$QSM_QC == 1 & QC.info$has_QSM == 1
QC.info$R1_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$R2s_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$PrD_QC <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
QC.info$MTsat_QC <- QC.info$MPM_QC == 1 & QC.info$has_MT == 1

#### adjust qMRI values ####

# initialize a data frame for the adjusted values
adj.qMRI.data.dems <- qMRI.data.dems
adj.qMRI.ACG.data.dems <- qMRI.ACG.data.dems

# get mean and diff from mean age and sex for all subjects
mn.qMRI.age <- mean(qMRI.data.dems$age)
diff.qMRI.age <- mn.qMRI.age - qMRI.data.dems$age
mn.qMRI.sex <- mean(as.numeric(qMRI.data.dems$sex))
diff.qMRI.sex <- mn.qMRI.sex - as.numeric(qMRI.data.dems$sex)

# adjust all ROIs for age and sex using the covariance method
for (ROI in qMRI.ROI.names) {
  # DLPFC
  mdl.qMRI.adj <- lm(qMRI.data.dems[,ROI] ~ age + sex, data = qMRI.data.dems)
  adj.qMRI.vals <- qMRI.data.dems[,ROI] +
    diff.qMRI.age * mdl.qMRI.adj$coefficients[2] +
    diff.qMRI.sex * mdl.qMRI.adj$coefficients[3]
  adj.qMRI.data.dems[,ROI] <- adj.qMRI.vals
  # ACG
  mdl.qMRI.ACG.adj <- lm(qMRI.ACG.data.dems[,ROI] ~ age + sex, data = qMRI.ACG.data.dems)
  adj.qMRI.ACG.vals <- qMRI.ACG.data.dems[,ROI] +
    diff.qMRI.age * mdl.qMRI.ACG.adj$coefficients[2] +
    diff.qMRI.sex * mdl.qMRI.ACG.adj$coefficients[3]
  adj.qMRI.ACG.data.dems[,ROI] <- adj.qMRI.ACG.vals
}


#### qMRI ACG boxplots ----

# subjects
subj <- adj.qMRI.ACG.data.dems$subject

# define spatial domains / subregions
dom <- c("SWM", "SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

# colours for spatial domains
SpD.cols <- c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677")

# include in boxplots?
adj.qMRI.ACG.data.dems$inc_QSM <- QC.info$QSM_QC == 1 & QC.info$has_QSM == 1
adj.qMRI.ACG.data.dems$inc_R1 <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
adj.qMRI.ACG.data.dems$inc_R2s <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
adj.qMRI.ACG.data.dems$inc_PrD <- QC.info$MPM_QC == 1 & QC.info$has_MPM == 1
adj.qMRI.ACG.data.dems$inc_MTsat <- QC.info$MPM_QC == 1 & QC.info$has_MT == 1

# do this for qMRI
adj.qMRI.ACG.data.long <- expand.grid(subject = subj, SpD = dom, region = reg)
adj.qMRI.ACG.data.dems.melt <- melt(adj.qMRI.ACG.data.dems,id.vars=c("subject","age","sex","diagnosis",
                                                             "inc_QSM","inc_R1","inc_R2s","inc_PrD","inc_MTsat"))

for (m in qMRI.metrics){
  metric.data <- adj.qMRI.ACG.data.dems.melt[grepl(m, adj.qMRI.ACG.data.dems.melt$variable),11]
  QC.data <- adj.qMRI.ACG.data.dems.melt[grepl(m, adj.qMRI.ACG.data.dems.melt$variable),paste("inc_",m,sep = "")]
  adj.qMRI.ACG.data.long[,paste(m, sep = "")] <- metric.data
  adj.qMRI.ACG.data.long[,paste("inc_", m, sep = "")] <- QC.data
}

adj.qMRI.ACG.data.long$diagnosis <- adj.qMRI.ACG.data.dems.melt$diagnosis[1:nrow(adj.qMRI.ACG.data.long)]

# plot for controls & PD
adj.qMRI.ACG.data.long.Ctrl <- subset(adj.qMRI.ACG.data.long, diagnosis == "Control")
adj.qMRI.ACG.data.long.PD <- subset(adj.qMRI.ACG.data.long, diagnosis == "PD")


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

# run again without SWM

adj.qMRI.ACG.data.long.no.SWM <- subset(adj.qMRI.ACG.data.long, SpD != "SWM")

# plot for controls & PD
adj.qMRI.ACG.data.long.Ctrl.no.SWM <- subset(adj.qMRI.ACG.data.long.no.SWM, diagnosis == "Control")
adj.qMRI.ACG.data.long.PD.no.SWM <- subset(adj.qMRI.ACG.data.long.no.SWM, diagnosis == "PD")

dom <- c("SpD4", "SpD3","SpD2","SpD1")
# colours for spatial domains
SpD.cols <- c("#66CCEE","#228833","#CCBB44","#EE6677")

# qMRI boxplots for Ctrl
qMRI.ACG.boxplotlist.Ctrl.no.SWM <- lapply(1:length(MOI), function(i) 
  qMRI.boxplot(metric = MOI[i], adj.qMRI.ACG.data.long.Ctrl.no.SWM))

png(filename = paste(results.dir.qMRIonly, "/Control_qMRI_ACG_boxplots_no_SWM.png", sep = ""),
    width = 1400, height = 600, units = "px", res = 120)
qMRI.ACG.multiboxplot.Ctrl.no.SWM <- ggarrange(plotlist = qMRI.ACG.boxplotlist.Ctrl.no.SWM, 
                                           nrow = 1, ncol = 5, 
                                           common.legend = TRUE, 
                                           legend = "right")
plot(qMRI.ACG.multiboxplot.Ctrl.no.SWM)
dev.off()

# qMRI boxplots for PD
qMRI.ACG.boxplotlist.PD.no.SWM <- lapply(1:length(MOI), function(i) 
  qMRI.boxplot(metric = MOI[i], adj.qMRI.ACG.data.long.PD.no.SWM))

png(filename = paste(results.dir.qMRIonly, "/PD_qMRI_ACG_boxplots_no_SWM.png", sep = ""),
    width = 1400, height = 600, units = "px", res = 120)
qMRI.ACG.multiboxplot.PD.no.SWM <- ggarrange(plotlist = qMRI.ACG.boxplotlist.PD.no.SWM, 
                                         nrow = 1, ncol = 5, 
                                         common.legend = TRUE, 
                                         legend = "right")
plot(qMRI.ACG.multiboxplot.PD.no.SWM)
dev.off()


# remove info about boxplot inclusion
adj.qMRI.ACG.data.dems<- subset(adj.qMRI.ACG.data.dems, 
                            select = -c(inc_QSM, 
                                        inc_R1,
                                        inc_R2s,
                                        inc_PrD,
                                        inc_MTsat))

#### load in & reorganise cell data ####

# load in ST patient dems
ST.dems <- read.csv(paste(dems.dir,"ST_patient_info.csv",sep = ""))
ST.dems <- ST.dems[order(ST.dems$sample_id),]
ST.dems$sex <- factor(ST.dems$sex)

# specify the cell types we're interested in
COI <- c("Astro", "Micro", "Oligo", "OPC", "Inhib")

# get names of 4 layer cell count files
all.cell.fnames <- dir(cell.data.dir,"count_(.+)_4SpD")

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
  cell.data <- cell.data[,c(paste(COI[i], "_WM",sep = ""),
                            paste(COI[i], "_SpD4",sep = ""),
                            paste(COI[i], "_SpD3",sep = ""),
                            paste(COI[i], "_SpD2",sep = ""),
                            paste(COI[i], "_SpD1",sep = ""))]
  cell.data.dems <- cbind(cell.data.dems,cell.data)
}
colnames(cell.data.dems) <- gsub("_WM", "_SWM", colnames(cell.data.dems))
#colnames(cell.data.dems) <- lapply(colnames(cell.data.dems), rename_SpDs)

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

#### calculate average across subjects

# making columns
dom <- c("SWM", "SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

adj.qMRI.ST.medians <- expand.grid(SpD = dom, region = reg)
adj.qMRI.ACG.medians <- expand.grid(SpD = dom, region = reg)

# go through each qMRI metric and remove QC-failed subjects / nan values before averaging
qMRI.metrics = c("QSM","R1","R2s","PrD","MTsat")

# do this for controls and PD separately
# (although only looking at Ctrl vs Ctrl for now)
adj.qMRI.data.dems.Ctrl <- subset(adj.qMRI.data.dems, diagnosis == "Control")
adj.qMRI.ACG.data.dems.Ctrl <- subset(adj.qMRI.ACG.data.dems, diagnosis == "Control")
QC.info.Ctrl <- subset(QC.info, adj.qMRI.data.dems$diagnosis == "Control" )
adj.qMRI.data.dems.PD <- subset(adj.qMRI.data.dems, diagnosis == "PD")
adj.qMRI.ACG.data.dems.PD <- subset(adj.qMRI.ACG.data.dems, diagnosis == "PD")
QC.info.PD <- subset(QC.info, adj.qMRI.data.dems$diagnosis == "PD" )

for (m in qMRI.metrics){
  # remove QC failed subs before calculating mean
  # DLPFC
  metric.data.Ctrl <- adj.qMRI.data.dems.Ctrl[QC.info.Ctrl[,paste(m,"_QC",sep = "")], 
                                              grepl(m, colnames(adj.qMRI.data.dems))]
  metric.medians.Ctrl <- colMedians(as.matrix(metric.data.Ctrl),na.rm = TRUE)
  metric.data.PD <- adj.qMRI.data.dems.PD[QC.info.PD[,paste(m,"_QC",sep = "")], 
                                          grepl(m, colnames(adj.qMRI.data.dems))]
  metric.medians.PD <- colMedians(as.matrix(metric.data.PD),na.rm = TRUE)
  adj.qMRI.ST.medians[,paste(m, "_Control", sep = "")] <- metric.medians.Ctrl
  adj.qMRI.ST.medians[,paste(m, "_PD", sep = "")] <- metric.medians.PD
  # ACG
  metric.ACG.data.Ctrl <- adj.qMRI.ACG.data.dems.Ctrl[QC.info.Ctrl[,paste(m,"_QC",sep = "")], 
                                              grepl(m, colnames(adj.qMRI.ACG.data.dems))]
  metric.ACG.medians.Ctrl <- colMedians(as.matrix(metric.ACG.data.Ctrl),na.rm = TRUE)
  metric.ACG.data.PD <- adj.qMRI.ACG.data.dems.PD[QC.info.PD[,paste(m,"_QC",sep = "")], 
                                          grepl(m, colnames(adj.qMRI.ACG.data.dems))]
  metric.ACG.medians.PD <- colMedians(as.matrix(metric.ACG.data.PD),na.rm = TRUE)
  adj.qMRI.ACG.medians[,paste(m, "_Control", sep = "")] <- metric.ACG.medians.Ctrl
  adj.qMRI.ACG.medians[,paste(m, "_PD", sep = "")] <- metric.ACG.medians.PD
}

# removing neuronal_death and mineral (for now) as these were not in the main analysis
adj.GO.data.dems <- adj.GO.data.dems[,
                                     !(grepl("neuron_death", colnames(adj.GO.data.dems)) | 
                                         grepl("mineral",colnames(adj.GO.data.dems)))]
GO.names <- GO.names[!(grepl("neuron_death", GO.names) | 
                         grepl("mineral",GO.names))]

# add in median cell values
for (c in COI){
  cell.data <- adj.cell.data.dems[,grepl(paste(c,"_",sep = ""), colnames(adj.cell.data.dems))]
  cell.medians <- colMedians(as.matrix(cell.data), na.rm = TRUE)
  adj.qMRI.ST.medians[,paste(c, "_Control", sep = "")] <- cell.medians
}

# add in median GO values
for (p in GO.names){
  GO.data <- adj.GO.data.dems[,grepl(paste(p,"_",sep = ""), colnames(adj.GO.data.dems))]
  GO.medians <- colMedians(as.matrix(GO.data), na.rm = TRUE)
  adj.qMRI.ST.medians[,paste(p, "_Control", sep = "")] <- GO.medians
}

write.csv(adj.qMRI.ST.medians, file = paste(
  results.dir2.stats, "/adj_qMRI_ST_medians_extended.csv", sep = ""), 
  row.names = FALSE)

#### predict PD ST scores with SWM ####

all.ST.names <- colnames(adj.qMRI.ST.medians)[13:ncol(adj.qMRI.ST.medians)]

# define a function to project the PD values
project.PD.ST <- function(dframe, ST) {
  
  output.dframe <- data.frame(matrix(NA, nrow = nrow(dframe), ncol = length(ST)))
  colnames(output.dframe) <- gsub("Control", "PD_predicted", ST)
  output.dframe <- cbind(dframe[,1:2], output.dframe)
  
  # loop though ST metrics
  for (ii in 1:length(ST)) {
    
    st <- ST[ii]
    
      
    # define formula
    control.fm <- as.formula(paste0(st, 
                                 "~ QSM_Control + R1_Control + 
                                 R2s_Control + PrD_Control + MTsat_Control"))
    # PD.data
    PD.data <- dframe[,c("QSM_PD","R1_PD", "R2s_PD","PrD_PD", "MTsat_PD")]
      
    # fit control
    control.fit <- lm(control.fm, data = dframe)
      
    # predict PD values
    # need to rename columns to get this to work
    colnames(PD.data) <- c("QSM_Control","R1_Control", "R2s_Control", "PrD_Control", "MTsat_Control")
    PD.pred <- predict(control.fit, newdata = PD.data)
    
    # add to data frame
    output.dframe[,gsub("Control", "PD_predicted", st)] <- PD.pred
    
  }
  
  return(output.dframe)
}

project.ACG.ST <- function(DLPFC, ACG, ST) {
  
  # want to output two dataframes, one for control ACG, and one for PD ACG
  # Control
  output.dframe1 <- data.frame(matrix(NA, nrow = nrow(ACG), ncol = length(ST)))
  colnames(output.dframe1) <- gsub("Control", "Control_predicted", ST)
  output.dframe1 <- cbind(ACG[,1:2], output.dframe1)
  # PD
  output.dframe2 <- data.frame(matrix(NA, nrow = nrow(ACG), ncol = length(ST)))
  colnames(output.dframe2) <- gsub("Control", "PD_predicted", ST)
  output.dframe2 <- cbind(ACG[,1:2], output.dframe2)
  
  # loop though ST metrics
  for (ii in 1:length(ST)) {
    
    st <- ST[ii]
    
    
    # define formula
    control.fm <- as.formula(paste0(st, 
                                    "~ QSM_Control + R1_Control + 
                                 R2s_Control + PrD_Control + MTsat_Control"))
    
    # control ACG data
    control.ACG.data <- ACG[,c("QSM_Control", "R1_Control", "R2s_Control", 
                               "PrD_Control", "MTsat_Control")]
    # PD ACG data
    PD.ACG.data <- ACG[,c("QSM_PD","R1_PD", "R2s_PD","PrD_PD", "MTsat_PD")]
    
    # fit control
    control.fit <- lm(control.fm, data = DLPFC)
    
    # predict PD values
    # need to temporarily rename columns to get this to work
    colnames(control.ACG.data) <- c("QSM_Control","R1_Control", "R2s_Control", "PrD_Control", "MTsat_Control")
    colnames(PD.ACG.data) <- c("QSM_Control","R1_Control", "R2s_Control", "PrD_Control", "MTsat_Control")
    control.ACG.pred <- predict(control.fit, newdata = control.ACG.data)
    PD.ACG.pred <- predict(control.fit, newdata = PD.ACG.data)
    
    # add to data frame
    output.dframe1[,gsub("Control", "Control_predicted", st)] <- control.ACG.pred
    output.dframe2[,gsub("Control", "PD_predicted", st)] <- PD.ACG.pred
    
  }
  
  return(list(
    control = output.dframe1,
    PD = output.dframe2
    ))
}

# PD DLPFC projections
PD.ST.proj <- project.PD.ST(adj.qMRI.ST.medians, all.ST.names)

PD.cell.proj <- PD.ST.proj[,1:6]
PD.GO.proj <- cbind(PD.ST.proj[,1:2], PD.ST.proj[,7:ncol(PD.ST.proj)])

write.csv(PD.cell.proj, file = paste(
  results.dir.stats, "/Projected_PD_cell_counts.csv", sep = ""), 
  row.names = FALSE)

write.csv(PD.GO.proj, file = paste(
  results.dir.stats, "/Projected_PD_pathway_scores.csv", sep = ""), 
  row.names = FALSE)

#### predict PD ST scores without SWM ####

adj.qMRI.ST.medians.no.SWM <- subset(adj.qMRI.ST.medians, SpD != "SWM")
adj.qMRI.ACG.medians.no.SWM <- subset(adj.qMRI.ACG.medians, SpD != "SWM")

PD.ST.proj.no.SWM <- project.PD.ST(adj.qMRI.ST.medians.no.SWM, all.ST.names)

PD.cell.proj.no.SWM <- PD.ST.proj.no.SWM[,1:6]
PD.GO.proj.no.SWM <- cbind(PD.ST.proj.no.SWM[,1:2], PD.ST.proj.no.SWM[,7:ncol(PD.ST.proj.no.SWM)])

write.csv(PD.cell.proj.no.SWM, file = paste(
  results.dir.stats, "/Projected_PD_cell_counts_no_SWM.csv", sep = ""), 
  row.names = FALSE)

write.csv(PD.GO.proj.no.SWM, file = paste(
  results.dir.stats, "/Projected_PD_pathway_scores_no_SWM.csv", sep = ""), 
  row.names = FALSE)

ACG.proj.no.SWM <- project.ACG.ST(adj.qMRI.ST.medians.no.SWM,
                                  adj.qMRI.ACG.medians.no.SWM, all.ST.names)

control.ACG.cell.proj.no.SWM <- ACG.proj.no.SWM$control[,1:6]
control.ACG.ST.proj.no.SWM <- cbind(ACG.proj.no.SWM$control[,1:2],
                                    ACG.proj.no.SWM$control[,7:ncol(ACG.proj.no.SWM$control)])
PD.ACG.cell.proj.no.SWM <- ACG.proj.no.SWM$PD[,1:6]
PD.ACG.ST.proj.no.SWM <- cbind(ACG.proj.no.SWM$PD[,1:2],
                               ACG.proj.no.SWM$PD[,7:ncol(ACG.proj.no.SWM$PD)])

write.csv(control.ACG.cell.proj.no.SWM, file = paste(
  results.dir.stats, "/Projected_control_ACG_cell_counts_no_SWM.csv", sep = ""), 
  row.names = FALSE)

write.csv(control.ACG.ST.proj.no.SWM, file = paste(
  results.dir.stats, "/Projected_control_ACG_pathway_scores_no_SWM.csv", sep = ""), 
  row.names = FALSE)

write.csv(PD.ACG.cell.proj.no.SWM, file = paste(
  results.dir.stats, "/Projected_PD_ACG_cell_counts_no_SWM.csv", sep = ""), 
  row.names = FALSE)

write.csv(PD.ACG.ST.proj.no.SWM, file = paste(
  results.dir.stats, "/Projected_PD_ACG_pathway_scores_no_SWM.csv", sep = ""), 
  row.names = FALSE)




