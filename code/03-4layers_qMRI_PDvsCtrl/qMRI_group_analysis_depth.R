# script to do some qMRI layer analysis

# load in some libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(emmeans)
library(lme4)

# define directories
d0 <- getwd();
depth.data.dir <- paste(d0, "/../../raw_data/4layers_qMRI/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir <- paste(d0, "/../../plots/03-4layers_PDvsCtrl/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/03-4layers_PDvsCtrl/", sep = "")

# load in the qMRI data and demographics for the combined DLPFC
qMRI.dems <- read.csv(paste(dems.dir,"anon_DLPFC_qMRI_demographics.csv",sep = ""))
depth.data <-  read.csv(paste(depth.data.dir,"anon_DLPFC_mean_depth_4_equivol_layers.csv",sep = ""))
QC.info <- read.csv(paste(dems.dir,"anon_7T_QC_info.csv",sep = ""))

# remove data from demographics and qMRI data using QC info
# if they failed MPRAGE QC or FS QC, then remove from everywhere
qMRI.dems <- subset(qMRI.dems,QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
depth.data <- subset(depth.data, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
QC.info <- subset(QC.info, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)

# merge into one data frame, this will make things easier
depth.data.dems <- cbind(qMRI.dems,depth.data)
depth.data.dems <- subset(depth.data.dems, select = -c(subjects))

# set the group and sex variables as factors
depth.data.dems$sex <- factor(depth.data.dems$sex)
depth.data.dems$diagnosis <- factor(depth.data.dems$diagnosis, levels=c("Control","PD"))

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

# rename columns for to match spatial domains
colnames(depth.data.dems) <- lapply(colnames(depth.data.dems), rename_layers)

# now do some group analyses

group.cols <- c("#C3B3A1","#4B8D86")

## PD AND CONTROLS (without SWM) ----

depth.data.dems.no.SWM <- depth.data.dems[ , !grepl("_SWM", names(depth.data.dems))]

# depth
print("#################### depth ####################")

depth.data.dems.no.SWM <- subset(depth.data.dems.no.SWM, select = c(grepl("subject|age|sex|diagnosis|DLPFC_CT|max_PFD|depth",colnames(depth.data.dems.no.SWM))))
ROI.names.depth.no.SWM <- colnames(depth.data.dems.no.SWM)[grepl("depth",colnames(depth.data.dems.no.SWM))]

depth.data.dems.no.SWM.melt <- melt(depth.data.dems.no.SWM,id.vars=c("subject","age","sex","diagnosis","DLPFC_CT","max_PFD"))
colnames(depth.data.dems.no.SWM.melt) <- c("subject","age","sex","diagnosis","DLPFC_CT","max_PFD","SpD", "depth")
depth.data.dems.no.SWM.melt$SpD <- gsub("depth_", "", depth.data.dems.no.SWM.melt$SpD)
depth.data.dems.no.SWM.melt$SpD <- factor(depth.data.dems.no.SWM.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.depth.no.SWM<- lmer(data = depth.data.dems.no.SWM.melt, depth ~ age + sex + DLPFC_CT + max_PFD + SpD*diagnosis + (1 | subject))
# calculate estimated marginal means (adjusting for age and sex)
mdl.depth.emns.no.SWM <- emmeans(mdl.depth.no.SWM, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.depth.emns.no.SWM))

# save some stats
write.csv(as.data.frame(summary(mdl.depth.emns.no.SWM)),
          paste(results.dir.stats,"/stats_depth_group_comparison_PDvsCtrl_no_SWM",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_depth_group_comparison_PDvsCtrl_no_SWM",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(depth.data.dems.no.SWM.melt, aes(x = SpD, y = depth, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("depth (mm)") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()