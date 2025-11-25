# script to do some qMRI layer analysis

# load in some libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(emmeans)

# define directories
d0 <- getwd();
qMRI.data.dir <- paste(d0, "/../../processed_data/02-4layers_qMRI/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir <- paste(d0, "/../../plots/03-4layers_PDvsCtrl/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/03-4layers_PDvsCtrl/", sep = "")

# load in the qMRI data and demographics for the combined DLPFC
qMRI.dems <- read.csv(paste(dems.dir,"anon_DLPFC_qMRI_demographics.csv",sep = ""))
qMRI.data <-  read.csv(paste(qMRI.data.dir,"anon_DLPFC_qMRI_medians_4_equivol_layers.csv",sep = ""))
QC.info <- read.csv(paste(dems.dir,"anon_7T_QC_info.csv",sep = ""))

# remove data from demographics and qMRI data using QC info
# if they failed MPRAGE QC or FS QC, then remove from everywhere
qMRI.dems <- subset(qMRI.dems,QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
qMRI.data <- subset(qMRI.data, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
QC.info <- subset(QC.info, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)

# merge into one data frame, this will make things easier
qMRI.data.dems <- cbind(qMRI.dems,qMRI.data)
qMRI.data.dems <- subset(qMRI.data.dems, select = -c(subjects))

# set the group and sex variables as factors
qMRI.data.dems$sex <- factor(qMRI.data.dems$sex)
qMRI.data.dems$diagnosis <- factor(qMRI.data.dems$diagnosis, levels=c("Control","PD","DLB"))

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

# rename columns
colnames(qMRI.data.dems) <- lapply(colnames(qMRI.data.dems), rename_layers)
# get all ROI names
qMRI.ROI.names <- colnames(qMRI.data.dems)[5:ncol(qMRI.data.dems)]

# now do some group analyses

## JUST PD AND CONTROLS ----

###############################################################  R1

# remove subjects that failed R1 QC or are missing their R1
qMRI.data.dems.R1 <- subset(qMRI.data.dems,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems$R1_SWM))
qMRI.data.dems.R1 <- subset(qMRI.data.dems.R1, select = c(grepl("subject|age|sex|diagnosis|R1",colnames(qMRI.data.dems.R1))))
ROI.names.R1 <- colnames(qMRI.data.dems.R1)[grepl("R1",colnames(qMRI.data.dems.R1))]

qMRI.data.dems.R1.melt <- melt(qMRI.data.dems.R1,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.R1.melt) <- c("subject","age","sex","diagnosis","SpD", "R1")
qMRI.data.dems.R1.melt$SpD <- gsub("R1_", "", qMRI.data.dems.R1.melt$SpD)
qMRI.data.dems.R1.melt$SpD <- factor(qMRI.data.dems.R1.melt$SpD, levels = c("SWM", "SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.R1 <- aov(data = qMRI.data.dems.R1.melt, R1 ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.R1.emns <- emmeans(mdl.R1, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "bonferroni")

print("#################### R1 ####################")

print(summary(mdl.R1.emns))

qMRI.data.dems.R1.melt = qMRI.data.dems.R1.melt[qMRI.data.dems.R1.melt$SpD!='SWM',]

# save some stats
write.csv(as.data.frame(summary(mdl.R1.emns)),
          paste(results.dir.stats,"/stats_R1_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_R1_group_comparison_PDvsCtrl_MGP",".png",sep=""), 
    width = 750, height = 375, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.R1.melt, aes(x = SpD, y = R1, fill = diagnosis)) + 
  scale_fill_manual(values = c("PD" = "#4B8D86", "Control" = "#C3B3A1")) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("R1") +
  theme_minimal()
print(plot.mdl)
dev.off()


############################################################### MT sat

# remove subjects that failed MTsat QC or are missing their MTsat
qMRI.data.dems.MTsat <- subset(qMRI.data.dems,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems$MTsat_SWM))
qMRI.data.dems.MTsat <- subset(qMRI.data.dems.MTsat, select = c(grepl("subject|age|sex|diagnosis|MTsat",colnames(qMRI.data.dems.MTsat))))
ROI.names.MTsat <- colnames(qMRI.data.dems.MTsat)[grepl("MTsat",colnames(qMRI.data.dems.MTsat))]

qMRI.data.dems.MTsat.melt <- melt(qMRI.data.dems.MTsat,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.MTsat.melt) <- c("subject","age","sex","diagnosis","SpD", "MTsat")
qMRI.data.dems.MTsat.melt$SpD <- gsub("MTsat_", "", qMRI.data.dems.MTsat.melt$SpD)
qMRI.data.dems.MTsat.melt$SpD <- factor(qMRI.data.dems.MTsat.melt$SpD, levels = c("SWM", "SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.MTsat <- aov(data = qMRI.data.dems.MTsat.melt, MTsat ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.MTsat.emns <- emmeans(mdl.MTsat, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "bonferroni")

print("#################### MTsat ####################")

print(summary(mdl.MTsat.emns))

qMRI.data.dems.MTsat.melt = qMRI.data.dems.MTsat.melt[qMRI.data.dems.MTsat.melt$SpD!='SWM',]

# save some stats
#write.csv(as.data.frame(summary(mdl.MTsat.emns)),
#          paste(results.dir.stats,"/stats_MTsat_group_comparison_PDvsCtrl",".csv",sep=""),
#          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_MTsat_group_comparison_PDvsCtrl_MGP",".png",sep=""), 
    width = 750, height = 375, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.MTsat.melt, aes(x = SpD, y = MTsat, fill = diagnosis)) + 
  scale_fill_manual(values = c("PD" = "#4B8D86", "Control" = "#C3B3A1")) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("MTsat") +
  theme_minimal()
print(plot.mdl)
dev.off()


###############################################################  QSM

# remove subjects that failed MTsat QC or are missing their MTsat
qMRI.data.dems.MTsat <- subset(qMRI.data.dems,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems$MTsat_SWM))
qMRI.data.dems.MTsat <- subset(qMRI.data.dems.MTsat, select = c(grepl("subject|age|sex|diagnosis|MTsat",colnames(qMRI.data.dems.MTsat))))
ROI.names.MTsat <- colnames(qMRI.data.dems.MTsat)[grepl("MTsat",colnames(qMRI.data.dems.MTsat))]

qMRI.data.dems.MTsat.melt <- melt(qMRI.data.dems.MTsat,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.MTsat.melt) <- c("subject","age","sex","diagnosis","SpD", "MTsat")
qMRI.data.dems.MTsat.melt$SpD <- gsub("MTsat_", "", qMRI.data.dems.MTsat.melt$SpD)
qMRI.data.dems.MTsat.melt$SpD <- factor(qMRI.data.dems.MTsat.melt$SpD, levels = c("SWM", "SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.MTsat <- aov(data = qMRI.data.dems.MTsat.melt, MTsat ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.MTsat.emns <- emmeans(mdl.MTsat, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "bonferroni")

print("#################### MTsat ####################")

print(summary(mdl.MTsat.emns))

qMRI.data.dems.MTsat = qMRI.data.dems.MTsat[qMRI.data.dems.MTsat$SpD!='SWM',]

# save some stats
#write.csv(as.data.frame(summary(mdl.MTsat.emns)),
#          paste(results.dir.stats,"/stats_MTsat_group_comparison_PDvsCtrl",".csv",sep=""),
#          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_MTsat_group_comparison_PDvsCtrl_MGP",".png",sep=""), 
    width = 750, height = 375, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.MTsat.melt, aes(x = SpD, y = MTsat, fill = diagnosis)) + 
  scale_fill_manual(values = c("PD" = "#4B8D86", "Control" = "#C3B3A1")) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("MTsat") +
  theme_minimal()
print(plot.mdl)
dev.off()

###############################################################  PrD

# remove subjects that failed MTsat QC or are missing their MTsat
qMRI.data.dems.MTsat <- subset(qMRI.data.dems,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems$MTsat_SWM))
qMRI.data.dems.MTsat <- subset(qMRI.data.dems.MTsat, select = c(grepl("subject|age|sex|diagnosis|MTsat",colnames(qMRI.data.dems.MTsat))))
ROI.names.MTsat <- colnames(qMRI.data.dems.MTsat)[grepl("MTsat",colnames(qMRI.data.dems.MTsat))]

qMRI.data.dems.MTsat.melt <- melt(qMRI.data.dems.MTsat,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.MTsat.melt) <- c("subject","age","sex","diagnosis","SpD", "MTsat")
qMRI.data.dems.MTsat.melt$SpD <- gsub("MTsat_", "", qMRI.data.dems.MTsat.melt$SpD)
qMRI.data.dems.MTsat.melt$SpD <- factor(qMRI.data.dems.MTsat.melt$SpD, levels = c("SWM", "SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.MTsat <- aov(data = qMRI.data.dems.MTsat.melt, MTsat ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.MTsat.emns <- emmeans(mdl.MTsat, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "bonferroni")

print("#################### MTsat ####################")

print(summary(mdl.MTsat.emns))

qMRI.data.dems.MTsat = qMRI.data.dems.MTsat[qMRI.data.dems.MTsat$SpD!='SWM',]

# save some stats
#write.csv(as.data.frame(summary(mdl.MTsat.emns)),
#          paste(results.dir.stats,"/stats_MTsat_group_comparison_PDvsCtrl",".csv",sep=""),
#          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_MTsat_group_comparison_PDvsCtrl_MGP",".png",sep=""), 
    width = 750, height = 375, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.MTsat.melt, aes(x = SpD, y = MTsat, fill = diagnosis)) + 
  scale_fill_manual(values = c("PD" = "#4B8D86", "Control" = "#C3B3A1")) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("MTsat") +
  theme_minimal()
print(plot.mdl)
dev.off()

###############################################################  R2*

# remove subjects that failed MTsat QC or are missing their MTsat
qMRI.data.dems.MTsat <- subset(qMRI.data.dems,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems$R2s_SWM))
qMRI.data.dems.MTsat <- subset(qMRI.data.dems.MTsat, select = c(grepl("subject|age|sex|diagnosis|MTsat",colnames(qMRI.data.dems.MTsat))))
ROI.names.MTsat <- colnames(qMRI.data.dems.MTsat)[grepl("MTsat",colnames(qMRI.data.dems.MTsat))]

qMRI.data.dems.MTsat.melt <- melt(qMRI.data.dems.MTsat,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.MTsat.melt) <- c("subject","age","sex","diagnosis","SpD", "MTsat")
qMRI.data.dems.MTsat.melt$SpD <- gsub("MTsat_", "", qMRI.data.dems.MTsat.melt$SpD)
qMRI.data.dems.MTsat.melt$SpD <- factor(qMRI.data.dems.MTsat.melt$SpD, levels = c("SWM", "SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.MTsat <- aov(data = qMRI.data.dems.MTsat.melt, MTsat ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.MTsat.emns <- emmeans(mdl.MTsat, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "bonferroni")

print("#################### MTsat ####################")

print(summary(mdl.MTsat.emns))

qMRI.data.dems.MTsat = qMRI.data.dems.MTsat[qMRI.data.dems.MTsat$SpD!='SWM',]

# save some stats
#write.csv(as.data.frame(summary(mdl.MTsat.emns)),
#          paste(results.dir.stats,"/stats_MTsat_group_comparison_PDvsCtrl",".csv",sep=""),
#          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_MTsat_group_comparison_PDvsCtrl_MGP",".png",sep=""), 
    width = 750, height = 375, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.MTsat.melt, aes(x = SpD, y = MTsat, fill = diagnosis)) + 
  scale_fill_manual(values = c("PD" = "#4B8D86", "Control" = "#C3B3A1")) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("MTsat") +
  theme_minimal()
print(plot.mdl)
dev.off()