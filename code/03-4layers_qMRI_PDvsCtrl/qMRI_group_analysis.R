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

group.cols <- c("#C3B3A1","#4B8D86")

## PD AND CONTROLS ----

# QSM
print("#################### QSM ####################")

# remove subjects that failed QSM QC or are missing their QSM
qMRI.data.dems.QSM <- subset(qMRI.data.dems,QC.info$QSM_QC == 1 & !is.nan(qMRI.data.dems$QSM_SWM))
qMRI.data.dems.QSM <- subset(qMRI.data.dems.QSM, select = c(grepl("subject|age|sex|diagnosis|QSM",colnames(qMRI.data.dems.QSM))))
ROI.names.QSM <- colnames(qMRI.data.dems.QSM)[grepl("QSM",colnames(qMRI.data.dems.QSM))]

qMRI.data.dems.QSM.melt <- melt(qMRI.data.dems.QSM,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.QSM.melt) <- c("subject","age","sex","diagnosis","SpD", "QSM")
qMRI.data.dems.QSM.melt$SpD <- gsub("QSM_", "", qMRI.data.dems.QSM.melt$SpD)
qMRI.data.dems.QSM.melt$SpD <- factor(qMRI.data.dems.QSM.melt$SpD, levels = c("SWM", "SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.QSM <- aov(data = qMRI.data.dems.QSM.melt, QSM ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.QSM.emns <- emmeans(mdl.QSM, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.QSM.emns))

# save some stats
write.csv(as.data.frame(summary(mdl.QSM.emns)),
          paste(results.dir.stats,"/stats_QSM_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_QSM_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.QSM.melt, aes(x = SpD, y = QSM, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("QSM") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()

# R1
print("#################### R1 ####################")

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
mdl.R1.emns <- emmeans(mdl.R1, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.R1.emns))

# save some stats
write.csv(as.data.frame(summary(mdl.R1.emns)),
          paste(results.dir.stats,"/stats_R1_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_R1_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.R1.melt, aes(x = SpD, y = R1, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("R1") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()

# R2*
print("#################### R2s ####################")

# remove subjects that failed R2s QC or are missing their R2s
qMRI.data.dems.R2s <- subset(qMRI.data.dems,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems$R2s_SWM))
qMRI.data.dems.R2s <- subset(qMRI.data.dems.R2s, select = c(grepl("subject|age|sex|diagnosis|R2s",colnames(qMRI.data.dems.R2s))))
ROI.names.R2s <- colnames(qMRI.data.dems.R2s)[grepl("R2s",colnames(qMRI.data.dems.R2s))]

qMRI.data.dems.R2s.melt <- melt(qMRI.data.dems.R2s,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.R2s.melt) <- c("subject","age","sex","diagnosis","SpD", "R2s")
qMRI.data.dems.R2s.melt$SpD <- gsub("R2s_", "", qMRI.data.dems.R2s.melt$SpD)
qMRI.data.dems.R2s.melt$SpD <- factor(qMRI.data.dems.R2s.melt$SpD, levels = c("SWM", "SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.R2s <- aov(data = qMRI.data.dems.R2s.melt, R2s ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.R2s.emns <- emmeans(mdl.R2s, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.R2s.emns))

# save some stats
write.csv(as.data.frame(summary(mdl.R2s.emns)),
          paste(results.dir.stats,"/stats_R2s_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_R2s_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.R2s.melt, aes(x = SpD, y = R2s, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("R2s") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()

# PrD
print("#################### PrD ####################")

# remove subjects that failed PrD QC or are missing their PrD
qMRI.data.dems.PrD <- subset(qMRI.data.dems,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems$PrD_SWM))
qMRI.data.dems.PrD <- subset(qMRI.data.dems.PrD, select = c(grepl("subject|age|sex|diagnosis|PrD",colnames(qMRI.data.dems.PrD))))
ROI.names.PrD <- colnames(qMRI.data.dems.PrD)[grepl("PrD",colnames(qMRI.data.dems.PrD))]

qMRI.data.dems.PrD.melt <- melt(qMRI.data.dems.PrD,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.PrD.melt) <- c("subject","age","sex","diagnosis","SpD", "PrD")
qMRI.data.dems.PrD.melt$SpD <- gsub("PrD_", "", qMRI.data.dems.PrD.melt$SpD)
qMRI.data.dems.PrD.melt$SpD <- factor(qMRI.data.dems.PrD.melt$SpD, levels = c("SWM", "SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.PrD <- aov(data = qMRI.data.dems.PrD.melt, PrD ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.PrD.emns <- emmeans(mdl.PrD, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.PrD.emns))

# save some stats
write.csv(as.data.frame(summary(mdl.PrD.emns)),
          paste(results.dir.stats,"/stats_PrD_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_PrD_group_comparison_PDvsCtrl",".png",sep=""),
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.PrD.melt, aes(x = SpD, y = PrD, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) +
  ylab("PrD") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()

# MT sat
print("#################### MTsat ####################")

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
mdl.MTsat.emns <- emmeans(mdl.MTsat, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.MTsat.emns))

# save some stats
write.csv(as.data.frame(summary(mdl.MTsat.emns)),
          paste(results.dir.stats,"/stats_MTsat_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_MTsat_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.MTsat.melt, aes(x = SpD, y = MTsat, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("MTsat") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()


## PD AND CONTROLS (without SWM) ----

qMRI.data.dems.no.SWM <- qMRI.data.dems[ , !grepl("_SWM", names(qMRI.data.dems))]

# QSM
print("#################### QSM ####################")

# remove subjects that failed QSM QC or are missing their QSM
qMRI.data.dems.QSM.no.SWM <- subset(qMRI.data.dems.no.SWM,QC.info$QSM_QC == 1 & !is.nan(qMRI.data.dems.no.SWM$QSM_SpD4))
qMRI.data.dems.QSM.no.SWM <- subset(qMRI.data.dems.QSM.no.SWM, select = c(grepl("subject|age|sex|diagnosis|QSM",colnames(qMRI.data.dems.QSM.no.SWM))))
ROI.names.QSM.no.SWM <- colnames(qMRI.data.dems.QSM.no.SWM)[grepl("QSM",colnames(qMRI.data.dems.QSM.no.SWM))]

qMRI.data.dems.QSM.no.SWM.melt <- melt(qMRI.data.dems.QSM.no.SWM,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.QSM.no.SWM.melt) <- c("subject","age","sex","diagnosis","SpD", "QSM")
qMRI.data.dems.QSM.no.SWM.melt$SpD <- gsub("QSM_", "", qMRI.data.dems.QSM.no.SWM.melt$SpD)
qMRI.data.dems.QSM.no.SWM.melt$SpD <- factor(qMRI.data.dems.QSM.no.SWM.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.QSM.no.SWM<- aov(data = qMRI.data.dems.QSM.no.SWM.melt, QSM ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.QSM.emns.no.SWM <- emmeans(mdl.QSM.no.SWM, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.QSM.emns.no.SWM))

# save some stats
write.csv(as.data.frame(summary(mdl.QSM.emns.no.SWM)),
          paste(results.dir.stats,"/stats_QSM_group_comparison_PDvsCtrl_no_SWM",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_QSM_group_comparison_PDvsCtrl_no_SWM",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.QSM.no.SWM.melt, aes(x = SpD, y = QSM, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("QSM") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()

# R1
print("#################### R1 ####################")

# remove subjects that failed R1 QC or are missing their R1
qMRI.data.dems.R1.no.SWM <- subset(qMRI.data.dems.no.SWM,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems.no.SWM$R1_SpD4))
qMRI.data.dems.R1.no.SWM <- subset(qMRI.data.dems.R1.no.SWM, select = c(grepl("subject|age|sex|diagnosis|R1",colnames(qMRI.data.dems.R1.no.SWM))))
ROI.names.R1.no.SWM <- colnames(qMRI.data.dems.R1.no.SWM)[grepl("R1",colnames(qMRI.data.dems.R1.no.SWM))]

qMRI.data.dems.R1.no.SWM.melt <- melt(qMRI.data.dems.R1.no.SWM,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.R1.no.SWM.melt) <- c("subject","age","sex","diagnosis","SpD", "R1")
qMRI.data.dems.R1.no.SWM.melt$SpD <- gsub("R1_", "", qMRI.data.dems.R1.no.SWM.melt$SpD)
qMRI.data.dems.R1.no.SWM.melt$SpD <- factor(qMRI.data.dems.R1.no.SWM.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.R1.no.SWM <- aov(data = qMRI.data.dems.R1.no.SWM.melt, R1 ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.R1.no.SWM.emns.no.SWM <- emmeans(mdl.R1.no.SWM, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.R1.no.SWM.emns.no.SWM))

# save some stats
write.csv(as.data.frame(summary(mdl.R1.no.SWM.emns.no.SWM)),
          paste(results.dir.stats,"/stats_R1_group_comparison_PDvsCtrl_no_SWM",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_R1_group_comparison_PDvsCtrl_no_SWM",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.R1.no.SWM.melt, aes(x = SpD, y = R1, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("R1") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()

# R2*
print("#################### R2s ####################")

# remove subjects that failed R2s QC or are missing their R2s
qMRI.data.dems.R2s.no.SWM <- subset(qMRI.data.dems.no.SWM,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems.no.SWM$R2s_SpD4))
qMRI.data.dems.R2s.no.SWM <- subset(qMRI.data.dems.R2s.no.SWM, select = c(grepl("subject|age|sex|diagnosis|R2s",colnames(qMRI.data.dems.R2s.no.SWM))))
ROI.names.R2s.no.SWM <- colnames(qMRI.data.dems.R2s.no.SWM)[grepl("R2s",colnames(qMRI.data.dems.R2s.no.SWM))]

qMRI.data.dems.R2s.no.SWM.melt <- melt(qMRI.data.dems.R2s.no.SWM,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.R2s.no.SWM.melt) <- c("subject","age","sex","diagnosis","SpD", "R2s")
qMRI.data.dems.R2s.no.SWM.melt$SpD <- gsub("R2s_", "", qMRI.data.dems.R2s.no.SWM.melt$SpD)
qMRI.data.dems.R2s.no.SWM.melt$SpD <- factor(qMRI.data.dems.R2s.no.SWM.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.R2s.no.SWM <- aov(data = qMRI.data.dems.R2s.no.SWM.melt, R2s ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.R2s.no.SWM.emns.no.SWM <- emmeans(mdl.R2s.no.SWM, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.R2s.no.SWM.emns.no.SWM))

# save some stats
write.csv(as.data.frame(summary(mdl.R2s.no.SWM.emns.no.SWM)),
          paste(results.dir.stats,"/stats_R2s_group_comparison_PDvsCtrl_no_SWM",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_R2s_group_comparison_PDvsCtrl_no_SWM",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.R2s.no.SWM.melt, aes(x = SpD, y = R2s, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("R2s") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()

# PrD
print("#################### PrD ####################")

# remove subjects that failed PrD QC or are missing their PrD
qMRI.data.dems.PrD.no.SWM <- subset(qMRI.data.dems.no.SWM,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems.no.SWM$PrD_SpD4))
qMRI.data.dems.PrD.no.SWM <- subset(qMRI.data.dems.PrD.no.SWM, select = c(grepl("subject|age|sex|diagnosis|PrD",colnames(qMRI.data.dems.PrD.no.SWM))))
ROI.names.PrD.no.SWM <- colnames(qMRI.data.dems.PrD.no.SWM)[grepl("PrD",colnames(qMRI.data.dems.PrD.no.SWM))]

qMRI.data.dems.PrD.no.SWM.melt <- melt(qMRI.data.dems.PrD.no.SWM,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.PrD.no.SWM.melt) <- c("subject","age","sex","diagnosis","SpD", "PrD")
qMRI.data.dems.PrD.no.SWM.melt$SpD <- gsub("PrD_", "", qMRI.data.dems.PrD.no.SWM.melt$SpD)
qMRI.data.dems.PrD.no.SWM.melt$SpD <- factor(qMRI.data.dems.PrD.no.SWM.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.PrD.no.SWM <- aov(data = qMRI.data.dems.PrD.no.SWM.melt, PrD ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.PrD.no.SWM.emns.no.SWM <- emmeans(mdl.PrD.no.SWM, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.PrD.no.SWM.emns.no.SWM))

# save some stats
write.csv(as.data.frame(summary(mdl.PrD.no.SWM.emns.no.SWM)),
          paste(results.dir.stats,"/stats_PrD_group_comparison_PDvsCtrl_no_SWM",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_PrD_group_comparison_PDvsCtrl_no_SWM",".png",sep=""),
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.PrD.no.SWM.melt, aes(x = SpD, y = PrD, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) +
  ylab("PrD") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()

# MT sat
print("#################### MTsat ####################")

# remove subjects that failed MTsat QC or are missing their MTsat
qMRI.data.dems.MTsat.no.SWM <- subset(qMRI.data.dems.no.SWM,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems.no.SWM$MTsat_SpD4))
qMRI.data.dems.MTsat.no.SWM <- subset(qMRI.data.dems.MTsat.no.SWM, select = c(grepl("subject|age|sex|diagnosis|MTsat",colnames(qMRI.data.dems.MTsat.no.SWM))))
ROI.names.MTsat.no.SWM <- colnames(qMRI.data.dems.MTsat.no.SWM)[grepl("MTsat",colnames(qMRI.data.dems.MTsat.no.SWM))]

qMRI.data.dems.MTsat.no.SWM.melt <- melt(qMRI.data.dems.MTsat.no.SWM,id.vars=c("subject","age","sex","diagnosis"))
colnames(qMRI.data.dems.MTsat.no.SWM.melt) <- c("subject","age","sex","diagnosis","SpD", "MTsat")
qMRI.data.dems.MTsat.no.SWM.melt$SpD <- gsub("MTsat_", "", qMRI.data.dems.MTsat.no.SWM.melt$SpD)
qMRI.data.dems.MTsat.no.SWM.melt$SpD <- factor(qMRI.data.dems.MTsat.no.SWM.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.MTsat.no.SWM <- aov(data = qMRI.data.dems.MTsat.no.SWM.melt, MTsat ~ age + sex + SpD*diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.MTsat.no.SWM.emns.no.SWM <- emmeans(mdl.MTsat.no.SWM, pairwise ~ diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.MTsat.no.SWM.emns.no.SWM))

# save some stats
write.csv(as.data.frame(summary(mdl.MTsat.no.SWM.emns.no.SWM)),
          paste(results.dir.stats,"/stats_MTsat_group_comparison_PDvsCtrl_no_SWM",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_MTsat_group_comparison_PDvsCtrl_no_SWM",".png",sep=""), 
    width = 750, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.MTsat.no.SWM.melt, aes(x = SpD, y = MTsat, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("MTsat") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal()
print(plot.mdl)
dev.off()
