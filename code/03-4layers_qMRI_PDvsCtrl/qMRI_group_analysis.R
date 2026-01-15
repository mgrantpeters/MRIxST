# script to do some qMRI layer analysis

# load in some libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(emmeans)

# define directories
d0 <- getwd();
qMRI.data.dir <- paste(d0, "/../../raw_data/4layers_qMRI/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir <- paste(d0, "/../../plots/03-4layers_PDvsCtrl/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/03-4layers_PDvsCtrl/", sep = "")

# load in the qMRI data and demographics for the combined DLPFC
qMRI.dems <- read.csv(paste(dems.dir,"anon_DLPFC_qMRI_demographics.csv",sep = ""))
qMRI.data <-  read.csv(paste(qMRI.data.dir,"anon_DLPFC_qMRI_medians_4_equivol_layers.csv",sep = ""))
qMRI.data.whole <- read.csv(paste(qMRI.data.dir,"anon_whole_DLPFC_qMRI_medians.csv",sep = ""))
QC.info <- read.csv(paste(dems.dir,"anon_7T_QC_info.csv",sep = ""))

# remove data from demographics and qMRI data using QC info
# if they failed MPRAGE QC or FS QC, then remove from everywhere
qMRI.dems <- subset(qMRI.dems,QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
qMRI.data <- subset(qMRI.data, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
qMRI.data.whole <- subset(qMRI.data.whole, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)
QC.info <- subset(QC.info, QC.info$MPRAGE_QC == 1 & QC.info$FS_QC == 1)

# merge into one data frame, this will make things easier
qMRI.data.dems <- cbind(qMRI.dems,qMRI.data)
qMRI.data.dems <- subset(qMRI.data.dems, select = -c(subjects))

qMRI.data.dems.whole <- cbind(qMRI.dems,qMRI.data.whole)
qMRI.data.dems.whole <- subset(qMRI.data.dems.whole, select = -c(subjects))

# set the group and sex variables as factors
qMRI.data.dems$sex <- factor(qMRI.data.dems$sex)
qMRI.data.dems$diagnosis <- factor(qMRI.data.dems$diagnosis, levels=c("Control","PD"))

qMRI.data.dems.whole$sex <- factor(qMRI.data.dems.whole$sex)
qMRI.data.dems.whole$diagnosis <- factor(qMRI.data.dems.whole$diagnosis, levels=c("Control","PD"))

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
colnames(qMRI.data.dems) <- lapply(colnames(qMRI.data.dems), rename_layers)

# now do some group analyses

group.cols <- c("#C3B3A1","#4B8D86")

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

# do for the whole DLPFC
qMRI.data.dems.whole.QSM <- subset(qMRI.data.dems.whole,QC.info$QSM_QC == 1 & !is.nan(qMRI.data.dems.whole$QSM))
qMRI.data.dems.whole.QSM <- subset(qMRI.data.dems.whole.QSM, select = c(grepl("subject|age|sex|diagnosis|QSM",colnames(qMRI.data.dems.whole.QSM))))

# simple anova
mdl.QSM.whole <- aov(data = qMRI.data.dems.whole.QSM, QSM ~ age + sex + diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
res.QSM.whole <- anova(mdl.QSM.whole)

print(res.QSM.whole)

# save some stats
write.csv(as.data.frame(res.QSM.whole),
          paste(results.dir.stats,"/stats_whole_DLPFC_QSM_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot
png(filename = paste(results.dir,"/plot_whole_DLPFC_QSM_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 250, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.whole.QSM, aes(x=diagnosis, y = QSM, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("QSM") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal() +
  theme(legend.position = "none") 
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

# do for the whole DLPFC
qMRI.data.dems.whole.R1 <- subset(qMRI.data.dems.whole,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems.whole$R1))
qMRI.data.dems.whole.R1 <- subset(qMRI.data.dems.whole.R1, select = c(grepl("subject|age|sex|diagnosis|R1",colnames(qMRI.data.dems.whole.R1))))

# simple anova
mdl.R1.whole <- aov(data = qMRI.data.dems.whole.R1, R1 ~ age + sex + diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
res.R1.whole <- anova(mdl.R1.whole)

print(res.R1.whole)

# save some stats
write.csv(as.data.frame(res.R1.whole),
          paste(results.dir.stats,"/stats_whole_DLPFC_R1_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot
png(filename = paste(results.dir,"/plot_whole_DLPFC_R1_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 250, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.whole.R1, aes(x=diagnosis, y = R1, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("R1") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal() +
  theme(legend.position = "none") 
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

# do for the whole DLPFC
qMRI.data.dems.whole.R2s <- subset(qMRI.data.dems.whole,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems.whole$R2s))
qMRI.data.dems.whole.R2s <- subset(qMRI.data.dems.whole.R2s, select = c(grepl("subject|age|sex|diagnosis|R2s",colnames(qMRI.data.dems.whole.R2s))))

# simple anova
mdl.R2s.whole <- aov(data = qMRI.data.dems.whole.R2s, R2s ~ age + sex + diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
res.R2s.whole <- anova(mdl.R2s.whole)

print(res.R2s.whole)

# save some stats
write.csv(as.data.frame(res.R2s.whole),
          paste(results.dir.stats,"/stats_whole_DLPFC_R2s_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot
png(filename = paste(results.dir,"/plot_whole_DLPFC_R2s_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 250, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.whole.R2s, aes(x=diagnosis, y = R2s, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("R2s") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal() +
  theme(legend.position = "none") 
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

# do for the whole DLPFC
qMRI.data.dems.whole.PrD <- subset(qMRI.data.dems.whole,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems.whole$PrD))
qMRI.data.dems.whole.PrD <- subset(qMRI.data.dems.whole.PrD, select = c(grepl("subject|age|sex|diagnosis|PrD",colnames(qMRI.data.dems.whole.PrD))))

# simple anova
mdl.PrD.whole <- aov(data = qMRI.data.dems.whole.PrD, PrD ~ age + sex + diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
res.PrD.whole <- anova(mdl.PrD.whole)

print(res.PrD.whole)

# save some stats
write.csv(as.data.frame(res.PrD.whole),
          paste(results.dir.stats,"/stats_whole_DLPFC_PrD_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot
png(filename = paste(results.dir,"/plot_whole_DLPFC_PrD_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 250, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.whole.PrD, aes(x=diagnosis, y = PrD, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("PrD") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal() +
  theme(legend.position = "none") 
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

# do for the whole DLPFC
qMRI.data.dems.whole.MTsat <- subset(qMRI.data.dems.whole,QC.info$MPM_QC == 1 & !is.nan(qMRI.data.dems.whole$MTsat))
qMRI.data.dems.whole.MTsat <- subset(qMRI.data.dems.whole.MTsat, select = c(grepl("subject|age|sex|diagnosis|MTsat",colnames(qMRI.data.dems.whole.MTsat))))

# simple anova
mdl.MTsat.whole <- aov(data = qMRI.data.dems.whole.MTsat, MTsat ~ age + sex + diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
res.MTsat.whole <- anova(mdl.MTsat.whole)

print(res.MTsat.whole)

# save some stats
write.csv(as.data.frame(res.MTsat.whole),
          paste(results.dir.stats,"/stats_whole_DLPFC_MTsat_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot
png(filename = paste(results.dir,"/plot_whole_DLPFC_MTsat_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 250, height = 500, units = "px", res = 150)
plot.mdl <- ggplot(qMRI.data.dems.whole.MTsat, aes(x=diagnosis, y = MTsat, fill = diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = diagnosis)) + 
  ylab("MTsat") +
  scale_fill_manual(values =(group.cols)) +
  theme_minimal() +
  theme(legend.position = "none") 
print(plot.mdl)
dev.off()
