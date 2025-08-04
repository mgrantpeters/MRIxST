# script to do some HALO layer analysis

# load in some libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(emmeans)

# define directories
d0 <- getwd();
HALO.data.dir <- paste(d0, "/../../processed_data/08_validation_oligos/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir <- paste(d0, "/../../plots/08_validation_oligos/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/08_validation_oligos/", sep = "")

# load in the HALO data and demographics
HALO.dems <- read.csv(paste(dems.dir,"oligo_validation_subject_info.csv",sep = ""))
HALO.data <-  read.csv(paste(HALO.data.dir,"SpD_quant_norm_equivol.csv",sep = ""))

# order both by subject
HALO.dems <- HALO.dems[order(HALO.dems$Subject),]
HALO.data <- HALO.data[order(HALO.data$Subject),]

# merge into one data frame, this will make things easier
HALO.data.dems <- cbind(HALO.dems,subset(HALO.data, select = -c(Subject)))

# set the group and sex variables as factors
HALO.data.dems$Sex <- factor(HALO.data.dems$Sex)
HALO.data.dems$Diagnosis <- factor(HALO.data.dems$Diagnosis, levels=c("HC","PD"))

# get all ROI names
HALO.ROI.names <- colnames(HALO.data.dems)[8:ncol(HALO.data.dems)]

# now do some group analyses

## PD vs HC ----

print("#################### Oligo count ####################")

# Oligo count

# select only oligo count columns
HALO.data.dems.OligCount <- subset(HALO.data.dems, select = c(grepl("Subject|AgeDeath|Sex|Diagnosis|OligCount",colnames(HALO.data.dems))))
ROI.names.OligCount <- colnames(HALO.data.dems.OligCount)[grepl("OligCount",colnames(HALO.data.dems.OligCount))]

HALO.data.dems.OligCount.melt <- melt(HALO.data.dems.OligCount,id.vars=c("Subject","AgeDeath","Sex","Diagnosis"))
colnames(HALO.data.dems.OligCount.melt) <- c("Subject","AgeDeath","Sex","Diagnosis","SpD", "OligCount")
HALO.data.dems.OligCount.melt$SpD <- gsub("OligCountNorm_", "", HALO.data.dems.OligCount.melt$SpD)
HALO.data.dems.OligCount.melt$SpD <- factor(HALO.data.dems.OligCount.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.OligCount <- aov(data = HALO.data.dems.OligCount.melt, OligCount ~ AgeDeath + Sex + SpD*Diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.OligCount.emns <- emmeans(mdl.OligCount, pairwise ~ Diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.OligCount.emns))

# save some stats
write.csv(as.data.frame(summary(mdl.OligCount.emns)),
          paste(results.dir.stats,"/stats_OligCount_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_OligCount_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 750, height = 750, units = "px", res = 150)
plot.mdl <- ggplot(HALO.data.dems.OligCount.melt, aes(x = SpD, y = OligCount, fill = Diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = Diagnosis)) + 
  ylab("Oligo. count per pixel (normalised)") +
  theme_minimal()
print(plot.mdl)
dev.off()

print("#################### MBP cytoplasm intensity ####################")

# MBP intensity in the cytoplasm

# select only oligo count columns
HALO.data.dems.MBPCyto <- subset(HALO.data.dems, select = c(grepl("Subject|AgeDeath|Sex|Diagnosis|MBPCyto",colnames(HALO.data.dems))))
ROI.names.MBPCyto <- colnames(HALO.data.dems.MBPCyto)[grepl("MBPCyto",colnames(HALO.data.dems.MBPCyto))]

HALO.data.dems.MBPCyto.melt <- melt(HALO.data.dems.MBPCyto,id.vars=c("Subject","AgeDeath","Sex","Diagnosis"))
colnames(HALO.data.dems.MBPCyto.melt) <- c("Subject","AgeDeath","Sex","Diagnosis","SpD", "MBPCyto")
HALO.data.dems.MBPCyto.melt$SpD <- gsub("MBPCytoNorm_", "", HALO.data.dems.MBPCyto.melt$SpD)
HALO.data.dems.MBPCyto.melt$SpD <- factor(HALO.data.dems.MBPCyto.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.MBPCyto <- aov(data = HALO.data.dems.MBPCyto.melt, MBPCyto ~ AgeDeath + Sex + SpD*Diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.MBPCyto.emns <- emmeans(mdl.MBPCyto, pairwise ~ Diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.MBPCyto.emns))

# save some stats
write.csv(as.data.frame(summary(mdl.MBPCyto.emns)),
          paste(results.dir.stats,"/stats_MBPCyto_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_MBPCyto_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 750, height = 750, units = "px", res = 150)
plot.mdl <- ggplot(HALO.data.dems.MBPCyto.melt, aes(x = SpD, y = MBPCyto, fill = Diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = Diagnosis)) + 
  ylab("MBP cytoplasm intensity (normalised)") +
  theme_minimal()
print(plot.mdl)
dev.off()

print("#################### MBP nucleus intensity ####################")

# MBP intensity in the nucleus

# select only MCP nucleus columns
HALO.data.dems.MBPNucl <- subset(HALO.data.dems, select = c(grepl("Subject|AgeDeath|Sex|Diagnosis|MBPNucl",colnames(HALO.data.dems))))
ROI.names.MBPNucl <- colnames(HALO.data.dems.MBPNucl)[grepl("MBPNucl",colnames(HALO.data.dems.MBPNucl))]

HALO.data.dems.MBPNucl.melt <- melt(HALO.data.dems.MBPNucl,id.vars=c("Subject","AgeDeath","Sex","Diagnosis"))
colnames(HALO.data.dems.MBPNucl.melt) <- c("Subject","AgeDeath","Sex","Diagnosis","SpD", "MBPNucl")
HALO.data.dems.MBPNucl.melt$SpD <- gsub("MBPNuclNorm_", "", HALO.data.dems.MBPNucl.melt$SpD)
HALO.data.dems.MBPNucl.melt$SpD <- factor(HALO.data.dems.MBPNucl.melt$SpD, levels = c("SpD4", "SpD3", "SpD2", "SpD1"), ordered = TRUE)

# putting everything into one model rather than running loads
mdl.MBPNucl <- aov(data = HALO.data.dems.MBPNucl.melt, MBPNucl ~ AgeDeath + Sex + SpD*Diagnosis)
# calculate estimated marginal means (adjusting for age and sex)
mdl.MBPNucl.emns <- emmeans(mdl.MBPNucl, pairwise ~ Diagnosis | SpD, cov.reduce = mean, adjust = "fdr")

print(summary(mdl.MBPNucl.emns))

# save some stats
write.csv(as.data.frame(summary(mdl.MBPNucl.emns)),
          paste(results.dir.stats,"/stats_MBPNucl_group_comparison_PDvsCtrl",".csv",sep=""),
          row.names = FALSE)

# do a boxplot across groups
png(filename = paste(results.dir,"/plot_MBPNucl_group_comparison_PDvsCtrl",".png",sep=""), 
    width = 750, height = 750, units = "px", res = 150)
plot.mdl <- ggplot(HALO.data.dems.MBPNucl.melt, aes(x = SpD, y = MBPNucl, fill = Diagnosis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, shape = 21, aes(fill = Diagnosis)) + 
  ylab("MBP nucleus intensity (normalised)") +
  theme_minimal()
print(plot.mdl)
dev.off()
