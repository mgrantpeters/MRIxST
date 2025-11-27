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
GO.data.dir <- paste(d0, "/../../processed_data/01b-4layers_withendo/05-pathway_scoring/", sep = "")
GO.data.dir <- paste(d0, "/../../processed_data/01b-4layers_withendo/05-pathway_scoring/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir <- paste(d0, "/../../plots/02c-4layers_qMRIvsST_withendo_with_PrD/", sep = "")
results.dir.GO.only <- paste(d0, "/../../plots/01b-4layers_withendo/05-pathway_scoring/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/02c-4layers_qMRIvsST_withendo_with_PrD/", sep = "")

adj.qMRI.ST.medians <- read.csv(paste(results.dir.stats, "adj_qMRI_ST_medians_extended.csv",sep = ""))
adj.qMRI.ST.medians.no.SWM <- subset(adj.qMRI.ST.medians, SpD != "SWM")


#### load in GO dems and info ----

# load in ST patient dems
ST.dems <- read.csv(paste(dems.dir,"ST_patient_info.csv",sep = ""))
ST.dems <- ST.dems[order(ST.dems$sample_id),]
ST.dems$sex <- factor(ST.dems$sex)

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
    factor(GO.data$region, levels = c("WM", "SpD4","SpD3","SpD2","SpD1"))),]
  GO.data <- pivot_wider(GO.data,names_from = region, values_from = POI.names[i],
                         names_prefix = paste(POI.names[i],"_", sep = ""))
  GO.data.dems <- cbind(GO.data.dems,GO.data[,2:ncol(GO.data)])
}
colnames(GO.data.dems) <- gsub("_WM", "_SWM", colnames(GO.data.dems))

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

#### adjust GO data ----

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

#### plot GO data ----

# subjects
subj <- adj.GO.data.dems$subject

# define spatial domains / subregions
dom <- c("SWM", "SpD4", "SpD3","SpD2","SpD1")
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
SpD.cols <- c("#EE6677","#CCBB44","#228833","#66CCEE")

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

#### define functions for correlation / LOO analyses ----

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
        list(
          ylab(ST),
          theme(axis.title.y = element_text(size = 10))  
        )
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
             size = 2.75,
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
        list(
          ylab(ST),
          theme(axis.title.y = element_text(size = 10)) 
        )
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


#### run correlation / LOO analyses ----

MOI <- c("MTsat","PrD", "QSM", "R1","R2s")

MOIPOI <- expand.grid(MOI,POI.names,stringsAsFactors = FALSE)

COI <- c("Astro","Micro", "Oligo", "OPC")

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
