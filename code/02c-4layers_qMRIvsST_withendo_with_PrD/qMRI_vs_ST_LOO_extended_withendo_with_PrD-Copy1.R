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
results.dir <- paste(d0, "/../../plots/02c-4layers_qMRIvsST_withendo_with_PrD/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/02c-4layers_qMRIvsST_withendo_with_PrD/", sep = "")

adj.qMRI.ST.medians <- read.csv(paste(results.dir.stats, "adj_qMRI_ST_medians_extended.csv",sep = ""))
adj.qMRI.ST.medians.no.SWM <- subset(adj.qMRI.ST.medians, SpD != "SWM")

##### UPDATETHIS
#POI.noSWM <- as.character(read.delim("pathways_of_interest_no_SWM.txt", sep = ",", header = FALSE, quote = "'")
POI.noSWM = colnames(adj.qMRI.ST.medians.no.SWM)[17:2333]
POI.noSWM <- sub("_Control$", "", POI.noSWM)

# colours for spatial domains
SpD.cols <- c("#EE6677","#CCBB44","#228833","#66CCEE")

# do LOO for these pathways

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


# define empty data frames for LOO results
GO.LOO.Ctrl.colnames <- paste(rep(POI.noSWM,each = 3),c("actual","predicted", "SE"),sep = "_")
GO.LOO.Ctrl.no.SWM <- data.frame(matrix(nrow = nrow(adj.qMRI.ST.medians.no.SWM), ncol = length(GO.LOO.Ctrl.colnames)))
colnames(GO.LOO.Ctrl.no.SWM) <- GO.LOO.Ctrl.colnames
GO.LOO.Ctrl.no.SWM <-cbind(adj.qMRI.ST.medians.no.SWM[,1:2],GO.LOO.Ctrl.no.SWM)

GO.LOO.stats.Ctrl.no.SWM <- data.frame(matrix(nrow = length(POI.noSWM), ncol = 6))
colnames(GO.LOO.stats.Ctrl.no.SWM) <- c("GO_type", "NRMSE", "T","rho", "R2", "p")

# run the LOO
GO.LOO.Ctrl.no.SWM.output <- run.qMRI.ST.LOO(dframe = adj.qMRI.ST.medians.no.SWM,
                                      ST = POI.noSWM,
                                      LOO = GO.LOO.Ctrl.no.SWM,
                                      stats = GO.LOO.stats.Ctrl.no.SWM,
                                      group.string = "Control")
GO.LOO.Ctrl.no.SWM <- GO.LOO.Ctrl.no.SWM.output$LOO
GO.LOO.stats.Ctrl.no.SWM <- GO.LOO.Ctrl.no.SWM.output$stats

# save some stats
GO.LOO.stats.Ctrl.no.SWM$corr_p <- p.adjust(GO.LOO.stats.Ctrl.no.SWM$p, method = "fdr")
write.csv(GO.LOO.stats.Ctrl.no.SWM, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_GO_pathway_LOO_no_SWM_extended_all.csv", sep = ""),
  row.names = FALSE)

