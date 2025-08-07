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
qMRI.data.dir <- paste(d0, "/../../processed_data/02-4layers_qMRI/", sep = "")
cell.data.dir <- paste(d0, "/../../processed_data/01b-4layers_withendo/02-top3_cells_counts/", sep = "")
GO.data.dir <- paste(d0, "/../../processed_data/01b-4layers_withendo/05-pathway_scoring/", sep = "")
dems.dir<- paste(d0, "/../../raw_data/participant_info/", sep = "")
results.dir <- paste(d0, "/../../plots/02c-4layers_qMRIvsST_withendo_with_PrD/", sep = "")
results.dir.stats <- paste(d0, "/../../stats/02c-4layers_qMRIvsST_withendo_with_PrD/", sep = "")
results.dir.cellonly <- paste(d0, "/../../plots/01b-4layers_withendo/02-top3_cells_counts/", sep = "")


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
qMRI.metrics <- c("QSM","R1","R2s","PrD","MTsat")
MOI <- c("QSM","R1","R2s","PrD","MTsat")

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

# adjust all ROIs for age and sex using the covariance method
for (ROI in qMRI.ROI.names) {
  mdl.qMRI.adj <- lm(qMRI.data.dems[,ROI] ~ age + sex, data = qMRI.data.dems)
  adj.qMRI.vals <- qMRI.data.dems[,ROI] +
    diff.qMRI.age * mdl.qMRI.adj$coefficients[2] +
    diff.qMRI.sex * mdl.qMRI.adj$coefficients[3]
  adj.qMRI.data.dems[,ROI] <- adj.qMRI.vals
}

#### load in & reorganise cell data ####

# load in ST patient dems
ST.dems <- read.csv(paste(dems.dir,"ST_patient_info.csv",sep = ""))
ST.dems <- ST.dems[order(ST.dems$sample_id),]
ST.dems$sex <- factor(ST.dems$sex)

# specify the cell types we're interested in
COI <- c("Excit_L2_3", "Excit_L3","Excit_L3_4_5","Excit_L4" ,"Excit_L5","Excit_L5_6","Excit_L6","Inhib")

# get names of 4 layer cell count files
all.cell.fnames <- dir(cell.data.dir,"count_(.+)_4SpD_withEndo")

# select only the interesting ones (preserve COI order)
COI.fnames <- unique(unlist(lapply(COI, function(coi) {
  grep(coi, all.cell.fnames, value = TRUE)
})))

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

#### adjust ST data ####

# initialize data frames for the adjusted values
adj.cell.data.dems <- cell.data.dems

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

#### cell type boxplots ####

# subjects
subj <- adj.cell.data.dems$subject

# define spatial domains / subregions
dom <- c("SWM", "SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

# colours for spatial domains
SpD.cols <- c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677")

# do this for cell types
adj.cell.data.long <- expand.grid(subject = subj, SpD = dom, region = reg)
adj.cell.data.dems.melt <- melt(adj.cell.data.dems,id.vars=c("subject","age","sex","diagnosis"))

for (c in COI){
  cell.data <- adj.cell.data.dems.melt[grepl(paste0(c,"_S"), adj.cell.data.dems.melt$variable),6]
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

# cell boxplots for Ctrl
cell.boxplotlist <- lapply(1:length(COI), function(i) 
  cell.boxplot(ST = COI[i], adj.cell.data.long))

png(filename = paste(results.dir.cellonly, "/Control_neuron_boxplots.png", sep = ""),
    width = 1400, height = 1000, units = "px", res = 120)
cell.multiboxplot <- ggarrange(plotlist = cell.boxplotlist, 
                               nrow = 2, ncol = 4, 
                               common.legend = TRUE, 
                               legend = "right")
plot(cell.multiboxplot)
dev.off()

# now run without SWM
dom <- c("SpD4", "SpD3","SpD2","SpD1")
SpD.cols <- c("#66CCEE","#228833","#CCBB44","#EE6677")
adj.cell.data.long.no.SWM <- subset(adj.cell.data.long, SpD != "SWM")

# cell boxplots for Ctrl
cell.boxplotlist.no.SWM <- lapply(1:length(COI), function(i) 
  cell.boxplot(ST = COI[i], adj.cell.data.long.no.SWM))

png(filename = paste(results.dir.cellonly, "/Control_neuron_boxplots_no_SWM.png", sep = ""),
    width = 1400, height = 1000, units = "px", res = 120)
cell.multiboxplot.no.SWM <- ggarrange(plotlist = cell.boxplotlist.no.SWM, 
                                      nrow = 2, ncol = 4, 
                                      common.legend = TRUE, 
                                      legend = "right")
plot(cell.multiboxplot.no.SWM)
dev.off()


#### qMRI vs ST correlation analysis set-up ####

# subjects
subj <- adj.qMRI.data.dems$subject

# define spatial domains / subregions
dom <- c("SWM", "SpD4", "SpD3","SpD2","SpD1")
reg <- c("anterior", "middle", "posterior")

# colours for spatial domains
SpD.cols <- c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677")

# set up table containing qMRI / ST median values

adj.qMRI.ST.medians <- expand.grid(SpD = dom,region = reg)

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
  cell.data <- adj.cell.data.dems[,grepl(paste(c,"_S",sep = ""), colnames(adj.cell.data.dems))]
  cell.medians <- colMedians(as.matrix(cell.data), na.rm = TRUE)
  adj.qMRI.ST.medians[,paste(c, "_Control", sep = "")] <- cell.medians
}

write.csv(adj.qMRI.ST.medians, file = paste(
  results.dir.stats, "/adj_qMRI_ST_neuron_medians.csv", sep = ""), 
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
      if (ST == COI[length(COI)])
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

# define a function for cooks distance
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
      if (ST == COI[length(COI)])
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

# define a fucntion for residual plots
plot.ST.qMRI.residuals <- function(dframe, MRIST, group.string) {
  MRI <- MRIST[1]
  ST <- MRIST[2]
  
  x.data <- dframe[[paste(MRI, "_", group.string, sep = "")]]
  y.data <- dframe[[paste(ST, "_Control", sep = "")]]
  SpD <- dframe$SpD
  region <- dframe$region
  
  model <- lm(y.data ~ x.data)
  
  df.resid <- data.frame(
    fitted = fitted(model),
    residuals = resid(model),
    SpD = SpD,
    region = region
  )
  
  ggplot(df.resid, aes(x = fitted, y = residuals)) +
    geom_point(aes(colour = SpD, shape = region)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_minimal() +
    scale_color_manual(values = SpD.cols) +
    {
      if (ST == COI[length(COI)])
        xlab(MRI)
      else
        theme(axis.title.x = element_blank())
    } +
    {
      if (MRI == MOI[1])
        ylab(ST)
      else
        theme(axis.title.y = element_blank())
    } 
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

#### qMRI vs ST correlation with SWM ####

# CELL COUNTS

# generate a stats table
cell.stats.Ctrl <- lapply(1:nrow(MOICOI), function(i) stats.ST.qMRI(
  dframe = adj.qMRI.ST.medians, MRIST = MOICOI[i,], group.string = "Control"))
cell.stats.Ctrl <- as.data.frame(do.call(rbind, cell.stats.Ctrl))
cell.stats.Ctrl <- cbind(MOICOI,cell.stats.Ctrl)
colnames(cell.stats.Ctrl) <- c("qMRI", "cell_type", "T", "rho", "R2","p")
cell.stats.Ctrl$corr_p <- p.adjust(cell.stats.Ctrl$p, method = "fdr")
write.csv(cell.stats.Ctrl, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_neuron_count.csv", sep = ""), 
  row.names = FALSE)

# lapply the plot function to create a multi plot
cell.plot.list.Ctrl <- lapply(1:nrow(MOICOI), function(i) plot.ST.qMRI(
  dframe = adj.qMRI.ST.medians, MRIST = MOICOI[i, ], group.string = "Control",sts = cell.stats.Ctrl))

# Ctrl
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count.png", sep = ""),
    width = 1200, height = 1800, units = "px", res = 120)
cell.multi.plot.Ctrl <- ggarrange(
  plotlist = cell.plot.list.Ctrl,
  nrow = 8,
  ncol = 5,
  common.legend = TRUE,
  legend = "bottom"
)
plot(cell.multi.plot.Ctrl)
annotate_figure(cell.multi.plot.Ctrl, top = text_grob("Control qMRI versus control cell count"))
dev.off()

# lapply for cook's distances
cell.plot.list.cooksd.Ctrl <- lapply(1:nrow(MOICOI), function(i) plot.ST.qMRI.cooksd(
  dframe = adj.qMRI.ST.medians, MRIST = MOICOI[i, ], group.string = "Control"))

# Ctrl
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_cooksd.png", sep = ""),
    width = 1200, height = 1800, units = "px", res = 120)
cell.multi.plot.cooksd.Ctrl <- ggarrange(
  plotlist = cell.plot.list.cooksd.Ctrl,
  nrow = 8,
  ncol = 5,
  common.legend = TRUE,
  legend = "right"
)
plot(cell.multi.plot.cooksd.Ctrl)
annotate_figure(cell.multi.plot.cooksd.Ctrl, 
                top = text_grob("Effect of regions on HC qMRI vs. HC cell count correlation"),
                left = text_grob("Cook's distance", rot = 90))
dev.off()

# lapply Residuals vs Fitted plots for Control group
cell.plot.list.resid.Ctrl <- lapply(1:nrow(MOICOI), function(i) plot.ST.qMRI.residuals(
  dframe = adj.qMRI.ST.medians, MRIST = MOICOI[i, ], group.string = "Control"))

# Save the combined plot for Control
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_residuals.png", sep = ""),
    width = 1200, height = 1800, units = "px", res = 120)

cell.multi.plot.resid.Ctrl <- ggarrange(
  plotlist = cell.plot.list.resid.Ctrl,
  nrow = 8,
  ncol = 5,
  common.legend = TRUE,
  legend = "bottom"
)

plot(cell.multi.plot.resid.Ctrl)
annotate_figure(cell.multi.plot.resid.Ctrl,
                top = text_grob("Residuals vs Fitted: Control qMRI vs. Neuron Count"),
                left = text_grob("Residuals", rot = 90),
                bottom = text_grob("Fitted Values"))
dev.off()


#### qMRI vs ST LOO with SWM ####

# CELL COUNTS

# define empty data frames for LOO results
cell.LOO.Ctrl.colnames <- paste(rep(COI,each = 3),c("actual","predicted", "SE"),sep = "_")
cell.LOO.Ctrl <- data.frame(matrix(nrow = nrow(adj.qMRI.ST.medians), ncol = length(cell.LOO.Ctrl.colnames)))
colnames(cell.LOO.Ctrl) <- cell.LOO.Ctrl.colnames
cell.LOO.Ctrl <-cbind(adj.qMRI.ST.medians[,1:2],cell.LOO.Ctrl)

cell.LOO.stats.Ctrl <- data.frame(matrix(nrow = length(COI), ncol = 6))
colnames(cell.LOO.stats.Ctrl) <- c("cell_type", "NRMSE", "T","rho", "R2", "p")

# run the LOO
cell.LOO.Ctrl.output <- run.qMRI.ST.LOO(dframe = adj.qMRI.ST.medians,
                                        ST = COI, 
                                        LOO = cell.LOO.Ctrl, 
                                        stats = cell.LOO.stats.Ctrl, 
                                        group.string = "Control")
cell.LOO.Ctrl <- cell.LOO.Ctrl.output$LOO
cell.LOO.stats.Ctrl <- cell.LOO.Ctrl.output$stats

# save some stats
cell.LOO.stats.Ctrl$corr_p <- p.adjust(cell.LOO.stats.Ctrl$p, method = "fdr")
write.csv(cell.LOO.stats.Ctrl, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_neuron_count_LOO.csv", sep = ""), 
  row.names = FALSE)

# plot the LOO
plot.list.cell.LOO.Ctrl <- lapply(1:length(COI), function(i)
  plot.LOO(ST = COI[i], LOO = cell.LOO.Ctrl, sts = cell.LOO.stats.Ctrl))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_LOO.png", sep = ""),
    width = 1200, height = 600, units = "px", res = 120)
multi.plot.cell.LOO.Ctrl <- ggarrange(
  plotlist = plot.list.cell.LOO.Ctrl,
  nrow = 2,
  ncol = 4,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.cell.LOO.Ctrl)
annotate_figure(multi.plot.cell.LOO.Ctrl, bottom = text_grob("Actual cell count"),
                left = text_grob("Predicted cell count", rot = 90),
                top = text_grob("Control qMRI medians predicting control cell count"))
dev.off()

# plot the LOO cook's distance
plot.list.cell.LOO.cooksd.Ctrl <- lapply(1:length(COI), function(i) plot.LOO.cooksd(ST = COI[i], cell.LOO.Ctrl))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_LOO_cooksd.png", sep = ""),
    width = 1200, height = 600, units = "px", res = 120)
multi.plot.cell.LOO.cooksd.Ctrl <- ggarrange(
  plotlist = plot.list.cell.LOO.cooksd.Ctrl,
  nrow = 2,
  ncol = 4,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.cell.LOO.cooksd.Ctrl)
annotate_figure(multi.plot.cell.LOO.cooksd.Ctrl,  
                left = text_grob("Cook's distance", rot = 90),
                top = text_grob("Influence of regions on cell LOO correlations"))
dev.off()


#### qMRI vs ST correlation no SWM ####

# remove the SWM ROIs
adj.qMRI.ST.medians.no.SWM <- subset(adj.qMRI.ST.medians, SpD != "SWM")

# define cols without SWM
SpD.cols <- c("#66CCEE","#228833","#CCBB44","#EE6677")

# CELL COUNTS

# generate a stats table
cell.stats.Ctrl.no.SWM <- lapply(1:nrow(MOICOI), function(i) stats.ST.qMRI(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOICOI[i,], group.string = "Control"))
cell.stats.Ctrl.no.SWM <- as.data.frame(do.call(rbind, cell.stats.Ctrl.no.SWM))
cell.stats.Ctrl.no.SWM <- cbind(MOICOI,cell.stats.Ctrl.no.SWM)
colnames(cell.stats.Ctrl.no.SWM) <- c("qMRI", "cell_type", "T", "rho", "R2","p")
cell.stats.Ctrl.no.SWM$corr_p <- p.adjust(cell.stats.Ctrl.no.SWM$p, method = "fdr")
write.csv(cell.stats.Ctrl.no.SWM, file = paste(
  results.dir.stats, "/Stats_control_qMRI_vs_neuron_count_no_SWM.csv", sep = ""), 
  row.names = FALSE)

# lapply the plot function to create a multi plot
cell.plot.list.Ctrl.no.SWM <- lapply(1:nrow(MOICOI), function(i) plot.ST.qMRI(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOICOI[i, ], group.string = "Control", sts = cell.stats.Ctrl.no.SWM))

# Ctrl
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_no_SWM.png", sep = ""),
    width = 1200, height = 1800, units = "px", res = 120)
cell.multi.plot.Ctrl.no.SWM <- ggarrange(
  plotlist = cell.plot.list.Ctrl.no.SWM,
  nrow = 8,
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
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_cooksd_no_SWM.png", sep = ""),
    width = 1200, height = 1800, units = "px", res = 120)
cell.multi.plot.cooksd.Ctrl.no.SWM <- ggarrange(
  plotlist = cell.plot.list.cooksd.Ctrl.no.SWM,
  nrow = 8,
  ncol = 5,
  common.legend = TRUE,
  legend = "right"
)
plot(cell.multi.plot.cooksd.Ctrl.no.SWM)
annotate_figure(cell.multi.plot.cooksd.Ctrl.no.SWM, 
                top = text_grob("Effect of regions on HC qMRI vs. HC cell count correlation"),
                left = text_grob("Cook's distance", rot = 90))
dev.off()

# lapply Residuals vs Fitted plots for Control group
cell.plot.list.resid.Ctrl.no.SWM <- lapply(1:nrow(MOICOI), function(i) plot.ST.qMRI.residuals(
  dframe = adj.qMRI.ST.medians.no.SWM, MRIST = MOICOI[i, ], group.string = "Control"))

# Save the combined plot for Control
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_residuals_no_SWM.png", sep = ""),
    width = 1200, height = 1800, units = "px", res = 120)

cell.multi.plot.resid.Ctrl.no.SWM <- ggarrange(
  plotlist = cell.plot.list.resid.Ctrl.no.SWM,
  nrow = 8,
  ncol = 5,
  common.legend = TRUE,
  legend = "bottom"
)

plot(cell.multi.plot.resid.Ctrl.no.SWM)
annotate_figure(cell.multi.plot.resid.Ctrl.no.SWM,
                top = text_grob("Residuals vs Fitted: Control qMRI vs. Neuron Count"),
                left = text_grob("Residuals", rot = 90),
                bottom = text_grob("Fitted Values"))
dev.off()

#### qMRI vs ST LOO no SWM ####

# CELL COUNTS

# define empty data frames for LOO results
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
  results.dir.stats, "/Stats_control_qMRI_vs_neuron_count_LOO_no_SWM.csv", sep = ""), 
  row.names = FALSE)

# plot the LOO
plot.list.cell.LOO.Ctrl.no.SWM <- lapply(1:length(COI), function(i) 
  plot.LOO(ST = COI[i], LOO = cell.LOO.Ctrl.no.SWM, sts = cell.LOO.stats.Ctrl.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_LOO_no_SWM.png", sep = ""),
    width = 1200, height = 500, units = "px", res = 120)
multi.plot.cell.LOO.Ctrl.no.SWM <- ggarrange(
  plotlist = plot.list.cell.LOO.Ctrl.no.SWM,
  nrow = 2,
  ncol = 4,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.cell.LOO.Ctrl.no.SWM)
annotate_figure(multi.plot.cell.LOO.Ctrl.no.SWM, bottom = text_grob("Actual cell count"),
                left = text_grob("Predicted cell count", rot = 90),
                top = text_grob("Control qMRI medians predicting control cell count"))
dev.off()

# plot the LOO cook's distance
plot.list.cell.LOO.cooksd.Ctrl.no.SWM <- lapply(1:length(COI), function(i) plot.LOO.cooksd(ST = COI[i], cell.LOO.Ctrl.no.SWM))
png(filename = paste(results.dir, "/Plot_control_qMRI_vs_neuron_count_LOO_cooksd_no_SWM.png", sep = ""),
    width = 1200, height = 500, units = "px", res = 120)
multi.plot.cell.LOO.cooksd.Ctrl.no.SWM <- ggarrange(
  plotlist = plot.list.cell.LOO.cooksd.Ctrl.no.SWM,
  nrow = 2,
  ncol = 4,
  common.legend = TRUE,
  legend = "right"
)
plot(multi.plot.cell.LOO.cooksd.Ctrl.no.SWM)
annotate_figure(multi.plot.cell.LOO.cooksd.Ctrl.no.SWM,  
                left = text_grob("Cook's distance", rot = 90),
                top = text_grob("Influence of regions on cell LOO correlations"))
dev.off()
