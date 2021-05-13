## K-means cluster analysis to assign spatial eco-regions based on species
## habitat importance scores

## Also produce ordination plots (nMDS and PCA) and fit environment vectors

## Ryan R Reisinger, Ben Raymond

setwd("~/RAATD_01/RAATD/additionalAnalysis/cluster")

# library(amap)
library(ggplot2)
library(viridis)
# library(ggridges)
library(raster)

#-------------------------------------
#-------------------------------------
# 1. Non-weighted

# Get the data
dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformed.RDS")

# Index to merge back later
dat$dx <- 1:nrow(dat)

# Only complete data
datS <- dat[complete.cases(dat[ , 3:18]), ]

#-------------------------------------
# K-means clustering
set.seed(20)

clusters <- Kmeans(datS[, 3:18], centers = 3, method = "manhattan", iter.max = 100)

# Save the cluster number in the data
datS$group <- clusters$cluster

# Plot spatially to check
ggplot(data = datS, aes(x = x, y = y, fill = as.factor(group))) + geom_raster() + coord_quickmap()

# ##### If smoothing is required
# 
# # # Must be applied on a raster
# # 
# # Merge back to complete dataframe
# dat <- merge(x = dat, y = datS[ , c("dx", "group")], all.x = T)
# dat$dx <- NULL
# 
# # Create the raster
# r <- rasterFromXYZ(dat[ , c("x", "y", "group")])
# 
# # Smooth
# library(rasterKernelEstimates)
# 
# W <- matrix(1,11,11)
# rSmooth <- rasterLocalCategoricalModes(r, W)
# 
# # And extract
# datS$group <- raster::extract(rSmooth, datS[ , c("x", "y")])
# 
# ##### End of smoothing

# Plot species compositions

# Go to long data.frame
library(tidyr)
datL <- gather(data = datS,
               ADPE, ANFS, ANPE, BBAL, CRAS, DMSA, EMPE, GHAL, HUWH, KIPE,
               LMSA, MAPE.ROPE, SOES, WAAL, WESE, WHCP,
               key = "Species", value = "Importance")

# Plot

# Histogram
library(gridExtra)

#### custom theme
theme_ryan <- function () { 
  theme_bw(base_size=9, base_family="") %+replace% 
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.line = element_line(colour = "black")
    )
}
####

# Colours
colrs <- c("#EE7733", "#CC3311", "#EE3377")

# Add foraging classifications
guilds <- read.csv("guilds.csv", stringsAsFactors = F)

datL <- merge(datL, guilds, by = "Species", all.x = T)

pdf("speciesComposition.pdf", width = ((75*0.0393701)*2)/0.77777, height = (70*0.0393701)/0.77777, pointsize = 9)

mn <- mean(datL[datL$group == 1, "Importance"])
p1 <- ggplot(datL[datL$group == 1, ], aes(x = reorder(Species, Importance, FUN = mean), y = Importance, fill = Guild, colour = Guild)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colrs) +
  scale_colour_manual(values = colrs) +
  guides(colour = "none", fill = "none") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) }) +
  facet_wrap(~group, ncol = 1) +
  geom_hline(yintercept = mn) +
  labs(x = "", y = "Habitat importance") +
  coord_flip() +
  theme_ryan()

mn <- mean(datL[datL$group == 2, "Importance"])
p2 <- ggplot(datL[datL$group == 2, ], aes(x = reorder(Species, Importance, FUN = mean), y = Importance, fill = Guild, colour = Guild)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colrs) +
  scale_colour_manual(values = colrs) +
  guides(colour = "none", fill = "none") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) }) +
  facet_wrap(~group, ncol = 1) +
  geom_hline(yintercept = mn) +
  labs(x = "", y = "Habitat importance") +
  coord_flip() +
  theme_ryan()

mn <- mean(datL[datL$group == 3, "Importance"])
p3 <- ggplot(datL[datL$group == 3, ], aes(x = reorder(Species, Importance, FUN = mean), y = Importance, fill = Guild, colour = Guild)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colrs) +
  scale_colour_manual(values = colrs) +
  guides(colour = "none", fill = "none") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) }) +
  facet_wrap(~group, ncol = 1) +
  geom_hline(yintercept = mn) +
  labs(x = "", y = "Habitat importance") +
  coord_flip() +
  theme_ryan()

grid.arrange(p1, p3, p2, nrow = 1, ncol = 3)
     
dev.off()


# Merge back to complete dataframe
datS <- datS[ , c("dx", "group")]
dat <- merge(x = dat, y = datS, all.x = T)
dat$dx <- NULL

# Write output
saveRDS(dat, "../../meanPredictions/griddedScoresEnvarsTransformedGroup.RDS")

#-------------------------------------
#-------------------------------------
# 2. Colony-weighted

rm(list = ls())

# Get the data
dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")

# Index to merge back later
dat$dx <- 1:nrow(dat)

# Only complete data
datS <- dat[complete.cases(dat[ , 3:18]), ]

#-------------------------------------
# K-means clustering
set.seed(20)

clusters <- Kmeans(datS[, 3:18], centers = 3, method = "manhattan", iter.max = 100)

# Match group assignment to non-weighted version
foo <- data.frame(oldcluster = clusters$cluster, newcluster = NA)

foo[foo$oldcluster == 1, "newcluster"] <- 2
foo[foo$oldcluster == 2, "newcluster"] <- 1
foo[foo$oldcluster == 3, "newcluster"] <- 3

# Save the cluster number in the data
datS$group <- foo$newcluster
rm(foo)

# Plot spatially to check
ggplot(data = datS, aes(x = x, y = y, fill = as.factor(group))) + geom_raster() + coord_quickmap()

# ##### If smoothing is required
# 
# # # Must be applied on a raster
# # 
# # Merge back to complete dataframe
# dat <- merge(x = dat, y = datS[ , c("dx", "group")], all.x = T)
# dat$dx <- NULL
# 
# # Create the raster
# r <- rasterFromXYZ(dat[ , c("x", "y", "group")])
# 
# # Smooth
# library(rasterKernelEstimates)
# 
# W <- matrix(1,11,11)
# rSmooth <- rasterLocalCategoricalModes(r, W)
# 
# # And extract
# datS$group <- raster::extract(rSmooth, datS[ , c("x", "y")])
# 
# ##### End of smoothing

# Plot species compositions

# Go to long data.frame
library(tidyr)
datL <- gather(data = datS,
               ADPE, ANFS, ANPE, BBAL, CRAS, DMSA, EMPE, GHAL, HUWH, KIPE,
               LMSA, MAPE.ROPE, SOES, WAAL, WESE, WHCP,
               key = "Species", value = "Importance")

# Plot

# Histogram
library(gridExtra)

#### custom theme
theme_ryan <- function () { 
  theme_bw(base_size=9, base_family="") %+replace% 
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.line = element_line(colour = "black")
    )
}
####

# Colours
colrs <- c("#EE7733", "#CC3311", "#EE3377")

# Add foraging classifications
guilds <- read.csv("guilds.csv", stringsAsFactors = F)

datL <- merge(datL, guilds, by = "Species", all.x = T)

pdf("speciesCompositionColonyWeighted.pdf", width = ((75*0.0393701)*2)/0.77777, height = (70*0.0393701)/0.77777, pointsize = 9)

mn <- mean(datL[datL$group == 1, "Importance"])
p1 <- ggplot(datL[datL$group == 1, ], aes(x = reorder(Species, Importance, FUN = mean), y = Importance, fill = Guild, colour = Guild)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colrs) +
  scale_colour_manual(values = colrs) +
  guides(colour = "none", fill = "none") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) }) +
  facet_wrap(~group, ncol = 1) +
  geom_hline(yintercept = mn) +
  labs(x = "", y = "Habitat importance") +
  coord_flip() +
  theme_ryan()

mn <- mean(datL[datL$group == 2, "Importance"])
p2 <- ggplot(datL[datL$group == 2, ], aes(x = reorder(Species, Importance, FUN = mean), y = Importance, fill = Guild, colour = Guild)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colrs) +
  scale_colour_manual(values = colrs) +
  guides(colour = "none", fill = "none") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) }) +
  facet_wrap(~group, ncol = 1) +
  geom_hline(yintercept = mn) +
  labs(x = "", y = "Habitat importance") +
  coord_flip() +
  theme_ryan()

mn <- mean(datL[datL$group == 3, "Importance"])
p3 <- ggplot(datL[datL$group == 3, ], aes(x = reorder(Species, Importance, FUN = mean), y = Importance, fill = Guild, colour = Guild)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colrs) +
  scale_colour_manual(values = colrs) +
  guides(colour = "none", fill = "none") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) }) +
  facet_wrap(~group, ncol = 1) +
  geom_hline(yintercept = mn) +
  labs(x = "", y = "Habitat importance") +
  coord_flip() +
  theme_ryan()

grid.arrange(p1, p3, p2, nrow = 1, ncol = 3)

dev.off()


# Merge back to complete dataframe
datS <- datS[ , c("dx", "group")]
dat <- merge(x = dat, y = datS, all.x = T)
dat$dx <- NULL

# Write output
saveRDS(dat, "../../meanPredictions/griddedScoresEnvarsTransformedColonyWeightedGroup.RDS")