## K-means cluster analysis to assign spatial eco-regions based on species
## habitat importance scores

## Also produce ordination plots (nMDS and PCA) and fit environment vectors

## Ryan R Reisinger, Ben Raymond

setwd("~/RAATD_01/RAATD/additionalAnalysis/hotspot")

# Plotting
library(ggplot2)
library(viridis)
# library(ggridges)
library(colorspace)
library(SOmap)
library(pals)

# Raster
library(raster)

# For hotspot
# library(devtools)
# devtools::install_github("biggis-project/soh", build_vignettes = TRUE)
library(soh)

# For clustering
library(cluster)
# library(factoextra)
library(vegan)

# Parallel processing
library(foreach)
library(parallel)
library(doParallel)

#-------------------------------------
## 1. Calculate hotspots

# Source the meanR function
source("~/RAATD_01/RAATD/Code/function_meanR.R")

# Get the data - weighted or non-weighted
if (TRUE) {
  dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
} else {
  dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformed.RDS")
}

# Calculate means
dat <- meanR(this.data = dat, these.groups = "ALL")

# Create a raster for hotspot
r <- rasterFromXYZ(dat[ , c("x", "y", "MEAN")])

## Hotspot
# Calculate Gi*

# Can calculate Gi* using a circular window with sigmod falloff
if (TRUE) {
w_sigma <- weight_matrix_circular_fade(7, 2)
g <- GetisOrd(r, w_sigma)
}

# Or circular window with sharp cutoffs (throws error)
if (FALSE) {
w_circular <- weight_matrix_circular(11)
g <- GetisOrdLocalStats(r, w_circular)
}

# Or a simple square matrix (throws error)
if (FALSE) {
w_square <- weight_matrix_squared(11)
g <- GetisOrdLocalStats(r, w_square)
}

# There is a strong edge effect
# Filter to try and deal with these

# Function to detecent any NAs
foo <- function(x) {
  if (anyNA(x)) {
    return(NA)
  } else {
    sum(x)
  }
}

# Weight matrix to prevent any smoothing
no_smooth <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,1,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0),
                    nrow=13)

g <- focal(g, w = no_smooth, fun = foo)

# Plot
plot(g)
plot_zscore(g, sigma_thresh = 5, hotspot_quantile_thresh = 0.10)

# Simple plot in raster
max(abs(minValue(g)), abs(maxValue(g)))

l <- length(seq(-20, 20, 0.5))
plot(g, col = ocean.balance(l),
     breaks = seq(-20, 20, 0.5),
     axis.args=list(at = c(-20, 0, 20),
                    labels = c(-20, 0, 20)))

# Create contours corresponding with the quantiles identified in
# 'plot_zscore'
crs(g) <- "+proj=longlat +datum=WGS84 +no_defs"
contour_lines <- rasterToContour(g, levels = c(-6.21, 7.16))

lines(contour_lines)

# SOmap

tiff("./plots/hotspot_polar.tiff", width = 8, height = 6, units = "in", res = 300)
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)

SOplot(g,
       col = ocean.balance(l),
       breaks = seq(-20, 20, 0.5),
       axis.args=list(at = c(-20, 0, 20),
                      labels = c(-20, 0, 20)))

SOplot(contour_lines)

dev.off()

# Extract the values
dat$z <- extract(g, dat[,c("x", "y")])

# Categorise
dat$hotspot <- NA
dat[dat$z > 0 & !is.na(dat$z), "hotspot"] <- "Y"
dat[dat$z < 0 & !is.na(dat$z), "hotspot"] <- "N"

#-------------------------------------
## 2. Clustering with CLARA

# Only complete data
dx <- which(complete.cases(dat[ , 3:18])) # Index to merge back later
datS <- dat[complete.cases(dat[ , 3:18]), ]

# Scale the habitat importance scores
datSc <- scale(datS[ , 3:18])

# set.seed(20)

# Go through different values of k to determine best number, based on silhoutte width
# following a one-dimensional example from Christian Hennig:

# Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. J. Comput. Appl. Math., 20, 53-65.


if (TRUE) {
  # ----------------
  # Parallel
  clust <- makeCluster(detectCores() - 1) # leave 1 core for OS
  registerDoParallel(clust)
  
  asw <- foreach (k = 2:40, .combine = 'c', .inorder = TRUE, .packages = "cluster") %dopar% {
    cluster::clara(x = datSc, k = k, metric = "manhattan", stand = FALSE, samples = 100)$silinfo$avg.width
  }
  
  stopCluster(clust)
  
  asw <- c(0, asw)
  
} else {
  
  # ----------------
  # Non-parallel
  asw <- numeric(40)
  
  for (k in 2:40) {
    print(k)
    asw[k] <- clara(x = datSc, k = k, metric = "manhattan", stand = FALSE, samples = 100)$silinfo$avg.width
  }
  
}

# ----------------
# Plot the profile
k1 <- which.max(asw)
k2 <-  which(asw == sort(asw, decreasing = TRUE)[2])
k3 <-  which(asw == sort(asw, decreasing = TRUE)[3])

cat("silhouette-optimal number of clusters:", k.best, "\n")
plot(1:40, asw, type= "h", main = "CLARA clustering assessment",
     xlab= "k  (number of clusters)", ylab = "Average silhouette width")
axis(1, k1, paste("",k1,sep="\n"), col = "red", col.axis = "red")
axis(1, k2, paste("",k2,sep="\n"), col = "red", col.axis = "red")
axis(1, k3, paste("",k3,sep="\n"), col = "red", col.axis = "red")

# ----------------
# Calculate the clusters again, using the best value of k
# k.best <- 17
# k.best <- 33
clusters <- clara(x = datSc, k = k3, metric = "manhattan", stand = FALSE, samples = 100)

# Silhouette width
clusters$silinfo$avg.width

# ----------------
# Save the cluster number in the data
datS$group <- clusters$cluster

# ----------------
# Plot spatially to check
ggplot(data = datS, aes(x = x, y = y, fill = as.factor(group))) + geom_raster() + coord_quickmap() +
  scale_fill_manual(values = rainbow(length(unique(datS$group))), name = "Cluster") +
  labs(x = "", y = "")

# Calculate if each cluster is mainly hot or cold
foo <- function(x) {
  x$cluster.type <- NA
  grps <- unique(x$group)
  for (i in grps) {
    print(i)
  dt <- x[x$group == i, ]
  hot <- nrow(dt[dt$hotspot == "Y", ])
  cold <- nrow(dt[dt$hotspot == "N", ])
  this.grp <- c("hot", "cold")[which(c(hot, cold) == max(hot, cold))]
  x[x$group == i, "cluster.type"] <- this.grp
  }
  return(x)
}

datS <- foo(datS)

k.hot <- length(unique(datS[datS$cluster.type == "hot", "group"]))
k.cold <- length(unique(datS[datS$cluster.type == "cold", "group"]))

# Plot only hotspots
ggplot(data = datS[datS$cluster.type == "hot", ], aes(x = x, y = y, fill = as.factor(group))) + geom_tile() + coord_quickmap() +
  scale_fill_manual(values = rainbow(k.hot))

## Plot combined

# New group names
datS$group2 <- paste0(datS$cluster.type, "_", datS$group)
datS$group2 <- paste(datS$cluster.type, formatC(datS$group, width=2, flag="0"), sep="_")

# Create a combined palette that groups cold clusters
pal <- c(sequential_hcl(k.cold, palette = "Light Grays"), rainbow(k.hot))

ggplot(data = datS, aes(x = x, y = y, fill = as.factor(group2))) + geom_tile() + coord_quickmap() +
  scale_fill_manual(values = pal, name = "Cluster") +
  labs(x = "", y = "")

## Plot with SOmap

# # Merge back to complete dataframe
dat$group2 <- NA
dat[dx, "group2"] <- datS$group2
dat[dx, "group"] <- datS$group

# # Create the raster
r <- rasterFromXYZ(dat[ , c("x", "y", "group")])
crs(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

tiff("plots/clusterMap.tiff", height = 8, width = 12, units = "in", res = 300)
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)
SOplot(r, col = pal, legend = F)
legend(x='right', legend = seq(1, 17, 1), fill = pal)
dev.off()

##### If smoothing is required

# Smooth
library(rasterKernelEstimates)

W <- matrix(1,11,11)
rSmooth <- rasterLocalCategoricalModes(r, W)

# And extract
datS$group <- raster::extract(rSmooth, datS[ , c("x", "y")])

# Plot smoothed
tiff("plots/clusterMapSmooth.tiff", height = 8, width = 12, units = "in", res = 300)
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)
SOplot(rSmooth, col = pal, legend = F)
legend(x='right', legend = seq(1, 17, 1), fill = pal)
dev.off()

##### End of smoothing

#-------------------------------------
## 3. Hierarchical clustering of these results
## By Ben Raymond

# Species
sp_names <- names(datS)[3:18]

# Desired number of groups for hierarchical clustering
# This would usually be less than k.best, but can be the same if all we want is a dendrogram
n_groups <- k3

# Mean values of each species
xc <- do.call(rbind, lapply(seq_len(k3),
                            function(z) colMeans(datS[datS$group == z, sp_names])))

## dissimilarities of these clusters
D <- vegdist(xc, method = "gower")

## hierarchical clustering
hcl <- hclust(D, method = "ave")

## Get the original cluster names
dorder <- order.dendrogram(as.dendrogram(hcl))
lbls <- unique(datS[, c("group", "group2")])
hcl$labels <- lbls[dorder, ]$group2

## And get the colours from the hot-cold palette made above
lbls$col <- NULL
lbls[order(lbls$group2), "col"] <- pal

## now extract the desired number of groups from the dendrogram
if (floor(n_groups) == n_groups) {
  ## we specified a number of groups directly
  cn_new <- cutree(hcl, k = n_groups)
  ## work out the dissimilarity level (height) that corresponds to this number of groups
  temph <- mean(c(hcl$height[length(hcl$height)+2-n_groups], hcl$height[length(hcl$height)+2-n_groups-1]))
} else {
  ## we specified a height at which to cut the dendrogram
  ## show on the dendrogram the height at which we are cutting
  temph <- n_groups
  cn_new <- cutree(hcl, h = n_groups)
  n_groups <- length(unique(cn_new))
}

## plot the dendrogram
plot(hcl, labels = hcl$labels, hang = -1)
if (k.best != n_groups) {
lines(c(1, k.best), c(temph, temph), lty = 2, col = 2)
}

## add markers for group labels
dorder <- order.dendrogram(as.dendrogram(hcl))

if (FALSE) {
  # Use the new groups for colours
  for (k in 1:n_groups) {
    temp <- which(cn_new[dorder] == k)
    points(temp, rep(0, length(temp)), col = rainbow(k.best)[k], bg = rainbow(k.best)[k], pch = 21, cex = 2)
  }
} else {
  # Or use predefined colours from the palette defined for the original groups
  for (k in 1:n_groups) {
    temp <- which(cn_new[dorder] == k)
    points(temp, rep(0, length(temp)), col = lbls$col[dorder][k], bg = lbls$col[dorder][k], pch = 21, cex = 2)
  }
  
}

## Map the k-means clusters to the final clusters from the hierarchical step
clst_final <- rep(NA_integer_, length(datS$group))
for (k in seq_len(k.best)) {
  clst_final[datS$group == k] <- cn_new[k]
}

## Add these new clusters to the data.frame
datS$cluster_h <- clst_final

## Spatial plot
## If we've defined fewer groups...
## otherwise they'll be the same as for k-means
ggplot(data = datS, aes(x = x, y = y, fill = as.factor(cluster_h))) + geom_tile() + coord_quickmap() +
  scale_fill_manual(values = rainbow(n_groups))


#-------------------------------------
## 3. Plot species compositions

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
colrs <- c("#0077BB", "#009988", "#EE3377")

# Add foraging classifications
guilds <- read.csv("../cluster/guilds.csv", stringsAsFactors = F)

datL <- merge(datL, guilds, by = "Species", all.x = T)

plist <- list()

for (i in unique(datL$group2)) {

p1 <- ggplot(datL[datL$group2 == i, ], aes(x = reorder(Species, Importance, FUN = mean), y = Importance, fill = Guild, colour = Guild)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colrs) +
  scale_colour_manual(values = colrs) +
  guides(colour = "none", fill = "none") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) }) +
  facet_wrap(~group2, ncol = 1, scales = "free_x") +
  labs(x = "", y = "Habitat importance") +
  coord_flip() +
  theme_ryan()

plist[i] <- list(p1)

}

library(gridExtra)

ggsave(file = "speciesCompositionClusters.pdf", arrangeGrob(grobs = plist, ncol = 4),
       width = 20, height = 29, units = "cm") 

# Merge back to complete dataframe
datS <- datS[ , c("dx", "group")]
dat <- merge(x = dat, y = datS, all.x = T)
dat$dx <- NULL

# Write output
# saveRDS(dat, "griddedGroup.RDS")

#-------------------------------------