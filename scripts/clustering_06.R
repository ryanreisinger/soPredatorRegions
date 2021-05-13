## K-means cluster analysis to assign spatial eco-regions based on species
## habitat importance scores

## Also produce ordination plots (nMDS and PCA) and fit environment vectors

## Ryan R Reisinger, Ben Raymond

setwd("~/RAATD_01/RAATDhotspot")

# Plotting
library(ggplot2)
library(viridis)
# library(ggridges)
library(colorspace)
library(SOmap)
library(pals)

# Raster
library(raster)

# For clustering
library(cluster)
# library(factoextra)
library(vegan)

# Parallel processing
library(foreach)
library(parallel)
library(doParallel)

#-------------------------------------
## 1. Calculate mean habitat importance upper and lower quantiles

# Source the meanR function
source("~/RAATD_01/RAATD/Code/function_meanR.R")

# Get the data - weighted or non-weighted
if (TRUE) {
  dat <- readRDS("../RAATD/meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
} else {
  dat <- readRDS("../RAATD/meanPredictions/griddedScoresEnvarsTransformed.RDS")
}

# Calculate means
dat <- meanR(this.data = dat, these.groups = "ALL")

# Create a raster for hotspot
r <- rasterFromXYZ(dat[ , c("x", "y", "MEAN")])
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"

# Simple plot in raster
plot(r, col = ocean.balance(125))

# Calculate top and bottom quantiles
lower_q <- quantile(r, probs = c(0.25, 0.50, 0.75))[1]
median_q <- quantile(r, probs = c(0.25, 0.50, 0.75))[2]
upper_q <- quantile(r, probs = c(0.25, 0.50, 0.75))[3]

# Create contours corresponding with the upper and lower quantiles
contour_lines <- rasterToContour(r, levels = c(lower_q, upper_q))
crs(contour_lines) <- "+proj=longlat +datum=WGS84 +no_defs"

# Plot
plot(r, col = ocean.balance(125))
lines(contour_lines)

# SOmap

tiff("./plots/hotspot_quartile_polar.tiff", width = 8, height = 6, units = "in", res = 300)
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.0001,
      straight = TRUE,
      graticules = TRUE)

SOplot(SOproj(r),
       col = ocean.balance(125))

SOplot(contour_lines)

dev.off()

# Categorise
dat$hotspot <- NA
dat[dat$MEAN > upper_q & !is.na(dat$MEAN), "hotspot"] <- "Y"
dat[dat$MEAN < lower_q & !is.na(dat$MEAN), "hotspot"] <- "N"

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
this_k <- 17

clusters <- clara(x = datSc, k = this_k, metric = "manhattan", stand = FALSE, samples = 100)

# Silhouette width
clusters$silinfo$avg.width

# ----------------
# Save the cluster number in the data
datS$group <- clusters$cluster

# ----------------
# Plot spatially to check
ggplot(data = datS, aes(x = x, y = y, fill = as.factor(group))) + geom_raster() + coord_quickmap() +
  scale_fill_manual(values = glasbey(length(unique(datS$group))), name = "Cluster") +
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
  scale_fill_manual(values = glasbey(k.hot))

## Plot combined

# New group names
datS$group2 <- paste0(datS$cluster.type, "_", datS$group)
datS$group2 <- paste(datS$cluster.type, formatC(datS$group, width=2, flag="0"), sep="_")

# Create a combined palette that groups cold clusters
pal <- c(sequential_hcl(k.cold, palette = "Light Grays"), glasbey(k.hot))

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

# Create a table to store names and colours
rat <- unique(dat[,c("group", "group2")])
rat <- rat[!is.na(rat$group), ]
rat <- rat[order(rat$group2), ]
rat$pal <- pal
rat$newgroup <- 1:nrow(rat)

r <- reclassify(r, rat[,c("group", "newgroup")])

# Need to reproject using nearest neighbour lookup, otherwise
# integer values become numerical and color palette is messed up
r_proj <- projectRaster(r, crs = SOcrs(), method = "ngb")

tiff("plots/clusterMap_quantile.tiff", height = 8, width = 12, units = "in", res = 300)
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)
SOplot(r_proj, col = rat$pal, legend = FALSE)
legend(x='right', legend = rat$group2, fill = rat$pal)
dev.off()

##### If smoothing is required

# # Smooth
# library(rasterKernelEstimates)
# 
# W <- matrix(1,11,11)
# rSmooth <- rasterLocalCategoricalModes(r, W)
# 
# # And extract
# datS$groupSmooth <- raster::extract(rSmooth, datS[ , c("x", "y")])
# 
# # New group names
# datS$group2Smooth <- paste0(datS$cluster.type, "_", datS$groupSmooth)
# datS$group2Smooth <- paste(datS$cluster.type, formatC(datS$groupSmooth, width=2, flag="0"), sep="_")
# 
# # Plot smoothed
# tiff("plots/clusterMapSmooth.tiff", height = 8, width = 12, units = "in", res = 300)
# SOmap(trim = -40,
#       bathy_legend = FALSE,
#       border_col = c("white", "white"),
#       border_width = 0.01,
#       straight = TRUE,
#       graticules = TRUE)
# SOplot(rSmooth, col = pal, legend = F)
# legend(x='right', legend = seq(1, this_k, 1), fill = pal)
# dev.off()

##### End of smoothing

#-------------------------------------
## 3. Hierarchical clustering of these results
## By Ben Raymond

# Species
sp_names <- names(datS)[3:18]

# Desired number of groups for hierarchical clustering
# This would usually be less than k.best, but can be the same if all we want is a dendrogram
n_groups <- this_k

# Mean values of each species
xc <- do.call(rbind, lapply(seq_len(this_k),
                            function(z) colMeans(datS[datS$group == z, sp_names])))

## dissimilarities of these clusters
D <- vegdist(xc, method = "gower")

## hierarchical clustering
hcl <- hclust(D, method = "ave")

## Get the original cluster names
dorder <- order.dendrogram(as.dendrogram(hcl))
lbls <- unique(datS[, c("group", "group2")])
# hcl$labels <- lbls[dorder, ]$group2 # Don't need to be ordered, plot seems to take care of the order
hcl$labels <- lbls$group2

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
pdf ("./plots/dendrogram_quartiles.pdf", height = 6, width = 6, useDingbats = F)
plot(hcl, hang = -1)
if (this_k != n_groups) {
lines(c(1, this_k), c(temph, temph), lty = 2, col = 2)
}

## add markers for group labels
dorder <- order.dendrogram(as.dendrogram(hcl))

if (FALSE) {
  # Use the new groups for colours
  for (k in 1:n_groups) {
    temp <- which(cn_new[dorder] == k)
    points(temp, rep(0, length(temp)), col = glasbey(this_k)[k], bg = glasbey(this_k)[k], pch = 21, cex = 2)
  }
} else {
  # Or use predefined colours from the palette defined for the original groups
  for (k in 1:n_groups) {
    temp <- which(cn_new[dorder] == k)
    points(temp, rep(0, length(temp)), col = lbls$col[k], bg = lbls$col[k], pch = 21, cex = 2)
  }
  
}

dev.off()

## Map the k-means clusters to the final clusters from the hierarchical step
clst_final <- rep(NA_integer_, length(datS$group))
for (k in seq_len(this_k)) {
  clst_final[datS$group == k] <- cn_new[k]
}

## Add these new clusters to the data.frame
datS$cluster_h <- clst_final

## Spatial plot
## If we've defined fewer groups...
## otherwise they'll be the same as for k-means
ggplot(data = datS, aes(x = x, y = y, fill = as.factor(cluster_h))) + geom_tile() + coord_quickmap() +
  scale_fill_manual(values = glasbey(n_groups))


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
colrs <- c("#0077BB", "#EE7733", "#EE3377")

# Add foraging classifications
guilds <- read.csv("./dat_in/guilds.csv", stringsAsFactors = F)

datL <- merge(datL, guilds, by = "Species", all.x = T)

# Add supercluster labels
library(dplyr)
datL <- datL %>%
  mutate(super_group = case_when(group2 == "hot_13" ~ "antarctic",
                                 group2 == "hot_15" ~ "antarctic",
                                 group2 == "hot_16" ~ "antarctic",
                                 group2 == "cold_14" ~ "antarctic",
                                 group2 == "hot_17" ~ "antarctic",
                                 group2 == "hot_09" ~ "scotia_arc",
                                 group2 == "hot_10" ~ "scotia_arc",
                                 group2 == "hot_03" ~ "subantarctic_distant",
                                 group2 == "hot_08" ~ "subantarctic_distant",
                                 group2 == "hot_04" ~ "subantarctic_distant",
                                 group2 == "cold_01" ~ "subantarctic_distant",
                                 group2 == "cold_02" ~ "subantarctic_distant",
                                 group2 == "hot_05" ~ "subantarctic",
                                 group2 == "hot_06" ~ "subantarctic",
                                 group2 == "hot_07" ~ "subantarctic",
                                 group2 == "hot_11" ~ "subantarctic",
                                 group2 == "hot_12" ~ "subantarctic"))

datL <- unite(data = datL, col = "super_group_cluster", c("super_group", "group2"), sep = "_", remove = FALSE)

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

ggsave(file = "./plots/speciesCompositionClusters_quartile.pdf", arrangeGrob(grobs = plist, ncol = 4),
       width = 20, height = 29, units = "cm") 

# Merge back to complete dataframe
dat[dx, "group"] <- datS$group

# Heatplot of mean habitat importance for each species X cluster
library(dplyr)
datsum <- group_by(datL, Species, group2)
datsum <- summarise_at(datsum, .vars = "Importance", .funs = mean)

lbls$order <- NA
lbls[dorder, ]$order <- 1:this_k
datsum <- merge(x = datsum, y = lbls[, c("group2", "order")], by = "group2", all.x = T)
datsum$group_order <- paste0(formatC(datsum$order, width=2, flag="0"), "_", datsum$group)

pdf("./plots/heatplotSpecies_quartile.pdf", width = 7, height = 6)
ggplot(data = datsum, aes(x = group_order, y = Species, fill = Importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean\nhabitat\nimportance") +
  labs(x = "Cluster", y = "Species") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(labels = lbls[dorder, ]$group2)
dev.off()

# Write output
saveRDS(dat, "./dat_out/gridded_with_clusters_quartile.RDS")

#-------------------------------------
## 4. Ordination

# Subsample per cluster
n_subs <- 500

dat_sub <- datS[complete.cases(datS[ , 19:37]), ] # Use only if envars are complete
dat_sub <- dat_sub %>% group_by(group2)
dat_sub <- sample_n(dat_sub, n_subs)

# Run an MDS
mds <- metaMDS(dat_sub[, 3:18],
               distance = "gower",
               k = 2)

saveRDS(mds, "./dat/mds.RDS")

# Check the stress
mds$stress

# Colours for groups
col <- data.frame("col" = pal,
                  "group2" = unique(datS$group2),
                  stringsAsFactors = FALSE)
grp <- data.frame("group2" = dat_sub$group2,
                  "dx" = 1:nrow(dat_sub),
                  stringsAsFactors = F)
grp <- merge(x = grp, y = col, by = "group2", all.x = T)
grp <- grp[order(grp$dx), ]

# Plot
ordiplot(mds, display = "sites")

plot(mds, type = "n", display = "sites")
points(mds, col = grp$col, bg = grp$col, pch = 21)
plot(envfit(mds, dat_sub[ , 19:37]))

# Better plot
env_fit <- envfit(mds, dat_sub[ , 19:37])

# Get the nMDS coordinates
dat_sub$MDS1 <- mds$points[ , 1]
dat_sub$MDS2 <- mds$points[ , 2]

# Get the environmental vectors
vec.sp.df <- as.data.frame(env_fit$vectors$arrows*sqrt(env_fit$vectors$r))
vec.sp.df$env <- rownames(vec.sp.df)

# Mutliplier for scaling the vector output
mult <- vegan:::ordiArrowMul(env_fit)

pdf("./plots/mds.pdf", width = 6, height = 6, useDingbats = FALSE)
ggplot(data = dat_sub, aes(x = MDS1, y = MDS2, colour = as.factor(group2))) +
  geom_point(alpha = 0.7, size = 0.7) +
  scale_color_manual(values = pal, name = "Cluster") +
  geom_segment(data = vec.sp.df, aes(x = 0, xend = NMDS1*mult, y = 0, yend = NMDS2*mult),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey30", inherit.aes = FALSE) + 
  geom_text(data = vec.sp.df, aes(x = NMDS1*mult, y = NMDS2*mult, label = env), size = 2.5, inherit.aes = FALSE) +
  coord_fixed(1) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  annotate(geom = "text", x = - 0.30, y = -0.27, label = paste0("Stress: ", round(mds$stress, 2)), size = 2.5) +
  theme_bw() +
  theme(
    text = element_text(colour = "black",
                        size = 9),
    axis.text = element_text(colour = "black",
                             size = 8),
    axis.title = element_text(colour = "black",
                              size = 9),
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.background = element_rect(fill="transparent", colour=NA),
    legend.key        = element_rect(fill="transparent", colour=NA),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )
dev.off()

# Check Rsquared values for the environmental covariates
env_fit$vectors

env_fit_results <- data.frame("covariate" = names(env_fit$vectors$r),
                              "r2" = round(env_fit$vectors$r, 2),
                              "p" = env_fit$vectors$pvals)
env_fit_results <- env_fit_results[order(env_fit_results$r2, decreasing = T), ]
write.csv(env_fit_results, "./dat/mds_envfit_results.csv", row.names = F)
