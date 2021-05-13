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

library(tidyr)
library(gridExtra)
library(dplyr)
library(gridExtra)

#-------------------------------------
## Get data
dat <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

#-------------------------------------
## 3. Plot species compositions

# Go to long data.frame
datL <- gather(data = dat,
               ADPE, ANFS, ANPE, BBAL, CRAS, DMSA, EMPE, GHAL, HUWH, KIPE,
               LMSA, MAPE.ROPE, SOES, WAAL, WESE, WHCP,
               key = "Species", value = "Importance")

# Plot

# Histogram

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

datL <- datL[!is.na(datL$group2), ]

ordered_clusters <- c("antarctic_hot_13",
                               "antarctic_hot_15",
                               "antarctic_hot_16",
                               "antarctic_cold_14",
                               "antarctic_hot_17",
                               "scotia_arc_hot_09",
                               "scotia_arc_hot_10",
                               "subantarctic_distant_hot_03",
                               "subantarctic_distant_hot_08",
                               "subantarctic_distant_hot_04",
                               "subantarctic_distant_cold_01",
                               "subantarctic_distant_cold_02",
                               "subantarctic_hot_05",
                               "subantarctic_hot_06",
                               "subantarctic_hot_07",
                               "subantarctic_hot_11",
                               "subantarctic_hot_12")

plist <- list()

for (i in ordered_clusters) {

p1 <- ggplot(datL[datL$super_group_cluster == i, ], aes(x = reorder(Species, Importance, FUN = mean), y = Importance, fill = Guild, colour = Guild)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colrs) +
  scale_colour_manual(values = colrs) +
  guides(colour = "none", fill = "none") +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) }) +
  facet_wrap(~ super_group + group2, ncol = 1, scales = "free_x") +
  labs(x = "", y = "Habitat importance") +
  coord_flip() +
  theme_ryan()

plist[i] <- list(p1)

}

ggsave(file = "./plots/speciesCompositionClusters_quartile_02.pdf", arrangeGrob(grobs = plist, ncol = 4),
       width = 20, height = 29, units = "cm") 


# Heatplot of mean habitat importance for each species X cluster
datsum <- group_by(datL, Species, group2)
datsum <- summarise_at(datsum, .vars = "Importance", .funs = mean)

lbls <- data.frame("group2" = ordered_clusters,
                   "order" = 1:17)
datsum <- merge(x = datsum, y = lbls[, c("group2", "order")], by = "group2", all.x = T)
datsum$group_order <- paste0(formatC(datsum$order, width=2, flag="0"), "_", datsum$group)

datsum$Species <- factor(datsum$Species, levels = c("ANFS",
                                                    "CRAS",
                                                    "SOES",
                                                    "WESE",
                                                    "HUWH",
                                                    "ADPE",
                                                    "EMPE",
                                                    "KIPE",
                                                    "MAPE.ROPE",
                                                    "ANPE",
                                                    "BBAL",
                                                    "DMSA",
                                                    "GHAL",
                                                    "LMSA",
                                                    "WAAL",
                                                    "WHCP"))

pdf("./plots/heatplotSpecies_quartile_02.pdf", width = 7, height = 6)
ggplot(data = datsum, aes(x = group_order, y = Species, fill = Importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean\nhabitat\nimportance") +
  labs(x = "Cluster", y = "Species") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(labels = lbls$group2)
dev.off()
