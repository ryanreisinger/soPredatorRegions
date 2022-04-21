## K-means cluster analysis to assign spatial eco-regions based on species
## habitat importance scores

## Also produce ordination plots (nMDS and PCA) and fit environment vectors

## Ryan R Reisinger, Ben Raymond

setwd("D:/soPredatorRegions")

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
guilds <- read.csv("./dat_in_public/guilds.csv", stringsAsFactors = F)
datL <- merge(datL, guilds, by = "Species", all.x = T)

# Add supercluster labels
datL <- datL %>%
  mutate(super_group = case_when(group2 == "13" ~ "antarctic",
                                 group2 == "15" ~ "antarctic",
                                 group2 == "16" ~ "antarctic",
                                 group2 == "14" ~ "antarctic",
                                 group2 == "17" ~ "antarctic",
                                 group2 == "09" ~ "scotia_arc",
                                 group2 == "10" ~ "scotia_arc",
                                 group2 == "03" ~ "subantarctic_distant",
                                 group2 == "08" ~ "subantarctic_distant",
                                 group2 == "04" ~ "subantarctic_distant",
                                 group2 == "01" ~ "subantarctic_distant",
                                 group2 == "02" ~ "subantarctic_distant",
                                 group2 == "05" ~ "subantarctic",
                                 group2 == "06" ~ "subantarctic",
                                 group2 == "07" ~ "subantarctic",
                                 group2 == "11" ~ "subantarctic",
                                 group2 == "12" ~ "subantarctic"))

datL <- unite(data = datL, col = "super_group_cluster", c("super_group", "group2"), sep = "_", remove = FALSE)

datL <- datL[!is.na(datL$group2), ]

ordered_clusters <- c("antarctic_13",
                               "antarctic_15",
                               "antarctic_16",
                               "antarctic_14",
                               "antarctic_17",
                               "scotia_arc_09",
                               "scotia_arc_10",
                               "subantarctic_distant_03",
                               "subantarctic_distant_08",
                               "subantarctic_distant_04",
                               "subantarctic_distant_01",
                               "subantarctic_distant_02",
                               "subantarctic_05",
                               "subantarctic_06",
                               "subantarctic_07",
                               "subantarctic_11",
                               "subantarctic_12")

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
datsum <- group_by(datL, Species, super_group_cluster)
datsum <- summarise_at(datsum, .vars = "Importance", .funs = mean)

datsum$super_group_cluster <- factor(datsum$super_group_cluster, levels = ordered_clusters)

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
ggplot(data = datsum, aes(x = super_group_cluster, y = Species, fill = Importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean\nhabitat\nimportance") +
  labs(x = "Cluster", y = "Species") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
