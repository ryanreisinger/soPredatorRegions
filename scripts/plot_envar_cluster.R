# Plot environmental covariates

library(ggplot2)
library(tidyr)
library(dplyr)
library(pals)
library(colorspace)
library(gridExtra)

setwd("D:/RAATD/hotspot/")

#-----------------------------------------
## Read in gridded data with groups classified
dat <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

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

# Covariate names
covars <- names(dat)[19:37]

# Rather get them in order of correlation in the MDS
covar_r <- read.csv("./dat_out/mds_envfit_results.csv", stringsAsFactors = F)
covars <- covar_r$covariate

# Data to long format
dat_long <- pivot_longer(dat,
                              cols = all_of(covars),
                              names_to = "Covariate",
                              values_to = "Value")

# Remove NAs
dat_long <- dat_long[!is.na(dat_long$group2), ]

# Colours
## Creat a colour palette
k.hot <- length(which(grepl("hot_", unique(dat$group2))))
k.cold <- length(which(grepl("cold_", unique(dat$group2))))
pal <- c(sequential_hcl(k.cold, palette = "Light Grays"), glasbey(k.hot))

# Reorder factors to match clustering
dat_long$group2 <- factor(dat_long$group2, levels = c("hot_13",
                                                      "hot_15",
                                                      "hot_16",
                                                      "cold_14",
                                                      "hot_17",
                                                      "hot_09",
                                                      "hot_10",
                                                      "hot_03",
                                                      "hot_08",
                                                      "hot_04",
                                                      "cold_01",
                                                      "cold_02",
                                                      "hot_05",
                                                      "hot_06",
                                                      "hot_07",
                                                      "hot_11",
                                                      "hot_12"))

# Reorder palette
# Right order
ordered_clusters <- c("hot_13",
                      "hot_15",
                      "hot_16",
                      "cold_14",
                      "hot_17",
                      "hot_09",
                      "hot_10",
                      "hot_03",
                      "hot_08",
                      "hot_04",
                      "cold_01",
                      "cold_02",
                      "hot_05",
                      "hot_06",
                      "hot_07",
                      "hot_11",
                      "hot_12")

# Match
dx <- vector(length = length(ordered_clusters))
for (i in 1:length(ordered_clusters)) {
dx[i] <- which(unique(as.character(dat_long$group2)) == ordered_clusters[i])
}

# Reorder the palette
pal <- pal[c(1,2, 4:14, 3, 15:17)] # reorder to match colours to original order 
pal <-pal[dx] # reorder by factor order

# Add supercluster
dat_long <- dat_long %>%
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

dat_long <- dat_long %>%
  mutate(super_group_ab = case_when(group2 == "hot_13" ~ "ANT",
                                 group2 == "hot_15" ~ "ANT",
                                 group2 == "hot_16" ~ "ANT",
                                 group2 == "cold_14" ~ "ANT",
                                 group2 == "hot_17" ~ "ANT",
                                 group2 == "hot_09" ~ "SCO",
                                 group2 == "hot_10" ~ "SCO",
                                 group2 == "hot_03" ~ "DIS",
                                 group2 == "hot_08" ~ "DIS",
                                 group2 == "hot_04" ~ "DIS",
                                 group2 == "cold_01" ~ "DIS",
                                 group2 == "cold_02" ~ "DIS",
                                 group2 == "hot_05" ~ "SUB",
                                 group2 == "hot_06" ~ "SUB",
                                 group2 == "hot_07" ~ "SUB",
                                 group2 == "hot_11" ~ "SUB",
                                 group2 == "hot_12" ~ "SUB"))

# Plot in a list
plist <- list()

for (i in covars) {
  
  ylim1 <- boxplot.stats(dat_long[dat_long$Covariate == i,]$Value)$stats[c(1, 5)] # Create better limits for boxplots
  
  p1 <- ggplot(dat_long[dat_long$Covariate == i, ], aes(y = Value, x = group2, col = group2, fill = group2)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_segment(x = 5.5, xend = 5.5, y = ylim1[1], yend = ylim1[2]*1.2) +
    # geom_segment(x = 7.5, xend = 7.5, y = ylim1[1], yend = ylim1[2]*1.2) +
    # geom_segment(x = 12.5, xend = 12.5, y = ylim1[1], yend = ylim1[2]*1.2) +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal) +
    guides(colour = "none", fill = "none") +
    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
    facet_grid(super_group_ab ~ Covariate, scales = "free") +
    labs(y = "", x = "Cluster") +
    # coord_cartesian(ylim = ylim1*1.05) +
    coord_flip(ylim = ylim1*1.05) +
    theme_ryan()
  
  plist[i] <- list(p1)
  
}

ggsave(file = "./plots/covariates_by_cluster.pdf", arrangeGrob(grobs = plist, ncol = 4),
       width = 20, height = 29, units = "cm") 
