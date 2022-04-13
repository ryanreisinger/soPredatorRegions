# Plot environmental covariates

library(ggplot2)
library(tidyr)
library(dplyr)
library(pals)
library(colorspace)
library(gridExtra)

setwd("D:/soPredatorRegions")

#-----------------------------------------
## Read in gridded data with groups classified
dat <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

#### custom theme
theme_ryan <- function () { 
  theme_bw(base_size=8, base_family="") %+replace% 
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
# k.hot <- length(which(grepl("hot_", unique(dat$group2))))
# k.cold <- length(which(grepl("cold_", unique(dat$group2))))
# pal <- c(sequential_hcl(k.cold, palette = "Light Grays"), glasbey(k.hot))
pal <- glasbey(17)

# Reorder factors to match clustering
dat_long$group2 <- factor(dat_long$group2, levels = c("13",
                                                      "15",
                                                      "16",
                                                      "14",
                                                      "17",
                                                      "09",
                                                      "10",
                                                      "03",
                                                      "08",
                                                      "04",
                                                      "01",
                                                      "02",
                                                      "05",
                                                      "06",
                                                      "07",
                                                      "11",
                                                      "12"))

# Reorder palette
# Right order
ordered_clusters <- c("13",
                      "15",
                      "16",
                      "14",
                      "17",
                      "09",
                      "10",
                      "03",
                      "08",
                      "04",
                      "01",
                      "02",
                      "05",
                      "06",
                      "07",
                      "11",
                      "12")

# Match
dx <- vector(length = length(ordered_clusters))
for (i in 1:length(ordered_clusters)) {
dx[i] <- which(unique(as.character(dat_long$group2)) == ordered_clusters[i])
}

# Reorder the palette
# pal <- pal[c(1,2, 4:14, 3, 15:17)] # reorder to match colours to original order 
pal <-pal[dx] # reorder by factor order

# Add supercluster
dat_long <- dat_long %>%
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

dat_long <- dat_long %>%
  mutate(super_group_ab = case_when(group2 == "13" ~ "ANT",
                                 group2 == "15" ~ "ANT",
                                 group2 == "16" ~ "ANT",
                                 group2 == "14" ~ "ANT",
                                 group2 == "17" ~ "ANT",
                                 group2 == "09" ~ "SCO",
                                 group2 == "10" ~ "SCO",
                                 group2 == "03" ~ "DIS",
                                 group2 == "08" ~ "DIS",
                                 group2 == "04" ~ "DIS",
                                 group2 == "01" ~ "DIS",
                                 group2 == "02" ~ "DIS",
                                 group2 == "05" ~ "SUB",
                                 group2 == "06" ~ "SUB",
                                 group2 == "07" ~ "SUB",
                                 group2 == "11" ~ "SUB",
                                 group2 == "12" ~ "SUB"))

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
