setwd("~/RAATD_01/RAATDhotspot")

library(indicspecies)
library(vegan)
library(dplyr)

# Get cluster data
clu <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

# Remove NAs
temp <- clu[complete.cases(clu[,c(3:17, 40)]), ]

# Sample
temp <- temp %>%
  group_by(., group2) %>%
  sample_n(., 500)

# With package indicspecies
multipatt(x = temp[,3:17], cluster = temp$group2)

#SIMPER from vegan
simper(comm = temp[,3:17], group = temp$group2)
