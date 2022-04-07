## Plot MDS envar correlations

library(ggplot2)

setwd("~/soPredatorRegions")

dat <- read.csv("./dat_out/mds_envfit_results.csv",
                stringsAsFactors = F)



dat$name <- c("Sea ice concentration",
              "Accessibility through sea ice",
              "Standard deviation of sea ice concentration",
              "Surface wind speed",
              "Sea surface temperature",
              "Sea surface temperature gradient",
              "Geostrophic current velocity",
              "Sea surface height standard deviation",
              "Eddy kinetic energy",
              "Clorophyll a concentration",
              "Standard deviation of vertical velocity",
              "Vertical velocity",
              "Sea surfache height anomaly",
              "Depth",
              "Distance to shelf",
              "Standard deviation of surface heat flux",
              "Surface heat flux",
              "Depth gradient",
              "Salinity difference")

dat$cov_name <- paste(dat$name, " | ", dat$order)

dat$covariate <- factor(dat$covariate, levels = rev(dat$covariate))
dat$cov_name <- factor(dat$cov_name, levels = rev(dat$cov_name))

pdf("./plots/mds_correlation.pdf", height = 4, width = 5)
ggplot(data = dat, aes(x = r2, y = cov_name)) +
  geom_bar(stat = "identity") +
  labs(x = "Correlation with nMDS (R^2)",
       y = "Covariate") +
  theme_bw()
dev.off()
