setwd("~/RAATD_01/RAATDhotspot")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(pals)
library(colorspace)

# Get cluster data
clu <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

clu$CHLA <- log10(clu$CHLA)

# Data to long format
dat_long <- pivot_longer(clu,
                         cols = c("CHLA",
                                  "CURR",
                                  "DEPTH",
                                  "DEPTHg",
                                  "dSHELF",
                                  "EKE",
                                  "ICE",
                                  "ICEA",
                                  "ICEsd",
                                  "SAL",
                                  "SHFLUX",
                                  "SHFLUXsd",
                                  "SSHa",
                                  "SSHsd",
                                  "SST",
                                  "SSTg",
                                  "VMIX",
                                  "VMIXsd",
                                  "WIND"),
                         names_to = "Variable",
                         values_to = "Value")



k.hot <- length(unique(clu[clu$hotspot == "Y", "group2"]))
k.cold <- length(unique(clu[clu$hotspot == "N", "group2"]))

k.hot <- 14
k.cold <- 3
pal <- c(sequential_hcl(k.cold, palette = "Light Grays"), glasbey(k.hot))

dat_long <- dat_long[!is.na(dat_long$group2), ]

if (FALSE) {
pdf("./plots/envar_joyplots.pdf", paper = "a4")
ggplot(data = dat_long, aes(x = Value, y = group2, colour = group2, fill = group2)) +
  geom_density_ridges2() +
  scale_y_discrete() +
  scale_fill_manual(values = pal, name = "Cluster") +
  scale_colour_manual(values = pal, name = "Cluster") +
  facet_wrap(~ Variable, scales = "free_x", ncol = 4) +
  labs(x = "Covariate value", y = "Cluster") +
  theme_bw() +
  theme_ridges()
dev.off()
}

if (TRUE) {
tiff("./plots/envar_scatterjitter.tiff", width = 16, height = 8, units = "in", res = 300)
ggplot(data = dat_long, aes(x = Value, y = group2, colour = group2, fill = group2)) +
  geom_jitter(alpha = 0.01) +
  scale_y_discrete() +
  scale_fill_manual(values = pal, name = "Cluster") +
  scale_colour_manual(values = pal, name = "Cluster") +
  facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
  labs(x = "Covariate value", y = "Cluster") +
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
}

if (FALSE) {
pdf("./plots/envar_violin.pdf", paper = "a4")
ggplot(data = dat_long, aes(x = Value, y = group2, fill = group2)) +
  geom_violin() +
  scale_y_discrete() +
  scale_fill_manual(values = pal, name = "Cluster") +
  facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
  labs(x = "Covariate value", y = "Cluster") +
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
}