## MPA calculations

library(raster)
library(sf)
library(fasterize)

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("D:/RAATD/hotspot/")

## Read in gridded data with groups classified
dat <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

# Reorder clusters
# Reorder factors to match clustering
dat$group2 <- factor(dat$group2, levels = c("hot_13",
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

# Add info on the add-hoc superclusters from hierarchical clustering
dat <- dat %>%
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

dat <- unite(data = dat, col = "super_group_cluster", c("super_group", "group2"), sep = "_", remove = FALSE)

# Stuff
wgs_84 <- "+proj=longlat +datum=WGS84 +no_defs"

# Raster template
raster_temp <- rasterFromXYZ(select(dat, x, y, MEAN), crs = wgs_84)

## Get pre-processed MPA polygon and join with the grid
mpa <- readRDS("./dat_out/mpa.RDS")

# Doesn't run...
# df <- st_as_sf(x = select(dat, x, y, MEAN),                         
#                coords = c("x", "y"),
#                crs = wgs_84)
# 
# cell_in_mpa <- st_join(df, mpa, join = st_within)

# This isn't working either...
# tst <- raster::extract(raster_temp, mpa)

## Calculate area per cluster
dat_raster <- rasterFromXYZ(select(dat, x, y, MEAN), crs = wgs_84)
dat$cell_area <- raster::extract(raster::area(dat_raster), select(dat, x, y))
dat[is.na(dat$MEAN), ]$cell_area <- NA # Set land to NA area

# Extract whether MPA
mpa_raster <- raster("./dat_out/mpa_raster.grd")
dat$is_mpa <- raster::extract(mpa_raster, select(dat, x, y)) # Extract whether MPA from raster
dat[!is.na(dat$is_mpa), ]$is_mpa <- dat[!is.na(dat$is_mpa), ]$cell_area # Transfer area value if MPA

# Designated?
mpa_designated_raster <- raster("./dat_out/mpa_designated_raster.grd")
dat$is_mpa_designated <- raster::extract(mpa_designated_raster, select(dat, x, y)) # Extract whether designated MPA from raster
dat[!is.na(dat$is_mpa_designated), ]$is_mpa_designated <- dat[!is.na(dat$is_mpa_designated), ]$cell_area # Transfer area value if designated MPA

# Proposed?
mpa_proposed_raster <- raster("./dat_out/mpa_proposed_raster.grd")
dat$is_mpa_proposed <- raster::extract(mpa_proposed_raster, select(dat, x, y)) # Extract whether proposed MPA from raster
dat[!is.na(dat$is_mpa_proposed), ]$is_mpa_proposed <- dat[!is.na(dat$is_mpa_proposed), ]$cell_area # Transfer area value if proposed MPA

# No take?
no_take_raster <- raster("./dat_out/no_take_mpa_raster.grd")
dat$is_no_take <- raster::extract(no_take_raster, select(dat, x, y)) # Extract whether no take MPA from raster
dat[!is.na(dat$is_no_take), ]$is_no_take <- dat[!is.na(dat$is_no_take), ]$cell_area # Transfer area value if no take MPA

# Only land
dat <- dat[dat$super_group_cluster != "NA_NA", ]

# Summarise
# Add up areas
dat_summary <- dat %>%
  group_by(., super_group_cluster) %>%
  summarise(., total_area = sum(cell_area, na.rm = T),
            mpa_area = sum(is_mpa, na.rm = T),
            mpa_designated_area = sum(is_mpa_designated, na.rm = T),
            mpa_proposed_area = sum(is_mpa_proposed, na.rm = T),
            no_take_area = sum(is_no_take, na.rm = T))

# Percentages

# No take
dat_summary$no_take_proportion <- dat_summary$no_take_area/dat_summary$total_area
dat_summary$no_take_proportion <- dat_summary$no_take_proportion * 100

# MPA
dat_summary$mpa_proportion <- dat_summary$mpa_area/dat_summary$total_area
dat_summary$mpa_proportion <- dat_summary$mpa_proportion * 100

# Designated MPA
dat_summary$mpa_designated_proportion <- dat_summary$mpa_designated_area/dat_summary$total_area
dat_summary$mpa_designated_proportion <- dat_summary$mpa_designated_proportion * 100

# Subtract no-take from MPA
dat_summary$mpa_designated_proportion <- dat_summary$mpa_designated_proportion - dat_summary$no_take_proportion

# Proposed MPA
dat_summary$mpa_proposed_proportion <- dat_summary$mpa_proposed_area/dat_summary$total_area
dat_summary$mpa_proposed_proportion <- dat_summary$mpa_proposed_proportion * 100

# Area and percentage outside mpa
dat_summary$outside_mpa_area <- dat_summary$total_area - dat_summary$mpa_area
dat_summary$outside_mpa_proportion <- dat_summary$outside_mpa_area/dat_summary$total_area
dat_summary$outside_mpa_proportion <- dat_summary$outside_mpa_proportion * 100

# Update designated MPA area to exclude no-take
dat_summary$mpa_designated_area <- dat_summary$mpa_designated_area - dat_summary$no_take_area

# Go to long and plot
# Area
dat_summary_long_area <- pivot_longer(dat_summary,
                                 c(outside_mpa_area,
                                   no_take_area,
                                 mpa_designated_area,
                                 mpa_proposed_area),
                                 names_to = "area_type",
                                 values_to = "area")

dat_summary_long_area$area <- dat_summary_long_area$area / 1000000

# Reorder factor
dat_summary_long_area$area_type <- factor(dat_summary_long_area$area_type, levels = c("no_take_area", "mpa_designated_area",
                                                                                      "mpa_proposed_area", "outside_mpa_area"))

# Reorder clusters
dat_summary_long_area$super_group_cluster <- factor(dat_summary_long_area$super_group_cluster,
                                                   levels = c("antarctic_hot_13",
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
                                                                 "subantarctic_hot_12"))

pdf("./plots/mpa_area_km.pdf", width = 6, height = 5, useDingbats = FALSE)
plot_a <- ggplot(data = dat_summary_long_area,
       aes(y = area, x = super_group_cluster, fill = area_type)) +
  geom_col(position = position_stack(reverse = T)) +
  scale_fill_manual(values = c("#0077BB", "#009988", "#EE7733", "darkgrey"),
                    name = "",
                    labels = c("Designated no-take MPA", "Designated mixed-use MPA", "Proposed MPA", "Outside MPA")) +
  labs(x = "Cluster", y = "Area (million km^2)") +
  coord_flip() +
  theme_bw()

print(plot_a)

dev.off()

# Percentages
dat_summary_long_proportion <- pivot_longer(dat_summary,
                                      c(outside_mpa_proportion,
                                        mpa_designated_proportion,
                                        no_take_proportion,
                                        mpa_proposed_proportion),
                                      names_to = "area_type",
                                      values_to = "percentage")

# Reorder factor
dat_summary_long_proportion$area_type <- factor(dat_summary_long_proportion$area_type, levels = c("no_take_proportion", "mpa_designated_proportion",
                                                                                      "mpa_proposed_proportion", "outside_mpa_proportion"))

# Reorder clusters
dat_summary_long_proportion$super_group_cluster <- factor(dat_summary_long_proportion$super_group_cluster,
                                                    levels = c("antarctic_hot_13",
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
                                                               "subantarctic_hot_12"))


pdf("./plots/mpa_area_percentage.pdf", width = 6, height = 5, useDingbats = FALSE)
plot_b <- ggplot(data = dat_summary_long_proportion,
       aes(y = percentage, x = super_group_cluster, fill = area_type)) +
  geom_col(position = position_stack(reverse = T)) +
  scale_fill_manual(values = c("#0077BB", "#009988", "#EE7733", "darkgrey"),
                    name = "",
                    labels = c("Designated no-take MPA", "Designated mixed-use MPA", "Proposed MPA", "Outside MPA")) +
  labs(x = "Cluster", y = "Area (%)") +
  geom_hline(yintercept = c(10, 30)) +
  coord_flip() +
  theme_bw()
print(plot_b)
dev.off()

#------------------------------------------------
## Cumulative by super-cluster
dat_super_summary <- dat %>%
  group_by(., super_group) %>%
  summarise(., total_area = sum(cell_area, na.rm = T),
            mpa_area = sum(is_mpa, na.rm = T),
            mpa_designated_area = sum(is_mpa_designated, na.rm = T),
            mpa_proposed_area = sum(is_mpa_proposed, na.rm = T),
            no_take_area = sum(is_no_take, na.rm = T))

# Percentages
dat_super_summary$mpa_proportion <- dat_super_summary$mpa_area/dat_super_summary$total_area
dat_super_summary$mpa_proportion <- dat_super_summary$mpa_proportion * 100

dat_super_summary$mpa_designated_proportion <- dat_super_summary$mpa_designated_area/dat_super_summary$total_area
dat_super_summary$mpa_designated_proportion <- dat_super_summary$mpa_designated_proportion * 100

dat_super_summary$mpa_proposed_proportion <- dat_super_summary$mpa_proposed_area/dat_super_summary$total_area
dat_super_summary$mpa_proposed_proportion <- dat_super_summary$mpa_proposed_proportion * 100

dat_super_summary$no_take_proportion <- dat_super_summary$no_take_area/dat_super_summary$total_area
dat_super_summary$no_take_proportion <- dat_super_summary$no_take_proportion * 100

dat_super_summary$mpa_designated_proportion <- dat_super_summary$mpa_designated_proportion - dat_super_summary$no_take_proportion

# Area and percentage outside mpa
dat_super_summary$outside_mpa_area <- dat_super_summary$total_area - dat_super_summary$mpa_area
dat_super_summary$outside_mpa_proportion <- dat_super_summary$outside_mpa_area/dat_super_summary$total_area
dat_super_summary$outside_mpa_proportion <- dat_super_summary$outside_mpa_proportion * 100

# Go to long and plot
dat_super_summary_long_proportion <- pivot_longer(dat_super_summary,
                                            c(outside_mpa_proportion,
                                              mpa_designated_proportion,
                                              mpa_proposed_proportion,
                                              no_take_proportion),
                                            names_to = "area_type",
                                            values_to = "percentage")

# Reorder factor
dat_super_summary_long_proportion$area_type <- factor(dat_super_summary_long_proportion$area_type, levels = c("no_take_proportion", "mpa_designated_proportion",
                                                                                                              "mpa_proposed_proportion", "outside_mpa_proportion"))


pdf("./plots/mpa_area_super_percentage.pdf", width = 6, height = 3, useDingbats = FALSE)
plot_c <- ggplot(data = dat_super_summary_long_proportion,
       aes(y = percentage, x = super_group, fill = area_type)) +
  geom_col(position = position_stack(reverse = T)) +
  scale_fill_manual(values = c("#0077BB", "#009988", "#EE7733", "darkgrey"),
                    name = "",
                    labels = c("Designated no-take MPA", "Designated mixed-use MPA", "Proposed MPA", "Outside MPA")) +
  labs(x = "Cluster", y = "Area (%)") +
  geom_hline(yintercept = c(10, 30)) +
  coord_flip() +
  theme_bw()
print(plot_c)
dev.off()

# Plot all three together
library(patchwork)
plots <- (plot_a | plot_b) / (plot_c | guide_area())
plots <- plots + plot_annotation(tag_levels = 'a') + plot_layout(guides = 'collect')

pdf("./plots/mpa_cover_combined.pdf", width = 12, height = 7, useDingbats = FALSE)
print(plots)
dev.off()
