# Prepare and extract prey layers

#--------------------------------------
setwd("~/soPredatorRegions")

library(raster)
library(ggplot2)
library(ggridges)
library(tidyr)
library(pals)
library(colorspace)
library(tidyselect)
library(vegan)
library(dplyr)

#--------------------------------------
# Get cluster data
clu <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

# Remove NAs
# temp <- clu[complete.cases(clu[,c(3:17, 40)]), ]

# Colours for plotting
k.hot <- length(unique(clu[clu$hotspot == "Y", "group2"]))
k.cold <- length(unique(clu[clu$hotspot == "N", "group2"]))

k.hot <- 14
k.cold <- 3
pal <- c(sequential_hcl(k.cold, palette = "Light Grays"), glasbey(k.hot))

#--------------------------------------
# Habitat importance function
suitability_to_percentiles <- function(x, N=21) {
  ## expect that x is a raster containing predicted usage probability (suitability multiplied by availability)
  cell_areas <- raster::values(area(x))
  vals <- raster::values(x)
  total_area <- sum(cell_areas, na.rm = T)
  tst <- seq(min(vals,na.rm=TRUE),max(vals,na.rm=TRUE),length.out=N)
  ## calculate the percentage of area corresponding to each of these tst values
  s2p <- function(z) sum(cell_areas[which(vals<=z)])/total_area*100
  arp <- vapply(tst,s2p,FUN.VALUE=1)
  values(x) <- approx(tst,arp,vals)$y
  x
}

#--------------------------------------
# Cephalopods
# tst <- raster("~/../shared/prey/cephalopods/species_richness.nc")
# plot(tst, axes = F)

# List of files
fls <- list.files("~/../shared/prey/cephalopods/", full.names = TRUE)
fls <- fls[-which(grepl("species_richness.nc", fls))] # Drop the species richness layer

# Corresponding species names
ceph_species <- list.files("~/../shared/prey/cephalopods/", full.names = FALSE)
ceph_species <- ceph_species[-which(grepl("species_richness.nc", ceph_species))] # Drop the species richness layer
ceph_species <- substr(ceph_species,1,nchar(ceph_species)-3) # Drop file extension
ceph_species <- paste0("CE_", ceph_species)

# Replace species names with updated taxonomy from Cherel 2020.
ceph_speces[ceph_species == "CE_moroteuthis_robsoni"] <- "CE_onykia_robsoni"
ceph_speces[ceph_species == "CE_kondakovia_longimana"] <- "CE_moroteuthopsis_longimana"
ceph_speces[ceph_species == "CE_moroteuthis_ingens"] <- "CE_moroteuthopsis_ingens"
ceph_speces[ceph_species == "CE_loligo_gahi"] <- "CE_doryteuthis_gahi"

cephalopods <- stack()

for (i in 1:length(fls)) {
# this_file <- raster("~/../shared/prey/cephalopods/alluroteuthis_antarcticus.nc")
  this_file <- raster(fls[i])
  # this_name <- species[i]
  this_file <- suitability_to_percentiles(this_file)
  cephalopods <- stack(cephalopods, this_file)
}

# Species names as layer names
names(cephalopods) <- ceph_species

# Extract
ceph_df <- raster::extract(cephalopods, clu[,c("x", "y")])

# Bind
clu_prey <- cbind(clu, ceph_df)

# Plot
# Data to long format
ceph_dat_long <- pivot_longer(clu_prey,
                         cols = all_of(ceph_species),
                         names_to = "Species",
                         values_to = "Suitability")

ceph_dat_long$taxon <- "Cephalopod"

ceph_dat_long <- ceph_dat_long[!is.na(ceph_dat_long$group2), ]

pdf("./plots/joyplots_ceph.pdf", paper = "a4r")
ggplot(data = ceph_dat_long, aes(x = Suitability, y = Species, colour = group2, fill = group2)) +
  geom_density_ridges2() +
  scale_y_discrete() +
  scale_fill_manual(values = pal, name = "Cluster") +
  scale_colour_manual(values = pal, name = "Cluster") +
  facet_wrap(~ group2, scales = "free_x", ncol = 6) +
  labs(x = "Habitat suitability", y = "Species") +
  theme_bw() +
  theme_ridges()
dev.off()

#--------------------------------------
# Mycthophids
fls <- list.files("./dat_in/prey/myctophids/", full.names = TRUE)
fls <- fls[-which(grepl("readme.txt", fls))] # Drop the readme

# tst <- raster(fls[1])
# plot(tst, axes = F)

# Corresponding species names
myct_species <- list.files("./dat_in/prey/myctophids/", full.names = FALSE)
myct_species <- myct_species[-which(grepl("readme.txt", myct_species))] # Drop the readme
myct_species <- substr(myct_species,1,nchar(myct_species)-8) # Drop last part and extension
myct_species <- paste0("MY_", myct_species)

myctophids <- stack()

for (i in 1:length(fls)) {
  this_file <- raster(fls[i])
  this_file <- suitability_to_percentiles(this_file)
  myctophids <- stack(myctophids, this_file)
}

# Species names as layer names
names(myctophids) <- myct_species

# Extract
myct_df <- raster::extract(myctophids, clu[,c("x", "y")])

# Bind
clu_prey <- cbind(clu, ceph_df, myct_df)

# Plot
# Data to long format
myct_dat_long <- pivot_longer(clu_prey,
                         cols = all_of(myct_species),
                         names_to = "Species",
                         values_to = "Suitability")

myct_dat_long$taxon <- "Myctophid"

myct_dat_long <- myct_dat_long[!is.na(myct_dat_long$group2), ]

pdf("./plots/joyplots_myct.pdf", paper = "a4r")
ggplot(data = myct_dat_long, aes(x = Suitability, y = Species, colour = group2, fill = group2)) +
  geom_density_ridges2() +
  scale_y_discrete() +
  scale_fill_manual(values = pal, name = "Cluster") +
  scale_colour_manual(values = pal, name = "Cluster") +
  facet_wrap(~ group2, scales = "free_x", ncol = 6) +
  labs(x = "Habitat suitability", y = "Species") +
  theme_bw() +
  theme_ridges()
dev.off()

#--------------------------------------
# Euphasiids
fls <- list.files("./dat_in/prey/euphasiids/", full.names = TRUE)
fls <- fls[which(grepl("vars-full-pred-boot", fls))] # Keep only the final predictions

# Corresponding species names
euph_species <- list.files("./dat_in/prey/euphasiids/", full.names = FALSE)
euph_species <- euph_species[which(grepl("vars-full-pred-boot", fls))]
euph_species <- substr(euph_species,34,nchar(euph_species)-4) # Drop extension
euph_species <- paste0("EU_", euph_species)

euphasiids <- stack()

for (i in 1:length(fls)) {
  print(i)
  this_file <- read.csv(fls[i], stringsAsFactors = F)
  this_raster <- rasterFromXYZ(this_file[,1:3])
  crs(this_raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  this_raster <- projectRaster(this_raster, cephalopods, method = "bilinear")
  this_species <- euph_species[i]
  names(this_raster) <- this_species
  this_raster <- suitability_to_percentiles(this_raster)
  euphasiids <- stack(euphasiids, this_raster)
}

# Extract
euph_df <- raster::extract(euphasiids, clu[,c("x", "y")])

# Bind
clu_prey <- cbind(clu, ceph_df, myct_df, euph_df)

# Plot
# Data to long format
euph_dat_long <- pivot_longer(clu_prey,
                              cols = all_of(euph_species),
                              names_to = "Species",
                              values_to = "Suitability")

euph_dat_long$taxon <- "Euphasiid"

euph_dat_long <- euph_dat_long[!is.na(euph_dat_long$group2), ]

pdf("./plots/joyplots_euph.pdf", paper = "a4r")
ggplot(data = euph_dat_long, aes(x = Suitability, y = Species, colour = group2, fill = group2)) +
  geom_density_ridges2() +
  scale_y_discrete() +
  scale_fill_manual(values = pal, name = "Cluster") +
  scale_colour_manual(values = pal, name = "Cluster") +
  facet_wrap(~ group2, scales = "free_x", ncol = 6) +
  labs(x = "Habitat suitability", y = "Species") +
  theme_bw() +
  theme_ridges()
dev.off()

#-------------------------------------
# Test differences in prey among clusters

# Subsample per cluster
n_subs <- 500

dat_sub <- clu_prey[complete.cases(clu_prey[ , c(40, 42:70)]), ] # Use only if envars are complete
dat_sub <- dat_sub %>% group_by(group2)
dat_sub <- sample_n(dat_sub, n_subs)

# PERMANOVA
permanova <- adonis(dat_sub[ , 42:70] ~ dat_sub$group2,
                    method = "gower")

# Check the dispersion - signficant differences may be due to
# different within-group variation (dispersion) rather than different
# mean values. See ?adonis
mod <- betadisper(d = vegdist(dat_sub[ , 42:70], method = "gower"),
                  group = dat_sub$group2)

# saveRDS(mod, "./dat_out/betadisper_output.RDS")

## Perform ANOVA on dispersion
anova(mod)

## Or Perform a permutation test for F
permutest(mod, pairwise = TRUE, permutations = 999)

## Or even Tukey's Honest Significant Differences
mod.HSD <- TukeyHSD(mod)
png("./plots/betadisper_tukey.png", width = 5, height = 5, units = "in", res = 300)
plot(mod.HSD)
dev.off()

#--------------------------------------

# Random forest classification
library(caret)
library(ranger)
library(MLmetrics)

rf_dat <- dplyr::select(clu_prey, "group2", all_of(c(ceph_species, myct_species, euph_species)))
rf_dat <- dplyr::filter(rf_dat, complete.cases(rf_dat))

# Select training data, preserving class distribution
trainIndex <- createDataPartition(as.factor(rf_dat$group2), p = 0.5, 
                                  list = FALSE, 
                                  times = 1)

rf_dat_train <- rf_dat[trainIndex, ]

## Fit in ranger
rf <- ranger(y = as.factor(rf_dat_train$group2),
             x = rf_dat_train[,-1],
             importance = "impurity_corrected",
             splitrule = "gini")

# Fit in caret
# method = "ranger"
# Tuning parameters:
#   
#   mtry (#Randomly Selected Predictors)
#   splitrule (Splitting Rule)
#   min.node.size (Minimal Node Size)

fitControl <- trainControl(method = "none",
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)
metric <- "Accuracy"

rf_caret <- train(group2 ~ .,
                  data = rf_dat_train, 
                  method = "ranger",
                  importance = "impurity_corrected", # ranger specific importance measure
                  trControl = fitControl, 
                  verbose = FALSE, 
                  ## Only a single model can be passed to the
                  ## function when no resampling is used:
                  tuneGrid = data.frame(mtry = 5,
                                        splitrule = "gini",
                                        min.node.size = 1),
                  metric = metric)

rf_caret

#--------------------------------------
# Heat map of prey mean habitat importance by cluster

library(dplyr)
library(patchwork)
datL <- pivot_longer(clu_prey,
                     cols = c(ceph_species, myct_species, euph_species),
                     names_to = "Species",
                     values_to = "Importance")

datsum <- group_by(datL, Species, group2)
datsum <- summarise_at(datsum, .vars = "Importance", .funs = mean, na.rm = T)

datsum <- datsum[!is.na(datsum$group2), ]

# Labels & ordering
# Order of clusters for matching dendrogram
dorder <- c(13, 15, 16, 14, 17, 9, 10, 3, 8, 4, 1, 2, 5, 6, 7, 11, 12)
lbls <- unique(clu[, c("group", "group2")])
lbls <- lbls[complete.cases(lbls), ]
lbls$order <- NA
lbls[dorder, ]$order <- 1:17

datsum <- merge(x = datsum, y = lbls[, c("group2", "order")], by = "group2", all.x = T)
datsum$group_order <- paste0(formatC(datsum$order, width=2, flag="0"), "_", datsum$group)

# Order prey species according to random forest variable importance
col_index <- varImp(rf_caret)$importance %>% 
  mutate(names=row.names(.)) %>%
  arrange(-Overall)

datsum$Species <- factor(datsum$Species, levels = rev(col_index$names))

# Create heatmap
p1 <- ggplot(data = datsum, aes(x = group_order, y = Species, fill = Importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean\nhabitat\nimportance") +
  scale_y_discrete(position = "right") +
  labs(x = "Cluster", y = "") +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        legend.position = "left") +
  scale_x_discrete(labels = lbls[dorder, ]$group2)

# Histogram to go next to the heatmap
col_index$names <- factor(col_index$names, levels = rev(col_index$names)) # Order factors

# Create plot
p2 <- ggplot(data = col_index, aes(x = Overall, y = names)) +
  geom_point() +
  labs(x = "Relative variable importance", y = "") +
  theme(axis.text.y = element_text(hjust = 0.5))

# Plot together
pdf("./plots/prey_heat.pdf", width = 10, height = 6)
p1 + p2
dev.off()