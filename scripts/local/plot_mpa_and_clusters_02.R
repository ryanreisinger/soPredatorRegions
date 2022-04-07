## Plot hotspots and MPAs together

setwd("D:/soPredatorRegions")

library(sf)
library(raster)
library(SOmap)
library(colorspace)
library(pals)

## Map stuff
wgs_84 <- "+proj=longlat +datum=WGS84 +no_defs"
study_area <- st_bbox(c(xmin = -180, xmax = +180, ymax = -40, ymin = -80), crs = st_crs(4326))
raster_temp <- raster::raster(res = 0.1,
                              xmn = -180, xmx = +180, ymx = -40, ymn = -80,
                              crs =  "+proj=longlat +datum=WGS84")
SOmap()

## Read in gridded data with groups classified
dat <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

# # Create the raster
r <- rasterFromXYZ(dat[ , c("x", "y", "group")])
crs(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

## Creat a colour palette
# k.hot <- length(which(grepl("hot_", unique(dat$group2))))
# k.cold <- length(which(grepl("cold_", unique(dat$group2))))
# pal <- c(sequential_hcl(k.cold, palette = "Light Grays"), glasbey(k.hot))
this_k <- length(unique(dat$group2))
pal <- glasbey(this_k)

# Create a table to store names and colours
rat <- unique(dat[,c("group", "group2")])
rat <- rat[!is.na(rat$group), ]
# rat <- rat[order(rat$group2), ]
rat$pal <- pal
rat$newgroup <- 1:nrow(rat)

r <- reclassify(r, rat[,c("group", "newgroup")])

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
   dx[i] <- which(rat$group2 == ordered_clusters[i])
}

full_cluster_names <- c("Antarctic 13",
                        "Antarctic 15",
                        "Antarctic 16",
                        "Antarctic 14",
                        "Antarctic 17",
                        "Scoti Arc 09",
                        "Scotia Arc 10",
                        "Distant Subantarctic 03",
                        "Distant Subantarctic 08",
                        "Distant Subantarctic 04",
                        "Distant Subantarctic 01",
                        "Distant Subantarctic 02",
                        "Subantarctic 05",
                        "Subantarctic 06",
                        "Subantarctic 07",
                        "Subantarctic 11",
                        "Subantarctic 12")

# Need to reproject using nearest neighbour lookup, otherwise
# integer values become numerical and color palette is messed up
SOmap() # Call SOmap to initialize projection
r_proj <- projectRaster(r, crs = SOcrs(), method = "ngb")


## Get MPA information
mpa_designated <- readRDS("./dat_out/mpa_designated.RDS")
mpa_proposed <- readRDS("./dat_out/mpa_proposed.RDS")
no_take_mpa <- readRDS("./dat_out/no_take_mpa.RDS")

# Get ATCS area
# ccamlr <- st_read(".\\dat_in\\ccamlr_area\\asd-shapefile-WGS84\\asd-shapefile-WGS84.shp")
# ccamlr <- st_union(ccamlr)

# Try projected to avoid geometry problems
ccamlr <- st_read(".\\dat_in\\ccamlr_area\\asd-shapefile-EPSG102020\\asd-shapefile-EPSG102020.shp")
ccamlr <- st_union(ccamlr)

# Maps
tiff("plots/mpa_cluster_map.tiff", height = 6, width = 10, units = "in", res = 300)

SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)

SOplot(r_proj, col = rat$pal, legend = FALSE)

# legend(x = 'right',
#        title = "Cluster",
#        title.adj = 0.2,
#        bty = "n",
#        legend = rat$group2,
#        fill = rat$pal)

legend(x = 6000000,
       y = 4800000,
       title = "Antarctic",
       title.adj = 0,
       inset = -0.08,
       bty = "n",
       # legend = full_cluster_names,
       legend = ordered_clusters[1:5],
       fill = rat$pal[dx][1:5])

legend(x = 6000000,
       y = 2000000,
       title = "Scotia Arc",
       title.adj = 0,
       inset = -0.08,
       bty = "n",
       # legend = full_cluster_names,
       legend = ordered_clusters[6:7],
       fill = rat$pal[dx][6:7])

legend(x = 6000000,
       y = -25000,
       title = "Distant\nSubantarctic",
       title.adj = 0,
       inset = -0.08,
       bty = "n",
       # legend = full_cluster_names,
       legend = ordered_clusters[8:12],
       fill = rat$pal[dx][8:12])

legend(x = 6000000,
       y = -2800000,
       title = "Subantarctic",
       title.adj = 0,
       inset = -0.08,
       bty = "n",
       # legend = full_cluster_names,
       legend = ordered_clusters[13:17],
       fill = rat$pal[dx][13:17])

SOplot(mpa_designated, col = NA, border = "white", lwd = 4, legend = FALSE)
SOplot(mpa_proposed, col = NA, border = "white", lwd = 4, legend = FALSE)

SOplot(mpa_designated, col = NA, border = "#009988", lwd = 2, legend = FALSE)
SOplot(mpa_proposed, col = NA, border = "#EE7733", lwd = 2, legend = FALSE)

SOplot(ccamlr, col = NA, border = "black", lwd = 1, legend = FALSE)

# legend(x = 'right', legend = c("Designated", "Proposed"), fill = c("#009988", "#EE7733"))
dev.off()