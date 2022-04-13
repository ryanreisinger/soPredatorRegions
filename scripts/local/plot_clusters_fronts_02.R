## Plot hotspots and MPAs together

setwd("D:/soPredatorRegions")

library(sf)
library(raster)
library(SOmap)
library(colorspace)
library(pals)
library(sp)
library(ncdf4)
library(dplyr)

#-----------------------------------------
## Map stuff
wgs_84 <- "+proj=longlat +datum=WGS84 +no_defs"
study_area <- st_bbox(c(xmin = -180, xmax = +180, ymax = -40, ymin = -80), crs = st_crs(4326))
raster_temp <- raster::raster(res = 0.1,
                              xmn = -180, xmx = +180, ymx = -40, ymn = -80,
                              crs =  "+proj=longlat +datum=WGS84")

# Initialize SOmap
SOmap()

#-----------------------------------------
# Colony locations
prj <-  SOcrs() # get current projection
cloc <- "./dat_in/colonies/colonies.csv"

# Get colony locations and re-project
cloc <- read.csv(cloc)
cloc <- filter(cloc, lat <= -40) %>%
   dplyr::select("lon", "lat") %>%
   SpatialPoints(., proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")) 
cloc <- sp::spTransform(cloc, CRS(prj))

#---------------------------------------------------
# Fronts from Park & Durand 2019
# Get Southern Ocean fronts from Park & Durand 2019
# https://doi.org/10.17882/59800

frnts <- nc_open("./dat_in/fronts/62985.nc")
# NB <- data.frame(
#   "lat" = ncvar_get(frnts, "LatNB"),
#   "lon" = ncvar_get(frnts, "LonNB"),
#   "name" = "NB"
# )
SAF <- data.frame(
   "lat" = ncvar_get(frnts, "LatSAF"),
   "lon" = ncvar_get(frnts, "LonSAF"),
   "name" = "SAF"
)
PF <- data.frame(
   "lat" = ncvar_get(frnts, "LatPF"),
   "lon" = ncvar_get(frnts, "LonPF"),
   "name" = "PF"
)
SACCF <- data.frame(
   "lat" = ncvar_get(frnts, "LatSACCF"),
   "lon" = ncvar_get(frnts, "LonSACCF"),
   "name" = "SACCF"
)
# SB <- data.frame(
#   "lat" = ncvar_get(frnts, "LatSB"),
#   "lon" = ncvar_get(frnts, "LonSB"),
#   "name" = "SB"
# )
nc_close(frnts)

# frnts <- rbind(NB, SAF, PF, SACCF, SB)
frnts <- rbind(SAF, PF, SACCF)
frnts <- frnts[complete.cases(frnts),]
frnts <- SpatialPoints(coords = frnts[,c("lon", "lat")],
                       proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
frnts.polar <- spTransform(frnts, CRS(prj))

#-----------------------------------------
## Get ice
## Median sea ice extent for 1981 - 2010

## Max in September:
maxice <- st_read(dsn = "./dat_in/ice/median_extent_S_09_1981-2010_polyline_v3.0.shp", layer = "median_extent_S_09_1981-2010_polyline_v3.0")

## Min in March
minice <- st_read(dsn = "./dat_in/ice/median_extent_S_03_1981-2010_polyline_v3.0.shp", layer = "median_extent_S_03_1981-2010_polyline_v3.0")

#-----------------------------------------
## Read in gridded data with groups classified
dat <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")

#-----------------------------------------
## Hot and coldspots
rq <- rasterFromXYZ(dat[ , c("x", "y", "MEAN")])
crs(rq) <- "+proj=longlat +datum=WGS84 +no_defs"
# Calculate top and bottom quantiles
lower_q <- quantile(rq, probs = c(0.25, 0.50, 0.75))[1]
median_q <- quantile(rq, probs = c(0.25, 0.50, 0.75))[2]
upper_q <- quantile(rq, probs = c(0.25, 0.50, 0.75))[3]
# Create contours corresponding with the upper and lower quantiles
contour_lines_lower <- rasterToContour(rq, levels = c(lower_q))
contour_lines_upper <- rasterToContour(rq, levels = c(upper_q))
crs(contour_lines_lower) <- "+proj=longlat +datum=WGS84 +no_defs"
crs(contour_lines_upper) <- "+proj=longlat +datum=WGS84 +no_defs"

SOmap() # Call SOmap to initialize projection
rq_proj <- projectRaster(rq, crs = SOcrs())

#-----------------------------------------
## Plot clusters
# # Create the raster
r <- rasterFromXYZ(dat[ , c("x", "y", "group")])
crs(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

## Creat a colour palette
# k.hot <- length(which(grepl("hot_", unique(dat$group2))))
# k.cold <- length(which(grepl("cold_", unique(dat$group2))))
# pal <- c(sequential_hcl(k.cold, palette = "Light Grays"), glasbey(k.hot))
this_k <- 17
pal <- glasbey(this_k)

# Create a table to store names and colours
rat <- unique(dat[,c("group", "group2")])
rat <- rat[!is.na(rat$group), ]
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

#rat <- rat[c(13, 15, 16, 14, 17, 9, 10, 3, 8, 4, 1, 2, 5, 6, 7, 11, 12),]

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
r_proj <- projectRaster(r, crs = SOcrs(), method = "ngb")

#--------------------------------
# Maps
#--------------------------------
# Clusters
tiff("plots/cluster_map_fronts.tiff", height = 150/0.66666666666, width = 190/0.66666666666, units = "mm", res = 600)

SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)

SOplot(r_proj, col = rat$pal, legend = FALSE)

legend(x = 6000000,
       y = 3000000,
       title = "Antarctic",
       title.adj = 0,
       inset = -0.08,
       bty = "n",
       # legend = full_cluster_names,
       legend = ordered_clusters[1:5],
       fill = rat$pal[dx][1:5])

legend(x = 6000000,
       y = 1000000,
       title = "Scotia Arc",
       title.adj = 0,
       inset = -0.08,
       bty = "n",
       # legend = full_cluster_names,
       legend = ordered_clusters[6:7],
       fill = rat$pal[dx][6:7])

legend(x = 6000000,
       y = -500000,
       title = "Distant\nSubantarctic",
       title.adj = 0,
       inset = -0.08,
       bty = "n",
       # legend = full_cluster_names,
       legend = ordered_clusters[8:12],
       fill = rat$pal[dx][8:12])

legend(x = 6000000,
       y = -2500000,
       title = "Subantarctic",
       title.adj = 0,
       inset = -0.08,
       bty = "n",
       # legend = full_cluster_names,
       legend = ordered_clusters[13:17],
       fill = rat$pal[dx][13:17])
# SOplot(mpa_designated, col = NA, border = "white", lwd = 3.5, legend = FALSE)
# SOplot(mpa_proposed, col = NA, border = "white", lwd = 3.5, legend = FALSE)
# 
# SOplot(mpa_designated, col = NA, border = "#009988", lwd = 1.0, legend = FALSE)
# SOplot(mpa_proposed, col = NA, border = "#EE7733", lwd = 1.0, legend = FALSE)

# fronts
SOplot(frnts.polar,
       col = "black",
       cex = 0.1,
       add = T)

# sea-ice
SOplot(minice,
       col = "white",
       lwd = 2,
       add = T)

SOplot(maxice,
       col = "white",
       lwd = 2,
       add = T)

SOplot(cloc,
       pch = 16,
       col = "white",
       cex = 1,
       add = T)

SOplot(cloc,
       pch = 16,
       col = "black",
       cex = 0.8,
       add = T)

# legend(x = 'right', legend = c("Designated", "Proposed"), fill = c("#009988", "#EE7733"))
dev.off()

#--------------------------------
# Hot and coldspots
# Clusters
tiff("plots/hotspot_map_fronts.tiff", height = 70/0.66666666666, width = 100/0.66666666666, units = "mm", res = 600)

SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)

SOplot(rq_proj, col = ocean.balance(125), legend = T,
       legend.args = list(text = 'Mean\nhabitat\nimportance', title.adj = 0))

# legend(x = 'right',
#        title = "Mean habitat\nimportance",
#        title.adj = 0.2,
#        legend = seq(0, 100, 20),
#        inset = -0.08,
#        bty = "n")

# SOplot(mpa_designated, col = NA, border = "white", lwd = 3.5, legend = FALSE)
# SOplot(mpa_proposed, col = NA, border = "white", lwd = 3.5, legend = FALSE)
# 
# SOplot(mpa_designated, col = NA, border = "#009988", lwd = 1.0, legend = FALSE)
# SOplot(mpa_proposed, col = NA, border = "#EE7733", lwd = 1.0, legend = FALSE)

# contours
# SOplot(contour_lines_lower,
#        lwd = 2,
#        col = ocean.balance(5)[2],
#        add = T)
# SOplot(contour_lines_upper,
#        lwd = 2,
#        col = ocean.balance(5)[4],
#        add = T)

# fronts
SOplot(frnts.polar,
       col = "black",
       cex = 0.1,
       add = T)

# sea-ice
SOplot(minice,
       col = "white",
       lwd = 2,
       add = T)

SOplot(maxice,
       col = "white",
       lwd = 2,
       add = T)

# colonies
SOplot(cloc,
       pch = 16,
       col = "white",
       cex = 1,
       add = T)

SOplot(cloc,
       pch = 16,
       col = "black",
       cex = 0.8,
       add = T)

# legend(x = 'right', legend = c("Designated", "Proposed"), fill = c("#009988", "#EE7733"))
dev.off()
