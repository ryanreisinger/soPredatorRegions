## Prepare World MPA Database

setwd("D:/soPredatorRegions")

library(sf)
library(fasterize)
library(raster)
library(SOmap)

## Received by email from Beth Pike,
## 15 April 2020
# Updated March 2021
## beth.pike@marine-conservation.org
dat <- st_read(dsn = "./dat_in/mpa_data/mpatlas_20201223_clean/mpatlas_20201223_clean.shp", layer = "mpatlas_20201223_clean")

# Reference area
study_area <- st_bbox(c(xmin = -180, xmax = +180, ymax = -40, ymin = -90), crs = st_crs(4326))

# Raster template
raster_temp <- raster::raster(res = 0.1,
                              xmn = -180, xmx = +180, ymx = -40, ymn = -80,
                              crs =  "+proj=longlat +datum=WGS84")

# Crop
# dat <- st_crop(st_buffer(dat, dist = 0), study_area) # 0 buffer to resolve geometry problem
dat <- st_crop(dat, study_area)

# Remove IWC MPA (mpa_id 9258) and 'Indian Ocean' mpa (mpa_id 68813327)
dat <- dat[dat$mpa_id != 9258, ]
dat <- dat[dat$mpa_id != 68813327, ]

# Plot the valid MPAs
plot(dat["is_mpa"])

# Plot status
plot(dat["status"])

# Plot fishing
plot(dat["no_take"])

# Subsets
# MPAs
mpa <- dat["is_mpa"]
# Designated and proposed MPAs
mpa_designated <- dat[dat$status == "Designated","status"]
mpa_proposed <- dat[dat$status == "Proposed","status"]
# No-take MPAs
no_take_mpa <- dat[dat$"no_take" == "All" & dat$status == "Designated", "no_take"]

# Write
saveRDS(mpa, "./dat_out/mpa.RDS")
writeRaster(fasterize(st_cast(mpa), raster_temp), "./dat_out/mpa_raster.grd", format = "raster", overwrite = T)

saveRDS(mpa_designated, "./dat_out/mpa_designated.RDS")
writeRaster(fasterize(st_cast(mpa_designated), raster_temp), "./dat_out/mpa_designated_raster.grd", format = "raster", overwrite = T)

saveRDS(mpa_proposed, "./dat_out/mpa_proposed.RDS")
writeRaster(fasterize(st_cast(mpa_proposed), raster_temp), "./dat_out/mpa_proposed_raster.grd", format = "raster", overwrite = T)

saveRDS(no_take_mpa, "./dat_out/no_take_mpa.RDS")
writeRaster(fasterize(st_cast(no_take_mpa), raster_temp), "./dat_out/no_take_mpa_raster.grd", format = "raster", overwrite = T)

# Get CCAMLR boundary
ccamlr <- st_read(".\\dat_in\\ccamlr_area\\asd-shapefile-EPSG102020\\asd-shapefile-EPSG102020.shp")
ccamlr <- st_union(ccamlr)

# Maps
par(mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
tiff("plots/mpaMap.tiff", height = 70/0.66666666666, width = 130/0.66666666666, units = "mm", res = 600)
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)
SOplot(mpa_designated, col = "#009988", legend = FALSE)
SOplot(mpa_proposed, col = "#EE7733", legend = FALSE)
SOplot(no_take_mpa, col = "#0077BB", legend = FALSE)
SOplot(ccamlr, col = NA, border = "black", lwd = 1, legend = FALSE)
legend(x='right',
       title = "MPA",
       title.adj = 0.2,
       legend = c("Designated\nmixed-use", "Designated\nno-take", "Proposed"),
       fill = c("#009988", "#0077BB", "#EE7733"),
       bty = "n",
       inset = 0,
       y.intersp = 1.25)
dev.off()