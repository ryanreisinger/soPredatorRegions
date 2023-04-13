library(rgdal)
library(rgeos)
library(raster)
library(smoothr)
library(sf)
library(sp)

foo <- readRDS("./dat_out/gridded_with_clusters_quartile.RDS")
foo <- rasterFromXYZ(xyz = foo[,c("x", "y", "group2")],
                               res = c(0.1, 0.1),
                               crs = "+proj=longlat +datum=WGS84 +no_defs")
foo_poly <- rasterToPolygons(foo, dissolve = T)

#plot(foo_poly)

foo_poly <- smoothr::drop_crumbs(foo_poly,
                                 threshold = 1)

foo_poly <- smoothr::smooth(foo_poly)

# Write raster
writeRaster(foo,
            "./dat_out/ecoregions/ecoregions.asc",
            format = "ascii")

# Write shape file
st_write(obj = st_as_sf(foo_poly),
         dsn = "./dat_out/ecoregions/ecoregions.shp")

#writeOGR(out.df, dsn = './dat_out/ecoregions/', layer = 'ecoregions', driver = "ESRI Shapefile")
