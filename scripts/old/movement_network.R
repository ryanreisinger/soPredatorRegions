# Network representation of animal movements

setwd("/mnt/home/ryan/RAATD_01/RAATDhotspot")

library(raster)
library(dplyr)
library(igraph)
library(colorspace) # sequential_hcl
library(pals) # glasbey
library(ggplot2)
library(SOmap)
library(sf)

##----------------------------------------
# Create graph
gr <- select(tst, cluster, cluster_to, n)
gr <- gr[complete.cases(gr), ]
gr <- rename(gr, from = cluster, to = cluster_to, weight = n)

gr <- graph.data.frame(gr)

# Add colours
V(gr)$color <- pal

plot(gr, layout = layout.fruchterman.reingold, edge.width = log(E(gr)$weight),
     vertex.label.family = "sans-serif")

# Some centrality measures
degree(gr) # Degree centrality
closeness(gr) # Closeness centrality
betweenness(gr) # Betweenness centrality
evcent(gr)$vector # Eigenvector centrality

##----------------------------------------
# Plotting in ggplot and soMap

# Calculate cluster geographic centers
foo <- function(cluster_name, frame) {
  prj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  hold <- st_sfc(st_multipoint(as.matrix(frame[ , c("x", "y")])))
  st_crs(hold) <- 4326
  hold <- st_transform(hold, prj)
  frame$x <- st_coordinates(hold)[,1]
  frame$y <- st_coordinates(hold)[,2]
  rm(hold)
  this <- filter(frame, group2 == cluster_name)
  out <- data.frame("group2" = cluster_name,
                    "x" = mean(this$x, na.rm = T),
                    "y" = mean(this$y, na.rm = T))
  return(out)
}

cluster_locs <- do.call("rbind", lapply(unique(clu$group2), FUN = foo, frame = clu))
cluster_locs <- cluster_locs[complete.cases(cluster_locs),]

# Initialize a SOmap
# Create basemap
basemap <- SOmap(trim = -40,
                 bathy_legend = FALSE,
                 border_col = c("white", "white"),
                 border_width = 0.01,
                 straight = TRUE,
                 graticules = TRUE)

raster::values(basemap$bathy[[1]]$plotargs$x)[which(raster::values(basemap$bathy[[1]]$plotargs$x) < 0)] <- 0
basemap$bathy_legend <- NULL

# Create ggplot object
gg <- plot(SOgg(basemap))

# Convert to polar coordinates
# prj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# hold <- st_sfc(st_multipoint(as.matrix(cluster_locs[ , c("x", "y")])))
# st_crs(hold) <- 4326
# hold <- st_transform(hold, prj)
# cluster_locs$x <- st_coordinates(hold)[,1]
# cluster_locs$y <- st_coordinates(hold)[,2]
# rm(hold)

cluster_locs_from <- rename(cluster_locs,
                            "cluster" = group2,
                            "x_from" = x,
                            "y_from" = y)
cluster_locs_to <- rename(cluster_locs,
                          "cluster_to" = group2,
                          "x_to" = x,
                          "y_to" = y)

geo_net <- left_join(tst, cluster_locs_from)
geo_net <- left_join(geo_net, cluster_locs_to)

# Remove incomplete rows
geo_net <- filter(geo_net, complete.cases(geo_net))

# Summarize for plotting
geo_net <- mutate(geo_net, "cluster_link" = paste0(cluster, "-", cluster_to)) %>% 
  group_by(., cluster_link) %>% 
  summarise(.,
            "cluster" = unique(cluster),
            "cluster_to" = unique(cluster_to),
            "weight" = sum(n),
            "x_from" = mean(x_from),
            "y_from" = mean(y_from),
            "x_to" = mean(x_to),
            "y_to" = mean(y_to)) %>%
  ungroup(.)

# Plot in ggplot
ggplot(data = geo_net[geo_net$cluster != geo_net$cluster_to, ], aes(x = x_from, y = y_from)) +
  geom_point() +
  # coord_map("ortho", orientation = c(-90, 0, 0)) +
  # coord_quickmap() +
  geom_curve(
             aes(x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to,
                 # alpha = log(n),
                 colour = log(n)),
             size = 1,
             # arrow = arrow(length = unit(0.2, "cm")),
             inherit.aes = F)

# In SOmap
gg +
  geom_point(data = geo_net[geo_net$cluster != geo_net$cluster_to, ], aes(x = x_from, y = y_from)) +
  geom_curve(data = geo_net[geo_net$cluster != geo_net$cluster_to, ],
             aes(x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to,
                 alpha = weight),
                 arrow = arrow(length = unit(0.2, "cm")),
                 inherit.aes = F) +
  geom_point(data = geo_net[geo_net$cluster != geo_net$cluster_to, ], aes(x = x_from, y = y_from, colour = cluster))
  
