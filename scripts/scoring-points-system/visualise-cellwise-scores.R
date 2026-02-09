# Visualiser to show the scores on a zoomable map

library(raster)  
library(mapview)

# load maps
maps = terra::rast("outputs/spatial-layers/scores_spatial_maps.rds")
score = terra::rast("outputs/spatial-layers/scores_spatial_total.rds")

# set global mapview options
mapview::mapviewOptions(na.color = "transparent",
                        raster.palette = viridis::turbo,basemaps = "CartoDB.Positron")

# create raster stack
maps_stack = stack(c(maps, score))
mapview(maps_stack)

# save as html
mapshot(mapview(maps_stack), url = "docs/scoremaps_inatamphibians.html")
