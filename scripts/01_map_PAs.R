library(tidyterra)
library(ggplot2)
library(colorspace)
library(patchwork)
library(terra)

# load maps
# can be downloaded from: https://mcgill.sharepoint.com/:f:/r/sites/QuantitativeBiodiversityLab_Group/Shared%20Documents/BC%20Parks/Spatial%20data/BC%20endemism%20and%20PA%20maps%20(Eckert%20et%20al%202023)?csf=1&web=1&e=qBwflX
m1 = terra::rast("data/isaac_prioritymaps/BC maps/endemism.bc.tif")
m2 = terra::rast("data/isaac_prioritymaps/BC maps/ranks.bc.tif")
m3 = terra::rast("data/isaac_prioritymaps/BC maps/pas.bc.tif")

# bc polygon
bc = terra::vect("data/polygons/bc_polygon.shp")
bc <- terra::project(bc, "EPSG:3005")
bc = sf::st_as_sf(bc)

# project
m1 <- terra::project(m1, "EPSG:3005")
# clip out the blank cells
m1[m1 == 0] <- NA
# trim whitespace around map
m1 = terra::trim(m1)
# plot!
(A = ggplot() +
  geom_sf(data = bc, col = "black", fill = "white") +
  geom_spatraster(data = m1) +
  colorspace::scale_fill_continuous_sequential(palette = "Batlow",
                                           na.value = "transparent", 
                                           breaks = c(5e-7, 1e-6, 1.5e-6, 2.5e-6),
                                           rev = FALSE, 
                                           trans = "sqrt",
                                           name = "Endemism") +
    theme(legend.position = "right"))

# mask the existing PAs
m2_maskpas = m2
m2_maskpas[m2 >= 0.807] <- NA
m2_maskpas <- m2_maskpas*100
m2_maskpas <- terra::project(m2_maskpas, "EPSG:3005")
terra::writeRaster(m2_maskpas, "outputs/spatial-layers/PA_priorityranks.tif")
m2_maskpas = terra::trim(m2_maskpas)

(B = ggplot() +
    geom_sf(data = bc, col = "black", fill = "white") +
  geom_spatraster(data = m2_maskpas) +
  colorspace::scale_fill_binned_sequential(palette = "Batlow",
                                           breaks = c(25, 50, 75),
                                           na.value = "transparent", 
                                           rev = FALSE,
                                           name = "Priority\nrank (%)") +
  theme(legend.position = "right"))
ggsave("figures/map_PA_priorityranks.png", width = 9.62, height = 7.8)

# mapview::mapview(m2_maskpas)
A / B + plot_annotation(tag_levels = "a")
ggsave("figures/map_BC_endemism_PApriorities.png", width = 6, height = 8)
