# Script to score each cell of the map according to spatial criteria
# Dataset: BC parks iNaturalist (subset for amphibians)
# Criteria: observation density, last sample date, climate frequency, 
# estimated completeness, SDM disagreement

library(terra)
library(tidyterra)
library(tidyverse)
library(rgbif)
library(scales)

#theme_set(hrbrthemes::theme_ipsum_rc())
theme_set(theme_minimal())

# Load data --------------------------------------------------------------------

# load base grid
base5k = terra::rast("data/for_scoring/canada-basegrids/canada.base.5k.tiff")

# load observations (this is from iNat downloads)
dat = readRDS("data/for_scoring/iNaturalist/BCParks_obs_unique_with_BigTeams_flagged.RDS")

# load bc parks polygons
parks = terra::vect("data/for_scoring/bc-parks-pol/bc_parks.shp")

# make base grid with 1 and 0 to outline BC Parks
base5k = project(base5k, crs(parks))
parks5k = rasterize(parks, base5k, touches = TRUE) |> trim()
# terra::writeRaster(parks5k, "outputs/spatial-layers/parks5k.tif")

# make a park = 0 layer to add scores to!
parks0 = parks5k
parks0[parks5k == 1] <- 0


# Data for specific scores =====================================================

## Climate ----
# import mean annual temperature and precipitation
clim = c(
  terra::rast("data/for_scoring/climate/wc2.1_10m_bio_1.tif"),
  terra::rast("data/for_scoring/climate/wc2.1_10m_bio_12.tif")
)
clim <- project(clim, crs(parks))
clim <- mask(clim, parks)
clim <- trim(clim)
clim <- resample(clim, parks5k, method = "bilinear")
plot(clim)


# Process data =================================================================

# convert data into points layer
obs = dat |>
  dplyr::filter(captive_cultivated == "false",
         iconic_taxon_name == "Amphibia") |> ## FILTER TO MAKE THIS LIGHT FOR AN EXAMPLE
  dplyr::select(c(latitude, longitude, year, month, day,
           coordinates_obscured, iconic_taxon_name, 
           scientific_name, dataset)) |>
  dplyr::mutate_at(vars("longitude", "latitude"), as.numeric) |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")

# project obs to match the grid
obs = project(obs, crs(parks5k))
saveRDS(obs, "outputs/inatobs_amphibians_bcparks.rds")



# Spatial coverage scores ======================================================

# Density map ------------------------------------------------------------------

# save this for the completeness calculation (without sqrt)
map_obsdens = rasterize(obs, parks5k, fun = "length") 
map_obsdens[is.na(map_obsdens)] <- 0
map_obsdens[is.na(parks5k)] <- NA
saveRDS(map_obsdens, "outputs/inatdensity_amphibians_bcparks.rds")

map_obsdens = rasterize(obs, parks5k, fun = "length") |> sqrt()
map_obsdens = sqrt(map_obsdens)
map_obsdens[is.na(map_obsdens)] <- 0
map_obsdens[is.na(parks5k)] <- NA
plot(map_obsdens)

# Last update map --------------------------------------------------------------

map_lastupdate = rasterize(obs, parks5k, field = "year", fun = "max")
map_lastupdate[is.na(map_lastupdate)] <- min(obs$year, na.rm = TRUE)
map_lastupdate[is.na(parks5k)] <- NA
plot(map_lastupdate)

# Climate uniqueness score -----------------------------------------------------

plot(clim[[1]], clim[[2]])
clim_scale = scale(clim)
clim_df = data.frame(
  "temp" = values(clim_scale[[1]]),
  "prec" = values(clim_scale[[2]])
)
# Compute equal-width bin breaks
n_bins = 15
x_breaks <- seq(minmax(clim_scale[[1]])[1], minmax(clim_scale[[1]])[2], length.out = n_bins + 1)
y_breaks <- seq(minmax(clim_scale[[2]])[1], minmax(clim_scale[[2]])[2], length.out = n_bins + 1)
plot(clim_scale[[1]], clim_scale[[2]])
abline(h = y_breaks, col="red")
abline(v = x_breaks, col="red")

# Cut into bins and count 
BINS_range <- data.frame(
    "x_bin" = cut(values(clim_scale[[1]]), breaks = x_breaks, include.lowest = TRUE),
    "y_bin"= cut(values(clim_scale[[2]]), breaks = y_breaks, include.lowest = TRUE)
  ) |> na.omit()
BINS_cellcount = count(BINS_range, x_bin, y_bin)
BINS = left_join(BINS_range, BINS_cellcount)

clim_freq = clim_scale[[1]]
clim_freq[!is.na(clim_freq)] <- BINS$n
clim_freq[is.na(parks5k)] <- NA
plot(clim_freq)
# save the map
map_climfrequency = clim_freq

# plot to demo the approach
ggplot() +
  geom_point(data = clim_df, 
             aes(x = wc2.1_10m_bio_1, y = wc2.1_10m_bio_12),
             alpha = .4, col = "steelblue3") +
  geom_hline(yintercept = y_breaks, linewidth = .1) +
  geom_vline(xintercept = x_breaks, linewidth = .1) +
  theme_classic() 
ggsave("figures/clim_grid_bins.png", width = 5, height = 5)
ggplot() +
  geom_bin2d(data = clim_df, 
             aes(x = wc2.1_10m_bio_1, y = wc2.1_10m_bio_12), 
             bins = n_bins, alpha = 1,) +
  scale_fill_viridis_c(option = "plasma") +
  theme_classic() +
  labs(fill = "Count")
ggsave("figures/clim_grid_count.png", width = 5, height = 5)


# Completeness score -----------------------------------------------------------

# calculated in estimate_completeness.R
map_comp = readRDS("outputs/spatial-layers/scores_map_coverage.rds")
# invert so high coverage becomes low score
map_comp = 1-map_comp



# Optional purpose-oriented scores =============================================

# SDM discrepancy ----

cv = terra::rast("data/for_scoring/sdm-cv/current_cv.tif")
cv = project(cv, crs(map_obsdens))
cv_crop = crop(cv, map_obsdens)
cv_crop = resample(cv_crop, map_obsdens)
plot(cv_crop)
cv_crop[is.na(map_obsdens)] <- NA
map_modeldisagree = cv_crop

# make list of maps
maps = c("score_spatial_density" = map_obsdens, 
         "score_spatial_lastupdate" = map_lastupdate, 
         "score_spatial_climfreq" = map_climfrequency, 
         "score_spatial_completeness" = map_comp, 
         "score_spatial_modeldisagree" = map_modeldisagree) |> 
  rast()

# function to rescale maps
rescale_map = function(map) {
  vals = values(map) |> as.vector()
  rescaled_vals = (vals-min(vals, na.rm = T))/diff(range(vals, na.rm=T))
  values(map) <- rescaled_vals
  return(map)
}
maps_rescaled = sapply(maps, rescale_map) |> rast()
saveRDS(maps_rescaled, "outputs/spatial-layers/scores_spatial_maps.rds")


# Make a scores table ----------------------------------------------------------

scores = values(maps_rescaled) |> 
  as.data.frame() |> na.omit() 

# scale columns between 0 and 1
scores_scaled = scores |>
  # reverse scores so all follow 1 = high priority and  0 = low priority
  mutate(across(-c(score_spatial_completeness,
                   score_spatial_modeldisagree), ~ 1 - .x))

# Weight importance of the scores
scores_weights = rep(1, ncol(scores_scaled)) / ncol(scores_scaled) # equal for now!

# Tally total score per cell
scores_weighted = scores_scaled 
for(i in 1:ncol(scores_weighted)){
  scores_weighted[,i] = scores_weighted[,i]*scores_weights[i]
}
scores_m = as.matrix(scores_weighted)
scores_total = rowSums(scores_m)

# make a clean table
scores_df = scores_scaled
scores_df$total = scores_total
scores_df$cellID = cells(maps_rescaled)
write.csv(scores_df, "outputs/scores_spatial.csv")

# save output
saveRDS(list("scaled" = scores_scaled,
             "weighted" = scores_weighted,
             "total" = scores_total), "outputs/spatial-layers/scores_spatial.rds")
# basic sum (no weights) of all the layers
total = sum(maps_rescaled) |> sapply(rescale_map) |> rast()

# map the scores 
map_scores = total
map_scores[!is.na(map_scores)] <- scores_total
saveRDS(map_scores, "outputs/spatial-layers/scores_spatial_total.rds")

# pretty maps ===================================================================

ggplot() +
  tidyterra::geom_spatraster(data = map_scores) +
  scale_fill_viridis_c(option = "turbo", limits = c(0,1),
                       na.value = "transparent",name = "VOI") +
  labs(title = "Total spatial priority score")
ggsave("figures/scores_spatialpriorities.png", width = 9, height = 8)

p = list()
for(i in 1:nlyr(maps_rescaled)){
  p[[i]] = ggplot() +
    tidyterra::geom_spatraster(data = maps_rescaled[[i]]) +
    scale_fill_viridis_c(option = "turbo", limits = c(0,1), direction = 1,
                         na.value = "transparent", name = "VOI") +
    labs(title = gsub("score_spatial_", "", names(maps_rescaled)[i]))
  ggsave(plot = p[[i]], filename = paste0("figures/scores_spatialpriorities_",names(maps_rescaled)[i],".png"), width = 9, height = 8)
  
}
patchwork::wrap_plots(p, nrow = 3, ncol = 2)
ggsave("figures/scores_spatialpriorities_eachscore.png", width = 9, height = 8)

