# Script to plot and explore BC parks data

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
  terra::rast("for_scoring/climate/wc2.1_10m/wc2.1_10m_bio_1.tif"),
  terra::rast("for_scoring/climate/wc2.1_10m/wc2.1_10m_bio_12.tif")
)
clim <- project(clim, crs(parks))
clim <- mask(clim, parks)
clim <- trim(clim)
clim <- resample(clim, parks5k, method = "bilinear")
plot(clim)


# Process data =================================================================

# convert data into points layer
obs = dat |>
  filter(captive_cultivated == "false",
         iconic_taxon_name == "Amphibia") |> ## FILTER TO MAKE THIS LIGHT FOR AN EXAMPLE
  select(c(latitude, longitude, year, month, day,
           coordinates_obscured, iconic_taxon_name, 
           scientific_name, dataset)) |>
  mutate_at(vars("longitude", "latitude"), as.numeric) |>
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
    "x_bin" = cut(values(clim_scale[[1]]), breaks = x_breaks, include.lowest = FALSE),
    "y_bin"= cut(values(clim_scale[[2]]), breaks = y_breaks, include.lowest = FALSE)
  ) |> na.omit()
BINS_cellcount = count(BINS_range, x_bin, y_bin)
BINS = left_join(BINS_range, BINS_cellcount)
clim_freq = clim_scale[[1]]
values(clim_freq) <- BINS$n
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

# Taxonomic coverage scores ====================================================



# Optional purpose-oriented scores =============================================

# SDM discrepancy ----

# sdm_cv = terra::rast("~/Documents/GitHub/BTG-analyse-the-gap/data/heavy/noah-uncertainty/cv_8meth.tif")



# Make a scores table ----------------------------------------------------------

scores = data.frame(
  "cellID" = cells(parks5k),
  "score_spatial_density" = values(map_obsdens) |> na.omit() |> as.vector(),
  "score_spatial_lastupdate" = values(map_lastupdate) |> na.omit() |> as.vector(),
  "score_spatial_climfreq" = values(map_climfrequency) |> na.omit() |> as.vector(),
  "score_spatial_completeness" = values(map_comp) |> na.omit() |> as.vector()
)

# scale columns between 0 and 1
scores_scaled = scores |>
  # scale between 0 and 1
  mutate(across(2:score_spatial_climfreq, rescale)) |>
  # reverse scores so 1 = high priority and  0 = low priority
  mutate(across(2:score_spatial_climfreq, ~ 1 - .x))

# Weight importance of the scores
scores_weights = rep(1, c(ncol(scores_scaled)-1)) / c(ncol(scores_scaled)-1)

# Tally total score per cell
scores_weighted = scores_scaled 
for(i in 2:ncol(scores_weighted)){
  scores_weighted[,i] = scores_weighted[,i]*scores_weights[i-1]
}
scores_m = as.matrix(scores_weighted[,-1])
scores_total = rowSums(scores_m)

# make a clean table
scores_df = scores_scaled
scores_df$total = scores_total
write.csv(scores_df, "outputs/scores_spatial.csv")

# save output
saveRDS(list("scaled" = scores_scaled,
             "weighted" = scores_weighted,
             "total" = scores_total), "outputs/spatial-layers/scores_spatial.rds")


# map the scores 
map_scores = parks5k
map_scores[parks5k==1] <- scores_total

# pretty maps ===================================================================

ggplot() +
  tidyterra::geom_spatraster(data = map_scores) +
  scale_fill_viridis_c(option = "turbo", limits = c(0,1),
                       na.value = "transparent",name = "VOI") +
  labs(title = "Total spatial priority score")
ggsave("figures/scores_spatialpriorities.png", width = 9, height = 8)

# save maps of the scaled scores
map_spatial = c("obsdensity" = map_scores, 
                "lastupdate" = map_scores, 
                "climfreq" = map_scores,
                "completeness" = map_scores)
map_spatial$obsdensity[parks5k==1] <- scores_scaled$score_spatial_density
map_spatial$lastupdate[parks5k==1] <- scores_scaled$score_spatial_lastupdate
map_spatial$climfreq[parks5k==1] <- scores_scaled$score_spatial_climfreq
map_spatial$completeness[parks5k==1] <- scores_scaled$score_spatial_completeness
map_spatial$total = map_scores
saveRDS(map_spatial, "outputs/spatial-layers/scores_map_spatial.rds")

p = list()
for(i in 1:length(map_spatial)){
  p[[i]] = ggplot() +
    tidyterra::geom_spatraster(data = map_spatial[[i]]) +
    scale_fill_viridis_c(option = "plasma", limits = c(0,1), direction = -1,
                         na.value = "transparent", name = "VOI") +
    labs(title = names(map_spatial)[i])
  ggsave(plot = p[[i]], filename = paste0("figures/scores_spatialpriorities_",names(map_spatial)[i],".png"), width = 9, height = 8)
  
}
patchwork::wrap_plots(p[1:4], nrow = 2, ncol = 2)
ggsave("figures/scores_spatialpriorities_eachscore.png", width = 9, height = 8)