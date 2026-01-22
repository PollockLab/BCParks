# Script to plot and explore BC parks data

library(terra)
library(tidyterra)
library(tidyverse)
library(rgbif)
library(scales)

# Load data --------------------------------------------------------------------

# load base grid
base5k = terra::rast("data/for_shiny/data/canada-basegrids/canada.base.5k.tiff")

# load observations (this is from iNat downloads)
dat = readRDS("data/for_shiny/data/BCParks_obs_unique_with_BigTeams_flagged.RDS")

# load bc parks polygons
parks = terra::vect("data/for_shiny/data/bc-parks-pol/bc_parks.shp")

# make base grid with 1 and 0 to outline BC Parks
base5k = project(base5k, crs(parks))
parks5k = rasterize(parks, base5k, touches = TRUE) |> trim()
terra::writeRaster(parks5k, "outputs/spatial-layers/parks5k.tif")

# make a park = 0 layer to add scores to!
parks0 = parks5k
parks0[parks5k == 1] <- 0


# Data for specific scores =====================================================

## Climate ----
# import mean annual temperature and precipitation
clim = c(
  terra::rast("~/Documents/GitHub/sampling-scenarios/sampling-scenarios/data-raw/climate/wc2.1_10m/wc2.1_10m_bio_1.tif"),
  terra::rast("~/Documents/GitHub/sampling-scenarios/sampling-scenarios/data-raw/climate/wc2.1_10m/wc2.1_10m_bio_12.tif")
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


# Spatial coverage scores ======================================================

# Density map ------------------------------------------------------------------

map_obsdens = rasterize(obs, parks5k, fun = "length")
map_obsdens[is.na(map_obsdens)] <- 0
map_obsdens[is.na(parks5k)] <- NA
plot(map_obsdens)

# Last update map --------------------------------------------------------------

map_lastupdate = rasterize(obs, parks5k, field = "year", fun = "max")
map_lastupdate[is.na(map_lastupdate)] <- 0
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


# Completeness score -----------------------------------------------------------



# Taxonomic coverage scores ====================================================



# Optional purpose-oriented scores =============================================

# SDM discrepancy ----

# sdm_cv = terra::rast("~/Documents/GitHub/BTG-analyse-the-gap/data/heavy/noah-uncertainty/cv_8meth.tif")



# Make a scores table ----------------------------------------------------------

scores = data.frame(
  "cellID" = cells(parks5k),
  "score_spatial_density" = values(map_obsdens) |> na.omit() |> as.vector(),
  "score_spatial_lastupdate" = values(map_lastupdate) |> na.omit() |> as.vector(),
  "score_spatial_climfreq" = values(map_climfrequency) |> na.omit() |> as.vector()
)

# scale columns between 0 and 1
scores_scaled = scores |>
  # scale between 0 and 1
  mutate(across(2:score_spatial_climfreq, rescale)) |>
  # reverse scores so 1 = high priority and  0 = low priority
  mutate(across(2:score_spatial_climfreq, ~ 1 - .x))

# Weight importance of the scores
scores_weights = rep(1, ncol) / c(ncol(scores_scaled)-1)

# Tally total score per cell
scores_weighted = scores_scaled 
for(i in 2:ncol(scores_weighted)){
  scores_weighted[,i] = scores_weighted[,i]*scores_weights[i-1]
}
scores_m = as.matrix(scores_weighted[,-1])
scores_total = rowSums(scores_m)

# map the scores 
map_scores = parks5k
map_scores[parks5k==1] <- scores_total
plot(map_scores)