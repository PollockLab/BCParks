# Script to prepare inat parquet file for the BC maps

# load libraries
library(arrow)
library(sf)
library(tidyverse)

# load data files
inat = open_dataset("data/for_reporting/private/BCiNat_PRIVATE.parquet")

# filter to research only and keep only essential columns
inat_sf = inat |>  
  filter(quality_grade == "research") |>
  select(lon, lat, scientific_name, iconic_taxon_name, observed_on, user_login) |>
  collect() |> 
  # spatialise
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# change species names with > 2 characters (subsp) to just 2 words
# and remove ones with only 1 word
temp = 
  inat_sf |>
  separate(scientific_name, 
           into = c("genus", "species"), 
           sep =  " ", extra = "drop") # drop more than 2 vals
# find species with nothing in the species column
to_remove = which(is.na(temp$species))
# make species name again
temp$scientific_name = paste(temp$genus, temp$species, sep = " ")
inat_sf$scientific_name <- temp$scientific_name
# remove species with only 1 value
inat_sf <- inat_sf[-to_remove,]
rm(temp)

# make a year column ----
inat_sf <- inat_sf |>
  separate(col = "observed_on", into = c("day", "month", "year"), sep = "/")
inat_sf$year = as.numeric(inat_sf$year)

# Transform to BC Albers (EPSG:3005) for accurate distance-based hexagons
inat_sf <- st_transform(inat_sf, 3005)

# save parquet 
# arrow::write_parquet(sf::st_drop_geometry(inat_sf), 
#                      "data/processed/inat_bc_PRIVATEDONOTSHARE_cleanedsp.parquet")

# save spatial data
saveRDS(inat_sf, "data/processed/inat_bc_PRIVATEDONOTSHARE_cleanedsp_sf.rds")
