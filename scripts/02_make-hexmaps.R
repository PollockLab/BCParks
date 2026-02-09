# Script to make hexagon map of species richness and observation density
# Resoutions: 100, 50, 20, 10, 5

# load libraries
library(arrow)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(pmtiles)

# load data files
inat = open_dataset("data/for_reporting/private/BCiNat_PRIVATE.parquet")
# species list
splist = read.csv("data/for_reporting/private/BC_S-listed_spp.csv")
# # Provincial = BC.List (filter to Blue & Red listed)
# # National = COSEWIC 
bc = st_read("data/polygons/bc_polygon.shp")
# bc = st_read("~/Downloads/BC_Boundary_Terrestrial.gpkg")
# bc = st_simplify(bc)
# buffer_bc <- st_buffer(bc, dist = 0.1)
#bc = st_cast(bc, to = "POLYGON")

# filter to research only and keep only essential columns
inat_sf = inat |>  
  filter(quality_grade == "research") |>
  select(lon, lat, scientific_name, iconic_taxon_name) |>
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

# Transform to BC Albers (EPSG:3005) for accurate distance-based hexagons
bc_albers <- st_transform(bc, 3005)
wt_albers <- st_transform(inat_sf, 3005)



# =============================================================================
# 4. CREATE HEXAGON GRIDS AT MULTIPLE RESOLUTIONS
# =============================================================================

resolutions <- c(5000, 10000, 25000, 50000, 100000)  # meters
hex_list <- list()

for (res in resolutions) {
  hex <- st_make_grid(bc_albers, cellsize = res, square = FALSE)
  hex <- st_sf(hex_id = 1:length(hex), geometry = hex)
  hex <- hex[st_intersects(hex, bc_albers, sparse = FALSE), ]
  hex$n_obs <- lengths(st_intersects(hex, wt_albers))
  hex_wgs <- st_transform(hex, 4326)
  hex_list[[paste0(res/1000, 'km')]] <- hex_wgs
}
names(hex_list) <- paste0("obsdens-", names(hex_list))

# save rds
saveRDS(hex_list, "outputs/report-pmtiles/obsdens-layers.rds")



# =============================================================================
# 5. CREATE HEXAGON GRIDS AT MULTIPLE RESOLUTIONS FOR SPECIES RICHNESS
# =============================================================================

resolutions <- c(5000, 10000, 25000, 50000, 100000)  # meters

for (res in resolutions) {
  
  # make the hex grid
  hex <- st_make_grid(bc_albers, cellsize = res, square = FALSE)
  hex <- st_sf(hex_id = 1:length(hex), geometry = hex)
  hex <- hex[st_intersects(hex, bc_albers, sparse = FALSE), ]
  hex_join = st_join(wt_albers, hex, join = st_intersects)
  hex_sr = hex_join |>
    group_by(hex_id) |>
    distinct(scientific_name) |>
    summarize("n_sp" = n())
  hex = left_join(hex, hex_sr, by = "hex_id") |> select("n_sp")
  hex_wgs <- st_transform(hex, 4326)
  hex_list[[paste0(res/1000, 'km')]] <- hex_wgs
}
names(hex_list) <- paste0("spdens-", names(hex_list))
names(hex_list) <- gsub("spdens-obsdens", "obsdens", names(hex_list))
names(hex_list) <- paste0("alltaxa-", names(hex_list))

# save rds
saveRDS(hex_list, "outputs/report-pmtiles/alltaxa-hex.rds")

# make a cloud-optimized version to make the app faster
pmtiles::pm_create(hex_list, "outputs/report-pmtiles/alltaxa-hex.pmtiles",
                   layer_name = names(hex_list))


## Species at risk (National - COSEWIC) ----------------------------------------


# get list of species listed in COSEWIC
sarN = splist |> 
  filter(COSEWIC %in% c("E", "T", "SC", 
                        "E/T", "XT", "T/SC", "XX"))

inat_sarN = inat_sf |>
  filter(scientific_name %in% sarN$Scientific.Name)

# =============================================================================
# 4. CREATE HEXAGON GRIDS AT MULTIPLE RESOLUTIONS
# =============================================================================

resolutions <- c(5000, 10000, 25000, 50000, 100000)  # meters
hex_list <- list()

for (res in resolutions) {
  hex <- st_make_grid(bc_albers, cellsize = res, square = FALSE)
  hex <- st_sf(hex_id = 1:length(hex), geometry = hex)
  hex <- hex[st_intersects(hex, bc_albers, sparse = FALSE), ]
  hex$n_obs <- lengths(st_intersects(hex, wt_albers))
  hex_wgs <- st_transform(hex, 4326)
  hex_list[[paste0(res/1000, 'km')]] <- hex_wgs
}
names(hex_list) <- paste0("sarN-obsdens-", names(hex_list))

# save rds
saveRDS(hex_list, "outputs/report-pmtiles/sarN-obsdens-layers.rds")

## RICHNESS --------------------------------------------------------------------

for (res in resolutions) {
  
  # make the hex grid
  hex <- st_make_grid(bc_albers, cellsize = res, square = FALSE)
  hex <- st_sf(hex_id = 1:length(hex), geometry = hex)
  hex <- hex[st_intersects(hex, bc_albers, sparse = FALSE), ]
  hex_join = st_join(wt_albers, hex, join = st_intersects)
  hex_sr = hex_join |>
    group_by(hex_id) |>
    distinct(scientific_name) |>
    summarize("n_sp" = n())
  hex = left_join(hex, hex_sr, by = "hex_id") |> select("n_sp")
  hex_wgs <- st_transform(hex, 4326)
  hex_list[[paste0(res/1000, 'km')]] <- hex_wgs
}
names(hex_list)[6:10] <- paste0("sarN-spdens-", names(hex_list)[6:10])
names(hex_list) <- gsub("spdens-sarN-obsdens-", "sarN-spdens-", names(hex_list))

# save rds
saveRDS(hex_list, "outputs/report-pmtiles/sarN-hex.rds")



hex_list = readRDS("outputs/report-pmtiles/sarN-hex.rds")

## Species at risk (Provincial)

sarBC = splist |> filter(BC.List %in% c("Blue", "Red"))
inat_sarBC = inat_sf |>
  filter(scientific_name %in% sarBC$Scientific.Name)
