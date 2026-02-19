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
invasives = read.csv("data/for_reporting/BC_invasive_sp_priorities.csv", 
                     header = FALSE)
bc = st_read("data/polygons/bc_polygon.shp")

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

## Uncomment this to make the hexmaps for subset species lists -----------------

# ## Species at risk (National - COSEWIC) --------------------------------------
# 
# # get list of species listed in COSEWIC
# sarN = splist |> 
#   filter(COSEWIC %in% c("E", "T", "SC", 
#                         "E/T", "XT", "T/SC", "XX"))
# 
# inat_sf = inat_sf |>
#   filter(scientific_name %in% sarN$Scientific.Name)


# ## Species at risk (Provincial) ------------------------------------------------
# 
# sarBC = splist |> 
#   filter(BC.List %in% c("Blue", "Red"))
# inat_sf = inat_sf |>
#   filter(scientific_name %in% sarBC$scientific_name_cut)

## Invasive species ------------------------------------------------------------

invasives$scientific_name = paste(invasives$V4, invasives$V5) |> 
  # clean up the extra spaces
  stringr::str_squish()

inat_sf = inat_sf |>
  filter(scientific_name %in% invasives$scientific_name)



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
hex_list = readRDS("outputs/report-pmtiles/obsdens-layers.rds")


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
# names(hex_list) <- paste0("alltaxa-", names(hex_list))
# if SAR:
#names(hex_list) <- paste0("SAR-", names(hex_list))
# if BC SR:
# names(hex_list) <- paste0("SARprov-", names(hex_list))
# if invasive:
names(hex_list) <- paste0("invasives-", names(hex_list))

# save rds
# saveRDS(hex_list, "outputs/report-pmtiles/alltaxa-hex.rds")
# saveRDS(hex_list, "outputs/report-pmtiles/SAR-hex.rds")
# saveRDS(hex_list, "outputs/report-pmtiles/SARprov-hex.rds")
saveRDS(hex_list, "outputs/report-pmtiles/invasives-hex.rds")

## PMTILES ---------------------------------------------------------------------

# read
# hex_list = readRDS("outputs/report-pmtiles/alltaxa-hex.rds")
# hex_list = readRDS("outputs/report-pmtiles/SAR-hex.rds")
# hex_list = readRDS("outputs/report-pmtiles/SARprov-hex.rds")
hex_list = readRDS("outputs/report-pmtiles/invasives-hex.rds")

# rename files
# names(hex_list) = gsub("SAR-obsdens-", "obsdens_all_", names(hex_list))
# names(hex_list) = gsub("SAR-spdens-", "spdens_all_", names(hex_list))
names(hex_list) = gsub("invasives-obsdens-", "obsdens_all_", names(hex_list))
names(hex_list) = gsub("invasives-spdens-", "spdens_all_", names(hex_list))

for(i in 1:length(hex_list)){
  colnames(hex_list[[i]]) <- gsub("n_obs", "obsdens", colnames(hex_list[[i]]))
  colnames(hex_list[[i]]) <- gsub("n_sp", "spdens", colnames(hex_list[[i]]))
  if("hex_id" %in% colnames(hex_list[[i]])){
  hex_list[[i]] = select(hex_list[[i]], -hex_id)
  }
}

# stick obsdens and spdens together
map_100km = st_join(hex_list$obsdens_all_100km, hex_list$spdens_all_100km)
map_50km = st_join(hex_list$obsdens_all_50km, hex_list$spdens_all_50km)
map_25km = st_join(hex_list$obsdens_all_25km, hex_list$spdens_all_25km)
map_10km = st_join(hex_list$obsdens_all_10km, hex_list$spdens_all_10km)
map_5km = st_join(hex_list$obsdens_all_5km, hex_list$spdens_all_5km)

# all_list = list(all_5km, all_10km, all_25km, all_50km, all_100km)
# names(all_list) = c("all_5km", "all_10km", "all_25km", "all_50km", "all_100km")
map_list = list(map_5km, map_10km, map_25km, map_50km, map_100km)
names(map_list) = c("map_5km", "map_10km", "map_25km", "map_50km", "map_100km")

for(i in 1:length(map_list)){
  # st_write(all_list[[i]], paste0("outputs/report-pmtiles/SF/alltaxa_", resolutions[i]/1000, "km.shp"))
  # st_write(map_list[[i]], paste0("outputs/report-pmtiles/SF/SAR_", resolutions[i]/1000, "km.shp"))
  # st_write(map_list[[i]], paste0("outputs/report-pmtiles/SF/SARprov_", resolutions[i]/1000, "km.shp"))
  st_write(map_list[[i]], paste0("outputs/report-pmtiles/SF/invasives_", resolutions[i]/1000, "km.shp"))
  
}

# make a cloud-optimized version 
for(i in 2:length(map_list)){
  # pmtiles::pm_create(all_list[[i]], 
  #                    paste0("outputs/report-pmtiles/", names(all_list)[i], ".pmtiles"), 
  #                    layer_name = names(all_list)[i])
  # pmtiles::pm_create(map_list[[i]],
  #                    paste0("outputs/report-pmtiles/SAR/", names(map_list)[i], ".pmtiles"),
  #                    layer_name = names(map_list)[i])
  # pmtiles::pm_create(map_list[[i]],
  #                    paste0("outputs/report-pmtiles/SARprov/", names(map_list)[i], ".pmtiles"),
  #                    layer_name = names(map_list)[i])
  pmtiles::pm_create(map_list[[i]],
                     paste0("outputs/report-pmtiles/invasives/", names(map_list)[i], ".pmtiles"),
                     layer_name = names(map_list)[i])
}