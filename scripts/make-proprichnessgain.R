# Script to make hexagon map of species richness and observation density
# Resoutions: 100, 50, 20, 10, 5

# load libraries
library(arrow)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(pmtiles)
library(patchwork)
theme_set(theme_void())

# load data files
inat_sf = readRDS("data/processed/inat_bc_PRIVATEDONOTSHARE_cleanedsp_sf.rds")
# species list
splist = read.csv("data/for_reporting/private/BC_S-listed_spp.csv")
# # Provincial = BC.List (filter to Blue & Red listed)
# # National = COSEWIC 
# user logins
logins = read.csv("outputs/bigteam_userlogins.csv", row.names = 1)
# bc polygon
bc = st_read("data/polygons/bc_polygon.shp")

# Transform to BC Albers (EPSG:3005) for accurate distance-based hexagons
bc_albers <- st_transform(bc, 3005)

## make base hexmaps -----------------------------------------------------------

resolutions <- c(25000, 50000, 100000)  # meters
n_maps = length(resolutions)

# make base hexgrids
hexbase = list()
for (res in resolutions) {
  hex <- st_make_grid(bc_albers, cellsize = res, square = FALSE)
  hex <- st_sf(hex_id = 1:length(hex), geometry = hex)
  hex <- hex[st_intersects(hex, bc_albers, sparse = FALSE), ]
  hexbase[[paste0(res/1000, 'km')]] <- hex 
}

# function to count species per hex --------------------

hex_sr = function(inat_data, hex_basemap, resolution){
  # join the inat data to the hexmap
  hex_join = st_join(inat_data, 
                     hex_basemap,                    
                     join = st_intersects)
  # count distinct species per hex
  hex_sr = hex_join |>
    group_by(hex_id) |>
    distinct(scientific_name) |>
    summarize("n_sp" = n())
  # add number of species to hexmap
  hex = left_join(hex_basemap, hex_sr, by = "hex_id") |> 
    select("n_sp")
  # project
  return(hex)
}


## Proportional richness map ---------------------------------------------------

for(iconic_taxa in unique(inat_sf$iconic_taxon_name)[-14]){
  
  # filter to the group
  inat_taxa = inat_sf |> filter(iconic_taxon_name == iconic_taxa)
  
  # split data -------------------------------------------------------------------
  
  # all users 
  pre = inat_taxa |> filter(year < 2019)
  post = inat_taxa 
  
  # teams vs. visitors
  post.visitors = inat_taxa |> filter(!user_login %in% team)
  pre.visitors = post.visitors |> filter(year < 2019)
  post.teams = inat_taxa |> filter(user_login %in% team)
  pre.teams = post.teams |> filter(year < 2019) 
  
  # run for the datasets
  df_list = list(pre, post, pre.teams, post.teams, pre.visitors, post.visitors)
  df_names = c("pre", "post", "pre.teams", "post.teams", "pre.visitors", "post.visitors")
  names(df_list) = df_names
  
  # run for all datasets
  map_list <- list()
  for (i in 1:length(df_list)) {
    temp = list()
    for(r in 1:length(resolutions)){
      temp[[r]] = hex_sr(inat_data = df_list[[i]], 
                         hex_basemap = hexbase[[r]],
                         resolution = resolutions[r])
      names(temp)[r] = paste0(df_names[i], "-", resolutions[r]/1000, "km")
    }
    map_list[[i]] <- temp
  }
  names(map_list) = df_names
  
  # change pre NAs to zeros to avoid NAs downstream
  for(i in grep("pre", names(df_list))){
    for(r in 1:length(resolutions)){
      map_list[[i]][[r]]$n_sp[which(is.na(map_list[[i]][[r]]$n_sp))] <- 0
    }
  }
  
  # calculate ratios
  map_allbc = list()
  for (r in 1:length(resolutions)) {
    map_allbc[[r]] = st_join(map_list$pre[[r]], 
                             map_list$post[[r]], 
                             suffix = c(".pre", ".post")) |>
      mutate(ratio = (n_sp.post - n_sp.pre)/(n_sp.pre + 1))
  }
  
  map_teams = list()
  for (r in 1:length(resolutions)) {
    map_teams[[r]] = st_join(map_list$pre.teams[[r]], 
                             map_list$post.teams[[r]], 
                             suffix = c(".pre", ".post")) |>
      mutate(ratio = (n_sp.post - n_sp.pre)/(n_sp.pre + 1))
  }
  
  map_visitors = list()
  for (r in 1:length(resolutions)) {
    map_visitors[[r]] = st_join(map_list$pre.visitors[[r]], 
                                map_list$post.visitors[[r]], 
                                suffix = c(".pre", ".post")) |>
      mutate(ratio = (n_sp.post - n_sp.pre + 1)/(n_sp.pre + 1))
  }
  
  # save all maps
  maps_taxa = list("All" = map_allbc,
                   "Teams" = map_teams,
                   "Visitors" = map_visitors)
  saveRDS(maps_taxa, paste0("outputs/spatial-layers/hexmaps_relativegain_",iconic_taxa,".rds"))
  
  
  # plot the ratio maps ----------------------------------------------------------
  
  # get common color limits
  temp1 = maps_taxa$All[[3]]$ratio |> range(na.rm = T)
  temp2 = maps_taxa$Teams[[3]]$ratio |> range(na.rm = T)
  temp3 = maps_taxa$Visitors[[3]]$ratio |> range(na.rm = T)
  temp = rbind(temp1, temp2, temp3)
  lims = c(min(temp[,1]), max(temp[,2]))
  
  A = ggplot() +
    geom_sf(data = maps_taxa$All[[3]], aes(fill = ratio), linewidth = .3) +
    labs(title = "All") +
    scale_fill_viridis_c(option = "plasma", 
                         name = "Relative\ngain (%)",
                         trans = "sqrt", 
                         na.value = "transparent", 
                         limits = lims) +
    theme(legend.position = "none")
  
  B = ggplot() +
    geom_sf(data = maps_taxa$Teams[[3]], aes(fill = ratio), linewidth = .3) +
    labs(title = "Teams") +
    scale_fill_viridis_c(option = "plasma", 
                         name = "Relative\ngain (%)",
                         trans = "sqrt", 
                         na.value = "transparent", 
                         limits = lims) +
    theme(legend.position = "none")
  
  C = ggplot() +
    geom_sf(data = maps_taxa$Visitors[[3]], aes(fill = ratio), linewidth = .3) +
    labs(title = "Visitors") +
    scale_fill_viridis_c(option = "plasma", 
                         name = "Relative\ngain (%)",
                         trans = "sqrt", 
                         na.value = "transparent", 
                         limits = lims) +
    theme(legend.position = "right")
  A + B + C 
  ggsave(paste0("figures/hexmap_relativegain_", iconic_taxa, ".png"), 
         width = 13.2, height = 4.05)
  
}
