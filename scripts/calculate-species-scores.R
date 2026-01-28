# Script to calculate species-level scores for priority sampling
# Scope = Canada

# load libraries
library(rgbif)
library(taxize)
library(tidyverse)

# load species lists -----------------------------------------------------------

# full species list from Lab Sharepoint, Blitz the Gap/species list organization/data/processed
full = read.csv("data/for_scoring/species-lists/completeList.csv")
mw = read.csv("data/full_list_canadas_most_wanted(in).csv", sep = ";")

## Invasive species ------------------------------------------------------------

sp_inv = mw |>
  filter(is_Invasive == 1) |> 
  select(Scientific.Name)

## At risk (national risk)
sp_atrisk = mw |>
  filter(is_National == 1) |> 
  select(Scientific.Name)

# Taxonomic representativeness score -------------------------------------------

taxrep = full |>
  group_by(iNatClass) |>
  distinct(species) |>
  summarise("n_sp" = n(),
            "prop_sp" = n()/nrow(full))


