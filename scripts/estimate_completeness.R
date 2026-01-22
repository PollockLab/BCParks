# Sampling priority index derived from applying a general curve from 
# filled cells in an ecoregion 

library(tidyverse)
library(sf)
library(terra)
library(iNEXT)
library(arrow)
library(betareg)

theme_set(ggpubr::theme_pubr())

# read density map
map_obsdens = readRDS("outputs/inatdensity_amphibians_bcparks.rds")

## Select well-sampled cells ===================================================

hist(map_obsdens)
quantile(values(map_obsdens), na.rm = T)

# define threshold to find "well sampled" cells
sample_threshold = 90

# threshold cells (to define more thoughtfully later)
sampled_well = map_obsdens
sampled_well[map_obsdens < sample_threshold] <- NA
plot(sampled_well)

# extract observations in these cells (NA = not selected as well sampled)
obs_sel = terra::extract(sampled_well, obs, cell = TRUE) 

# clean up the data 
obs_df = obs
obs_df$selected = !is.na(obs_sel$V1_length)
obs_df$group = obs_sel$cell |> as.character()

# convert to sf
obs_sf = st_as_sf(obs_df)


# Rarefy in the cells -----------------------------------------------------------

# subset data
df = obs_sf |> 
  filter(selected == TRUE) |>
  select(c(group, scientific_name)) |>
  group_by(group, scientific_name) |>
  summarise(n = n()) |>
  sf::st_drop_geometry()

# prep data for iNEXT
groups = unique(df$group)

# there are nicer ways to do this but whatever! for now :)
occ = list()
for(i in 1:length(groups)){
  df2 = df |> filter(group == groups[i])
  occ[[i]] = df2$n |> sort(decreasing = TRUE)
}
names(occ) = groups

# rarefaction
out_cell = iNEXT(occ, q = 0, datatype = "abundance")
ggiNEXT(out_cell)
ggsave("figures/accumulationcurve_SR_highestsampledcells_amphibians.png")
ggiNEXT(out_cell, type = 2) + scale_x_sqrt()
ggsave("figures/accumulationcurve_coverage_highestsampledcells_amphibians.png")
ggiNEXT(out_cell, type = 3)
ggsave("figures/accumulationcurve_SRcoverage_highestsampledcells_amphibians.png")


## Model the curves to extend them to other cells ------------------------------

# make a dataframe with the outputs we need
mdf = data.frame(
  "coverage" = out_cell$iNextEst$coverage_based$SC, # standardized sample coverage
  "n_individuals" = out_cell$iNextEst$coverage_based$m, # sample size 
  "diversity" = out_cell$iNextEst$coverage_based$qD,
  "group" = out_cell$iNextEst$coverage_based$Assemblage
)

# plot sample coverage ~ effort (number of individuals)
plot(mdf$coverage ~ mdf$n_individuals)

# plot in a prettier style for saving
ggplot(data = mdf) +
  geom_line(aes(x = n_individuals, y = coverage, col = group)) +
  theme(legend.position = "right") +
  labs(x = "Number of individuals", y = "Sample coverage", col = "Cell") +
  hrbrthemes::theme_ipsum_es()
ggsave("figures/accumulationcurve_coverageVSeffort_highestsampledcells_amphibians.png")

# model the curves together
m = betareg::betareg(coverage ~ n_individuals + 0, 
                     data = mdf, link = "logit")
pred = predict(m, type = "response")
predvar = predict(m, type = "variance")

# plot the model predictions
plot(pred ~ mdf$n_individuals, col = "red", ylim = c(0,1))
points(pred-predvar ~ mdf$n_individuals, col = "blue")
points(pred+predvar ~ mdf$n_individuals, col = "blue")
points(mdf$coverage ~ mdf$n_individuals)
# we'll need a better model here... but let's use it for now.

# model the variance too just to see
plot(predvar ~ mdf$n_individuals, col = "blue")


# Predict coverage for the remaining cells -------------------------------------

# predicted coverage = completeness metric (I think?)
completeness = data.frame(
  "cellID" = cells(map_obsdens),
  "density" = na.omit(values(map_obsdens)) |> as.vector(),
  "coverage" = NA
)

# predict coverage per cell
preds = predict(m, type = "response", 
                data.frame("n_individuals" = completeness$density))
# join the the table
completeness$coverage = preds

# cells that were empty (0 density) should be 0 coverage too
completeness$coverage[which(completeness$density == 0)] <- 0

# map the completeness on the map
map_coverage = map_obsdens
map_coverage[parks5k==1] <- completeness$coverage
plot(map_coverage)
saveRDS(map_coverage, "outputs/spatial-layers/scores_map_coverage.rds")

ggplot() +
  tidyterra::geom_spatraster(data = map_coverage) +
  scale_fill_viridis_c(option = "plasma", 
                       na.value = "transparent", limits = c(0,1)) +
  labs(fill = "Predicted\nSample coverage") +
  theme_void() +
  theme(legend.position = "right")
ggsave("figures/map_predictedsamplecoverage.png", width = 7.36, height = 6.76)
