# Script to calculate spatial and climate distances between neighbouring observations

library(terra)
library(sf)
library(nngeo)
library(ggplot2)

df = readRDS("outputs/inatobs_amphibians_bcparks.rds") |> 
  sf::st_as_sf()

# 2. Calculate the 3 nearest neighbors' distances
# st_nn returns a list of indices; 'returnDist = TRUE' provides the distances
nn_data <- st_nn(df, df, k = 4, returnDist = TRUE, progress = TRUE)

# Note: k = 4 because the nearest neighbor to a point is always itself (distance 0).
# We will exclude the first column to get the true 3 nearest external neighbors.

# 3. Extract distances and format into a dataframe
dist_matrix <- do.call(rbind, nn_data$dist)

# Remove the first column (self-distance) and keep the 3 nearest
dist_only <- dist_matrix[, 2:4]
colnames(dist_only) <- c("dist_1", "dist_2", "dist_3")
dist_avg = apply(dist_only, 1, mean, na.rm = T)

# 4. Bind back to the original object
my_points_result <- cbind(df, dist_only)
my_points_result <- cbind(my_points_result, dist_avg)

# View results
print(my_points_result)

ggplot(data = my_points_result) +
  geom_sf(aes(color = dist_avg, size = dist_avg)) +
  scale_color_viridis_c(option = "plasma", name = "Avg Distance\n(3-NN)") +
  labs(
    title = "Average Distance to 3 Nearest Neighbors",
    subtitle = "Brighter/larger points are more isolated",
    size = "Distance"
  ) +
  theme_minimal()
