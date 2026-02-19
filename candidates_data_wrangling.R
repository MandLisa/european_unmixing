library(sf)

gpkg_path <- "/mnt/eo/EU_unmixing/endmembers_vis_filtered_new.gpkg"
data <- st_read(gpkg_path)

data <- data[, -c(2:6)]
data <- data[, -c(9:11)]
data <- data[, -c(12:20)]
data <- data[, -c(9:11)]


coords <- st_coordinates(data)

data$X <- coords[, 1]
data$Y <- coords[, 2]

st_crs(data)

# 1) extract projected coordinates (EPSG:3035) from the original geometry
xy_3035 <- st_coordinates(data)
data$X_3035 <- xy_3035[, 1]
data$Y_3035 <- xy_3035[, 2]

# 2) transform to lon/lat (EPSG:4326) and extract geographic coordinates
data_ll <- st_transform(data, 4326)
xy_4326 <- st_coordinates(data_ll)

data$lon <- xy_4326[, 1]
data$lat <- xy_4326[, 2]

data <- data[, -c(18:19)]

library(sf)

out_gpkg <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/candidates_vis_filtered.gpkg"
out_gpkg <- "/mnt/eo/EU_unmixing/candidates_vis_filtered_clean.gpkg"

st_write(data, out_gpkg, layer = "candidates", delete_layer = TRUE)


plot_df <- plot_df %>%
  group_by(ecoregion) %>%
  mutate(prop = n / sum(n))


ggplot(plot_df, aes(x = letter_group, y = n)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ ecoregion) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Letter group",
    y = "Number of points"
  )


library(dplyr)

summary_table <- data %>%
  st_drop_geometry() %>%              # remove geometry if sf object
  count(letter_group, ecoregion) %>%  # count combinations
  tidyr::pivot_wider(
    names_from  = ecoregion,
    values_from = n,
    values_fill = 0
  )

summary_table

summary_table <- summary_table %>%
  mutate(
    land_cover = case_when(
      letter_group == "A" ~ "Artificial surfaces",
      letter_group == "C" ~ "Woodland",
      letter_group == "D" ~ "Grassland",
      letter_group == "E" ~ "Shrubland",
      letter_group == "F" ~ "Bare land",
      TRUE ~ NA_character_
    )
  ) %>%
  relocate(land_cover, .after = letter_group)

summary_table





