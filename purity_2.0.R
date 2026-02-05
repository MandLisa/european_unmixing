# ============================================================
# Spectral purity filtering (time-series level) in NDVI–SWIR space
# - Input: df_wide (must contain: point_id, year, lc1_label, ndvi, b5, b6)
# - Output:
#   * purity_scores: purity metrics per point_id
#   * pure_ids: selected pure point_ids per class
#   * df_pure: filtered observations (all years) for pure point_ids
# - Notes:
#   * Shrubland is excluded.
#   * Works with dplyr >= 1.1.0 (uses cross_join()).
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---------------------------
# USER SETTINGS
# ---------------------------
target_classes <- c("Broadleaved woodland", "Coniferous woodland", "Grassland")  # adjust case/spelling if needed
exclude_classes <- c("Shrubland")                                               # excluded explicitly
min_years_required <- 20                                                        # require at least this many annual records per point_id
selection_mode <- c("top_n", "top_quantile")[1]                                  # "top_n" or "top_quantile"
top_n_per_class <- 300                                                          # used if selection_mode == "top_n"
top_quantile <- 0.90                                                            # used if selection_mode == "top_quantile" (e.g., 0.90 = top 10%)
use_yearwise_scaling <- FALSE                                                   # TRUE: z-score within year; FALSE: z-score globally
swir_ratio_upper_cap <- Inf                                                     # optional: set e.g. 5 to remove extreme ratios

# ---------------------------
# 1) Preprocess + filter classes (exclude shrubland)
# ---------------------------
df_fs <- df_wide %>%
  mutate(
    swir1 = b5,
    swir2 = b6,
    swir_ratio = swir2 / swir1
  ) %>%
  filter(
    lc1_label %in% target_classes,
    is.finite(ndvi),
    is.finite(swir_ratio),
    swir1 > 0, swir2 > 0,
    swir_ratio > 0,
    swir_ratio <= swir_ratio_upper_cap
  ) %>%
  (if (use_yearwise_scaling) {
    . %>%
      group_by(year) %>%
      mutate(
        ndvi_z = as.numeric(scale(ndvi)),
        swir_z = as.numeric(scale(swir_ratio))
      ) %>%
      ungroup()
  } else {
    . %>%
      mutate(
        ndvi_z = as.numeric(scale(ndvi)),
        swir_z = as.numeric(scale(swir_ratio))
      )
  })

# ---------------------------
# 2) Time-series (point_id) mean position in feature space
# ---------------------------
point_means <- df_fs %>%
  group_by(point_id, lc1_label) %>%
  summarise(
    n_years = n(),
    px = mean(ndvi_z, na.rm = TRUE),
    py = mean(swir_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_years >= min_years_required)

if (nrow(point_means) == 0) {
  stop("No point_id trajectories left after filtering. Check class names, min_years_required, and input columns.")
}

# ---------------------------
# 3) Class centroids (in feature space)
# ---------------------------
class_centroids <- point_means %>%
  group_by(lc1_label) %>%
  summarise(
    cx = mean(px, na.rm = TRUE),
    cy = mean(py, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------
# 4) Distance to own class centroid
# ---------------------------
purity_dist_own <- point_means %>%
  left_join(class_centroids, by = "lc1_label") %>%
  mutate(
    d_own = sqrt((px - cx)^2 + (py - cy)^2)
  ) %>%
  select(point_id, lc1_label, n_years, px, py, d_own)

# ---------------------------
# 5) Minimum distance to competing class centroids (margin component)
# ---------------------------
purity_margin <- point_means %>%
  cross_join(
    class_centroids %>% rename(other_label = lc1_label, ox = cx, oy = cy)
  ) %>%
  filter(lc1_label != other_label) %>%
  mutate(
    d_other = sqrt((px - ox)^2 + (py - oy)^2)
  ) %>%
  group_by(point_id, lc1_label) %>%
  summarise(
    d_other_min = min(d_other, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------
# 6) Combine into purity score
#    margin = (nearest competing centroid distance) - (own centroid distance)
#    larger margin => purer and more class-distinct
# ---------------------------
purity_scores <- purity_dist_own %>%
  inner_join(purity_margin, by = c("point_id", "lc1_label")) %>%
  mutate(
    margin = d_other_min - d_own,
    purity_score = as.numeric(scale(margin))
  ) %>%
  arrange(desc(purity_score))

# ---------------------------
# 7) Select pure trajectories per class
# ---------------------------
pure_ids <- purity_scores %>%
  group_by(lc1_label) %>%
  filter(
    case_when(
      lc1_label == "Coniferous woodland"  ~
        d_own  <= quantile(d_own,  0.30, na.rm = TRUE) &
        margin >= quantile(margin, 0.85, na.rm = TRUE),
      
      lc1_label == "Broadleaved woodland" ~
        d_own  <= quantile(d_own,  0.30, na.rm = TRUE) &
        margin >= quantile(margin, 0.85, na.rm = TRUE),
      
      lc1_label == "Grassland"            ~
        d_own  <= quantile(d_own,  0.40, na.rm = TRUE) &
        margin >= quantile(margin, 0.70, na.rm = TRUE),
      
      TRUE ~ FALSE
    )
  ) %>%
  ungroup()

# ---------------------------
# 8) Filter the full dataset to retain only selected pure point_ids
# ---------------------------
df_pure <- df_fs %>%
  semi_join(pure_ids %>% select(point_id), by = "point_id")

# ---------------------------
# 9) Quick diagnostics (optional)
# ---------------------------
message("Selected pure trajectories per class:")
print(pure_ids %>% count(lc1_label, name = "n_selected"))

message("Purity margin summary by class (selected only):")
print(
  pure_ids %>%
    group_by(lc1_label) %>%
    summarise(
      min = min(margin, na.rm = TRUE),
      q25 = quantile(margin, 0.25, na.rm = TRUE),
      med = median(margin, na.rm = TRUE),
      q75 = quantile(margin, 0.75, na.rm = TRUE),
      max = max(margin, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(lc1_label)
)

# Optional: visual check of selected vs all (in z-space)
# (comment out if not needed)
p_check <- ggplot(df_fs, aes(ndvi_z, swir_z)) +
  geom_point(alpha = 0.03, size = 0.35) +
  geom_point(
    data = df_pure,
    aes(color = lc1_label),
    alpha = 0.25,
    size = 0.45
  ) +
  #facet_wrap(~ lc1_label) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "NDVI (z)", y = "SWIR ratio (z)", color = "Class")

print

# not in z-scores
p_check_raw <- ggplot(df_fs, aes(ndvi, swir_ratio)) +
  geom_point(alpha = 0.02, size = 0.35, color = "grey30") +
  geom_point(
    data = df_pure,
    aes(color = lc1_label),
    alpha = 0.25,
    size = 0.45
  ) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank()) +
  labs(
    x = "NDVI",
    y = "SWIR ratio (SWIR2 / SWIR1)",
    color = "Class"
  )

plot(p_check_raw)

ggsave(
  filename = "/mnt/eo/EU_unmixing/figs/purity_diagnostic_ndvi_swir.png",
  plot = p_check_raw,
  width = 10,
  height = 7,
  dpi = 300
)

# Objects created:
# - df_fs         : filtered + feature-engineered data (shrubland excluded)
# - point_means   : point_id mean position in feature space
# - class_centroids
# - purity_scores : purity metrics per point_id
# - pure_ids      : selected pure time series IDs per class
# - df_pure       : observations for selected pure time series



# ============================================================
# Annual NDVI–SWIR feature-space figures (3 classes)
# - Grey: all stable observations in that year
# - Colour: globally selected "pure" point_ids (subset present in that year)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---------------------------
# USER SETTINGS
# ---------------------------
target_classes <- c("Broadleaved woodland", "Coniferous woodland", "Grassland")
use_z_space <- FALSE                 # TRUE: plot ndvi_z vs swir_z; FALSE: plot ndvi vs swir_ratio
ncol_years <- 6                     # columns in facet grid for the multi-year figure
alpha_all <- 0.06                   # transparency for grey background
alpha_pure <- 0.35                  # transparency for coloured points
size_all <- 0.35
size_pure <- 0.45
save_one_big <- TRUE
save_per_year <- FALSE              # set TRUE to save one file per year
out_dir <- "/mnt/eo/EU_unmixing/figs"                      # output directory
dpi_out <- 300

# ---------------------------
# 0) Ensure df_fs exists and contains required columns
#    df_fs must have: point_id, year, lc1_label, ndvi, b5, b6 (or swir_ratio), and optionally ndvi_z/swir_z
# ---------------------------

# If df_fs does not exist yet, compute it from df_wide:
if (!exists("df_fs")) {
  df_fs <- df_wide %>%
    mutate(
      swir1 = b5,
      swir2 = b6,
      swir_ratio = swir2 / swir1
    ) %>%
    filter(
      lc1_label %in% target_classes,
      is.finite(ndvi),
      is.finite(swir_ratio),
      swir1 > 0, swir2 > 0,
      swir_ratio > 0
    ) %>%
    mutate(
      ndvi_z = as.numeric(scale(ndvi)),
      swir_z = as.numeric(scale(swir_ratio))
    )
}

# If pure_ids does not exist, stop (you said you already have it)
if (!exists("pure_ids")) {
  stop("pure_ids not found. Please create pure_ids first (must contain point_id and lc1_label).")
}

# Keep only the 3 classes and mark pure vs not
df_plot <- df_fs %>%
  filter(lc1_label %in% target_classes) %>%
  mutate(is_pure = point_id %in% pure_ids$point_id)

# Choose axes (raw vs z-space)
xvar <- if (use_z_space) "ndvi_z" else "ndvi"
yvar <- if (use_z_space) "swir_z" else "swir_ratio"
xlabel <- if (use_z_space) "NDVI (z)" else "NDVI"
ylabel <- if (use_z_space) "SWIR ratio (z)" else "SWIR ratio (SWIR2 / SWIR1)"

# Consistent shapes across classes (as in your earlier figure style)
shape_map <- c(
  "Broadleaved woodland" = 21,
  "Coniferous woodland"  = 24,
  "Grassland"            = 22
)

# ---------------------------
# 1) One multi-panel figure: facets = year
# ---------------------------
p_annual <- ggplot(df_plot, aes(x = .data[[xvar]], y = .data[[yvar]])) +
  # all points in grey (background)
  geom_point(
    data = df_plot %>% filter(!is_pure),
    color = "grey20",
    alpha = alpha_all,
    size = size_all
  ) +
  # pure points in colour, shaped by class
  geom_point(
    data = df_plot %>% filter(is_pure),
    aes(color = lc1_label, shape = lc1_label),
    alpha = alpha_pure,
    size = size_pure
  ) +
  facet_wrap(~ year, ncol = ncol_years) +
  scale_shape_manual(values = shape_map) +
  labs(x = xlabel, y = ylabel, color = "Class", shape = "Class") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 9),
    legend.position = "right"
  )

print(p_annual)

if (save_one_big) {
  ggsave(
    filename = file.path(out_dir, paste0("ndvi_swir_feature_space_annual_", ifelse(use_z_space, "z", "raw"), ".png")),
    plot = p_annual,
    width = 15, height = 9, dpi = dpi_out
  )
}

# ---------------------------
# 2) Optional: save one figure per year (more readable if many years)
# ---------------------------
if (save_per_year) {
  years <- sort(unique(df_plot$year))
  for (yy in years) {
    p_y <- ggplot(df_plot %>% filter(year == yy), aes(x = .data[[xvar]], y = .data[[yvar]])) +
      geom_point(
        data = df_plot %>% filter(year == yy, !is_pure),
        color = "grey20", alpha = alpha_all, size = size_all
      ) +
      geom_point(
        data = df_plot %>% filter(year == yy, is_pure),
        aes(color = lc1_label, shape = lc1_label),
        alpha = alpha_pure, size = size_pure
      ) +
      scale_shape_manual(values = shape_map) +
      labs(
        title = paste0("NDVI–SWIR feature space (", yy, ")"),
        x = xlabel, y = ylabel, color = "Class", shape = "Class"
      ) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "right")
    
    ggsave(
      filename = file.path(out_dir, paste0("ndvi_swir_feature_space_", yy, "_", ifelse(use_z_space, "z", "raw"), ".png")),
      plot = p_y,
      width = 8, height = 6, dpi = dpi_out
    )
  }
}


# uniwue per coord_id
library(dplyr)
library(tidyr)

coords_unique <- LUCAS_multiyear_filtered %>%
  filter(!is.na(coord_id), coord_id != "") %>%
  count(point_id, coord_id, name = "n") %>%
  group_by(point_id) %>%
  slice_max(n, with_ties = FALSE) %>%
  ungroup() %>%
  separate(coord_id,
           into = c("gps_long", "gps_lat"),
           sep = "_",
           convert = TRUE,
           remove = FALSE) %>%   # <- key change
  select(point_id, coord_id, gps_long, gps_lat)


# df plot unique
df_plot_unique <- df_plot %>%
  arrange(point_id, desc(year)) %>%
  distinct(point_id, .keep_all = TRUE)




df_plot2 <- df_plot_unique %>%
  left_join(coords_unique, by = "point_id")

df_plot2_pure <- df_plot2 %>%
  filter(is_pure == TRUE)



df_plot_sf <- df_plot2_pure %>%
  filter(!is.na(coord_id), coord_id != "") %>%
  separate(
    coord_id,
    into = c("lon_coord", "lat_coord"),
    sep = "_",
    convert = TRUE,
    remove = FALSE
  ) %>%
  st_as_sf(
    coords = c("lon_coord", "lat_coord"),
    crs = 4326,      # WGS84
    remove = FALSE
  )


# Geometrie vorhanden?
st_geometry_type(df_plot_sf)

# CRS korrekt?
st_crs(df_plot_sf)

# ungültige Geometrien?
all(st_is_valid(df_plot_sf))

# crs to etrs
df_plot_sf <- st_transform(df_plot_sf, 3035)

# wrtir
st_write(
  df_plot_sf,
  "/mnt/eo/EU_unmixing/data/LUCAS/candidates.gpkg",
  layer = "pure_points",
  delete_layer = TRUE
)



