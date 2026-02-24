library(terra)
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)

# ---------------------------
# Inputs
# ---------------------------
gpkg_path   <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/candidates_vis_filtered_addmanual_cleaned.gpkg"
gpkg_layer  <- NULL  # set e.g. "layername" if needed; NULL reads first layer

bap_dir     <- "/mnt/dss_europe/mosaics_eu/mosaics_eu_baps"
out_csv     <- "/mnt/eo/EU_unmixing/data/LUCAS/candidates_BAPs.csv"

# ---------------------------
# Read points from GPKG
# ---------------------------
library(terra)
library(dplyr)
library(data.table)

# ===========================
# Paths / settings
# ===========================
gpkg_path  <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/candidates_vis_filtered_addmanual_cleaned.gpkg"
gpkg_layer <- NULL  # set e.g. "layername" if needed; NULL reads first layer

bap_dir <- "/mnt/dss_europe/mosaics_eu/mosaics_eu_baps"
out_csv <- "/mnt/eo/EU_unmixing/data/LUCAS/candidates_BAPs.csv"

# ===========================
# 1) Read points from GPKG
# ===========================
pts_v <- if (is.null(gpkg_layer)) {
  terra::vect(gpkg_path)
} else {
  terra::vect(gpkg_path, layer = gpkg_layer)
}

stopifnot(tolower(terra::geomtype(pts_v)) %in% c("points", "multipoints"))

# Attribute table
pts_df <- as.data.frame(pts_v)

# ===========================
# 2) Ensure there is a point_id attribute (prefer existing; else create)
#    IMPORTANT: keep point_id as CHARACTER for safe joins/extraction output
# ===========================
nms <- names(pts_df)
id_candidates <- c("point_id", "id", "ID", "fid", "FID", "objectid", "OBJECTID")
id_field <- id_candidates[id_candidates %in% nms][1]

if (is.na(id_field) || length(id_field) == 0) {
  pts_df$point_id <- as.character(seq_len(nrow(pts_df)))
} else {
  pts_df$point_id <- as.character(pts_df[[id_field]])
}

# If there are missing point_id values, fill ONLY those with new IDs continuing after max numeric ID
# (keeps existing IDs unchanged; uses numeric max if possible; falls back to 0 if none numeric)
pts_df <- pts_df %>%
  mutate(point_id_chr = point_id) %>%  # keep original
  mutate(point_id_num = suppressWarnings(as.numeric(point_id_chr)))

next_id <- max(pts_df$point_id_num, na.rm = TRUE)
if (!is.finite(next_id)) next_id <- 0

pts_df <- pts_df %>%
  mutate(
    point_id = if_else(
      is.na(point_id_chr) | point_id_chr == "" | point_id_chr == "NA",
      as.character(next_id + cumsum(is.na(point_id_chr) | point_id_chr == "" | point_id_chr == "NA")),
      point_id_chr
    )
  ) %>%
  select(-point_id_chr, -point_id_num)

# ===========================
# 3) Fill missing coordx/coordy (native CRS) + lon/lat (EPSG:4326)
#    using geometry from pts_v (row-aligned with pts_df)
# ===========================
coords_xy <- terra::geom(pts_v)[, c("x", "y")]
pts_v_ll  <- terra::project(pts_v, "EPSG:4326")
coords_ll <- terra::geom(pts_v_ll)[, c("x", "y")]

stopifnot(nrow(coords_xy) == nrow(pts_df), nrow(coords_ll) == nrow(pts_df))

# Ensure columns exist
if (!("coordx" %in% names(pts_df))) pts_df$coordx <- NA_real_
if (!("coordy" %in% names(pts_df))) pts_df$coordy <- NA_real_
if (!("lon"    %in% names(pts_df))) pts_df$lon    <- NA_real_
if (!("lat"    %in% names(pts_df))) pts_df$lat    <- NA_real_

# Fill only where NA
pts_df <- pts_df %>%
  mutate(
    coordx = if_else(is.na(coordx), coords_xy[, 1], coordx),
    coordy = if_else(is.na(coordy), coords_xy[, 2], coordy),
    lon    = if_else(is.na(lon),    coords_ll[, 1], lon),
    lat    = if_else(is.na(lat),    coords_ll[, 2], lat)
  )

# ===========================
# 4) Build a clean extraction vector from pts_df
#    (so that geometry + point_id are guaranteed consistent with your edited table)
#    Use coordx/coordy in native CRS.
# ===========================
pts_use_base <- terra::vect(
  pts_df,
  geom = c("coordx", "coordy"),
  crs  = terra::crs(pts_v)
)

# ===========================
# 5) List yearly VRTs
# ===========================
vrt_files <- list.files(
  bap_dir,
  pattern = "^[0-9]{4}_mosaic_eu_cog_v2\\.vrt$",
  full.names = TRUE
)
years <- as.integer(sub("^([0-9]{4}).*$", "\\1", basename(vrt_files)))

keep <- years >= 1984 & years <= 2024
vrt_files <- vrt_files[keep]
years     <- years[keep]
o <- order(years)
vrt_files <- vrt_files[o]
years     <- years[o]
stopifnot(length(vrt_files) > 0)

# ===========================
# 6) Project points once to CRS of first raster
# ===========================
r0 <- terra::rast(vrt_files[1])
pts_use0 <- pts_use_base
if (!terra::same.crs(pts_use0, r0)) pts_use0 <- terra::project(pts_use_base, terra::crs(r0))

# ===========================
# 7) Write header (overwrite existing)
# ===========================
data.table::fwrite(
  data.table(point_id = character(), year = integer(), band = integer(), BAP = numeric()),
  out_csv
)

# ===========================
# 8) Loop years and extract
# ===========================
for (i in seq_along(vrt_files)) {
  
  r <- terra::rast(vrt_files[i])
  
  # Safety: if CRS differs across years, reproject on the fly
  pts_use <- pts_use0
  if (!terra::same.crs(pts_use, r)) pts_use <- terra::project(pts_use_base, terra::crs(r))
  
  # Extract: returns a data.frame with 6 band columns (plus optionally an ID col, but ID=FALSE suppresses)
  vals <- terra::extract(r, pts_use, ID = FALSE)  # n_points x 6
  
  # Wide -> long
  out_i <- data.table(point_id = pts_use$point_id, vals)
  
  out_long <- data.table::melt(
    out_i,
    id.vars = "point_id",
    variable.name = "band_name",
    value.name = "BAP"
  )
  
  # Robust band index: parse trailing integer if present, otherwise use column order
  out_long[, band := suppressWarnings(as.integer(sub(".*_([0-9]+)$", "\\1", band_name)))]
  if (any(is.na(out_long$band))) {
    out_long[, band := as.integer(factor(band_name, levels = unique(band_name)))]
  }
  
  out_long[, band_name := NULL]
  out_long[, year := years[i]]
  
  data.table::setcolorder(out_long, c("point_id", "year", "band", "BAP"))
  data.table::fwrite(out_long, out_csv, append = TRUE)
  
  rm(r, vals, out_i, out_long)
  gc(FALSE)
}


library(data.table)

bap <- fread("/mnt/eo/EU_unmixing/data/LUCAS/candidates_BAPs.csv")


# --- sample 10 random points ---
set.seed(1)
pts10 <- sample(unique(bap$point_id), 10)

bap10 <- bap[point_id %in% pts10]

# --- (optional) map band -> wavelength for a more physical x-axis ---
band_cwl <- data.table(
  band = 1:6,
  cwl  = c(482, 561, 655, 865, 1610, 2200)
)

bap10 <- merge(bap10, band_cwl, by = "band", all.x = TRUE)

# --- plot: thin lines = yearly spectra; thick line = mean spectrum per point ---
ggplot(bap10, aes(x = cwl, y = BAP, group = interaction(point_id, year))) +
  geom_line(alpha = 0.15) +
  stat_summary(
    aes(group = point_id),
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  facet_wrap(~ point_id) +
  labs(
    x = "Wavelength (nm)",
    y = "BAP",
    title = "Spectral signatures for 10 random points (all years + mean)"
  ) +
  theme_bw()




library(terra)
library(dplyr)

gpkg_path <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/candidates_vis_filtered_addmanual_cleaned.gpkg"

# Read gpkg
pts_v <- vect(gpkg_path)

# Attribute table
pts_df <- as.data.frame(pts_v)

# Extract one coordinate pair per feature
coords <- geom(pts_v)[, c("x","y")]

# Ensure numeric IDs (or convert safely if stored as character)
pts_df$point_id <- as.numeric(pts_df$point_id)

# Find next available ID
next_id <- max(pts_df$point_id, na.rm = TRUE)

# Assign IDs only where missing
pts_df <- pts_df %>%
  mutate(
    point_id = if_else(
      is.na(point_id),
      next_id + cumsum(is.na(point_id)),
      point_id
    )
  )

sum(is.na(pts_df$point_id))


# --- choose columns to add from gpkg ---
pts_sel <- pts_df %>%
  select(
    point_id,
    coord_id,
    lc1,
    lc1_label,
    letter_group,
    lon,
    lat,
    re_checked,
    coordx,
    coordy,
    ECO_NAME,
    BIOME,
    biome_cor
  )

# IMPORTANT: remove missing keys and enforce one row per point_id on the RIGHT table
pts_sel_keyed <- pts_sel %>%
  filter(!is.na(point_id) & point_id != "") %>%
  distinct(point_id, .keep_all = TRUE)

# Join (keeps all LUCAS_BAP rows; adds attributes where keys match)

LUCAS_BAP <- bap
LUCAS_BAP <- LUCAS_BAP %>%
  mutate(point_id = as.character(point_id))

pts_sel_keyed <- pts_sel_keyed %>%
  mutate(point_id = as.character(point_id))

LUCAS_BAP_attr <- LUCAS_BAP %>%
  left_join(pts_sel_keyed, by = "point_id", relationship = "many-to-one")

# Sanity check: row count must remain unchanged
stopifnot(nrow(LUCAS_BAP_attr) == nrow(LUCAS_BAP))


### import soil data
library(readr)
library(dplyr)
library(tidyr)

LUCAS_topsoil_reflectance <- read_csv("/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/2015/final_datasets/topsoil_wide_subset_varaware.csv")

LUCAS_topsoil_reflectance_long <- LUCAS_topsoil_reflectance %>%
  pivot_longer(
    cols = starts_with("L_B"),
    names_to = "band",
    values_to = "BAP"
  )






library(dplyr)

LUCAS_topsoil_reflectance_long <- LUCAS_topsoil_reflectance_long %>%
  mutate(BAP = BAP * 10000)

LUCAS_topsoil_reflectance_long <- LUCAS_topsoil_reflectance_long %>%
  filter(!band %in% c("L_B1", "L_B8", "L_B9")) %>%              # remove unwanted bands
  mutate(
    band = as.integer(sub("L_B", "", band)) - 1                 # L_B2->1, ..., L_B7->6
  )

LUCAS_topsoil_reflectance <- LUCAS_topsoil_reflectance %>%
  mutate(letter_group = "F")

LUCAS_topsoil_reflectance <- LUCAS_topsoil_reflectance %>%
  mutate(biome_cor = NA_character_)

biome_pts <- st_read("/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/terrestrial_ecoregions_olson_europe.gpkg")

lucas_pts <- st_as_sf(
  LUCAS_topsoil_reflectance_long,
  coords = c("X", "Y"),
  remove = FALSE,
  crs = 4326
) %>%
  st_transform(st_crs(biome_pts))

lucas_with_biome <- st_join(
  lucas_pts,
  biome_pts %>% select(BIOME),
  join = st_within,
  left = TRUE
)

lucas_with_biome <- lucas_with_biome %>%
  mutate(
    biome_cor = case_when(
      BIOME %in% c(4, 5) ~ "Temperate",
      BIOME == 6 ~ "Boreal",
      BIOME == 12 ~ "Mediterranean",
      TRUE ~ NA_character_
    )
  )
lucas_with_biome <- lucas_with_biome %>%
  mutate(
    biome_cor = factor(
      biome_cor,
      levels = c("Boreal", "Temperate", "Mediterranean")
    )
  )

lucas_with_biome <- lucas_with_biome %>%
  rename(
    lon = X,
    lat = Y
  )


lucas_with_biome <- lucas_with_biome[, -c(2, 5)]

lucas_sf <- st_as_sf(
  lucas_with_biome,
  coords = c("lon", "lat"),
  remove = FALSE,
  crs = 4326
)

# Transform to meters (European LAEA projection)
lucas_sf_m <- st_transform(lucas_sf, 3035)

# Extract projected coordinates
coords_m <- st_coordinates(lucas_sf_m)

# Add to table
lucas_with_biome$xcoord <- coords_m[,1]
lucas_with_biome$ycoord <- coords_m[,2]

lucas_with_biome <- lucas_with_biome %>%
  mutate(letter_group = "F")


LUCAS_combined <- bind_rows(
  LUCAS_BAP_attr,
  lucas_with_biome
)


# Compute projected coordinates in meters (ETRS89 / LAEA Europe)
LUCAS_combined_sf <- st_as_sf(
  LUCAS_combined,
  coords = c("lon", "lat"),
  remove = FALSE,
  crs = 4326
) %>%
  st_transform(3035)

xy <- st_coordinates(LUCAS_combined_sf)

LUCAS_combined <- LUCAS_combined_sf %>%
  mutate(
    coordx = xy[, 1],
    coordy = xy[, 2]
  )

LUCAS_combined <- LUCAS_combined[, -c(23, 24)]

# Convert to sf
lucas_sf <- st_as_sf(
  LUCAS_combined,
  coords = c("lon", "lat"),
  remove = FALSE,
  crs = 4326
)

lucas_sf <- lucas_sf[, names(lucas_sf) != "LC1"]

# Write gpkg
st_write(
  lucas_sf,
  "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/LUCAS_time_series.gpkg",
  layer = "LUCAS_combined_full",
  delete_layer = TRUE
)


lucas_unique <- LUCAS_combined %>%
  distinct(point_id, .keep_all = TRUE)

lucas_unique_sf <- st_as_sf(
  lucas_unique,
  coords = c("lon", "lat"),
  remove = FALSE,
  crs = 4326
)

lucas_unique_sf <- lucas_unique_sf[, names(lucas_unique_sf) != "LC1"]

st_write(
  lucas_unique_sf,
  "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/LUCAS_unique.gpkg",
  layer = "LUCAS_combined_unique",
  delete_layer = TRUE
)


### plot
library(dplyr)
library(ggplot2)

set.seed(1)

# ---- 1) Filter to group F and years 1984–2024 ----
library(dplyr)
library(ggplot2)

set.seed(1)

# 1) filter
dat <- LUCAS_combined %>%
  mutate(
    year = as.integer(year),
    band = as.integer(band)
  ) %>%
  filter(
    letter_group == "F",
    year >= 1984, year <= 2024,
    !is.na(point_id), !is.na(year), !is.na(band), !is.na(BAP)
  )

# 2) sample 10 distinct point_ids
ids_10 <- dat %>%
  distinct(point_id) %>%
  slice_sample(n = 10) %>%
  pull(point_id)

dat10 <- dat %>%
  filter(point_id %in% ids_10)

# 3) plot: each year is one spectral signature (line across bands)
ggplot(dat10, aes(x = band, y = BAP, group = year)) +
  geom_line(alpha = 0.25) +
  # optional: show the band observations
  geom_point(alpha = 0.25, size = 0.6) +
  facet_wrap(~ point_id) +
  labs(
    x = "Band",
    y = "BAP",
    title = "Spectral signatures by year for 10 randomly sampled points (letter_group = F, 1984–2024)"
  ) +
  theme_bw()
