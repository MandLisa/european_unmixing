# ============================================================
# LUCAS candidate workflow: extract annual BAP spectra, enrich
# with point attributes, merge with topsoil spectra, write GPKG,
# and create quick-look spectral signature plots.
#
# Notes
# - Functionality is preserved; code is reorganized + commented.
# - Extraction is done with {terra}; vector joins/outputs use {sf}.
# - point_id is always treated as CHARACTER during joins/extraction
#   to avoid accidental numeric coercion issues.
# ============================================================

# ----------------------------
# Packages
# ----------------------------
suppressPackageStartupMessages({
  library(terra)      # raster + vector IO/extraction
  library(data.table) # fast IO + melt
  library(readr)      # read_csv
  library(dplyr)      # data wrangling
  library(tidyr)      # pivot_longer
  library(ggplot2)    # plotting
  library(sf)         # spatial joins + writing GPKG
})

# ----------------------------
# User inputs
# ----------------------------
gpkg_path   <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/candidates_vis_filtered_addmanual_cleaned.gpkg"
gpkg_layer  <- NULL  # if needed: set a layer name; NULL reads first layer

bap_dir     <- "/mnt/dss_europe/mosaics_eu/mosaics_eu_baps"
out_csv     <- "/mnt/eo/EU_unmixing/data/LUCAS/candidates_BAPs.csv"

topsoil_csv <- "/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/2015/final_datasets/topsoil_wide_subset_varaware.csv"
biome_gpkg  <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/terrestrial_ecoregions_olson_europe.gpkg"

out_gpkg_ts <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/LUCAS_time_series.gpkg"
out_gpkg_uq <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/1902/LUCAS_unique.gpkg"

years_min   <- 1984
years_max   <- 2024

# Sentinel-2-like BAP band center wavelengths (for nicer plots)
band_cwl <- data.table(
  band = 1:6,
  cwl  = c(482, 561, 655, 865, 1610, 2200)
)

# ============================================================
# 1) Read candidate points (terra) + create a robust point table
# ============================================================
pts_v <- if (is.null(gpkg_layer)) {
  terra::vect(gpkg_path)
} else {
  terra::vect(gpkg_path, layer = gpkg_layer)
}

stopifnot(tolower(terra::geomtype(pts_v)) %in% c("points", "multipoints"))

# Attribute table (row order is aligned with pts_v features)
pts_df <- as.data.frame(pts_v)

# ------------------------------------------------------------
# Ensure point_id exists and is usable for joins/extraction
# - Prefer an existing identifier column if present
# - Otherwise create sequential IDs
# - Fill missing/blank IDs without overwriting existing ones
# ------------------------------------------------------------
id_candidates <- c("point_id", "id", "ID", "fid", "FID", "objectid", "OBJECTID")
id_field <- id_candidates[id_candidates %in% names(pts_df)][1]

if (is.na(id_field) || length(id_field) == 0) {
  pts_df$point_id <- as.character(seq_len(nrow(pts_df)))
} else {
  pts_df$point_id <- as.character(pts_df[[id_field]])
}

# Fill only missing/empty point_id values with incremental IDs
pts_df <- pts_df %>%
  mutate(
    point_id_chr = point_id,
    point_id_num = suppressWarnings(as.numeric(point_id_chr))
  )

next_id <- max(pts_df$point_id_num, na.rm = TRUE)
if (!is.finite(next_id)) next_id <- 0

missing_key <- is.na(pts_df$point_id_chr) | pts_df$point_id_chr == "" | pts_df$point_id_chr == "NA"

pts_df <- pts_df %>%
  mutate(
    point_id = if_else(
      missing_key,
      as.character(next_id + cumsum(missing_key)),
      point_id_chr
    )
  ) %>%
  select(-point_id_chr, -point_id_num)

# ============================================================
# 2) Ensure coordinates exist (native CRS + lon/lat) using geometry
# ============================================================
coords_xy <- terra::geom(pts_v)[, c("x", "y")]

pts_v_ll  <- terra::project(pts_v, "EPSG:4326")
coords_ll <- terra::geom(pts_v_ll)[, c("x", "y")]

stopifnot(nrow(coords_xy) == nrow(pts_df), nrow(coords_ll) == nrow(pts_df))

# Create coord columns if absent; fill only missing values
if (!("coordx" %in% names(pts_df))) pts_df$coordx <- NA_real_
if (!("coordy" %in% names(pts_df))) pts_df$coordy <- NA_real_
if (!("lon"    %in% names(pts_df))) pts_df$lon    <- NA_real_
if (!("lat"    %in% names(pts_df))) pts_df$lat    <- NA_real_

pts_df <- pts_df %>%
  mutate(
    coordx = if_else(is.na(coordx), coords_xy[, 1], coordx),
    coordy = if_else(is.na(coordy), coords_xy[, 2], coordy),
    lon    = if_else(is.na(lon),    coords_ll[, 1], lon),
    lat    = if_else(is.na(lat),    coords_ll[, 2], lat)
  )

# Build a clean vector for extraction from the edited attribute table
# (uses coordx/coordy in the original CRS of pts_v)
pts_use_base <- terra::vect(
  pts_df,
  geom = c("coordx", "coordy"),
  crs  = terra::crs(pts_v)
)

# ============================================================
# 3) List annual BAP VRTs and extract point spectra year-by-year
# ============================================================
vrt_files <- list.files(
  bap_dir,
  pattern = "^[0-9]{4}_mosaic_eu_cog_v2\\.vrt$",
  full.names = TRUE
)

years <- as.integer(sub("^([0-9]{4}).*$", "\\1", basename(vrt_files)))

keep <- years >= years_min & years <= years_max
vrt_files <- vrt_files[keep]
years     <- years[keep]

o <- order(years)
vrt_files <- vrt_files[o]
years     <- years[o]

stopifnot(length(vrt_files) > 0)

# Project points once to CRS of the first raster (fast path)
r0 <- terra::rast(vrt_files[1])
pts_use0 <- pts_use_base
if (!terra::same.crs(pts_use0, r0)) {
  pts_use0 <- terra::project(pts_use_base, terra::crs(r0))
}

# Write header (overwrites any existing output file)
data.table::fwrite(
  data.table(point_id = character(), year = integer(), band = integer(), BAP = numeric()),
  out_csv
)

# Extraction loop (stream to CSV to keep memory low)
for (i in seq_along(vrt_files)) {
  
  r <- terra::rast(vrt_files[i])
  
  # Safety: if CRS differs across years, project on the fly
  pts_use <- pts_use0
  if (!terra::same.crs(pts_use, r)) {
    pts_use <- terra::project(pts_use_base, terra::crs(r))
  }
  
  # Extract at points; returns one row per point with 6 band columns
  vals <- terra::extract(r, pts_use, ID = FALSE)
  
  # Convert to long format: one row per (point_id, band, year)
  out_i <- data.table(point_id = pts_use$point_id, vals)
  
  out_long <- data.table::melt(
    out_i,
    id.vars       = "point_id",
    variable.name = "band_name",
    value.name    = "BAP"
  )
  
  # Robust band index:
  # - If column names end with _<int>, parse that
  # - Else fall back to their column order
  out_long[, band := suppressWarnings(as.integer(sub(".*_([0-9]+)$", "\\1", band_name)))]
  if (any(is.na(out_long$band))) {
    out_long[, band := as.integer(factor(band_name, levels = unique(band_name)))]
  }
  
  out_long[, band_name := NULL]
  out_long[, year := years[i]]
  data.table::setcolorder(out_long, c("point_id", "year", "band", "BAP"))
  
  # Append to disk
  data.table::fwrite(out_long, out_csv, append = TRUE)
  
  # Housekeeping
  rm(r, vals, out_i, out_long)
  gc(FALSE)
}

# ============================================================
# 4) Quick-look plot: 10 random points, all years + mean spectrum
# ============================================================
bap <- data.table::fread(out_csv)

set.seed(1)
pts10 <- sample(unique(bap$point_id), 10)
bap10 <- bap[point_id %in% pts10]

# Optional: wavelength axis
bap10 <- merge(bap10, band_cwl, by = "band", all.x = TRUE)

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

# ============================================================
# 5) Enrich extracted time series with selected GPKG attributes
#    (many-to-one: many BAP rows per point_id, one attribute row)
# ============================================================
# Select attribute columns (keep exactly one row per point_id)
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

pts_sel_keyed <- pts_sel %>%
  filter(!is.na(point_id) & point_id != "") %>%
  distinct(point_id, .keep_all = TRUE)

# Join into long BAP table; keep row count unchanged
LUCAS_BAP <- as_tibble(bap) %>% mutate(point_id = as.character(point_id))
pts_sel_keyed <- pts_sel_keyed %>% mutate(point_id = as.character(point_id))

LUCAS_BAP_attr <- LUCAS_BAP %>%
  left_join(pts_sel_keyed, by = "point_id", relationship = "many-to-one")

stopifnot(nrow(LUCAS_BAP_attr) == nrow(LUCAS_BAP))

# ============================================================
# 6) Import topsoil spectra and convert to compatible long format
# ============================================================
LUCAS_topsoil_reflectance <- readr::read_csv(topsoil_csv, show_col_types = FALSE)

LUCAS_topsoil_reflectance_long <- LUCAS_topsoil_reflectance %>%
  pivot_longer(
    cols      = starts_with("L_B"),
    names_to  = "band",
    values_to = "BAP"
  ) %>%
  mutate(
    # Match scaling used elsewhere
    BAP = BAP * 10000
  ) %>%
  # Remove unwanted bands and map to 1..6
  filter(!band %in% c("L_B1", "L_B8", "L_B9")) %>%
  mutate(
    band = as.integer(sub("L_B", "", band)) - 1
  )

# Add fields required downstream (consistent schema with LUCAS_BAP_attr)
LUCAS_topsoil_reflectance_long <- LUCAS_topsoil_reflectance_long %>%
  mutate(
    letter_group = "F",
    biome_cor    = NA_character_
  )

# ============================================================
# 7) Assign biome class to topsoil points via spatial join
# ============================================================
biome_pts <- sf::st_read(biome_gpkg, quiet = TRUE)

lucas_pts <- sf::st_as_sf(
  LUCAS_topsoil_reflectance_long,
  coords = c("X", "Y"),
  remove = FALSE,
  crs = 4326
) %>%
  sf::st_transform(sf::st_crs(biome_pts))

lucas_with_biome <- sf::st_join(
  lucas_pts,
  biome_pts %>% dplyr::select(BIOME),
  join = sf::st_within,
  left = TRUE
) %>%
  mutate(
    biome_cor = case_when(
      BIOME %in% c(4, 5) ~ "Temperate",
      BIOME == 6         ~ "Boreal",
      BIOME == 12        ~ "Mediterranean",
      TRUE               ~ NA_character_
    ),
    biome_cor = factor(biome_cor, levels = c("Boreal", "Temperate", "Mediterranean"))
  ) %>%
  rename(lon = X, lat = Y)

# Drop a couple of columns by position (kept as in original to preserve behavior)
# If you want this safer: replace with explicit column names.
lucas_with_biome <- lucas_with_biome[, -c(2, 5)]

# Add projected coordinates (ETRS89 / LAEA Europe, EPSG:3035)
lucas_sf <- sf::st_as_sf(
  lucas_with_biome,
  coords = c("lon", "lat"),
  remove = FALSE,
  crs = 4326
)

lucas_sf_m <- sf::st_transform(lucas_sf, 3035)
coords_m <- sf::st_coordinates(lucas_sf_m)

lucas_with_biome$xcoord <- coords_m[, 1]
lucas_with_biome$ycoord <- coords_m[, 2]

lucas_with_biome <- lucas_with_biome %>%
  mutate(letter_group = "F")

# ============================================================
# 8) Combine time series (candidates) + topsoil spectra
# ============================================================
LUCAS_combined <- bind_rows(
  LUCAS_BAP_attr,
  lucas_with_biome
)

# Ensure coordx/coordy exist in EPSG:3035 derived from lon/lat
LUCAS_combined_sf <- sf::st_as_sf(
  LUCAS_combined,
  coords = c("lon", "lat"),
  remove = FALSE,
  crs = 4326
) %>%
  sf::st_transform(3035)

xy <- sf::st_coordinates(LUCAS_combined_sf)

LUCAS_combined <- LUCAS_combined_sf %>%
  mutate(
    coordx = xy[, 1],
    coordy = xy[, 2]
  )

# Drop columns by position (as in original; consider making explicit later)
LUCAS_combined <- LUCAS_combined[, -c(23, 24)]

# ============================================================
# 9) Write outputs: full time series + unique points
# ============================================================
lucas_sf_out <- sf::st_as_sf(
  LUCAS_combined,
  coords = c("lon", "lat"),
  remove = FALSE,
  crs = 4326
)

# Remove LC1 column if present (case-sensitive, per original)
lucas_sf_out <- lucas_sf_out[, names(lucas_sf_out) != "LC1"]

sf::st_write(
  lucas_sf_out,
  out_gpkg_ts,
  layer = "LUCAS_combined_full",
  delete_layer = TRUE,
  quiet = TRUE
)

# One row per point_id (keep first occurrence)
lucas_unique <- LUCAS_combined %>%
  distinct(point_id, .keep_all = TRUE)

lucas_unique_sf <- sf::st_as_sf(
  lucas_unique,
  coords = c("lon", "lat"),
  remove = FALSE,
  crs = 4326
)

lucas_unique_sf <- lucas_unique_sf[, names(lucas_unique_sf) != "LC1"]

sf::st_write(
  lucas_unique_sf,
  out_gpkg_uq,
  layer = "LUCAS_combined_unique",
  delete_layer = TRUE,
  quiet = TRUE
)

# ============================================================
# 10) Plot: spectral signatures by year for 10 random points
# ============================================================
set.seed(1)

dat <- LUCAS_combined %>%
  mutate(
    year = as.integer(year),
    band = as.integer(band)
  ) %>%
  filter(
    letter_group == "F",
    year >= years_min, year <= years_max,
    !is.na(point_id), !is.na(year), !is.na(band), !is.na(BAP)
  )

ids_10 <- dat %>%
  distinct(point_id) %>%
  slice_sample(n = 10) %>%
  pull(point_id)

dat10 <- dat %>%
  filter(point_id %in% ids_10)

ggplot(dat10, aes(x = band, y = BAP, group = year)) +
  geom_line(alpha = 0.25) +
  geom_point(alpha = 0.25, size = 0.6) +
  facet_wrap(~ point_id) +
  labs(
    x = "Band",
    y = "BAP",
    title = sprintf(
      "Spectral signatures by year for 10 randomly sampled points (letter_group = F, %d–%d)",
      years_min, years_max
    )
  ) +
  theme_bw()