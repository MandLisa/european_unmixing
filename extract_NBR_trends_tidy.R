# ============================================================
# LUCAS workflow: NBR Sen slope extraction + join + plotting
# + spectral-signature visualizations (by group/biome/year)
# ============================================================

# -----------------------
# Libraries
# -----------------------
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# A) NBR TREND (Sen slope) from yearly rasters (1984–2024)
# ============================================================

# -----------------------
# User inputs
# -----------------------
nbr_dir     <- "/mnt/eo/eu_mosaics/NBR_comp"
year_min    <- 1984L
year_max    <- 2024L
chunk_size  <- 20000L   # RAM knob: 5000 if tight, 50000 if plenty
min_years   <- 5L       # require at least this many valid annual values
nodata_thr  <- -9990    # treat values <= this as NA (adjust if needed)

# -----------------------
# 0) Checks
# -----------------------
stopifnot(exists("lucas_unique_sf"))
stopifnot(inherits(lucas_unique_sf, "sf"))

# -----------------------
# 1) sf -> terra points + stable id
# -----------------------
pts <- terra::vect(lucas_unique_sf)

if (!("point_id" %in% names(pts))) {
  pts$point_id <- as.character(seq_len(nrow(pts)))
} else {
  pts$point_id <- as.character(pts$point_id)
}

# -----------------------
# 2) List & sort yearly NBR rasters
# -----------------------
nbr_files <- list.files(nbr_dir, pattern = "^NBR_[0-9]{4}\\.tif$", full.names = TRUE)
stopifnot(length(nbr_files) > 0)

yrs <- as.integer(sub("^.*NBR_([0-9]{4})\\.tif$", "\\1", nbr_files))
o   <- order(yrs)
nbr_files <- nbr_files[o]
yrs       <- yrs[o]

keep <- yrs >= year_min & yrs <= year_max
nbr_files <- nbr_files[keep]
yrs       <- yrs[keep]
stopifnot(length(nbr_files) >= min_years)

# -----------------------
# 3) Project points once to raster CRS
# -----------------------
r0   <- rast(nbr_files[1])
ptsR <- project(pts, crs(r0))

# -----------------------
# 4) Robust Sen slope (median of pairwise slopes; slope only)
# -----------------------
sen_slope_only <- function(y, x, min_years = 5L) {
  ok <- is.finite(y) & is.finite(x)
  y  <- y[ok]
  x  <- x[ok]
  n  <- length(y)
  
  if (n < min_years) return(NA_real_)
  if (all(y == y[1])) return(0)  # constant series
  
  # Pairwise slopes
  s <- numeric(0)
  for (i in 1:(n - 1)) {
    dx <- x[(i + 1):n] - x[i]
    s  <- c(s, (y[(i + 1):n] - y[i]) / dx)
  }
  stats::median(s, na.rm = TRUE)
}

# -----------------------
# 5) Chunked extraction + slope computation
# -----------------------
n_pts   <- nrow(ptsR)
n_years <- length(nbr_files)

out_slope <- rep(NA_real_, n_pts)
out_nobs  <- rep(0L,       n_pts)

starts <- seq.int(1L, n_pts, by = chunk_size)

for (s0 in starts) {
  s1 <- min(s0 + chunk_size - 1L, n_pts)
  pts_chunk <- ptsR[s0:s1]
  
  # time series matrix for this chunk
  Y <- matrix(NA_real_, nrow = nrow(pts_chunk), ncol = n_years)
  
  for (j in seq_along(nbr_files)) {
    r <- rast(nbr_files[j])
    v <- terra::extract(r, pts_chunk)[, 2]
    
    v[!is.finite(v)] <- NA_real_
    v[v <= nodata_thr] <- NA_real_
    
    Y[, j] <- v
  }
  
  n_ok <- rowSums(is.finite(Y))
  out_nobs[s0:s1] <- n_ok
  
  out_slope[s0:s1] <- apply(
    Y, 1, sen_slope_only,
    x = yrs, min_years = min_years
  )
  
  message(sprintf("Done %d-%d / %d points", s0, s1, n_pts))
}

# -----------------------
# 6) Attach results back to sf object
# -----------------------
lucas_unique_sf$nbr_sen_slope <- out_slope
lucas_unique_sf$nbr_n_years   <- out_nobs

# -----------------------
# 7) Optional: rescale if stored as scaled int (e.g., *10000)
# -----------------------
# lucas_unique_sf$nbr_sen_slope <- lucas_unique_sf$nbr_sen_slope / 10000

# -----------------------
# 8) Quick diagnostics
# -----------------------
print(summary(lucas_unique_sf$nbr_n_years))
print(summary(lucas_unique_sf$nbr_sen_slope))

# ============================================================
# B) WRITE NBR-TREND POINTS + JOIN TREND TO FULL TIME SERIES
# ============================================================

# -----------------------
# Write point layer with trend
# -----------------------
st_write(
  lucas_unique_sf,
  "/mnt/dss_project/lmandl/_unmixing/_spectral_library/2502/candidates_with_NBRtrend.gpkg",
  layer = "lucas_unique_sf",
  delete_layer = TRUE
)

# -----------------------
# Join trend back to your full LUCAS time-series sf
# Requires: LUCAS_combined_sf in env
# -----------------------
stopifnot(exists("LUCAS_combined_sf"))
stopifnot(inherits(LUCAS_combined_sf, "sf"))

lucas_unique_sf <- lucas_unique_sf %>%
  mutate(point_id = as.character(point_id))

LUCAS_combined_sf <- LUCAS_combined_sf %>%
  mutate(point_id = as.character(point_id))

trend_lut <- lucas_unique_sf %>%
  st_drop_geometry() %>%
  select(point_id, nbr_sen_slope) %>%
  distinct(point_id, .keep_all = TRUE)

LUCAS_combined_sf_NBR <- LUCAS_combined_sf %>%
  left_join(trend_lut, by = "point_id") %>%
  select(-LC1)

st_write(
  LUCAS_combined_sf_NBR,
  "/mnt/dss_project/lmandl/_unmixing/_spectral_library/2502/lucas_time_series_NBRTrend.gpkg",
  layer = "lucas_full",
  delete_layer = TRUE
)

# ============================================================
# C) QUICK MAPS (letter groups; NBR slope)
# ============================================================

p_letter <- ggplot(lucas_unique_sf) +
  geom_sf(aes(color = letter_group), size = 1) +
  theme_bw() +
  labs(
    title = "Distribution of LUCAS points by letter group",
    color = "Letter group"
  )

p_trend <- ggplot(lucas_unique_sf) +
  geom_sf(
    aes(fill = nbr_sen_slope),
    shape  = 21,
    color  = "black",
    stroke = 0.2,
    size   = 2
  ) +
  scale_fill_gradient2(
    low = "orange",
    mid = "white",
    high = "#135861",
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "NBR trend (Sen slope, 1984–2024)",
    fill  = "NBR trend"
  )

print(p_letter)
print(p_trend)

# ============================================================
# D) SPECTRAL SIGNATURES (BAP) — SUMMARIES & SAMPLES
# ============================================================
# Assumes your long table has (at least):
# point_id, band, BAP, year, letter_group, biome_cor
# Here we use: LUCAS_combined_sf_NBR (sf) as main source.

# -----------------------
# Helper: enforce clean plotting columns
# -----------------------
prep_ts_all <- function(x_sf) {
  x_sf %>%
    st_drop_geometry() %>%
    mutate(
      point_id     = as.character(point_id),
      year         = as.integer(year),
      band         = as.integer(band),
      letter_group = as.character(letter_group),
      biome_cor    = as.character(biome_cor),
      biome_cor    = ifelse(biome_cor == "Alpine", "Temperate", biome_cor)
    ) %>%
    filter(!is.na(biome_cor)) %>%
    mutate(
      letter_group = factor(letter_group, levels = c("A", "C", "D", "E", "F")),
      biome_cor    = factor(biome_cor, levels = c("Boreal", "Mediterranean", "Temperate"))
    )
}

ts_all <- prep_ts_all(LUCAS_combined_sf_NBR)

# -----------------------
# 1) Simple class summary: median + (0.1, 0.9) ribbons by letter_group
# -----------------------
spectra_summary <- ts_all %>%
  group_by(letter_group, band) %>%
  summarise(
    median_BAP = median(BAP, na.rm = TRUE),
    q10        = quantile(BAP, 0.1, na.rm = TRUE),
    q90        = quantile(BAP, 0.9, na.rm = TRUE),
    .groups    = "drop"
  )

p_summary <- ggplot(spectra_summary, aes(x = band, y = median_BAP, color = letter_group)) +
  geom_ribbon(
    aes(ymin = q10, ymax = q90, fill = letter_group),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    x = "Band",
    y = "BAP reflectance",
    color = "Letter group",
    fill  = "Letter group"
  )

print(p_summary)

# -----------------------
# 2) Thin point spectra + thick class median (facet by letter_group)
#    Collapses years: one spectrum per point_id (median across years)
# -----------------------
spectra_point <- ts_all %>%
  group_by(letter_group, point_id, band) %>%
  summarise(BAP = median(BAP, na.rm = TRUE), .groups = "drop")

# Sample N point_ids per class
set.seed(1)
n_per_class <- 100

sample_ids <- spectra_point %>%
  distinct(letter_group, point_id) %>%
  group_by(letter_group) %>%
  group_modify(~ slice_sample(.x, n = min(n_per_class, nrow(.x)))) %>%
  ungroup()

spectra_sample <- spectra_point %>%
  inner_join(sample_ids, by = c("letter_group", "point_id"))

spectra_class <- spectra_point %>%
  group_by(letter_group, band) %>%
  summarise(med = median(BAP, na.rm = TRUE), .groups = "drop")

p_facet_groups <- ggplot() +
  geom_line(
    data = spectra_sample,
    aes(x = band, y = BAP, group = point_id),
    linewidth = 0.4,
    alpha = 0.2,
    color = "#135861"
  ) +
  geom_line(
    data = spectra_class,
    aes(x = band, y = med, group = letter_group),
    linewidth = 1.0,
    color = "black"
  ) +
  facet_wrap(~ letter_group) +
  theme_bw() +
  labs(x = "Band", y = "Reflectance", title = NULL)

print(p_facet_groups)

# -----------------------
# 3) By biome × letter group: sample thin spectra + thick median
#    Ensure complete 6-band spectra per point (n_distinct(band) == 6)
# -----------------------
spectra_point_biome <- ts_all %>%
  group_by(biome_cor, letter_group, point_id, band) %>%
  summarise(BAP = median(BAP, na.rm = TRUE), .groups = "drop")

spectra_point_complete <- spectra_point_biome %>%
  group_by(biome_cor, letter_group, point_id) %>%
  filter(n_distinct(band) == 6) %>%
  ungroup()

set.seed(1)
n_per_class_biome <- 50

sample_ids_biome <- spectra_point_complete %>%
  distinct(biome_cor, letter_group, point_id) %>%
  group_by(biome_cor, letter_group) %>%
  group_modify(~ slice_sample(.x, n = min(n_per_class_biome, nrow(.x)))) %>%
  ungroup()

spectra_sample_biome <- spectra_point_complete %>%
  inner_join(sample_ids_biome, by = c("biome_cor", "letter_group", "point_id"))

spectra_class_biome <- spectra_point_complete %>%
  group_by(biome_cor, letter_group, band) %>%
  summarise(med = median(BAP, na.rm = TRUE), .groups = "drop")

p_biome_grid <- ggplot() +
  geom_line(
    data = spectra_sample_biome,
    aes(band, BAP, group = interaction(biome_cor, letter_group, point_id)),
    linewidth = 0.4,
    alpha = 0.25,
    color = "#135861"
  ) +
  geom_line(
    data = spectra_class_biome,
    aes(band, med, group = interaction(biome_cor, letter_group)),
    linewidth = 1,
    color = "black"
  ) +
  facet_grid(biome_cor ~ letter_group) +
  theme_bw() +
  labs(
    x = "Band",
    y = "BAP",
    title = "Spectral signatures by biome and letter group"
  )

print(p_biome_grid)

# ============================================================
# E) YEAR-SPECIFIC SIGNATURES (1990/2000/2010/2020)
#    Option: keep F fixed across years (your “inject F” logic)
# ============================================================

years_sel <- c(1990L, 2000L, 2010L, 2020L)

# Precompute "fixed F": median across years per point for group F
F_fixed <- ts_all %>%
  filter(letter_group == "F") %>%
  group_by(biome_cor, letter_group, point_id, band) %>%
  summarise(BAP = median(BAP, na.rm = TRUE), .groups = "drop")

plot_spectra_year_with_fixedF <- function(ts_all_df, year_label, F_fixed_df, n_per_class = 50, seed = 1) {
  
  df_year <- ts_all_df %>%
    filter(year == year_label, letter_group != "F") %>%
    select(biome_cor, letter_group, point_id, band, BAP)
  
  df_F <- F_fixed_df %>%
    mutate(year = year_label) %>%
    select(biome_cor, letter_group, point_id, band, BAP)
  
  df_plot <- bind_rows(df_year, df_F)
  
  # keep complete 6-band spectra
  df_complete <- df_plot %>%
    group_by(biome_cor, letter_group, point_id) %>%
    filter(n_distinct(band) == 6) %>%
    ungroup()
  
  set.seed(seed)
  sample_ids <- df_complete %>%
    distinct(biome_cor, letter_group, point_id) %>%
    group_by(biome_cor, letter_group) %>%
    group_modify(~ slice_sample(.x, n = min(n_per_class, nrow(.x)))) %>%
    ungroup()
  
  df_sample <- df_complete %>%
    semi_join(sample_ids, by = c("biome_cor", "letter_group", "point_id"))
  
  df_med <- df_complete %>%
    group_by(biome_cor, letter_group, band) %>%
    summarise(med = median(BAP, na.rm = TRUE), .groups = "drop")
  
  ggplot() +
    geom_line(
      data = df_sample,
      aes(band, BAP, group = interaction(biome_cor, letter_group, point_id)),
      linewidth = 0.35,
      alpha = 0.20,
      color = "#135861"
    ) +
    geom_line(
      data = df_med,
      aes(band, med, group = interaction(biome_cor, letter_group)),
      linewidth = 1.0,
      color = "black"
    ) +
    facet_grid(biome_cor ~ letter_group) +
    theme_bw() +
    labs(
      x = "Band",
      y = "BAP",
      title = paste0("Spectral signatures by biome and letter group — ", year_label)
    )
}

plots_years <- lapply(
  years_sel,
  function(y) plot_spectra_year_with_fixedF(ts_all, y, F_fixed, n_per_class = 50, seed = 1)
)

# Print (or save) year plots
print(plots_years[[1]])
print(plots_years[[2]])
print(plots_years[[3]])
print(plots_years[[4]])