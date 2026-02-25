# ============================================================
# Extract per-point NBR Sen slope (1984–2024) from yearly rasters
# - Uses lucas_unique_sf (sf POINT object) already in your environment
# - Memory-efficient: reads one raster at a time, processes points in chunks
# - Robust: handles nodata, constant series, and avoids qnorm/CI warnings
# ============================================================

library(sf)
library(terra)

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
# 1) sf -> terra points
# -----------------------
pts <- terra::vect(lucas_unique_sf)

# stable id (for joining/debugging)
if (!("point_id" %in% names(pts))) {
  pts$point_id <- as.character(seq_len(nrow(pts)))
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

# restrict to target period
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
#    Avoids CI/p-value computations that trigger qnorm warnings
# -----------------------
sen_slope_only <- function(y, x, min_years = 5L) {
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]; x <- x[ok]
  n <- length(y)
  if (n < min_years) return(NA_real_)
  if (all(y == y[1])) return(0)  # constant time series => slope 0
  
  # pairwise slopes (n ~ 41 max here; fine)
  s <- c()
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
out_nobs  <- rep(0L,       n_pts)  # number of valid annual obs used

starts <- seq(1L, n_pts, by = chunk_size)

for (s0 in starts) {
  s1 <- min(s0 + chunk_size - 1L, n_pts)
  pts_chunk <- ptsR[s0:s1]
  
  # time series matrix for this chunk only
  Y <- matrix(NA_real_, nrow = nrow(pts_chunk), ncol = n_years)
  
  for (j in seq_along(nbr_files)) {
    r <- rast(nbr_files[j])
    
    v <- terra::extract(r, pts_chunk)[, 2]
    
    # nodata handling (adjust threshold to your product if needed)
    v[!is.finite(v)] <- NA_real_
    v[v <= nodata_thr] <- NA_real_
    
    Y[, j] <- v
  }
  
  # count valid years per point (for diagnostics / filtering)
  n_ok <- rowSums(is.finite(Y))
  out_nobs[s0:s1] <- n_ok
  
  # compute Sen slope per point
  out_slope[s0:s1] <- apply(Y, 1, sen_slope_only, x = yrs, min_years = min_years)
  
  message(sprintf("Done %d-%d / %d points", s0, s1, n_pts))
}

# -----------------------
# 6) Attach results back to your sf object
# -----------------------
lucas_unique_sf$nbr_sen_slope <- out_slope
lucas_unique_sf$nbr_n_years   <- out_nobs

# -----------------------
# 7) Optional: rescale if your NBR is stored as scaled int (e.g., *10000)
#    Uncomment if needed.
# -----------------------
# lucas_unique_sf$nbr_sen_slope <- lucas_unique_sf$nbr_sen_slope / 10000

# -----------------------
# 8) Optional quick diagnostics
# -----------------------
summary(lucas_unique_sf$nbr_n_years)
summary(lucas_unique_sf$nbr_sen_slope)

# If many points have low nbr_n_years, your nodata_thr may be wrong.
# To inspect value range in a sample raster:
# global(r0, range, na.rm=TRUE)

library(sf)

st_write(
  lucas_unique_sf,
  "/mnt/dss_project/lmandl/_unmixing/_spectral_library/2502/candidates_with_NBRtrend.gpkg",
  layer = "lucas_unique_sf",
  delete_layer = TRUE
)


library(sf)
library(dplyr)

lucas_unique_sf <- lucas_unique_sf %>%
  mutate(point_id = as.character(point_id))

LUCAS_combined_sf <- LUCAS_combined_sf %>%
  mutate(point_id = as.character(point_id))


trend_lut <- lucas_unique_sf %>%
  st_drop_geometry() %>%
  select(point_id, nbr_sen_slope) %>%
  distinct(point_id, .keep_all = TRUE)


LUCAS_combined_sf_NBR <- LUCAS_combined_sf %>%
  left_join(trend_lut, by = "point_id")



LUCAS_combined_sf_NBR <- LUCAS_combined_sf_NBR %>%
  select(-LC1)



st_write(
  LUCAS_combined_sf_NBR,
  "/mnt/dss_project/lmandl/_unmixing/_spectral_library/2502/lucas_time_series_NBRTrend.gpkg",
  layer = "lucas_full",
  delete_layer = TRUE
)



library(ggplot2)
library(sf)

ggplot(lucas_unique_sf) +
  geom_sf(aes(color = letter_group), size = 1) +
  theme_bw() +
  labs(
    title = "Distribution of LUCAS Points by Letter Group",
    color = "Letter group"
  )



ggplot(lucas_unique_sf) +
  geom_sf(
    aes(fill = nbr_sen_slope),
    shape = 21,          # filled circle with outline
    color = "black",     # outline color
    stroke = 0.2,        # outline thickness (slim)
    size = 2
  ) +
  scale_fill_gradient2(
    low = "orange",
    mid = "white",
    high = "#135861",
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "NBR Trend (Sen Slope 1984–2024)",
    fill = "NBR trend"
  )



### visualise signatures

library(dplyr)

spectra_summary <- LUCAS_combined %>%
  group_by(letter_group, band) %>%
  summarise(
    median_BAP = median(BAP, na.rm=TRUE),
    q10 = quantile(BAP, 0.1, na.rm=TRUE),
    q90 = quantile(BAP, 0.9, na.rm=TRUE),
    .groups="drop"
  )


library(ggplot2)

ggplot(spectra_summary,
       aes(x = band, y = median_BAP, color = letter_group)) +
  
  geom_ribbon(
    aes(ymin = q10, ymax = q90, fill = letter_group),
    alpha = 0.2,
    color = NA
  ) +
  
  geom_line(size = 1) +
  
  theme_bw() +
  
  labs(
    x = "Band",
    y = "BAP reflectance",
    color = "Letter group",
    fill = "Letter group"
  )


library(dplyr)
library(ggplot2)

library(dplyr)
library(ggplot2)

# lucas_BAP must be long format with: point_id, letter_group, band, BAP, year (or similar)
# 1) Make sure band is ordered numeric (important for correct line drawing)
LUCAS_combined_sf_NBRd <- LUCAS_combined_sf_NBR %>%
  mutate(band = as.integer(band))

# 2) Collapse repeated years: one spectrum per point_id (median across years)
spectra_point <- LUCAS_combined_sf_NBR %>%
  group_by(letter_group, point_id, band) %>%
  summarise(BAP = median(BAP, na.rm = TRUE), .groups = "drop")

# 3) Sample 5–10 point spectra per class
n_per_class <- 100
set.seed(1)

sample_ids <- spectra_point %>%
  distinct(letter_group,point_id) %>%
  group_by(letter_group) %>%
  group_modify(~ slice_sample(.x, n = min(n_per_class, nrow(.x)))) %>%
  ungroup()

spectra_sample <- spectra_point %>%
  inner_join(sample_ids, by = c("letter_group", "point_id"))

# 4) Class-level summary (median across points) for a thick reference line
spectra_class <- spectra_point %>%
  group_by(letter_group, band) %>%
  summarise(med = median(BAP, na.rm = TRUE), .groups = "drop")

# 5) Plot: thin raw point spectra + thick class median, facetted
ggplot() +
  geom_line(
    data = spectra_sample,
    aes(x = band, y = BAP, group =point_id),
    linewidth = 0.4,
    alpha = 0.2,
    color = "#135861"
  ) +
  geom_line(
    data = spectra_class,
    aes(x = band, y = med),
    linewidth = 1.0,
    color = "black"
  ) +
  facet_wrap(~letter_group) +
  theme_bw() +
  labs(x = "Band", y = "Reflectance", title = "")




library(dplyr)

spectra_point <- LUCAS_combined_sf_NBR %>%
  mutate(band = as.integer(band)) %>%
  group_by(biome_cor, letter_group, point_id, band) %>%
  summarise(
    BAP = median(BAP, na.rm = TRUE),
    .groups = "drop"
  )


spectra_point_complete <- spectra_point %>%
  group_by(biome_cor, letter_group, point_id) %>%
  filter(n_distinct(band) == 6) %>%
  ungroup()


set.seed(1)
n_per_class <- 50

sample_ids <- spectra_point_complete %>%
  distinct(biome_cor, letter_group, point_id) %>%
  group_by(biome_cor, letter_group) %>%
  group_modify(~ slice_sample(.x, n = min(n_per_class, nrow(.x)))) %>%
  ungroup()

spectra_sample <- spectra_point_complete %>%
  inner_join(sample_ids,
             by = c("biome_cor","letter_group","point_id"))


spectra_class <- spectra_point_complete %>%
  group_by(biome_cor, letter_group, band) %>%
  summarise(
    med = median(BAP, na.rm = TRUE),
    .groups = "drop"
  )

library(dplyr)

spectra_sample_clean <- spectra_sample %>%
  mutate(
    biome_cor = ifelse(biome_cor == "Alpine", "Temperate", biome_cor)
  ) %>%
  filter(!is.na(biome_cor))


spectra_sample_clean$biome_cor <- factor(
  spectra_sample_clean$biome_cor,
  levels = c("Boreal","Mediterranean","Temperate")
)

spectra_class_clean$biome_cor <- factor(
  spectra_class_clean$biome_cor,
  levels = c("Boreal","Mediterranean","Temperate")
)

ggplot() +
  
  geom_line(
    data = spectra_sample_clean %>%
      filter(!is.na(biome_cor)),
    aes(band, BAP,
        group = interaction(biome_cor, letter_group, point_id)),
    linewidth = 0.4,
    alpha = 0.25,
    color = "#135861"
  ) +
  
  geom_line(
    data = spectra_class_clean %>%
      filter(!is.na(biome_cor)),
    aes(band, med),
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



library(dplyr)
library(ggplot2)
library(sf)

years_sel <- c(1990, 2000, 2010, 2020)

ts_df <- LUCAS_combined_sf_NBR %>%
  # if sf, drop geometry for plotting spectra
  st_drop_geometry() %>%
  mutate(
    point_id    = as.character(point_id),
    year        = as.integer(year),
    band        = as.integer(band),
    biome_cor   = ifelse(biome_cor == "Alpine", "Temperate", biome_cor),
    biome_cor   = ifelse(is.na(biome_cor), NA_character_, biome_cor),
    letter_group = as.character(letter_group)
  ) %>%
  filter(year %in% years_sel, !is.na(biome_cor))


ts_df <- ts_df %>%
  mutate(
    letter_group = factor(letter_group, levels = c("A","C","D","E","F")),
    biome_cor    = factor(biome_cor, levels = c("Boreal","Mediterranean","Temperate"))
  )


plot_spectra_year <- function(df_year, year_label, n_per_class = 50) {
  
  # keep only complete 6-band spectra per point within biome×group
  df_complete <- df_year %>%
    group_by(biome_cor, letter_group, point_id) %>%
    filter(n_distinct(band) == 6) %>%
    ungroup()
  
  # sample point_ids per biome×group (robust if small groups)
  set.seed(1)
  sample_ids <- df_complete %>%
    distinct(biome_cor, letter_group, point_id) %>%
    group_by(biome_cor, letter_group) %>%
    group_modify(~ slice_sample(.x, n = min(n_per_class, nrow(.x)))) %>%
    ungroup()
  
  df_sample <- df_complete %>%
    inner_join(sample_ids, by = c("biome_cor","letter_group","point_id"))
  
  # class median per biome×group×band
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
      aes(band, med),
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


p1990 <- plot_spectra_year(filter(ts_df, year == 1990), 1990)
p2000 <- plot_spectra_year(filter(ts_df, year == 2000), 2000)
p2010 <- plot_spectra_year(filter(ts_df, year == 2010), 2010)
p2020 <- plot_spectra_year(filter(ts_df, year == 2020), 2020)

p1990
p2000
p2010
p2020




library(dplyr)
library(ggplot2)
library(sf)

years_sel <- c(1990, 2000, 2010, 2020)

ts_all <- LUCAS_combined_sf_NBR %>%
  st_drop_geometry() %>%
  mutate(
    point_id     = as.character(point_id),
    year         = as.integer(year),
    band         = as.integer(band),
    biome_cor    = ifelse(biome_cor == "Alpine", "Temperate", as.character(biome_cor)),
    letter_group = as.character(letter_group)
  ) %>%
  filter(!is.na(biome_cor)) %>%
  mutate(
    letter_group = factor(letter_group, levels = c("A","C","D","E","F")),
    biome_cor    = factor(biome_cor, levels = c("Boreal","Mediterranean","Temperate"))
  )


F_fixed <- ts_all %>%
  filter(letter_group == "F") %>%
  group_by(biome_cor, letter_group, point_id, band) %>%
  summarise(BAP = median(BAP, na.rm = TRUE), .groups = "drop")



plot_spectra_year <- function(year_label, n_per_class = 50) {
  
  # year-specific spectra for non-F groups only
  df_year <- ts_all %>%
    filter(year == year_label, letter_group != "F") %>%
    select(biome_cor, letter_group, point_id, band, BAP)
  
  # inject fixed F and assign current year (so it shows up in that year's plot)
  df_F <- F_fixed %>%
    mutate(year = year_label) %>%
    select(biome_cor, letter_group, point_id, band, BAP)
  
  df_plot <- bind_rows(df_year, df_F)
  
  # keep only complete 6-band spectra per point within biome×group
  df_complete <- df_plot %>%
    group_by(biome_cor, letter_group, point_id) %>%
    filter(n_distinct(band) == 6) %>%
    ungroup()
  
  # sample point_ids per biome×group
  set.seed(1)
  sample_ids <- df_complete %>%
    distinct(biome_cor, letter_group, point_id) %>%
    group_by(biome_cor, letter_group) %>%
    group_modify(~ slice_sample(.x, n = min(n_per_class, nrow(.x)))) %>%
    ungroup()
  
  df_sample <- df_complete %>%
    semi_join(sample_ids, by = c("biome_cor","letter_group","point_id"))  # preserves rows safely
  
  # class median per biome×group×band
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
      aes(band, med),
      linewidth = 1.0,
      color = "black"
    ) +
    facet_grid(biome_cor ~ letter_group) +
    theme_bw() +
    labs(
      x = "Band",
      y = "BAP",
      title = paste0(" — ", year_label,
                     " ")
    )
}



p1990 <- plot_spectra_year(1990)
p2000 <- plot_spectra_year(2000)
p2010 <- plot_spectra_year(2010)
p2020 <- plot_spectra_year(2020)

p1990; p2000; p2010; p2020


