library(dplyr)
library(tidyr)
library(readr)

df_full <- read_csv("/mnt/eo/EU_unmixing/data/LUCAS/LUCAS_full_timeseries.csv")


# ------------------------------------------------------------
# Endmember candidate selection using PPI (Option B: baseline per LC)
# Goal:
# 1) Compute PPI per observation within each (lc1_label, year)
# 2) Derive an LC-wide purity baseline (PPI_q) by pooling years within each lc1_label
# 3) Select the purest observations per lc1_label while keeping enough samples per LC
# ------------------------------------------------------------

# ----------------------------
# 0) Reshape long -> wide
# df must contain: point_id, year, lc1, lc1_label, ndvi, band, BAP
# band 1..6 correspond to Blue, Green, Red, NIR, SWIR1, SWIR2 (in that order)
# ----------------------------
df_wide <- df_full %>%
  select(point_id, year, lc1, lc1_label, ndvi, band, BAP) %>%
  pivot_wider(
    names_from  = band,
    values_from = BAP,
    names_prefix = "b"
  )

spec_cols <- paste0("b", 1:6)

# Optional sanity check: ensure all band columns exist
stopifnot(all(spec_cols %in% names(df_wide)))

# ----------------------------
# 1) PPI function (classic random projection extremes)
# ----------------------------
compute_ppi <- function(X, n_proj = 2000, n_extreme = 1) {
  X <- as.matrix(X)
  n_pixels <- nrow(X)
  n_bands  <- ncol(X)
  
  ppi <- numeric(n_pixels)
  
  for (i in seq_len(n_proj)) {
    v <- rnorm(n_bands)
    v <- v / sqrt(sum(v^2))     # unit vector
    proj <- drop(X %*% v)
    
    hi <- order(proj, decreasing = TRUE)[seq_len(n_extreme)]
    lo <- order(proj, decreasing = FALSE)[seq_len(n_extreme)]
    
    ppi[hi] <- ppi[hi] + 1
    ppi[lo] <- ppi[lo] + 1
  }
  
  # normalize to [0,1] (share of times selected as extreme)
  ppi / (2 * n_proj * n_extreme)
}

# ----------------------------
# 2) Compute PPI per observation within each (lc1_label, year)
#    (one PPI value per row in df_wide)
# ----------------------------
set.seed(42)

min_group_n <- 30     # minimum obs per (lc, year) group to compute PPI
n_proj      <- 2000   # adjust for runtime vs stability
n_extreme   <- 1

df_ppi <- df_wide %>%
  group_by(lc1_label, year) %>%
  group_modify(\(d, keys) {
    
    ok <- complete.cases(d[, spec_cols])
    d$PPI <- NA_real_
    
    if (sum(ok) >= min_group_n) {
      X <- scale(as.matrix(d[ok, spec_cols]))  # scale within group
      d$PPI[ok] <- compute_ppi(X, n_proj = n_proj, n_extreme = n_extreme)
    }
    
    d
  }) %>%
  ungroup()

# ----------------------------
# 3) Option B baseline: compute LC-wide purity rank by pooling years
#    PPI_q is the percentile rank of PPI within each lc1_label (across all years)
# ----------------------------
df_ppi <- df_ppi %>%
  group_by(lc1_label) %>%
  mutate(PPI_q = percent_rank(PPI)) %>%
  ungroup()

# ----------------------------
# 4) Select purest endmember candidates per LC while keeping enough samples per LC
#    Choose ONE of the following selection approaches:
#    A) fixed n per LC
#    B) top percentile per LC
#    C) hybrid: purity floor + cap
# ----------------------------


# ---- B) Top percentile per LC (scale-free; class sizes will differ)
q_cut <- 0.995  # top 0.5% per LC (try 0.99 or 0.995)

pure_obs <- df_ppi %>%
  filter(!is.na(PPI_q), !is.na(PPI)) %>%
  group_by(lc1_label) %>%
  filter(PPI_q >= q_cut) %>%
  arrange(desc(PPI)) %>%      # within LC, order by raw PPI
  ungroup()

pure_points <- pure_obs %>%
  group_by(point_id, lc1_label) %>%
  summarise(
    n_pure_years = n(),
    first_year = min(year),
    last_year  = max(year),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pure_years))

pure_obs %>% count(lc1_label) %>% arrange(n)


### plot in feature space
spec_cols <- paste0("b", 1:6)

# Fit PCA on all observations (or per LC if you prefer)
X <- scale(df_ppi[, spec_cols])
pca <- prcomp(X, center = FALSE, scale. = FALSE)

df_ppi$pca1 <- pca$x[, 1]
df_ppi$pca2 <- pca$x[, 2]

# mark purest observations
df_ppi <- df_ppi %>%
  mutate(
    pure_flag = PPI_q >= 0.995   # or whatever threshold you used
  )

ggplot(df_ppi, aes(pca1, pca2)) +
  geom_point(color = "grey80", size = 0.6) +
  geom_point(
    data = df_ppi %>% filter(pure_flag),
    aes(color = lc1_label),
    size = 2
  ) +
  labs(
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()


### per year
# choose your purity definition
q_cut <- 0.995
df_ppi2 <- df_ppi %>%
  mutate(pure_flag = !is.na(PPI_q) & PPI_q >= q_cut)

ggplot(df_ppi2, aes(pca1, pca2)) +
  geom_point(color = "grey85", size = 0.25) +
  geom_point(
    data = df_ppi2 %>% filter(pure_flag),
    aes(color = lc1_label),
    size = 0.8,
    alpha = 0.9
  ) +
  facet_wrap(~ year, ncol = 7) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "right"
  ) +
  labs(x = "PC1 (spectral)", y = "PC2 (spectral)", color = "LC (pure)")




### compute time-series purity score per lc
pixel_score <- df_ppi %>%
  filter(!is.na(PPI)) %>%
  group_by(point_id, lc1_label) %>%
  summarise(
    n_years   = n(),
    ppi_mean  = mean(PPI, na.rm = TRUE),
    ppi_sum   = sum(PPI, na.rm = TRUE),
    ppi_median = median(PPI, na.rm = TRUE),
    .groups = "drop"
  )


### select top % pixel per lc class
top_prop <- 0.15
min_n <- 40

endmember_pixels <- pixel_score %>%
  group_by(lc1_label) %>%
  arrange(desc(ppi_median), .by_group = TRUE) %>%
  group_modify(\(d, key) {
    n_keep <- max(min_n, ceiling(top_prop * nrow(d)))
    d[seq_len(min(n_keep, nrow(d))), , drop = FALSE]
  }) %>%
  ungroup()

endmember_pixels %>%
  count(lc1_label) %>%
  arrange(n)


endmember_pixels <- endmember_pixels %>%
  left_join(
    LUCAS_multiyear_filtered %>%
      select(point_id, gps_lat, gps_long),
    by = "point_id"
  )

write.csv(endmember_pixels,
          "/mnt/eo/EU_unmixing/data/endmembers_filtered.csv",
          row.names = FALSE)


library(sf)

endmember_pixels_unique <- endmember_pixels %>%
  distinct(point_id, .keep_all = TRUE)



sf_df <- st_as_sf(
  endmember_pixels_unique,
  coords = c("gps_long", "gps_lat"),
  crs = 4326,          # WGS84
  remove = FALSE      # keep original lon/lat columns
)

st_write(
  sf_df,
  "/mnt/eo/EU_unmixing/data/LUCAS/endmembers.gpkg",
  layer = "endmember_pixels",
  delete_layer = TRUE
)






             
