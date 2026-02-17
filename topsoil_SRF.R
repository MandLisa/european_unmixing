library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# ============================================================
# A) INPUTS
# ============================================================

# 1) Your raw wide LUCAS spectra data frame:
#    It must have metadata columns + thousands of wavelength columns named like "400", "400.5", ..., "2499.5"
#    Here I assume it is called `raw_wide`. If yours is called differently, rename accordingly.
# raw_wide <- read_csv("path/to/lucas.csv", show_col_types = FALSE)

# 2) Landsat SRF directory (ch_res_1..9.txt or NWP-SAF rtcoef... files)
srf_dir <- "/mnt/dss_project/lmandl/_unmixing/SRF"

# 3) Metadata columns you want to retain (adjust to your table)
meta_cols <- c("SampleID", "PointID", "CLIMA_Cod", "Elevation", "Slope", "Aspect",
               "BioGeo", "Soil_Group", "Soil_Code")

# ============================================================
# B) 1) WIDE -> LONG SPECTRA (ROBUST, KEEPS SWIR)
# ============================================================

# Identify wavelength columns purely by column NAME being numeric (allows decimals)
get_wl_cols <- function(df) {
  wl_cols <- names(df)[str_detect(names(df), "^\\d+(?:\\.\\d+)?$")]
  if (length(wl_cols) == 0) stop("No wavelength columns detected by name. Check column names.")
  wl_cols
}

# Safely coerce wavelength columns to numeric without dropping columns
coerce_wl_cols_numeric <- function(df, wl_cols) {
  df %>%
    mutate(
      across(all_of(wl_cols), ~ suppressWarnings(as.numeric(.)))
    )
}






library(dplyr)
library(ggplot2)

id <- "18211_1"  # replace with any sample from your plot

raw_one <- ds %>%
  filter(SampleID == id) %>%
  filter(!is.na(reflectance))

range(raw_one$reflectance, na.rm = TRUE)

ggplot(raw_one, aes(wl, reflectance)) +
  geom_line() +
  labs(title = paste("Raw hyperspectral curve:", id),
       x = "Wavelength (nm)", y = "reflectance (as stored)")


# Report columns that still contain many NAs after coercion (often parsing problems)
report_bad_wl_cols <- function(df, wl_cols, na_frac_threshold = 0.5) {
  na_frac <- sapply(df[wl_cols], function(v) mean(is.na(v)))
  bad <- names(na_frac)[na_frac > na_frac_threshold]
  list(na_fraction = na_frac, bad_cols = bad)
}

# ---- MAIN: create ds_long ----
# Replace `raw_wide` with your wide data frame object name
wl_cols <- get_wl_cols(spectra)

# Confirm the wavelength range encoded in the column names
wl_numeric <- as.numeric(wl_cols)
cat("Wavelength columns detected:", length(wl_cols), "\n")
cat("Wavelength range (from column names):", range(wl_numeric, na.rm = TRUE), "\n")

# Coerce all wavelength columns to numeric
raw_wide_num <- coerce_wl_cols_numeric(spectra, wl_cols)

# Optional: detect problematic wavelength columns (helps debugging)
report_bad_wl_cols <- function(df, wl_cols, na_frac_threshold = 0.5) {
  na_frac <- sapply(wl_cols, function(col) mean(is.na(df[[col]])))
  bad <- names(na_frac)[na_frac > na_frac_threshold]
  list(na_fraction = na_frac, bad_cols = bad)
}


# Keep only metadata columns that exist
meta_present <- meta_cols[meta_cols %in% names(raw_wide_num)]
if (!("SampleID" %in% meta_present)) stop("SampleID must be present in raw_wide_num.")

# Pivot long (this keeps full VNIR+SWIR because wl_cols were detected by name)
ds <- raw_wide_num %>%
  select(all_of(meta_present), all_of(wl_cols)) %>%
  pivot_longer(cols = all_of(wl_cols),
               names_to = "wl",
               values_to = "reflectance") %>%
  mutate(
    wl = as.numeric(wl),
    reflectance = as.numeric(reflectance)
  )

# Sanity check: you should now see SWIR coverage
cat("After pivot_longer(): wl range =", paste(range(ds$wl, na.rm = TRUE), collapse = " - "), "\n")

# ============================================================
# C) 2) PREPARE SPECTRA FOR SRF JOIN
# ============================================================

round2 <- function(x, n) {
  posneg <- sign(x)
  z <- abs(x) * 10^n
  z <- z + 0.5 + sqrt(.Machine$double.eps)
  z <- trunc(z)
  z <- z / 10^n
  z * posneg
}

ds_spec <- ds %>%
  mutate(wl_r = round2(wl, 0)) %>%
  group_by(across(all_of(c("SampleID", setdiff(meta_present, "SampleID"), "wl_r")))) %>%
  summarise(mean_refl = mean(reflectance, na.rm = TRUE), .groups = "drop")

cat("Prepared spectra wl_r range =", paste(range(ds_spec$wl_r, na.rm = TRUE), collapse = " - "), "\n")

# ============================================================
# D) 3) READ SRFs (supports your ch_res_*.txt format)
# ============================================================

read_chres_srf <- function(path) {
  x <- readLines(path, warn = FALSE)
  
  two_num_pat <- "^\\s*[-+]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?\\s+[-+]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?\\s*$"
  first_data <- which(str_detect(x, two_num_pat))[1]
  if (is.na(first_data)) stop("No numeric SRF block found in: ", path)
  
  dat <- read_table(
    paste(x[first_data:length(x)], collapse = "\n"),
    col_names = c("wn_cm1", "response"),
    col_types = cols(wn_cm1 = col_double(), response = col_double())
  )
  
  dat %>%
    mutate(wl_nm = 1e7 / wn_cm1) %>%
    select(wl_nm, response)
}

wmean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & (w > 0)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# Locate SRF files: handle BOTH naming schemes
files_chres <- list.files(srf_dir, pattern = "^ch_res_\\d+\\.txt$", full.names = TRUE)
files_rtcoef <- list.files(srf_dir, pattern = "^rtcoef_landsat_9_oli_srf_ch\\d{2}\\.txt$", full.names = TRUE)

files <- if (length(files_chres) > 0) files_chres else files_rtcoef
if (length(files) != 9) stop("Expected 9 SRF files in ", srf_dir, " but found ", length(files))

# Channel -> band mapping
# For ch_res_1..9 : "1".."9"
# For rtcoef ... ch01..ch09 : "01".."09"
chan_to_band_1 <- c("1"="B1","2"="B2","3"="B3","4"="B4","5"="B5","6"="B6","7"="B7","8"="B8","9"="B9")
chan_to_band_2 <- c("01"="B1","02"="B2","03"="B3","04"="B4","05"="B5","06"="B6","07"="B7","08"="B8","09"="B9")

srf_long <- bind_rows(lapply(files, function(f) {
  
  bn <- basename(f)
  
  if (str_detect(bn, "^ch_res_\\d+\\.txt$")) {
    ch <- str_match(bn, "^ch_res_(\\d+)\\.txt$")[,2]
    band <- chan_to_band_1[[ch]]
  } else {
    ch <- str_match(bn, "ch(\\d{2})\\.txt$")[,2]
    band <- chan_to_band_2[[ch]]
  }
  
  if (is.null(band)) stop("Could not map channel for file: ", f)
  
  read_chres_srf(f) %>%
    mutate(band = band,
           wl_r = round(wl_nm)) %>%
    group_by(band, wl_r) %>%
    summarise(response = mean(response, na.rm = TRUE), .groups = "drop")
}))

srf_wide <- srf_long %>%
  pivot_wider(names_from = band, values_from = response, values_fill = 0) %>%
  rename(SR_WL = wl_r)

cat("SRF wl range =", paste(range(srf_wide$SR_WL, na.rm = TRUE), collapse = " - "), "\n")

# ============================================================
# E) 4) JOIN + CONVOLVE TO LANDSAT BANDS (B1..B9)
# ============================================================

# Restrict SRF to spectral range present in spectra (faster, cleaner)
wl_min <- floor(min(ds_spec$wl_r, na.rm = TRUE))
wl_max <- ceiling(max(ds_spec$wl_r, na.rm = TRUE))
srf_wide_sub <- srf_wide %>% filter(SR_WL >= wl_min, SR_WL <= wl_max)

ds_join <- ds_spec %>% left_join(srf_wide_sub, by = c("wl_r" = "SR_WL"))

band_cols <- intersect(c("B1","B2","B3","B4","B5","B6","B7","B8","B9"), names(ds_join))
bands_overlap <- band_cols[colSums(ds_join[, band_cols] > 0, na.rm = TRUE) > 0]

cat("Bands with overlap (computable):", paste(bands_overlap, collapse = ", "), "\n")

ds_l9 <- ds_join %>%
  group_by(across(all_of(c("SampleID", setdiff(meta_present, "SampleID"))))) %>%
  summarise(
    across(all_of(bands_overlap), ~ wmean_safe(mean_refl, .x), .names = "L_{.col}"),
    .groups = "drop"
  )

ds_l9_long <- ds_l9 %>%
  pivot_longer(cols = starts_with("L_"), names_to = "band", values_to = "value")

# ============================================================
# F) OPTIONAL: central wavelengths (only for bands computed)
# ============================================================
cwl_l89_full <- c(L_B1=443, L_B2=482, L_B3=561, L_B4=655, L_B5=865, L_B6=1610, L_B7=2200, L_B8=590, L_B9=1375)

ds_l9_long <- ds_l9_long %>%
  mutate(cwl = unname(cwl_l89_full[band]))

# Done:
# - ds      : long spectra (full VNIR+SWIR if present in raw_wide)
# - srf_wide: SRF table
# - ds_l9   : simulated Landsat bands per sample
# - ds_l9_long: long version for plotting

out_dir <- "/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/2015"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv_wide <- file.path(out_dir, "LUCAS_topsoil_reflectance_wide.csv")
out_csv_long <- file.path(out_dir, "LUCAS_topsoil_reflectance_long.csv")

write.csv(ds_l9, out_csv_wide, row.names = FALSE)
write.csv(ds_l9_long, out_csv_long, row.names = FALSE)



set.seed(1)
n_samples <- 24
ids <- sample(unique(ds_l9_long$SampleID), n_samples)

p1 <- ds_l9_long %>%
  filter(SampleID %in% ids) %>%
  filter(!is.na(value)) %>%
  # keep reflective bands only (optional; drop pan/cirrus if you want)
  filter(band %in% c("L_B1","L_B2","L_B3","L_B4","L_B5","L_B6","L_B7")) %>%
  ggplot(aes(x = cwl, y = value, group = SampleID)) +
  geom_line() +
  geom_point(size = 1.6) +
  facet_wrap(~SampleID) +
  labs(x = "Wavelength (nm)", y = "Simulated reflectance",
       title = "Landsat-simulated soil spectral signatures (sample subset)")

print(p1)
