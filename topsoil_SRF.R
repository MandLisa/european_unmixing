# ============================================================
# LUCAS 2015 soil spectra -> Landsat 8/9 (OLI/OLI-2) simulation
# Full end-to-end script (data.table-safe) + plots + write files
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(data.table)

# ------------------------------------------------------------
# 0) INPUTS (ADAPT THESE)
# ------------------------------------------------------------

# Your wide LUCAS table already loaded in memory:
# It must contain meta columns + wavelength columns named like "400", "400.5", ..., "2499.5"
# You said this object is called:
file_path <- "/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/2015/spectra_EU28_merged.csv"

spectra <- fread(file_path)

raw_wide <- spectra

# folder containing SRF files (either ch_res_1..9.txt or rtcoef...ch01..ch09.txt)
srf_dir <- "/mnt/dss_project/lmandl/_unmixing/SRF"

# output folder
out_dir <- "/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/2015"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# metadata columns to keep (only those present will be kept)
meta_cols <- c("source", "SampleID", "PointID", "NUTS_0", "SampleN",
               "CLIMA_Cod", "Elevation", "Slope", "Aspect",
               "BioGeo", "Soil_Group", "Soil_Code")

# IMPORTANT: LUCAS spectra value domain
# Colleague script uses: reflectance = 1/(10^value)  (== 10^(-value))
# Set this TRUE if your wavelength columns are in that "value" domain (common for LUCAS exports).
apply_lucas_transform <- TRUE

# If you need additional filtering (e.g. wl >= 426 as in colleague script)
min_wl <- 426

# Reflective bands to use for "spectral signature" plots (exclude pan/cirrus by default)
bands_plot <- c("L_B1","L_B2","L_B3","L_B4","L_B5","L_B6","L_B7")

# ------------------------------------------------------------
# 1) WIDE -> LONG (KEEP FULL VNIR+SWIR) + LUCAS TRANSFORM
# ------------------------------------------------------------

# Detect wavelength columns by name (numeric, allow decimals)
wl_cols <- names(raw_wide)[str_detect(names(raw_wide), "^\\d+(?:\\.\\d+)?$")]
if (length(wl_cols) == 0) stop("No wavelength columns detected. Check your column names.")

wl_numeric <- as.numeric(wl_cols)
cat("Detected", length(wl_cols), "wavelength columns.\n")
cat("Wavelength range from column names:", paste(range(wl_numeric, na.rm=TRUE), collapse=" - "), "\n")

# Keep only meta columns that exist
meta_present <- meta_cols[meta_cols %in% names(raw_wide)]
if (!("SampleID" %in% meta_present)) stop("SampleID must exist in your table.")

# Coerce wavelength columns to numeric robustly (data.table-safe)
raw_wide_num <- raw_wide %>%
  mutate(across(all_of(wl_cols), ~ suppressWarnings(as.numeric(.))))

# Pivot longer into "value" first
ds <- raw_wide_num %>%
  select(all_of(meta_present), all_of(wl_cols)) %>%
  pivot_longer(cols = all_of(wl_cols),
               names_to = "wl",
               values_to = "value") %>%
  mutate(
    wl = as.numeric(wl),
    value = as.numeric(value)
  ) %>%
  filter(!is.na(value)) %>%
  filter(wl >= min_wl)

# Apply LUCAS transform if needed (as in your colleague script)
if (apply_lucas_transform) {
  ds <- ds %>% mutate(reflectance = 1/(10^value))  # == 10^(-value)
} else {
  ds <- ds %>% mutate(reflectance = value)
}

cat("After pivot_longer(): wl range =", paste(range(ds$wl, na.rm=TRUE), collapse=" - "), "\n")
cat("Reflectance summary (after transform setting):\n")
print(ds %>% summarise(min=min(reflectance,na.rm=TRUE),
                       p01=quantile(reflectance,0.01,na.rm=TRUE),
                       med=median(reflectance,na.rm=TRUE),
                       p99=quantile(reflectance,0.99,na.rm=TRUE),
                       max=max(reflectance,na.rm=TRUE)))

# ------------------------------------------------------------
# 2) PREP SPECTRA FOR SRF JOIN (mean reflectance per integer nm)
# ------------------------------------------------------------

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

cat("Prepared ds_spec wl_r range =", paste(range(ds_spec$wl_r, na.rm=TRUE), collapse=" - "), "\n")

# ------------------------------------------------------------
# 3) READ SRFs (ch_res_*.txt or rtcoef_...txt) -> srf_wide
# ------------------------------------------------------------

read_chres_srf <- function(path) {
  x <- readLines(path, warn = FALSE)
  
  # first line that contains exactly two numbers (wn + response)
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

# locate SRF files
files_chres <- list.files(srf_dir, pattern = "^ch_res_\\d+\\.txt$", full.names = TRUE)
files_rtcoef <- list.files(srf_dir, pattern = "^rtcoef_landsat_9_oli_srf_ch\\d{2}\\.txt$", full.names = TRUE)
files <- if (length(files_chres) > 0) files_chres else files_rtcoef
if (length(files) != 9) stop("Expected 9 SRF files in ", srf_dir, " but found ", length(files))

# channel -> band mapping
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

cat("SRF wl range =", paste(range(srf_wide$SR_WL, na.rm=TRUE), collapse=" - "), "\n")

# ------------------------------------------------------------
# 4) JOIN + CONVOLVE TO LANDSAT BANDS
# ------------------------------------------------------------

wmean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & is.finite(x) & (w > 0)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# Restrict SRF to available wavelengths in ds_spec
wl_min <- floor(min(ds_spec$wl_r, na.rm = TRUE))
wl_max <- ceiling(max(ds_spec$wl_r, na.rm = TRUE))
srf_wide_sub <- srf_wide %>% filter(SR_WL >= wl_min, SR_WL <= wl_max)

ds_join <- ds_spec %>% left_join(srf_wide_sub, by = c("wl_r" = "SR_WL"))

band_cols <- intersect(c("B1","B2","B3","B4","B5","B6","B7","B8","B9"), names(ds_join))
bands_overlap <- band_cols[colSums(ds_join[, band_cols] > 0, na.rm = TRUE) > 0]
cat("Bands with overlap (computable):", paste(bands_overlap, collapse=", "), "\n")

ds_l9 <- ds_join %>%
  group_by(across(all_of(c("SampleID", setdiff(meta_present, "SampleID"))))) %>%
  summarise(
    across(all_of(bands_overlap), ~ wmean_safe(mean_refl, .x), .names = "L_{.col}"),
    .groups = "drop"
  )

ds_l9_long <- ds_l9 %>%
  pivot_longer(cols = starts_with("L_"), names_to = "band", values_to = "value")

# Attach approximate CWLs (nm) (used only for plotting x-axis)
cwl_l89_full <- c(L_B1=443, L_B2=482, L_B3=561, L_B4=655, L_B5=865, L_B6=1610, L_B7=2200, L_B8=590, L_B9=1375)
ds_l9_long <- ds_l9_long %>%
  mutate(cwl = unname(cwl_l89_full[band]))

# ------------------------------------------------------------
# 5) WRITE OUTPUTS
# ------------------------------------------------------------

out_csv_wide <- file.path(out_dir, "LUCAS_topsoil_reflectance_wide.csv")
out_csv_long <- file.path(out_dir, "LUCAS_topsoil_reflectance_long.csv")

write.csv(ds_l9, out_csv_wide, row.names = FALSE)
write.csv(ds_l9_long, out_csv_long, row.names = FALSE)

cat("Wrote:\n", out_csv_wide, "\n", out_csv_long, "\n")

# ------------------------------------------------------------
# 6) PLOTS: spectral signatures from Landsat-simulated bands (ds_l9_long)
# ------------------------------------------------------------

# (A) Subset of SampleIDs, faceted
set.seed(1)
n_samples <- 24
ids <- sample(unique(ds_l9_long$SampleID), min(n_samples, length(unique(ds_l9_long$SampleID))))

p_sig_subset <- ds_l9_long %>%
  filter(SampleID %in% ids) %>%
  filter(!is.na(value)) %>%
  filter(band %in% bands_plot) %>%
  ggplot(aes(x = cwl, y = value, group = SampleID)) +
  geom_line() +
  geom_point(size = 1.6) +
  facet_wrap(~SampleID) +
  theme_minimal() +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = "")

plot(p_sig_subset)

ggsave(file.path(out_dir, "soil_signatures_Landsat_subset.png"),
       plot = p_sig_subset, width = 14, height = 8, dpi = 300)

# (B) Median ± IQR by Soil_Group (if present)
if ("Soil_Group" %in% names(ds_l9_long)) {
  sig_group <- ds_l9_long %>%
    filter(!is.na(value)) %>%
    filter(band %in% bands_plot) %>%
    group_by(Soil_Group, band, cwl) %>%
    summarise(
      med = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # limit to top groups for readability
  top_k <- 9
  top_groups <- ds_l9 %>%
    count(Soil_Group, sort = TRUE) %>%
    slice_head(n = top_k) %>%
    pull(Soil_Group)
  
  sig_group <- sig_group %>% filter(Soil_Group %in% top_groups)
  
  p_sig_group <- ggplot(sig_group, aes(x = cwl, y = med, group = Soil_Group)) +
    geom_line() +
    geom_point(size = 1.6) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0) +
    facet_wrap(~Soil_Group) +
    theme_minimal() +
    labs(x = "Wavelength (nm)", y = "Reflectance (median ± IQR)",
         title = "Landsat-simulated soil signatures by Soil_Group (top groups)")
  
  ggsave(file.path(out_dir, "soil_signatures_Landsat_bySoilGroup.png"),
         plot = p_sig_group, width = 14, height = 8, dpi = 300)
}

# (C) Overlay for one sample: raw hyperspectral vs Landsat points (QC)
# pick one sample
one_id <- ids[1]

raw_one <- ds %>%
  filter(SampleID == one_id) %>%
  filter(!is.na(reflectance))

ls_one <- ds_l9_long %>%
  filter(SampleID == one_id) %>%
  filter(band %in% bands_plot)

p_overlay <- ggplot() +
  geom_line(data = raw_one, aes(x = wl, y = reflectance), linewidth = 0.7) +
  geom_line(data = ls_one,  aes(x = cwl, y = value), linewidth = 0.6) +
  geom_point(data = ls_one, aes(x = cwl, y = value), size = 2.5) +
  theme_minimal() +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = paste("QC overlay: hyperspectral vs Landsat-simulated:", one_id))

ggsave(file.path(out_dir, "QC_overlay_hyperspectral_vs_Landsat.png"),
       plot = p_overlay, width = 10, height = 6, dpi = 300)

cat("Plots saved in: ", out_dir, "\n")


### Re-import long spectra file
library(readr)

# file with Landsat-adjusted reflectance
topsoil_reflectance_long <- read_csv("/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/2015/LUCAS_topsoil_reflectance_long.csv")

# file with land cover types
topsoil_LC <- read_csv("/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/2015/LUCAS2015_topsoildata_20200323/LUCAS_Topsoil_2015_20200323.csv")

# ad lc info to topsoil_reflectance_long
topsoil_reflectance_long <- topsoil_reflectance_long %>%
  left_join(
    topsoil_LC %>% select(Point_ID, LC1),
    by = c("PointID" = "Point_ID")
  )

# filter only for woodland and bare land samples
topsoil_reflectance_long <- topsoil_reflectance_long %>%
  filter(grepl("^C", LC1) | LC1 == "F10")

# convert to long
topsoil_reflectance_wide <- topsoil_reflectance_long %>%
  pivot_wider(
    id_cols     = c(SampleID, source, PointID, NUTS_0, SampleN, LC1),
    names_from  = band,
    values_from = value
  )

# average sample pairs
topsoil_mean_long <- topsoil_reflectance_long %>%
  mutate(SampleID_base = str_remove(SampleID, "_[0-9]+$")) %>%  # 11001_1 -> 11001
  group_by(SampleID_base, band) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    # keep metadata (assuming constant within base ID)
    source  = first(source),
    PointID = first(PointID),
    NUTS_0  = first(NUTS_0),
    SampleN = first(SampleN),
    cwl     = first(cwl),
    LC1     = first(LC1),
    .groups = "drop"
  )

