############################################################
## FULL script (stable lc assignment + no year dropping)
## - Years within a trajectory are NEVER filtered out by ndvi != NA
## - Optional: remove point_ids that do NOT have full year coverage 1984–2024
## - Removes entire trajectories only (by point_id) via QC
## - Uses a STABLE point_id -> lc1 mapping computed ONCE from the raw table
##   (prevents lc migration / strange n per class effects)
## - lc1 recoding:
##     A11/A12/A13 -> A10
##     A21/A22     -> A20
## - F00 used
## - lc1_label added
## - Merge A10/A20 labels into "Artificial surfaces"
## - Median trajectories + random single trajectories (labels in facets)
############################################################

library(data.table)
library(ggplot2)

############################################################
## 0) Inputs + parameters
############################################################
setDT(LUCAS_BAP_6bands_long1)
LUCAS_BAP_6bands_long1[, lc1 := trimws(as.character(lc1))]
LUCAS_BAP_6bands_long1[, year := as.integer(year)]

years_full <- 1984:2024
n_full <- length(years_full)

# target lc after recoding
target_lc <- c("A10","A20","C10","C20","D20","E20","F00")

# --- completeness filter (your new idea) ---
# TRUE  -> keep only point_ids with all years 1984–2024 present in point-year table
# FALSE -> allow missing years in the table (but still never drop years by ndvi NA)
filter_incomplete_years <- TRUE

# QC settings
jump_step_cut <- 0.10
require_full_years <- FALSE  # QC requirement on table-level coverage; if filter_incomplete_years=TRUE, can stay FALSE
min_obs_keep  <- 10          # minimum observed NDVI values
max_na_rate   <- 0.30        # max fraction of NA NDVI values across years_full

robust_k <- 2.0
jump_abs_cut     <- 0.18
p_jump_rate_cut  <- 0.25
outlier_rate_cut <- 0.12

############################################################
## 0.1) Recode/merge lc1 codes ONCE (raw table)
############################################################
LUCAS_BAP_6bands_long1[
  ,
  lc1 := fcase(
    lc1 %in% c("A11","A12","A13"), "A10",
    lc1 %in% c("A21","A22"),      "A20",
    default = lc1
  )
]

############################################################
## 0.2) Build a STABLE point_id -> lc1 mapping computed ONCE
##      (prevents "n per lc increases" due to lc reassignment)
############################################################
lc_mode_global <- LUCAS_BAP_6bands_long1[
  lc1 %in% target_lc,
  .N,
  by = .(point_id, lc1)
][order(point_id, -N)][
  ,
  .SD[1],
  by = point_id
][, .(point_id, lc1)]

# sanity: how many points got an lc?
cat("\nStable lc mapping: n point_ids =", nrow(lc_mode_global), "\n")

############################################################
## 1) Build per-(point_id, year) NDVI table (KEEP NA years)
############################################################
ndvi_py <- unique(
  LUCAS_BAP_6bands_long1[
    lc1 %in% target_lc & year %in% years_full,
    .(point_id, year, ndvi)
  ],
  by = c("point_id", "year")
)

# attach stable lc
ndvi_py <- lc_mode_global[ndvi_py, on = "point_id", nomatch = 0L]
setorder(ndvi_py, point_id, year)

cat("\n--- BEFORE QC: points/point-years/NA% per lc1 (stable mapping) ---\n")
print(ndvi_py[, .(
  n_points    = uniqueN(point_id),
  n_pointyear = .N,
  ndvi_na_pct = 100 * mean(is.na(ndvi))
), by = lc1][order(lc1)])

############################################################
## 1.1) OPTIONAL: remove point_ids with incomplete year axis
############################################################
if (filter_incomplete_years) {
  full_points <- ndvi_py[
    ,
    .(n_years_total = uniqueN(year)),
    by = point_id
  ][n_years_total == n_full, point_id]
  
  cat("\nComplete-year points (1984–2024):", length(full_points),
      "of", uniqueN(ndvi_py$point_id), "\n")
  
  ndvi_py <- ndvi_py[point_id %in% full_points]
}

############################################################
## 2) Compute point-level trajectory metrics (per point_id)
##    - Full year axis is preserved in ndvi_py (if present)
##    - Variability computed on observed NDVI only
############################################################
ts_metrics <- ndvi_py[
  ,
  {
    y_all <- ndvi
    
    n_years_total <- .N
    na_rate <- mean(is.na(y_all))
    n_years_obs <- sum(!is.na(y_all))
    
    y_obs <- y_all[!is.na(y_all)]
    dy <- if (length(y_obs) >= 2) diff(y_obs) else numeric(0)
    
    base  <- if (length(y_obs) > 0) median(y_obs) else NA_real_
    res   <- if (length(y_obs) > 0) y_obs - base else numeric(0)
    s_res <- if (length(y_obs) > 0) mad(res, constant = 1, na.rm = TRUE) else NA_real_
    
    .(
      n_years_total = n_years_total,
      n_years_obs   = n_years_obs,
      na_rate       = na_rate,
      sd_ndvi       = if (length(y_obs) >= 2) sd(y_obs, na.rm = TRUE) else NA_real_,
      max_jump      = if (length(y_obs) >= 2) max(abs(dy), na.rm = TRUE) else NA_real_,
      sd_diff       = if (length(y_obs) >= 3) sd(dy, na.rm = TRUE) else NA_real_,
      mean_abs_diff = if (length(y_obs) >= 2) mean(abs(dy), na.rm = TRUE) else NA_real_,
      p_jump_step   = if (length(y_obs) >= 2) mean(abs(dy) > jump_step_cut, na.rm = TRUE) else NA_real_,
      outlier_rate  = if (length(y_obs) >= 5 && is.finite(s_res) && s_res > 0)
        mean(abs(res) > 3 * s_res, na.rm = TRUE) else 0
    )
  },
  by = point_id
]

############################################################
## 3) Universal QC thresholds (robust)
############################################################
robust_thr <- function(x, k = 2.0) {
  median(x, na.rm = TRUE) + k * mad(x, constant = 1, na.rm = TRUE)
}

sd_thr     <- robust_thr(ts_metrics$sd_ndvi,       k = robust_k)
jump_thr   <- robust_thr(ts_metrics$max_jump,      k = robust_k)
rough_thr  <- robust_thr(ts_metrics$sd_diff,       k = robust_k)
madiff_thr <- robust_thr(ts_metrics$mean_abs_diff, k = robust_k)

############################################################
## 4) Flag bad point_ids (remove whole trajectories only)
############################################################
bad_points <- ts_metrics[
  (require_full_years & n_years_total < n_full) |
    n_years_obs   < min_obs_keep |
    na_rate       > max_na_rate |
    sd_ndvi       > sd_thr |
    max_jump      > jump_thr |
    sd_diff       > rough_thr |
    mean_abs_diff > madiff_thr |
    max_jump      > jump_abs_cut |
    p_jump_step   > p_jump_rate_cut |
    outlier_rate  > outlier_rate_cut,
  point_id
]

cat("\nQC removed fraction:",
    round(length(bad_points) / uniqueN(ts_metrics$point_id), 3), "\n")

############################################################
## 5) Apply filtering to FULL long table (remove point_ids)
############################################################
LUCAS_BAP_6bands_long1_clean <- LUCAS_BAP_6bands_long1[!point_id %in% bad_points]
setDT(LUCAS_BAP_6bands_long1_clean)

############################################################
## 6) Rebuild cleaned point-year NDVI table (KEEP NA years)
############################################################
ndvi_py_clean <- unique(
  LUCAS_BAP_6bands_long1_clean[
    lc1 %in% target_lc & year %in% years_full,
    .(point_id, year, ndvi)
  ],
  by = c("point_id", "year")
)

# attach the SAME stable mapping (recomputed mapping is NOT used)
ndvi_py_clean <- lc_mode_global[ndvi_py_clean, on = "point_id", nomatch = 0L]
setorder(ndvi_py_clean, point_id, year)

# OPTIONAL: enforce complete-year again after QC
if (filter_incomplete_years) {
  full_points_clean <- ndvi_py_clean[
    ,
    .(n_years_total = uniqueN(year)),
    by = point_id
  ][n_years_total == n_full, point_id]
  
  cat("\nAfter QC, complete-year points:", length(full_points_clean),
      "of", uniqueN(ndvi_py_clean$point_id), "\n")
  
  ndvi_py_clean <- ndvi_py_clean[point_id %in% full_points_clean]
  LUCAS_BAP_6bands_long1_clean <- LUCAS_BAP_6bands_long1_clean[point_id %in% full_points_clean]
}

cat("\n--- AFTER QC: points/point-years/NA% per lc1 (stable mapping) ---\n")
print(ndvi_py_clean[, .(
  n_points    = uniqueN(point_id),
  n_pointyear = .N,
  ndvi_na_pct = 100 * mean(is.na(ndvi))
), by = lc1][order(lc1)])

############################################################
## 7) Labels + merge A10/A20 -> Artificial surfaces
############################################################
lc_lookup <- data.table(
  lc1 = c("A10","A20","C10","C20","D20","E20","F00"),
  lc1_label = c(
    "Roofed built-up",
    "Artificial non built-up areas",
    "Broadleaved woodland",
    "Coniferous woodland",
    "Shrubland",
    "Grassland",
    "Rocks, stones and sand"
  )
)

ndvi_py_clean_lab <- lc_lookup[ndvi_py_clean, on = "lc1"]
stopifnot(all(!is.na(ndvi_py_clean_lab$lc1_label)))

ndvi_py_clean_lab[
  lc1_label %in% c("Roofed built-up", "Artificial non built-up areas"),
  lc1_label := "Artificial surfaces"
]

############################################################
## 8) Aggregate to median/IQR per (label, year)
##    (years remain present; NA handled by na.rm=TRUE)
############################################################
ndvi_lab_ts <- ndvi_py_clean_lab[
  ,
  .(
    ndvi_median = median(ndvi, na.rm = TRUE),
    ndvi_q25    = quantile(ndvi, 0.25, na.rm = TRUE),
    ndvi_q75    = quantile(ndvi, 0.75, na.rm = TRUE),
    n_points    = uniqueN(point_id)
  ),
  by = .(lc1_label, year)
]

############################################################
## 9) Order labels for plotting
############################################################
label_order <- c(
  "Artificial surfaces",
  "Broadleaved woodland",
  "Coniferous woodland",
  "Shrubland",
  "Grassland",
  "Rocks, stones and sand"
)
label_order <- intersect(label_order, unique(ndvi_lab_ts$lc1_label))
ndvi_lab_ts[, lc1_label := factor(lc1_label, levels = label_order)]

############################################################
## 9.1) Add n=... to facet titles (TOTAL per label)
############################################################
lc_n <- ndvi_py_clean_lab[
  ,
  .(n_points_total = uniqueN(point_id)),
  by = lc1_label
]
lc_n[, lc1_label := factor(lc1_label, levels = label_order)]
lc_n[, facet_label := paste0(as.character(lc1_label), "\n(n = ", n_points_total, ")")]

ndvi_lab_ts <- lc_n[ndvi_lab_ts, on = "lc1_label"]
stopifnot("facet_label" %in% names(ndvi_lab_ts))

facet_levels <- lc_n[order(lc1_label), facet_label]
ndvi_lab_ts[, facet_label := factor(facet_label, levels = facet_levels)]

############################################################
## 10) Plot: median trajectories
############################################################
p_facet <- ggplot(ndvi_lab_ts, aes(year, ndvi_median)) +
  geom_ribbon(aes(ymin = ndvi_q25, ymax = ndvi_q75), alpha = 0.35) +
  geom_line(linewidth = 1) +
  facet_wrap(~ facet_label, ncol = 3) +
  theme_minimal(base_size = 16) +
  ylim(0, 1) +
  labs(x = "Year", y = "Median NDVI", title = "")

print(p_facet)


############################################################
## 10.1) Save plot
############################################################
out1 <- "/mnt/eo/EU_unmixing/figs/trajectories_all.png"
dir.create(dirname(out1), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = out1, plot = p_facet, width = 12, height = 7, dpi = 300)

############################################################
## 11) Random single trajectory per class label (different each run)
############################################################
setorder(ndvi_py_clean_lab, lc1_label, point_id, year)

# keep factor order
ndvi_py_clean_lab[, lc1_label := factor(lc1_label, levels = label_order)]

sampled_points <- ndvi_py_clean_lab[
  ,
  .(point_id = sample(unique(point_id), 1)),
  by = lc1_label
]

ndvi_single_ts <- ndvi_py_clean_lab[sampled_points, on = .(lc1_label, point_id)]

p_single <- ggplot(ndvi_single_ts, aes(year, ndvi)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~ lc1_label, ncol = 3, scales = "free_y") +
  theme_minimal(base_size = 12) +
  ylim(0, 1) +
  labs(x = "Year", y = "NDVI", title = "")

print(p_single)

out2 <- "/mnt/dss_europe/temp_lm/random_trajectories.png"
ggsave(filename = out2, plot = p_single, width = 12, height = 7, dpi = 300)

############################################################
## 12) Sanity checks: totals cannot increase after filtering
############################################################
cat("\nTotal unique point_ids BEFORE QC (post completeness filter):",
    uniqueN(ndvi_py$point_id), "\n")
cat("Total unique point_ids AFTER  QC:", uniqueN(ndvi_py_clean_lab$point_id), "\n")

cat("\nPoints per lc1 AFTER QC:\n")
print(ndvi_py_clean_lab[, .(n_points = uniqueN(point_id)), by = lc1][order(lc1)])



#### create df
library(data.table)

############################################################
## 0) Setup
############################################################
setDT(LUCAS_BAP_6bands_long1_clean)

# expected year axis
years_full <- 1984:2024
n_full <- length(years_full)

# lc1 recoding (apply again to be safe)
LUCAS_BAP_6bands_long1_clean[, lc1 := trimws(as.character(lc1))]
LUCAS_BAP_6bands_long1_clean[, year := as.integer(year)]

LUCAS_BAP_6bands_long1_clean[
  ,
  lc1 := fcase(
    lc1 %in% c("A11","A12","A13"), "A10",
    lc1 %in% c("A21","A22"),      "A20",
    default = lc1
  )
]

# target codes
target_lc <- c("A10","A20","C10","C20","D20","E20","F00")

############################################################
## 1) Stable point_id -> lc1 mapping (computed ONCE)
##    (use mode across all rows; stable trajectories)
############################################################
lc_mode_global <- LUCAS_BAP_6bands_long1_clean[
  lc1 %in% target_lc,
  .N,
  by = .(point_id, lc1)
][order(point_id, -N)][
  ,
  .SD[1],
  by = point_id
][, .(point_id, lc1)]

############################################################
## 2) Add labels + merge artificial classes
############################################################
lc_lookup <- data.table(
  lc1 = c("A10","A20","C10","C20","D20","E20","F00"),
  lc1_label = c(
    "Roofed built-up",
    "Artificial non built-up areas",
    "Broadleaved woodland",
    "Coniferous woodland",
    "Shrubland",
    "Grassland",
    "Rocks, stones and sand"
  )
)

# join labels to stable lc table
lc_mode_global <- lc_lookup[lc_mode_global, on = "lc1"]
stopifnot(all(!is.na(lc_mode_global$lc1_label)))

# merge the two artificial labels into one
lc_mode_global[
  lc1_label %in% c("Roofed built-up", "Artificial non built-up areas"),
  lc1_label := "Artificial surfaces"
]

############################################################
## 3) Keep only points with COMPLETE year coverage 1984–2024
##    (based on existence of rows per year in the cleaned table)
############################################################
# Here we check year coverage at point-year level (ignoring band multiplicity)
year_cov <- unique(
  LUCAS_BAP_6bands_long1_clean[
    point_id %in% lc_mode_global$point_id & year %in% years_full,
    .(point_id, year)
  ]
)

full_points <- year_cov[
  ,
  .(n_years = .N),
  by = point_id
][n_years == n_full, point_id]

cat("Points with full 1984–2024 coverage:", length(full_points), "\n")

# restrict stable mapping to full points
lc_mode_global <- lc_mode_global[point_id %in% full_points]

############################################################
## 4) Build a FULL grid (point_id x year x band) and merge BAP
##    -> ensures NO missing years within retained trajectories
############################################################
# Determine which bands exist (or set explicitly if you know: 1:6)
bands <- sort(unique(LUCAS_BAP_6bands_long1_clean$band))
# If you *know* it must be 1:6, you can enforce:
# bands <- 1:6

full_grid <- CJ(
  point_id = lc_mode_global$point_id,
  year     = years_full,
  band     = bands,
  unique   = TRUE
)

# bring in BAP from your long table
bap_long <- LUCAS_BAP_6bands_long1_clean[
  point_id %in% lc_mode_global$point_id & year %in% years_full & band %in% bands,
  .(point_id, year, band, BAP)
]

# If duplicates exist per point-year-band, collapse (median)
bap_long <- bap_long[
  ,
  .(BAP = median(BAP, na.rm = TRUE)),
  by = .(point_id, year, band)
]

df_full <- bap_long[full_grid, on = .(point_id, year, band)]

############################################################
## 5) Add stable lc1 + lc1_label to every row
############################################################
df_full <- lc_mode_global[df_full, on = "point_id"]

############################################################
## 6) Add NDVI as point-year attribute (replicated across bands)
##    Prefer: compute from df_full (band 3 red, band 4 nir) to guarantee consistency
############################################################
# pull band 3 and 4 per point-year from df_full
ndvi_wide <- dcast(
  df_full[band %in% c(3,4)],
  point_id + year ~ band,
  value.var = "BAP"
)
setnames(ndvi_wide, old = c("3","4"), new = c("red","nir"), skip_absent = TRUE)

ndvi_wide[
  ,
  ndvi := fifelse(
    is.na(red) | is.na(nir) | (red + nir) == 0,
    NA_real_,
    (nir - red) / (nir + red)
  )
]

# join ndvi back (replicate across bands)
df_full <- ndvi_wide[df_full, on = .(point_id, year)]

############################################################
## 7) Final column order + sanity checks
############################################################
setcolorder(df_full, c("point_id","year","band","BAP","lc1","lc1_label","ndvi"))

# check: every retained point has all years and all bands
stopifnot(df_full[, all(uniqueN(year) == n_full), by = point_id]$V1)
stopifnot(df_full[, all(uniqueN(band) == length(bands)), by = point_id]$V1)

# optional: check stable lc1 per point_id
stopifnot(df_full[, uniqueN(lc1), by = point_id][, all(V1 == 1)])

df_full

fwrite(df_full, "/mnt/eo/EU_unmixing/data/LUCAS/LUCAS_full_timeseries.csv")







