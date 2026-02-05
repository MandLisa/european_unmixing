############################################################
## LUCAS NDVI trajectory QC – THREE MEANINGFUL SCENARIOS
##
## Goals
## 1) Trajectory-level QC only: remove entire point_id (never single years)
## 2) Stable lc1 assignment: compute ONCE (mode per point_id)
## 3) Three *conceptually* different QC scenarios (not just k-MAD tweaks)
## 4) For each scenario:
##    - write full table (all original columns)
##    - add point-level linear NDVI trend (slope per year + per decade)
##    - add latitude-based ecoregion group
## 5) Visualise (balanced): NDVI trend distributions by ecoregion and LC facets
##
## Outputs (in base_dir):
## - LUCAS_BAP_comp_QC_<scenario>.csv
## - LUCAS_BAP_comp_QC_<scenario>_with_trend_ecogroup.csv
## - fig_balanced_trends_by_ecoregion_and_lc.png
############################################################

library(data.table)
library(ggplot2)

# ============================================================
# 0) Setup
# ============================================================
# --- Option A: load from .csv (fast & simple) ---
LUCAS_BAP_comp <- fread("/mnt/eo/EU_unmixing/data/LUCAS/LUCAS_BAP_full.csv")

setDT(LUCAS_BAP_comp)

base_dir <- "/mnt/eo/EU_unmixing/data/LUCAS"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

years_full <- 1984:2024
n_full     <- length(years_full)

target_lc <- c("A10","A20","C10","C20","C30","D20","E20","F00")

# Enforce complete year axis BEFORE QC? (recommended)
filter_incomplete_years <- TRUE

# Shared QC constraints across scenarios
min_obs_keep  <- 10
jump_step_cut <- 0.10

# ============================================================
# 1) Hygiene + lc recode ONCE
# ============================================================
LUCAS_BAP_comp[, lc1  := trimws(as.character(lc1))]
LUCAS_BAP_comp[, year := as.integer(year)]

# collapse subclasses ONCE
LUCAS_BAP_comp[
  ,
  lc1 := fcase(
    lc1 %in% c("A11","A12","A13"), "A10",
    lc1 %in% c("A21","A22"),      "A20",
    default = lc1
  )
]

# ============================================================
# 2) Stable point_id -> lc1 mapping (MODE)
# ============================================================
lc_mode <- LUCAS_BAP_comp[
  lc1 %in% target_lc,
  .N,
  by = .(point_id, lc1)
][order(point_id, -N)][
  , .SD[1], by = point_id
][, .(point_id, lc1)]

cat("Stable lc mapping for", nrow(lc_mode), "points\n")

# ============================================================
# 3) Build point–year NDVI table (unique per point_id-year)
#    Keep NA years as NA values (trajectory QC only)
# ============================================================
ndvi_py <- unique(
  LUCAS_BAP_comp[
    point_id %in% lc_mode$point_id & year %in% years_full,
    .(point_id, year, ndvi)
  ],
  by = c("point_id","year")
)

# attach stable lc1
ndvi_py <- lc_mode[ndvi_py, on = "point_id", nomatch = 0L]
setorder(ndvi_py, point_id, year)

# Optional: enforce full year axis (1984–2024 present for each point_id)
if (filter_incomplete_years) {
  full_points <- ndvi_py[, .(n_years = uniqueN(year)), by = point_id][n_years == n_full, point_id]
  cat("Points with full 1984–2024 coverage:",
      length(full_points), "of", uniqueN(ndvi_py$point_id), "\n")
  ndvi_py <- ndvi_py[point_id %in% full_points]
}

# ============================================================
# 4) Compute trajectory metrics ONCE (reused across scenarios)
# ============================================================
ts_metrics <- ndvi_py[
  ,
  {
    y_all <- ndvi
    y_obs <- y_all[!is.na(y_all)]
    dy    <- if (length(y_obs) >= 2) diff(y_obs) else numeric(0)
    
    base  <- if (length(y_obs) > 0) median(y_obs) else NA_real_
    res   <- y_obs - base
    s_res <- if (length(res) > 0) mad(res, constant = 1) else NA_real_
    
    .(
      n_years_obs   = length(y_obs),
      na_rate       = mean(is.na(y_all)),
      sd_ndvi       = if (length(y_obs) >= 2) sd(y_obs) else NA_real_,
      max_jump      = if (length(dy) > 0) max(abs(dy)) else NA_real_,
      sd_diff       = if (length(dy) >= 2) sd(dy) else NA_real_,
      mean_abs_diff = if (length(dy) > 0) mean(abs(dy)) else NA_real_,
      p_jump_step   = if (length(dy) > 0) mean(abs(dy) > jump_step_cut) else 0,
      outlier_rate  = if (length(res) >= 5 && !is.na(s_res) && s_res > 0)
        mean(abs(res) > 3 * s_res) else 0
    )
  },
  by = point_id
]

# ============================================================
# 5) QC scenarios (meaningfully different constraints)
# ============================================================
qc_scenarios <- data.table(
  scenario = c("permissive", "balanced", "strict"),
  k = c(3, 2, 1),
  max_na_rate      = c(0.40, 0.30, 0.20),
  jump_abs_cut     = c(0.22, 0.18, 0.12),
  p_jump_rate_cut  = c(0.35, 0.25, 0.15),
  outlier_rate_cut = c(0.20, 0.12, 0.08)
)

robust_thr <- function(x, k) median(x, na.rm = TRUE) + k * mad(x, constant = 1, na.rm = TRUE)

# ============================================================
# 6) QC runner: returns keep_points for a given scenario row
# ============================================================
get_keep_points <- function(sc_row, ts_metrics) {
  sd_thr     <- robust_thr(ts_metrics$sd_ndvi,       sc_row$k)
  jump_thr   <- robust_thr(ts_metrics$max_jump,      sc_row$k)
  rough_thr  <- robust_thr(ts_metrics$sd_diff,       sc_row$k)
  madiff_thr <- robust_thr(ts_metrics$mean_abs_diff, sc_row$k)
  
  bad_points <- ts_metrics[
    n_years_obs   < min_obs_keep |
      na_rate       > sc_row$max_na_rate |
      sd_ndvi       > sd_thr |
      max_jump      > jump_thr |
      sd_diff       > rough_thr |
      mean_abs_diff > madiff_thr |
      max_jump      > sc_row$jump_abs_cut |
      p_jump_step   > sc_row$p_jump_rate_cut |
      outlier_rate  > sc_row$outlier_rate_cut,
    point_id
  ]
  
  keep_points <- setdiff(ts_metrics$point_id, bad_points)
  list(keep_points = keep_points, bad_points = bad_points)
}

# ============================================================
# 7) Per-scenario processing:
#    - write QC’d full table
#    - add per-point trend
#    - add latitude-based ecoregion
# ============================================================
add_trend_per_point <- function(dt_full) {
  # unique point-year (avoid band duplication)
  ndvi_py <- unique(dt_full[, .(point_id, year, ndvi)], by = c("point_id","year"))
  
  trend_dt <- ndvi_py[
    ,
    {
      y <- ndvi
      x <- year
      if (sum(!is.na(y)) < 2) {
        .(ndvi_slope = NA_real_)
      } else {
        fit <- lm(y ~ x)
        .(ndvi_slope = coef(fit)[["x"]])
      }
    },
    by = point_id
  ]
  trend_dt[, ndvi_slope_decade := ndvi_slope * 10]
  trend_dt
}

add_ecoregion_from_lat <- function(dt) {
  stopifnot("gps_lat" %in% names(dt))
  dt[, ecoregion := fcase(
    gps_lat >= 55,                 "boreal",
    gps_lat >= 40 & gps_lat < 55,  "temperate",
    gps_lat < 40,                  "mediterranean",
    default = NA_character_
  )]
  dt
}

scenario_outputs <- qc_scenarios[, .(
  scenario,
  qc_file = file.path(base_dir, paste0("LUCAS_BAP_comp_QC_", scenario, ".csv")),
  out_file = file.path(base_dir, paste0("LUCAS_BAP_comp_QC_", scenario, "_with_trend_ecogroup.csv"))
)]

for (i in seq_len(nrow(qc_scenarios))) {
  sc <- qc_scenarios[i]
  
  res <- get_keep_points(sc, ts_metrics)
  keep_points <- res$keep_points
  bad_points  <- res$bad_points
  
  cat("\n--- Scenario:", sc$scenario, "---\n")
  cat("Points before QC:", uniqueN(ts_metrics$point_id), "\n")
  cat("Points removed   :", length(bad_points), "\n")
  cat("Points kept      :", length(keep_points), "\n")
  
  # Apply keep list to FULL raw table (keeps all metadata/columns)
  dt_qc <- LUCAS_BAP_comp[point_id %in% keep_points]
  
  # Write QC-only table
  qc_file <- scenario_outputs[scenario == sc$scenario, qc_file]
  fwrite(dt_qc, qc_file)
  cat("Wrote:", qc_file, "\n")
  
  # Add trends (point-level) back to full table
  dt_qc[, year := as.integer(year)]
  trend_dt <- add_trend_per_point(dt_qc)
  
  # Remove any old trend columns to avoid duplication
  if ("ndvi_slope" %in% names(dt_qc)) dt_qc[, c("ndvi_slope","ndvi_slope_decade") := NULL]
  dt_qc <- trend_dt[dt_qc, on = "point_id"]
  
  # Add ecoregion and write final output
  dt_qc <- add_ecoregion_from_lat(dt_qc)
  
  out_file <- scenario_outputs[scenario == sc$scenario, out_file]
  fwrite(dt_qc, out_file)
  cat("Wrote:", out_file, "\n")
}

# ============================================================
# 8) Plot (balanced only): NDVI trend distributions by LC facets
#    - point-level table (one row per point_id)
#    - fixed LC facet order + n in facet strip
# ============================================================
balanced_file <- scenario_outputs[scenario == "balanced", out_file]
dt_bal <- fread(balanced_file)
setDT(dt_bal)

# point-level unique table
pts_bal <- unique(
  dt_bal[, .(point_id, lc1_label, ecoregion, ndvi_slope_decade)],
  by = "point_id"
)

# merge selected lc labels to "Artificial surfaces" (optional; as in your script)
merge_to_artificial <- c(
  "Non built-up linear features",
  "Non built-up area features",
  "Buildings with 1 to 3 floors",
  "Buildings with more than 3 floors",
  "Greenhouses", "Bare land"
)
pts_bal[lc1_label %in% merge_to_artificial, lc1_label := "Artificial surfaces"]

# enforce desired facet order
lc_facet_order <- c(
  "Coniferous woodland",
  "Broadleaved woodland",
  "Mixed woodland",
  "Shrubland without tree cover",
  "Grassland without tree/shrub cover",
  "Artificial surfaces"
)
pts_bal[, lc1_label := factor(lc1_label, levels = lc_facet_order)]

# ecoregion order (match your plots)
pts_bal[, ecoregion := factor(ecoregion, levels = c("boreal", "temperate", "mediterranean"))]

# facet labels with n
lc_n <- pts_bal[!is.na(lc1_label), .(n = uniqueN(point_id)), by = lc1_label]
strip_labs <- setNames(
  paste0(as.character(lc_n$lc1_label), "\n(n = ", lc_n$n, ")"),
  as.character(lc_n$lc1_label)
)

p_bal <- ggplot(
  pts_bal[!is.na(ecoregion) & !is.na(lc1_label)],
  aes(x = ecoregion, y = ndvi_slope_decade, fill = ecoregion)
) +
  geom_violin(trim = TRUE, alpha = 0.80, color = NA) +
  geom_boxplot(width = 0.16, outlier.size = 0.6, alpha = 0.95, color = "grey15") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  facet_wrap(
    ~ lc1_label, ncol = 3,
    labeller = labeller(lc1_label = strip_labs)
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.position = "none"
  ) +
  labs(x = "Ecoregion", y = "NDVI trend (per decade)", title = "")

print(p_bal)

out_plot <- file.path(base_dir, "fig_balanced_trends_by_ecoregion_and_lc.png")
ggsave(out_plot, p_bal, width = 12, height = 8, dpi = 300)
cat("\nSaved plot:", out_plot, "\n")
