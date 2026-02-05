############################################################
## LUCAS NDVI trajectory QC – THREE MEANINGFUL SCENARIOS
## + Scenario-wise coord_id-based maps + GPKG export (NO Moran's I)
############################################################

library(data.table)
library(ggplot2)
library(sf)

# ============================================================
# 0) Setup + Load
# ============================================================
base_dir <- "/mnt/eo/EU_unmixing/data/LUCAS"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

LUCAS_BAP_comp <- fread(file.path(base_dir, "LUCAS_BAP_full.csv"))
setDT(LUCAS_BAP_comp)

years_full <- 1984:2024
target_lc  <- c("A10","A20","C10","C20","C30","D20","E20","F00")

filter_incomplete_years <- TRUE
min_obs_keep  <- 10
jump_step_cut <- 0.10

# ============================================================
# 1) Hygiene + lc recode
# ============================================================
LUCAS_BAP_comp[, lc1 := trimws(as.character(lc1))]
LUCAS_BAP_comp[, year := as.integer(year)]

LUCAS_BAP_comp[
  ,
  lc1 := fcase(
    lc1 %in% c("A11","A12","A13"), "A10",
    lc1 %in% c("A21","A22"),      "A20",
    default = lc1
  )
]

# ============================================================
# 2) Stable lc mapping (MODE)
# ============================================================
lc_mode <- LUCAS_BAP_comp[
  lc1 %in% target_lc, .N, by = .(point_id, lc1)
][order(point_id, -N)][, .SD[1], by = point_id][, .(point_id, lc1)]

# ============================================================
# 3) NDVI point–year table (trajectory QC basis)
# ============================================================
ndvi_py <- unique(
  LUCAS_BAP_comp[
    point_id %in% lc_mode$point_id & year %in% years_full,
    .(point_id, year, ndvi)
  ],
  by = c("point_id","year")
)

ndvi_py <- lc_mode[ndvi_py, on = "point_id", nomatch = 0L]
setorder(ndvi_py, point_id, year)

if (filter_incomplete_years) {
  full_points <- ndvi_py[, .(n = uniqueN(year)), by = point_id][n == length(years_full), point_id]
  cat("Points with full 1984–2024 coverage:", length(full_points), "\n")
  ndvi_py <- ndvi_py[point_id %in% full_points]
}

# ============================================================
# 4) Trajectory QC metrics
# ============================================================
ts_metrics <- ndvi_py[
  ,
  {
    y  <- ndvi
    yo <- y[!is.na(y)]
    dy <- if (length(yo) > 1) diff(yo) else numeric(0)
    
    .(
      n_years_obs   = length(yo),
      na_rate       = mean(is.na(y)),
      sd_ndvi       = if (length(yo) > 1) sd(yo) else NA_real_,
      max_jump      = if (length(dy)) max(abs(dy)) else NA_real_,
      sd_diff       = if (length(dy) > 1) sd(dy) else NA_real_,
      mean_abs_diff = if (length(dy)) mean(abs(dy)) else NA_real_,
      p_jump_step   = mean(abs(dy) > jump_step_cut),
      outlier_rate  = if (length(yo) > 4)
        mean(abs(yo - median(yo)) > 3 * mad(yo, constant = 1)) else 0
    )
  },
  by = point_id
]

# ============================================================
# 5) QC scenarios
# ============================================================
qc_scenarios <- data.table(
  scenario = c("permissive","balanced","strict"),
  k = c(3,2,1),
  max_na_rate      = c(0.40,0.30,0.20),
  jump_abs_cut     = c(0.22,0.18,0.12),
  p_jump_rate_cut  = c(0.35,0.25,0.15),
  outlier_rate_cut = c(0.20,0.12,0.08)
)

robust_thr <- function(x,k) median(x,na.rm=TRUE) + k*mad(x,constant=1,na.rm=TRUE)

get_keep_points <- function(sc, tm) {
  bad <- tm[
    n_years_obs < min_obs_keep |
      na_rate > sc$max_na_rate |
      sd_ndvi > robust_thr(tm$sd_ndvi, sc$k) |
      max_jump > robust_thr(tm$max_jump, sc$k) |
      sd_diff > robust_thr(tm$sd_diff, sc$k) |
      mean_abs_diff > robust_thr(tm$mean_abs_diff, sc$k) |
      max_jump > sc$jump_abs_cut |
      p_jump_step > sc$p_jump_rate_cut |
      outlier_rate > sc$outlier_rate_cut,
    point_id
  ]
  setdiff(tm$point_id, bad)
}

# ============================================================
# 6) Helpers: trend + ecoregion
# ============================================================
add_trend <- function(dt) {
  t <- unique(dt[, .(point_id, year, ndvi)], by=c("point_id","year"))[
    , .(
      ndvi_slope = if (sum(!is.na(ndvi)) > 1) coef(lm(ndvi ~ year))[2] else NA_real_
    ),
    by = point_id
  ]
  t[, ndvi_slope_decade := ndvi_slope * 10]
  t
}

add_ecoregion <- function(dt) {
  stopifnot("gps_lat" %in% names(dt))
  dt[, ecoregion := fcase(
    gps_lat >= 55, "boreal",
    gps_lat >= 40 & gps_lat < 55, "temperate",
    gps_lat < 40, "mediterranean",
    default = NA_character_
  )]
  dt
}

# ============================================================
# 7) Run QC for each scenario + write outputs
# ============================================================
scenario_files <- qc_scenarios[, .(
  scenario,
  out_csv = file.path(base_dir, paste0("LUCAS_BAP_comp_QC_", scenario, "_with_trend_ecogroup.csv"))
)]

for (i in seq_len(nrow(qc_scenarios))) {
  sc <- qc_scenarios[i]
  
  keep <- get_keep_points(sc, ts_metrics)
  cat("\nScenario:", sc$scenario, "-> keep points:", length(keep), "\n")
  
  dt_qc <- LUCAS_BAP_comp[point_id %in% keep]
  
  # add trends
  tr <- add_trend(dt_qc)
  dt_qc <- tr[dt_qc, on="point_id"]
  
  # add ecoregion
  dt_qc <- add_ecoregion(dt_qc)
  
  # write
  out_csv <- scenario_files[scenario == sc$scenario, out_csv]
  fwrite(dt_qc, out_csv)
  cat("Wrote:", out_csv, "\n")
}

# ============================================================
# 8) Scenario-wise MAPS + GPKG export using coord_id (lon_lat)
# ============================================================
make_sf_from_coord_id <- function(dt) {
  stopifnot("coord_id" %in% names(dt))
  dt[, coord_ok := grepl("^-?[0-9.]+_-?[0-9.]+$", coord_id)]
  dt[coord_ok == TRUE, c("lon","lat") := tstrsplit(coord_id, "_", fixed=TRUE)]
  dt[, `:=`(lon = as.numeric(lon), lat = as.numeric(lat))]
  
  # point-level unique table
  pts <- unique(
    dt[!is.na(lon) & !is.na(lat),
       .(point_id, coord_id, lon, lat, ndvi_slope_decade, ecoregion, lc1_label)],
    by = "point_id"
  )
  
  st_as_sf(pts, coords=c("lon","lat"), crs=4326, remove=FALSE)
}

plot_map <- function(pts_sf, scenario) {
  ggplot(pts_sf) +
    geom_sf(aes(color = ndvi_slope_decade), size=0.6) +
    scale_color_gradient2(low="darkblue", mid="white", high="brown", midpoint=0) +
    coord_sf(datum=NA) +
    theme_minimal(base_size=13) +
    theme(
      panel.grid.major = element_line(color="grey85", linewidth=0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_blank()
    ) +
    labs(color="NDVI trend (per decade)",
         title=paste0("NDVI trend map – ", scenario))
}

for (i in seq_len(nrow(scenario_files))) {
  sc <- scenario_files[i]
  
  cat("\n--- Mapping scenario:", sc$scenario, "---\n")
  dt_sc <- fread(sc$out_csv)
  setDT(dt_sc)
  
  # ensure ndvi_slope_decade exists
  if (!("ndvi_slope_decade" %in% names(dt_sc)) && ("ndvi_slope" %in% names(dt_sc))) {
    dt_sc[, ndvi_slope_decade := ndvi_slope * 10]
  }
  
  pts_sf <- make_sf_from_coord_id(dt_sc)
  
  # export gpkg
  out_gpkg <- file.path(base_dir, paste0("ndvi_trend_points_", sc$scenario, "_from_coord_id.gpkg"))
  st_write(pts_sf, out_gpkg, layer=paste0("ndvi_trends_", sc$scenario), delete_layer=TRUE)
  cat("Wrote:", out_gpkg, "\n")
  
  # plot + save
  p <- plot_map(pts_sf, sc$scenario)
  print(p)
  
  out_png <- file.path(base_dir, paste0("fig_map_ndvi_trend_", sc$scenario, "_coord_id.png"))
  ggsave(out_png, p, width=10, height=6, dpi=300)
  cat("Wrote:", out_png, "\n")
}

cat("\nDone. You now have one GPKG + one PNG map per QC scenario.\n")
