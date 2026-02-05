############################################################
## Add reliable coordinates from coord_id / coord_info
## to the THREE scenario files:
## - adds numeric columns: x (lon), y (lat)
## - extracted from string like "-9.360576_38.688748"
## - keeps existing gps_long/gps_lat untouched
## - writes UPDATED files (either overwrite or new suffix)
############################################################

library(data.table)

base_dir <- "/mnt/eo/EU_unmixing/data/LUCAS"

# the three scenario files you already produced
scenarios <- c("permissive", "balanced", "strict")

in_files <- file.path(
  base_dir,
  paste0("LUCAS_BAP_comp_QC_", scenarios, "_with_trend_ecogroup.csv")
)

# Choose output style:
# - overwrite = TRUE  -> overwrite the original files
# - overwrite = FALSE -> write new files with "_xy" suffix
overwrite <- FALSE

out_files <- if (overwrite) {
  in_files
} else {
  file.path(
    base_dir,
    paste0("LUCAS_BAP_comp_QC_", scenarios, "_with_trend_ecogroup_xy.csv")
  )
}

for (i in seq_along(in_files)) {
  
  cat("\n--- Adding x/y to:", basename(in_files[i]), "---\n")
  
  dt <- fread(in_files[i])
  setDT(dt)
  
  # support either column name: coord_id (your earlier name) or coord_info (as you wrote now)
  coord_col <- NULL
  if ("coord_id" %in% names(dt)) coord_col <- "coord_id"
  if (is.null(coord_col) && "coord_info" %in% names(dt)) coord_col <- "coord_info"
  
  stopifnot(!is.null(coord_col))  # fail fast if neither exists
  
  # Parse lon/lat from coord string
  dt[, coord_ok := grepl("^-?[0-9.]+_-?[0-9.]+$", get(coord_col))]
  
  # initialise x/y as NA
  dt[, `:=`(x = NA_real_, y = NA_real_)]
  
  # split only where format is OK
  dt[coord_ok == TRUE, c("x", "y") := {
    tmp <- tstrsplit(get(coord_col), "_", fixed = TRUE)
    list(as.numeric(tmp[[1]]), as.numeric(tmp[[2]]))
  }]
  
  # optional: sanity check
  dt[, xy_ok := !is.na(x) & !is.na(y) & x >= -180 & x <= 180 & y >= -90 & y <= 90]
  cat("Rows:", nrow(dt), "\n")
  cat("coord_ok fraction:", round(mean(dt$coord_ok, na.rm = TRUE), 4), "\n")
  cat("xy_ok fraction   :", round(mean(dt$xy_ok, na.rm = TRUE), 4), "\n")
  
  # cleanup helper columns if you prefer (keep if you want diagnostics)
  dt[, c("coord_ok", "xy_ok") := NULL]
  
  fwrite(dt, out_files[i])
  cat("Wrote:", out_files[i], "\n")
}

cat("\nDone.\n")
