library(terra)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
pts_path <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/candidates_strict.csv"
esa_path <- "/mnt/dss_project/lmandl/_unmixing/copernicus_cover_data/00_mosaics/mosaic_ESAlandcover2020_cog.tif"

radius <- 60   # meters (~2 pixels)
p_match <- 1.0 # 1.0 = strict, 0.9 = softer

out_csv  <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/candidates_buffer60.csv"
out_gpkg <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/candidates_buffer60.gpkg"

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
pts_df <- read.csv(pts_path)
pts_wgs84 <- vect(pts_df, geom=c("x","y"), crs="EPSG:4326")

lc <- rast(esa_path)[[1]]

# Project points to raster CRS
pts_r <- project(pts_wgs84, crs(lc))
pts_3035 <- project(pts_wgs84, "EPSG:3035")

# ------------------------------------------------------------
# Center pixel class
# ------------------------------------------------------------
center_vals <- extract(lc, pts_r)[,2]

# ------------------------------------------------------------
# Buffer extraction
# ------------------------------------------------------------
buf <- buffer(pts_r, width=radius)

vals <- extract(lc, buf)

# vals is a long table: ID + raster values
split_vals <- split(vals[,2], vals[,1])

keep <- logical(length(center_vals))

for (i in seq_along(split_vals)) {
  v0 <- center_vals[i]
  v <- split_vals[[i]]
  
  if (is.na(v0)) {
    keep[i] <- FALSE
    next
  }
  
  v <- v[!is.na(v)]
  
  if (length(v) == 0) {
    keep[i] <- FALSE
    next
  }
  
  keep[i] <- mean(v == v0) >= p_match
}

# ------------------------------------------------------------
# Filter points
# ------------------------------------------------------------
pts_3035$keep_buffer60 <- keep
pts_filt <- pts_3035[keep, ]

message("Kept points: ", nrow(pts_filt), " / ", nrow(pts_3035))

write.csv(as.data.frame(pts_filt), out_csv, row.names=FALSE)
writeVector(pts_filt, out_gpkg, overwrite=TRUE)



pts_filt <- pts_3035[keep, ]

# count unique points instead of rows
n_points_kept  <- length(unique(pts_filt$point_id))
n_points_total <- length(unique(pts_3035$point_id))

message("Kept point_ids: ", n_points_kept, " / ", n_points_total)

write.csv(as.data.frame(pts_filt), out_csv, row.names=FALSE)
writeVector(pts_filt, out_gpkg, overwrite=TRUE)



library(data.table)

dt <- as.data.frame(pts_filt)
dt <- as.data.table(dt)

# unique points per ecoregion x lc1
tab <- dt[, .(n_point_id = uniqueN(point_id)),
          by = .(ecoregion, lc1)][order(ecoregion, lc1)]

print(tab)

# (optional) totals per ecoregion
tab_eco <- dt[, .(n_point_id = uniqueN(point_id)), by = ecoregion][order(ecoregion)]
print(tab_eco)

# (optional) totals per lc1
tab_lc1 <- dt[, .(n_point_id = uniqueN(point_id)), by = lc1][order(lc1)]
print(tab_lc1)



library(data.table)
library(terra)

dt <- as.data.table(as.data.frame(pts_filt))
dt[, row_id := .I]
dt[, point_id := trimws(as.character(point_id))]

# one row per point_id
dt_unique <- dt[, .SD[1], by = point_id]

# SUBSET BY ROW INDEX (not %in%)
pts_unique <- pts_filt[dt_unique$row_id, ]

# sanity checks (by point_id)
stopifnot(length(unique(pts_unique$point_id)) == nrow(pts_unique))

# write
out_csv_unique  <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/candidates_buffer60_unique.csv"
out_gpkg_unique <- "/mnt/dss_project/lmandl/_unmixing/_spectral_library/candidates_buffer60_unique.gpkg"

write.csv(as.data.frame(pts_unique), out_csv_unique, row.names = FALSE)
writeVector(pts_unique, out_gpkg_unique, overwrite = TRUE)



#-------------------------------------------------------------------------------
### do the same for topsoil points
#-------------------------------------------------------------------------------

library(terra)

# ============================================================
# Radius-based land-cover homogeneity filtering (60 m)
# for ESDAC topsoil point dataset (GPKG)
#
# Steps:
#  1) Load point GPKG
#  2) Keep only LC0_Des in {Woodland, Bareland}
#  3) Load ESA land-cover raster
#  4) For each point, buffer by 60 m and extract raster values
#  5) Keep point if >= p_match of buffered values equals center pixel value
#  6) Write filtered points to CSV + GPKG
# ============================================================

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
library(data.table)
pck_list = c("lubridate", "ggplot2", "viridis","tidyr", "dplyr","purrr",
             "scales")
lapply(pck_list, require, character.only = TRUE)

file_path <- "/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/2015/spectra_EU28_merged.csv"

spectra <- fread(file_path)

str(spectra)

ds <- spectra
ds <- ds %>%
  mutate(across(-c(source, SampleID, PointID, NUTS_0, SampleN), as.numeric))

ds <- pivot_longer(ds, -c(source, SampleID, PointID, NUTS_0, SampleN))
ds$name <- substr(ds$name, 2,nchar(ds$name))
ds$wl <- as.numeric(ds$name)

ds$reflectance <- 1/(10^ds$value)
ds <- subset(ds, ds$wl >= 426)







gpkg_path <- "/mnt/dss_project/lmandl/_unmixing/esdac_topsoil/bd2018_fef_landcover_joined.gpkg"
esa_path  <- "/mnt/dss_project/lmandl/_unmixing/copernicus_cover_data/00_mosaics/mosaic_ESAlandcover2020_cog.tif"

radius  <- 60    # meters (~2 pixels for 30 m raster)
p_match <- 1.0   # 1.0 strict; try 0.9 if too harsh

# Filter classes in the point dataset BEFORE buffering
keep_lc0 <- c("Woodland", "Bareland")

# Outputs
out_csv  <- "/mnt/dss_project/lmandl/_unmixing/candidates_buffer60_topsoil.csv"
out_gpkg <- "/mnt/dss_project/lmandl/_unmixing/candidates_buffer60_topsoil.gpkg"

# ------------------------------------------------------------
# 1) Load point GPKG
# ------------------------------------------------------------
pts <- vect(gpkg_path)

# Basic checks
if (!("LC0_Desc" %in% names(pts))) {
  stop("Column 'LC0_Desc' not found in the GPKG attributes. Available columns: ",
       paste(names(pts), collapse = ", "))
}

# ------------------------------------------------------------
# 2) Filter points to LC0_Des == Woodland or Bareland
# ------------------------------------------------------------
pts$LC0_Desc <- trimws(as.character(pts$LC0_Desc))
pts_sub <- pts[pts$LC0_Desc %in% keep_lc0, ]

message("Points after LC0_Desc filter: ", nrow(pts_sub), " / ", nrow(pts))

if (nrow(pts_sub) == 0) {
  stop("After filtering for LC0_Desc in {", paste(keep_lc0, collapse = ", "),
       "}, no points remain. Check spelling/case in LC0_Des.")
}

# ------------------------------------------------------------
# 3) Load ESA land-cover raster (single band)
# ------------------------------------------------------------
lc <- rast(esa_path)[[1]]

# ------------------------------------------------------------
# 4) Project points to raster CRS for extraction/buffering
# ------------------------------------------------------------
pts_r <- project(pts_sub, crs(lc))

# ------------------------------------------------------------
# 5) Center pixel class at each point
# ------------------------------------------------------------
center_vals <- extract(lc, pts_r)[, 2]

# ------------------------------------------------------------
# 6) Buffer extraction within radius
# ------------------------------------------------------------
buf <- buffer(pts_r, width = radius)

vals <- extract(lc, buf)  # long table: ID + value
split_vals <- split(vals[, 2], vals[, 1])

# ------------------------------------------------------------
# 7) Decide keep/drop per point
# ------------------------------------------------------------
keep <- logical(length(center_vals))

for (i in seq_along(split_vals)) {
  v0 <- center_vals[i]
  v  <- split_vals[[i]]
  
  if (is.na(v0)) { keep[i] <- FALSE; next }
  
  v <- v[!is.na(v)]
  if (length(v) == 0) { keep[i] <- FALSE; next }
  
  keep[i] <- mean(v == v0) >= p_match
}

# Attach + filter (keep in original CRS of the GPKG for writing unless you prefer otherwise)
pts_sub$keep_buffer60 <- keep
pts_filt <- pts_sub[keep, ]

message("Kept rows: ", nrow(pts_filt), " / ", nrow(pts_sub))

# Optional: if there is a point_id column, also report unique IDs
if ("point_id" %in% names(pts_filt)) {
  n_points_kept  <- length(unique(pts_filt$point_id))
  n_points_total <- length(unique(pts_sub$point_id))
  message("Kept point_ids: ", n_points_kept, " / ", n_points_total)
}

# ------------------------------------------------------------
# 8) Write outputs
# ------------------------------------------------------------
write.csv(as.data.frame(pts_filt), out_csv, row.names = FALSE)
writeVector(pts_filt, out_gpkg, overwrite = TRUE)

message("Wrote CSV:  ", out_csv)
message("Wrote GPKG: ", out_gpkg)







