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


