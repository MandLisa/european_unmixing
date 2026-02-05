library(terra)
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)


points_csv <- "/mnt/eo/EU_unmixing/data/LUCAS/LUCAS_multiyear_filtered_unique.csv"
bap_dir    <- "/mnt/dss_europe/mosaics_eu/mosaics_eu_baps"
out_csv    <- "/mnt/eo/EU_unmixing/data/LUCAS/LUCAS_BAP_6bands_long.csv"



# --- read points ---
pts_dt <- fread(points_csv)
stopifnot(all(c("gps_lat", "gps_long") %in% names(pts_dt)))
if (!("point_id" %in% names(pts_dt))) pts_dt[, point_id := .I]

pts_v_wgs <- vect(pts_dt, geom = c("gps_long", "gps_lat"), crs = "EPSG:4326")

# --- list yearly VRTs ---
vrt_files <- list.files(
  bap_dir,
  pattern = "^[0-9]{4}_mosaic_eu_cog_v2\\.vrt$",
  full.names = TRUE
)
years <- as.integer(sub("^([0-9]{4}).*$", "\\1", basename(vrt_files)))

keep <- years >= 1984 & years <= 2024
vrt_files <- vrt_files[keep]
years     <- years[keep]
o <- order(years)
vrt_files <- vrt_files[o]
years     <- years[o]
stopifnot(length(vrt_files) > 0)

# --- project points once to CRS of first raster (assumes consistent CRS across years) ---
r0 <- rast(vrt_files[1])
pts_v <- pts_v_wgs
if (!same.crs(pts_v, r0)) pts_v <- project(pts_v, crs(r0))

# --- write header ---
fwrite(data.table(point_id=integer(), year=integer(), band=integer(), BAP=numeric()),
       out_csv)

# --- loop years ---
for (i in seq_along(vrt_files)) {
  
  r <- rast(vrt_files[i])
  
  # Safety: if some year differs in CRS (unlikely), reproject on the fly
  pts_use <- pts_v
  if (!same.crs(pts_use, r)) pts_use <- project(pts_v_wgs, crs(r))
  
  vals <- terra::extract(r, pts_use, ID = FALSE)  # n_points x 6
  
  # Convert to long: one row per point-year-band
  out_i <- data.table(point_id = pts_dt$point_id, vals)
  out_long <- melt(
    out_i,
    id.vars = "point_id",
    variable.name = "band_name",
    value.name = "BAP"
  )
  
  # turn band name into band index 1..6 (robust given your naming)
  out_long[, band := as.integer(sub(".*_([0-9]+)$", "\\1", band_name))]
  out_long[, c("band_name") := NULL]
  out_long[, year := years[i]]
  
  setcolorder(out_long, c("point_id", "year", "band", "BAP"))
  
  fwrite(out_long, out_csv, append = TRUE)
  
  rm(r, vals, out_i, out_long)
  gc(FALSE)
}

### import
LUCAS_BAP <- read_csv("/mnt/eo/EU_unmixing/data/LUCAS/LUCAS_BAP_6bands_long.csv")

df <- LUCAS_BAP_6bands_long %>%
  mutate(
    year = as.integer(year),
    band = as.integer(band),
    point_id = as.character(point_id)
  )

set.seed(1)
band_to_plot <- 3
n_points <- 12

ids <- df %>%
  filter(band == band_to_plot) %>%
  distinct(point_id) %>%
  slice_sample(n = n_points) %>%   # <- now n is constant
  pull(point_id)


ids <- df %>%
  filter(band == band_to_plot) %>%
  distinct(point_id) %>%
  slice_sample(n = min(n_points, n())) %>%
  pull(point_id)

df %>%
  filter(band == band_to_plot, point_id %in% ids) %>%
  ggplot(aes(x = year, y = BAP, group = point_id)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 0.9, alpha = 0.8) +
  facet_wrap(~ point_id, scales = "free_y") +
  labs(x = "Year", y = paste0("BAP (band ", band_to_plot, ")"))



