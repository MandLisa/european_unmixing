library(sf)
library(terra)
library(trend)

# -----------------------
# Inputs
# -----------------------
nbr_dir <- "/mnt/eo/eu_mosaics/NBR_comp"

# -----------------------
# 1) sf -> terra points
# -----------------------
stopifnot(inherits(lucas_unique_sf, "sf"))

pts <- terra::vect(lucas_unique_sf)

# stable id (so you can join later if needed)
if (!("point_id" %in% names(pts))) {
  pts$point_id <- as.character(seq_len(nrow(pts)))
}

# -----------------------
# 2) List yearly NBR rasters
# -----------------------
nbr_files <- list.files(nbr_dir, pattern = "^NBR_[0-9]{4}\\.tif$", full.names = TRUE)
stopifnot(length(nbr_files) > 0)

yrs <- as.integer(sub("^.*NBR_([0-9]{4})\\.tif$", "\\1", nbr_files))
o   <- order(yrs)
nbr_files <- nbr_files[o]
yrs       <- yrs[o]

# Restrict to 1984–2024
keep <- yrs >= 1984 & yrs <= 2024
nbr_files <- nbr_files[keep]
yrs       <- yrs[keep]

stopifnot(length(nbr_files) >= 5)

# -----------------------
# 3) Project points once to raster CRS
# -----------------------
r0   <- rast(nbr_files[1])
ptsR <- project(pts, crs(r0))

# -----------------------
# 4) Chunked extraction + Sen slope
# -----------------------
chunk_size <- 20000L  # lower if RAM is tight (e.g., 5000), higher if you have room

n_pts   <- nrow(ptsR)
n_years <- length(nbr_files)

out_slope <- rep(NA_real_, n_pts)
out_p     <- rep(NA_real_, n_pts)

sen1 <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 5) return(c(NA_real_, NA_real_))
  s <- trend::sens.slope(y[ok], x[ok])
  c(as.numeric(s$estimates), as.numeric(s$p.value))
}

starts <- seq(1L, n_pts, by = chunk_size)

for (s0 in starts) {
  s1 <- min(s0 + chunk_size - 1L, n_pts)
  pts_chunk <- ptsR[s0:s1]
  
  # matrix exists only for this chunk
  Y <- matrix(NA_real_, nrow = nrow(pts_chunk), ncol = n_years)
  
  for (j in seq_along(nbr_files)) {
    r <- rast(nbr_files[j])
    # extract returns ID + value; take value col
    Y[, j] <- terra::extract(r, pts_chunk)[, 2]
  }
  
  res <- t(apply(Y, 1, sen1, x = yrs))
  out_slope[s0:s1] <- res[, 1]
  out_p[s0:s1]     <- res[, 2]
  
  message(sprintf("Done %d-%d / %d points", s0, s1, n_pts))
}

# -----------------------
# 5) Attach results back to sf
# -----------------------
lucas_unique_sf$nbr_sen_slope <- out_slope
lucas_unique_sf$nbr_sen_p     <- out_p