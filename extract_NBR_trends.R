# ============================================================
# Extract per-point NBR Sen slope (1984–2024) from yearly rasters
# - Uses lucas_unique_sf (sf POINT object) already in your environment
# - Memory-efficient: reads one raster at a time, processes points in chunks
# - Robust: handles nodata, constant series, and avoids qnorm/CI warnings
# ============================================================

library(sf)
library(terra)

# -----------------------
# User inputs
# -----------------------
nbr_dir     <- "/mnt/eo/eu_mosaics/NBR_comp"
year_min    <- 1984L
year_max    <- 2024L
chunk_size  <- 20000L   # RAM knob: 5000 if tight, 50000 if plenty
min_years   <- 5L       # require at least this many valid annual values
nodata_thr  <- -9990    # treat values <= this as NA (adjust if needed)

# -----------------------
# 0) Checks
# -----------------------
stopifnot(exists("lucas_unique_sf"))
stopifnot(inherits(lucas_unique_sf, "sf"))

# -----------------------
# 1) sf -> terra points
# -----------------------
pts <- terra::vect(lucas_unique_sf)

# stable id (for joining/debugging)
if (!("point_id" %in% names(pts))) {
  pts$point_id <- as.character(seq_len(nrow(pts)))
}

# -----------------------
# 2) List & sort yearly NBR rasters
# -----------------------
nbr_files <- list.files(nbr_dir, pattern = "^NBR_[0-9]{4}\\.tif$", full.names = TRUE)
stopifnot(length(nbr_files) > 0)

yrs <- as.integer(sub("^.*NBR_([0-9]{4})\\.tif$", "\\1", nbr_files))
o   <- order(yrs)
nbr_files <- nbr_files[o]
yrs       <- yrs[o]

# restrict to target period
keep <- yrs >= year_min & yrs <= year_max
nbr_files <- nbr_files[keep]
yrs       <- yrs[keep]

stopifnot(length(nbr_files) >= min_years)

# -----------------------
# 3) Project points once to raster CRS
# -----------------------
r0   <- rast(nbr_files[1])
ptsR <- project(pts, crs(r0))

# -----------------------
# 4) Robust Sen slope (median of pairwise slopes; slope only)
#    Avoids CI/p-value computations that trigger qnorm warnings
# -----------------------
sen_slope_only <- function(y, x, min_years = 5L) {
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]; x <- x[ok]
  n <- length(y)
  if (n < min_years) return(NA_real_)
  if (all(y == y[1])) return(0)  # constant time series => slope 0
  
  # pairwise slopes (n ~ 41 max here; fine)
  s <- c()
  for (i in 1:(n - 1)) {
    dx <- x[(i + 1):n] - x[i]
    s  <- c(s, (y[(i + 1):n] - y[i]) / dx)
  }
  stats::median(s, na.rm = TRUE)
}

# -----------------------
# 5) Chunked extraction + slope computation
# -----------------------
n_pts   <- nrow(ptsR)
n_years <- length(nbr_files)

out_slope <- rep(NA_real_, n_pts)
out_nobs  <- rep(0L,       n_pts)  # number of valid annual obs used

starts <- seq(1L, n_pts, by = chunk_size)

for (s0 in starts) {
  s1 <- min(s0 + chunk_size - 1L, n_pts)
  pts_chunk <- ptsR[s0:s1]
  
  # time series matrix for this chunk only
  Y <- matrix(NA_real_, nrow = nrow(pts_chunk), ncol = n_years)
  
  for (j in seq_along(nbr_files)) {
    r <- rast(nbr_files[j])
    
    v <- terra::extract(r, pts_chunk)[, 2]
    
    # nodata handling (adjust threshold to your product if needed)
    v[!is.finite(v)] <- NA_real_
    v[v <= nodata_thr] <- NA_real_
    
    Y[, j] <- v
  }
  
  # count valid years per point (for diagnostics / filtering)
  n_ok <- rowSums(is.finite(Y))
  out_nobs[s0:s1] <- n_ok
  
  # compute Sen slope per point
  out_slope[s0:s1] <- apply(Y, 1, sen_slope_only, x = yrs, min_years = min_years)
  
  message(sprintf("Done %d-%d / %d points", s0, s1, n_pts))
}

# -----------------------
# 6) Attach results back to your sf object
# -----------------------
lucas_unique_sf$nbr_sen_slope <- out_slope
lucas_unique_sf$nbr_n_years   <- out_nobs

# -----------------------
# 7) Optional: rescale if your NBR is stored as scaled int (e.g., *10000)
#    Uncomment if needed.
# -----------------------
# lucas_unique_sf$nbr_sen_slope <- lucas_unique_sf$nbr_sen_slope / 10000

# -----------------------
# 8) Optional quick diagnostics
# -----------------------
summary(lucas_unique_sf$nbr_n_years)
summary(lucas_unique_sf$nbr_sen_slope)

# If many points have low nbr_n_years, your nodata_thr may be wrong.
# To inspect value range in a sample raster:
# global(r0, range, na.rm=TRUE)