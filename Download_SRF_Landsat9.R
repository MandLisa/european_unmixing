library(dplyr)
library(tidyr)
library(stringr)

# ------------------------------------------------------------
# 1) Download & unpack Landsat-9 OLI SRFs (works similarly for L8 OLI)
# ------------------------------------------------------------
out_dir <- "/mnt/dss_project/lmandl/_unmixing/SRF"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tar_url <- "https://nwp-saf.eumetsat.int/downloads/rtcoef_info/visir_srf/rtcoef_landsat_9_oli_srf/rtcoef_landsat_9_oli_srf.tar.gz"
tar_path <- file.path(out_dir, "rtcoef_landsat_9_oli_srf.tar.gz")

download.file(tar_url, tar_path, mode = "wb")
untar(tar_path, exdir = out_dir)

# After untar you should have files like:
# rtcoef_landsat_9_oli_srf_ch01.txt ... ch09.txt

# ------------------------------------------------------------
# 2) Helper: read one SRF TXT (it's one long line with header + numeric pairs)
# ------------------------------------------------------------
read_nwpsaf_srf_txt <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- paste(x, collapse = " ")
  
  # extract all numeric tokens (includes integers, decimals)
  nums <- str_extract_all(x, "-?\\d+\\.?\\d*(?:[eE][+-]?\\d+)?")[[1]]
  nums <- as.numeric(nums)
  
  # The file contains: ... "Number of data points: N" then 2*N numeric values (wn, resp)
  # N is typically the first "reasonable" integer after that phrase, but easiest:
  # take the last even-length block and split into pairs.
  if (length(nums) < 4) stop("No numeric content parsed from SRF file: ", path)
  
  # If the file includes extra numbers in the header, the last 2*k values are the data pairs.
  # Here we assume data pairs are the last even count of numbers.
  if (length(nums) %% 2 == 1) nums <- nums[-1]  # drop leading token if odd count
  
  wn   <- nums[seq(1, length(nums), by = 2)]     # wavenumber cm-1
  resp <- nums[seq(2, length(nums), by = 2)]     # relative response
  
  # Convert wavenumber -> wavelength (nm)
  wl_nm <- 1e7 / wn
  
  tibble(wl_nm = wl_nm, response = resp)
}

# ------------------------------------------------------------
# 3) Read all channels and create a *wide* SRF table on integer wavelengths
# ------------------------------------------------------------
srf_files <- list.files(out_dir, pattern = "rtcoef_landsat_9_oli_srf_ch\\d+\\.txt$", full.names = TRUE)

# Map channel -> OLI band naming (based on the NWP-SAF band listing)
band_names <- c(
  ch01 = "B1",  # Coastal/Aerosol
  ch02 = "B2",  # Blue
  ch03 = "B3",  # Green
  ch04 = "B4",  # Red
  ch05 = "B5",  # NIR
  ch06 = "B6",  # SWIR1
  ch07 = "B7",  # SWIR2
  ch08 = "B8",  # Pan
  ch09 = "B9"   # Cirrus
)

srf_long <- bind_rows(lapply(srf_files, function(f) {
  ch <- str_extract(basename(f), "ch\\d+")
  ch <- paste0("ch", str_pad(str_extract(ch, "\\d+"), 2, pad="0"))
  bn <- band_names[[ch]]
  
  df <- read_nwpsaf_srf_txt(f) %>%
    mutate(band = bn,
           wl_r = round(wl_nm)) %>%  # integer nm for joining to your wl_r
    group_by(band, wl_r) %>%         # collapse if multiple points round to same nm
    summarise(response = mean(response, na.rm = TRUE), .groups = "drop")
  
  df
}))

# Wide SRF table (one row per integer wavelength)
srf_wide <- srf_long %>%
  pivot_wider(names_from = band, values_from = response, values_fill = 0) %>%
  rename(SR_WL = wl_r)

# srf_wide now has:
# SR_WL, B1, B2, ..., B9  (weights)
