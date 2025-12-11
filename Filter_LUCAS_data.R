################################################################################
# Script:  Filter_LUCAS_data.R
# Purpose: 
#   1) Load multi-year LUCAS data from CSV
#   2) Remove unused columns
#   3) Summarise how many observations have revisit ∈ {4, 5}
#   4) Restrict data to observations with revisit ∈ {4, 5}
#   5) Keep only point_ids where all lc1_perc values are exactly "> 75 %"
#   6) Exclude point_ids whose letter_group changes over time
#   7) Summarise number and share of points per LUCAS Level-1 class
#
# Input (needs to be adapted individually:  
#   /mnt/eo/EU_unmixing/data/LUCAS/LUCAS_multiyear.csv
#
# Output objects (in R session):
#   - perc_45:        share of observations with revisit ∈ {4, 5} in full dataset
#   - dt_45:          subset with revisit ∈ {4, 5}
#   - dt_filtered:    dt_45 restricted to point_ids with lc1_perc always "> 75 %"
#   - dt_stable:      dt_filtered restricted to point_ids with stable letter_group
#   - points_per_lg:  table of counts and percentages per letter_group
################################################################################

library(data.table)

################################################################################
# 1) Load CSV and remove unused columns
################################################################################

dt <- fread(
  "/mnt/eo/EU_unmixing/data/LUCAS/LUCAS_multiyear.csv",
  nThread      = parallel::detectCores(),
  showProgress = TRUE
)

# Column indices to remove (integer positions in the original dt)
cols_to_remove <- c(
  1:5,
  10:36,
  43:45,
  47:54,
  65,
  67:68,
  88:94,
  97:107,
  113:116
)

setDT(dt)  # ensure data.table

# Drop unwanted columns by index
dt[, (cols_to_remove) := NULL]

################################################################################
# 2) Share of observations with revisit ∈ {4, 5} in the full dataset
################################################################################

n_all <- nrow(dt)
n_45  <- dt[revisit %in% c(4, 5), .N]
perc_45 <- 100 * n_45 / n_all

message(sprintf(
  "Share of observations with revisit = 4 or 5: %.2f%% (%d out of %d)",
  perc_45, n_45, n_all
))

################################################################################
# 3) Restrict to observations with revisit ∈ {4, 5}
################################################################################

dt_45 <- dt[revisit %in% c(4, 5)]

################################################################################
# 4) Keep only point_ids where ALL lc1_perc values are exactly '> 75 %'
#    (across the full time series in dt_45)
################################################################################

# For each point_id: TRUE if every lc1_perc is '> 75 %' (no other level, no NA)
good_ids <- dt_45[
  , .(
    all_gt75 = all(trimws(lc1_perc) == "> 75 %")
  ),
  by = point_id
][all_gt75 == TRUE, point_id]

# Subset dt_45 to these "high tree-cover" point_ids
dt_filtered <- dt_45[
  point_id %in% good_ids
]

################################################################################
# 5) Exclude point_ids with changing letter_group over time (based on dt_filtered)
################################################################################

# For each point_id, compute:
#   - n_rows:  number of rows
#   - n_years: number of distinct years
#   - n_lg:    number of distinct letter_groups
check_ts <- dt_filtered[
  , .(
    n_rows   = .N,
    n_years  = uniqueN(year),
    n_lg     = uniqueN(letter_group)
  ),
  by = point_id
]

# Keep only point_ids with a constant letter_group (n_lg == 1)
stable_ids <- check_ts[n_lg == 1, point_id]

# Subset dt_filtered to point_ids with stable letter_group
dt_stable <- dt_filtered[point_id %in% stable_ids]

################################################################################
# 6) Summarise number and percentage of points per letter_group
################################################################################

# Treat each point_id as one point
points_unique <- unique(dt_stable[, .(point_id, letter_group)])

# Count points per letter_group
points_per_lg <- points_unique[
  , .(n_points = .N),
  by = letter_group
][order(letter_group)]

# Add percentage of points per letter_group
points_per_lg[
  , perc_points := 100 * n_points / sum(n_points)
]

# Optional: decode letter_group to LUCAS Level-1 class labels
lucas_labels <- c(
  A = "Artificial land",
  B = "Cropland",
  C = "Woodland",
  D = "Shrubland",
  E = "Grassland",
  F = "Bare land",
  G = "Water",
  H = "Wetlands"
)

points_per_lg[
  , lucas_class := lucas_labels[letter_group]
]

# Final summary table:
#   letter_group, lucas_class, n_points, perc_points
points_per_lg

# write
# Write as CSV
fwrite(
  dt_stable,
  file = "/mnt/eo/EU_unmixing/data/LUCAS/LUCAS_multiyear_filtered.csv"
)
