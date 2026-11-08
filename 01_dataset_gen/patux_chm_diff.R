##### ----- Patuxent CHM Creation Workflow ----- #####
#### Author: Leah E. Staub
#### Creation Date: 04/25/2023
#### Update Date: 10/05/2026
#### Purpose: For two years of lidar data (2013, 2018), this script:
####          - Loads and inspects point cloud spatial info & acquisition metadata
####          - Loads and prepares a shared DTM
####          - Normalizes point clouds (Z - ground) using the DTM
####          - Cleans outliers from normalized clouds
####          - Assesses flightline vertical consistency ----- eh ---and applies strip corrections if needed
####          - Creates a 10.8 ft Canopy Height Model (CHM) for each year
####          - Produces a CHM difference raster (2018 - 2013)
####          - Saves diagnostic figures

library(pacman)
p_load(lidR, terra, ggplot2, sf, tidyterra, patchwork, dplyr)

## ---- Helper functions ----------------------------------------------------------------------------

# List present PointSourceIDs ---
# counts points per flight strip (PointSourceID)
pointsource_summary <- function(las) {
  df <- as.data.frame(las@data[, "PointSourceID", drop = FALSE])
  df %>%
    count(PointSourceID, name = "n_pts") %>%
    arrange(desc(n_pts))
}

# Basic acquisition stats ---
# computes density, pulse rate, Z range, flight line count per dataset

basic_point_stats <- function(las, label) {
  # Extent-based area (map units are US survey feet in EPSG:3582)
  e <- ext(las)
  area_ft2 <- (e$xmax - e$xmin) * (e$ymax - e$ymin)
  
  # Point counts
  n_points <- nrow(las@data)
  
  # Returns & pulses (approximate pulse density from first returns if present)
  has_return <- "ReturnNumber" %in% names(las@data)
  n_first_returns <- if (has_return) sum(las@data$ReturnNumber == 1, na.rm = TRUE) else NA_real_
  
  # Flight lines
  has_psid <- "PointSourceID" %in% names(las@data)
  n_flight_lines <- if (has_psid) length(unique(las@data$PointSourceID)) else NA_integer_
  
  # Z stats
  z <- las@data$Z
  z_min  <- suppressWarnings(min(z, na.rm = TRUE))
  z_mean <- suppressWarnings(mean(z, na.rm = TRUE))
  z_max  <- suppressWarnings(max(z, na.rm = TRUE))
  
  data.frame(
    year                 = label,
    n_points             = n_points,
    area_ft2             = area_ft2,
    density_pts_per_ft2  = n_points / area_ft2,
    density_pulses_per_ft2 = if (!is.na(n_first_returns)) n_first_returns / area_ft2 else NA_real_,
    n_flight_lines       = n_flight_lines,
    z_min_ft             = z_min,
    z_mean_ft            = z_mean,
    z_max_ft             = z_max,
    stringsAsFactors     = FALSE
  )
}


# Cleaning workflow --- 
# Statistical Outlier Removal (SOR)
# SOR removes points far from any neighboring points.

clean_las <- function(las_n, min_height = -0.5, max_height = 200) {
  las_f <- filter_poi(las_n, Z >= min_height & Z <= max_height)
  las_f <- filter_poi(las_f, !is.na(Z))
  las_f
}


# Flightline quality assessment on normalized heights #MAYBE REMOVE THIS ENTIRELY?
# Computes near-ground residual by strip; if offsets exceed tolerance, apply corrections.
flightline_quality_and_correction <- function(
    las_raw, las_norm_clean,
    good_strips = NULL,        # optionally provide known-good reference strips (e.g., c(30214, 30215))
    tol_ft = 0.5,              # apply correction only if |offset| > tol_ft
    ground_window = c(-1, 2),  # normalized Z window treated as near-ground residuals
    apply_to = c("normalized", "raw"),  # where to apply corrections
    dtm = NULL                 # required only if apply_to == "raw"
) {
  apply_to <- match.arg(apply_to)
  stopifnot(is.numeric(ground_window), length(ground_window) == 2)
  
  # Per-strip near-ground residual stats (median of normalized Z in ground window)
  df <- as.data.frame(las_norm_clean@data[, c("PointSourceID", "Z")])
  near <- df %>%
    filter(Z >= ground_window[1], Z <= ground_window[2]) %>%
    group_by(PointSourceID) %>%
    summarise(
      median_groundZ = median(Z, na.rm = TRUE),
      mean_groundZ   = mean(Z, na.rm = TRUE),
      sd_groundZ     = sd(Z, na.rm = TRUE),
      n_ground       = n(),
      .groups        = "drop"
    )
  
  if (nrow(near) == 0) {
    warning("No near-ground points found within specified window; skipping QA.")
    return(list(
      corrected_norm = las_norm_clean,
      offsets = NULL,
      reference = NA_real_,
      applied = FALSE,
      stats = near
    ))
  }
  
  # Choose reference: user-specified good strips or top-2 by n_ground
  if (!is.null(good_strips)) {
    ref_val <- near %>%
      filter(PointSourceID %in% good_strips) %>%
      summarise(ref = median(median_groundZ, na.rm = TRUE)) %>%
      pull(ref)
  } else {
    ref_val <- near %>%
      arrange(desc(n_ground)) %>%
      slice_head(n = min(2, n())) %>%
      summarise(ref = median(median_groundZ, na.rm = TRUE)) %>%
      pull(ref)
  }
  
  offsets <- near %>%
    mutate(offset_ft = ref_val - median_groundZ) %>%
    arrange(PointSourceID)
  
  # Decide whether correction is needed
  need_correction <- any(abs(offsets$offset_ft) > tol_ft, na.rm = TRUE)
  
  if (!need_correction) {
    cat("Flightline QA: All strips within tolerance (", tol_ft, " ft). No correction applied.\n", sep = "")
    return(list(
      corrected_norm = las_norm_clean,
      offsets = offsets,
      reference = ref_val,
      applied = FALSE,
      stats = near
    ))
  }
  
  cat("Flightline QA: Applying strip corrections (tolerance =", tol_ft, "ft).\n")
  
  if (apply_to == "normalized") {
    # Apply directly to normalized heights
    las_corr <- las_norm_clean
    for (sid in offsets$PointSourceID) {
      ov <- offsets$offset_ft[offsets$PointSourceID == sid]
      if (abs(ov) > tol_ft) {
        idx <- las_corr@data$PointSourceID == sid
        las_corr@data$Z[idx] <- las_corr@data$Z[idx] + ov
        cat("  Strip", sid, ": adjusted normalized Z by", round(ov, 3), "ft\n")
      }
    }
    return(list(
      corrected_norm = las_corr,
      offsets = offsets,
      reference = ref_val,
      applied = TRUE,
      stats = near
    ))
  } else {
    # Apply to raw Z, then renormalize using provided DTM
    if (is.null(dtm))
      stop("DTM is required when apply_to = 'raw' to re-normalize after correction.")
    las_raw_corr <- las_raw
    for (sid in offsets$PointSourceID) {
      ov <- offsets$offset_ft[offsets$PointSourceID == sid]
      if (abs(ov) > tol_ft) {
        idx <- las_raw_corr@data$PointSourceID == sid
        las_raw_corr@data$Z[idx] <- las_raw_corr@data$Z[idx] + ov
        cat("  Strip", sid, ": adjusted RAW Z by", round(ov, 3), "ft\n")
      }
    }
    las_norm_corr <- normalize_height(las_raw_corr, dtm)
    return(list(
      corrected_norm = las_norm_corr,
      offsets = offsets,
      reference = ref_val,
      applied = TRUE,
      stats = near
    ))
  }
}

# Diagnostic plot comparing flightline map vs. CHM difference
plot_flightline_diagnosis <- function(las, diff_raster, year_label, out_path) {
  pts <- as.data.frame(las@data[, c("X", "Y", "PointSourceID")])
  set.seed(42)
  if (nrow(pts) > 2e6) pts <- pts[sample(nrow(pts), 2e6), ]
  pts_sf <- st_as_sf(pts, coords = c("X", "Y"), crs = 3582)
  
  p_strips <- ggplot(pts_sf) +
    geom_sf(aes(color = factor(PointSourceID)), size = 0.1, alpha = 0.4) +
    scale_color_discrete(name = "Flight Line\n(PointSourceID)") +
    labs(title = paste(year_label, "Flight Lines by PointSourceID")) +
    theme_minimal()
  
  p_diff <- ggplot() +
    geom_spatraster(data = diff_raster) +
    scale_fill_viridis_c(
      name = "Height\nChange (ft)",
      option = "plasma",
      na.value = "transparent",
      limits = c(-20, 20),
      oob = scales::squish
    ) +
    labs(title = "CHM Difference (2018 - 2013)") +
    theme_minimal()
  
  combined <- p_strips + p_diff +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Flight Line Banding Diagnosis",
      subtitle = "If bands align with strip boundaries, banding is the cause"
    )
  
  ggsave(out_path, combined, width = 16, height = 8, dpi = 200)
  cat("Saved diagnosis plot to:", out_path, "\n")
  combined
}

## ==== SECTION 1  Load LAS files (2013 & 2018) ==========================================

las13 <- readLAS("D:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2013a.laz")
las18 <- readLAS("D:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Howard2018a.laz")

# 2013 Z is in meters — convert to US survey feet to match DTM and 2018
# 1 m = 3.280833333 US survey ft
#las13@data$Z <- las13@data$Z * 3.280833333
#cat("2013 Z converted from meters to US survey feet.\n")

cat("===== 2013 LAS =====\n")
crs(las13)
ext(las13)

cat("===== 2018 LAS =====\n")
crs(las18)
ext(las18)

# Assign CRS - Maryland State Plane NAD83(NSRS2007), US survey feet
# EPSG:3582 confirmed via QGIS/metadata
projection(las13) <- "EPSG:3582"
projection(las18) <- "EPSG:3582"

## ==== SECTION 2  Inspect spatial information ===========================================
##     - CRS
##     - Extent
##     - PointSourceIDs present

cat("===== Spatial Info: 2013 =====\n")
print(crs(las13))
print(ext(las13))
ps13 <- pointsource_summary(las13)
cat("PointSourceIDs present (2013):", nrow(ps13), "strips\n")
print(head(ps13, 20))

cat("\n===== Spatial Info: 2018 =====\n")
print(crs(las18))
print(ext(las18))
ps18 <- pointsource_summary(las18)
cat("PointSourceIDs present (2018):", nrow(ps18), "strips\n")
print(head(ps18, 20))

## ==== SECTION 3  Inspect acquisition metadata ================================================
##     - Point density
##     - Number of flight lines
##     - Basic point statistics

cat("\n===== Acquisition Metadata =====\n")
meta13 <- basic_point_stats(las13, "2013")
meta18 <- basic_point_stats(las18, "2018")
print(meta13)
print(meta18)
cat("Number of flight lines (2013):", nrow(ps13), "\n")
cat("Number of flight lines (2018):", nrow(ps18), "\n")

## ==== SECTION 4  Create DEMs ==========================================

# Shared DTM for normalizing both years
# CRS: EPSG:3582 (Maryland State Plane NAD83 NSRS2007, US survey feet)
dtm <- rast("D:/MASTERS/THESIS/2026_fresh/DTM/PatuxentTerrain.tif")

cat("===== DTM Info =====\n")
print(crs(dtm))
print(ext(dtm))
print(summary(values(dtm)))

## =============================================================================
## SECTION 4b  Build year-specific DEMs from raw point clouds
## =============================================================================

# Use ground points only for terrain interpolation
# TIN (triangulated irregular network) is the standard for lidar DEMs
dem13 <- rasterize_terrain(las13, res = 10.8, algorithm = tin())
dem18 <- rasterize_terrain(las18, res = 10.8, algorithm = tin())

writeRaster(dem13, "D:/MASTERS/THESIS/2026_fresh/Data/Patuxent_dem13_10.8ft.tif", overwrite = TRUE)
writeRaster(dem18, "D:/MASTERS/THESIS/2026_fresh/Data/Patuxent_dem18_10.8ft.tif", overwrite = TRUE)

cat("DEM extents:\n")
print(ext(dem13))
print(ext(dem18))
cat("DEM value ranges:\n")
print(global(dem13, c("min", "mean", "max"), na.rm = TRUE))
print(global(dem18, c("min", "mean", "max"), na.rm = TRUE))


## =============================================================================
## SECTION 6  Normalize point clouds
##     las13_n
##     las18_n
## =============================================================================

cat("===== Normalizing 2013 =====\n")
las13_n <- normalize_height(las13, dem13)

cat("\n===== Acquisition Metadata =====\n")
meta13 <- basic_point_stats(las13_n, "2013")
print(meta13)


writeLAS(las13_n, "D:/MASTERS/THESIS/2026_fresh/Data/las13_n.las")

cat("===== Normalizing 2018 =====\n")
las18_n <- normalize_height(las18, dem18)

cat("\n===== Acquisition Metadata =====\n")
meta18 <- basic_point_stats(las18_n, "2018")
print(meta18)

writeLAS(las18_n, "D:/MASTERS/THESIS/2026_fresh/Data/las18_n.las")

# Quick summaries
cat("\n--- 2013 normalized Z summary ---\n")
print(summary(las13_n@data$Z))
print(quantile(las13_n$Z, seq(0, 1, 0.1), na.rm = TRUE))

cat("\n--- 2018 normalized Z summary ---\n")
print(summary(las18_n@data$Z))
print(quantile(las18_n$Z, seq(0, 1, 0.1), na.rm = TRUE))

## =============================================================================
## SECTION 7  Clean normalized point clouds
##     - Remove negative heights
##     - Remove extreme outliers
##     - Statistical outlier removal
## =============================================================================

cat("\n===== Cleaning 2013 =====\n")
las13_clean <- clean_las(las13_n)
cat("2013 points remaining:", nrow(las13_clean@data), "\n")
print(summary(las13_clean@data$Z))

cat("\n===== Acquisition Metadata =====\n")
meta13 <- basic_point_stats(las13_clean, "2013")
print(meta13)

writeLAS(las13_clean,"D:/MASTERS/THESIS/2026_fresh/Data/las13_clean.las")


cat("\n===== Cleaning 2018 =====\n")
las18_clean <- clean_las(las18_n)
cat("2018 points remaining:", nrow(las18_clean@data), "\n")
print(summary(las18_clean@data$Z))

writeLAS(las18_clean,"D:/MASTERS/THESIS/2026_fresh/Data/las18_clean.las")

# Quick summaries
cat("\n--- 2013 cleaned Z summary ---\n")
print(summary(las13_clean@data$Z))
print(quantile(las13_clean$Z, seq(0, 1, 0.1), na.rm = TRUE))

cat("\n--- 2018 cleaned Z summary ---\n")
print(summary(las18_clean@data$Z))
print(quantile(las18_clean$Z, seq(0, 1, 0.1), na.rm = TRUE))


## =============================================================================
## Checking removed points
## =============================================================================
#Spatially mapped the removed points after las cleaning
# Get row indices of points remaining after cleaning
# Match on XYZ since there's no persistent point ID in lidR
kept <- paste(las13_clean@data$X, las13_clean@data$Y, las13_clean@data$Z)
all  <- paste(las13_n@data$X,     las13_n@data$Y,     las13_n@data$Z)

removed_idx <- which(!all %in% kept)
cat("Points removed:", length(removed_idx), "\n")

df_removed <- as.data.frame(las13_n@data[removed_idx, c("X", "Y", "Z")])

# Sample if large
set.seed(42)
if (nrow(df_removed) > 200000) df_removed <- df_removed[sample(nrow(df_removed), 200000), ]

ggsave(
  "D:/MASTERS/THESIS/2026_fresh/Data/removedpts13_hist.png",
  ggplot(df_removed, aes(x = X, y = Y, color = Z)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_viridis_c(name = "Height (ft)", option = "plasma") +
    labs(
      title = "Spatially Removed Points (las13_n → las13_clean)",
      subtitle = paste0("n = ", nrow(df_removed), " removed points shown"),
      x = "Easting (ft)", y = "Northing (ft)"
    ) +
    coord_equal() +
    theme_minimal()
)


kept <- paste(las18_clean@data$X, las18_clean@data$Y, las18_clean@data$Z)
all  <- paste(las18_n@data$X,     las18_n@data$Y,     las18_n@data$Z)

removed_idx <- which(!all %in% kept)
cat("Points removed:", length(removed_idx), "\n")

df_removed <- as.data.frame(las18_n@data[removed_idx, c("X", "Y", "Z")])

# Sample if large
set.seed(42)
if (nrow(df_removed) > 200000) df_removed <- df_removed[sample(nrow(df_removed), 200000), ]

ggsave(
  "D:/MASTERS/THESIS/2026_fresh/Data/removedpts18_hist.png",
  ggplot(df_removed, aes(x = X, y = Y, color = Z)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_viridis_c(name = "Height (ft)", option = "plasma") +
    labs(
      title = "Spatially Removed Points (las18_n → las18_clean)",
      subtitle = paste0("n = ", nrow(df_removed), " removed points shown"),
      x = "Easting (ft)", y = "Northing (ft)"
    ) +
    coord_equal() +
    theme_minimal()
)

# Compare raw Z ranges
summary(las13@data$Z)
summary(las18@data$Z)


# If 2013 Z range is ~3.28x larger than 2018, it's a units issue
# If 2013 Z range is shifted by ~100-115 ft, it's likely a datum issue

#looks good!

##=============================================================================
## SECTION 8  Pulse Density Raster
##    - Create pulse density raster
##    - Identify strips with poor pulse density
##    - Plot flightlines on top of density raster
##    - Pulse density raster can be used later on as a mask for the chm differences
##==============================================================================

# Pulse density raster (first returns only = one per pulse)
density13 <- rasterize_density(las13_clean, res = 10.8)
density18 <- rasterize_density(las18_clean, res = 10.8)

writeRaster(density13, "D:/MASTERS/THESIS/2026_fresh/Data/density13_10.8ft.tif", overwrite = TRUE)
writeRaster(density18, "D:/MASTERS/THESIS/2026_fresh/Data/density18_10.8ft.tif", overwrite = TRUE)

#understand distribution of densities
vals13 <- values(density13, na.rm = TRUE)
summary(vals13)
ggsave(
  "D:/MASTERS/THESIS/2026_fresh/Data/density13_hist.png",
  ggplot(data.frame(density = vals13), aes(x = density)) +
    geom_histogram(bins = 50, fill = "#377eb8", color = "white") +
    labs(title = "2013 Pulse Density Distribution", x = "pulses/ft²", y = "Count") +
    theme_minimal(),
  width = 8, height = 5, dpi = 150
)

vals18 <- values(density18, na.rm = TRUE)
summary(vals18)
ggsave(
  "D:/MASTERS/THESIS/2026_fresh/Data/density18_hist.png",
  ggplot(data.frame(density = vals18), aes(x = density)) +
    geom_histogram(bins = 50, fill = "#377eb8", color = "white") +
    labs(title = "2018 Pulse Density Distribution", x = "pulses/ft²", y = "Count") +
    theme_minimal(),
  width = 8, height = 5, dpi = 150
)

#understand density distribution based on flight lines
# --- Summarize mean pulse density per strip ---
density13_r <- rast(density13)  # ensure it's a SpatRaster

strip_density_summary <- function(las_clean, density_rast, label) {
  pts <- as.data.frame(las_clean@data[, c("X", "Y", "PointSourceID")])
  pts_sf <- st_as_sf(pts, coords = c("X", "Y"), crs = 3582)
  
  # Extract density value at each point location
  pts_sf$pixel_density <- terra::extract(density_rast, pts_sf)[, 2]
  
  pts_sf %>%
    st_drop_geometry() %>%
    group_by(PointSourceID) %>%
    summarise(
      mean_density = mean(pixel_density, na.rm = TRUE),
      min_density  = min(pixel_density, na.rm = TRUE),
      n_pts        = n(),
      .groups = "drop"
    ) %>%
    arrange(mean_density) %>%
    mutate(year = label)
}

strip13 <- strip_density_summary(las13_clean, density13, "2013")
strip18 <- strip_density_summary(las18_clean, density18,   "2018")

print(strip13)
print(strip18)

#all of the minimums are 0.00857, this must be the equivalent of 0. Let's get a better idea of what the floor looks like
thresholds <- c(0.00857, 0.02, 0.05, 0.10)
for (t in thresholds) {
  cat(sprintf("Below %.5f — 2013: %.1f%%  |  2018: %.1f%%\n",
      t,
      100 * mean(vals13 < t, na.rm = TRUE),
      100 * mean(vals18 < t, na.rm = TRUE)))
}

## =============================================================================
## SECTION 9  Build CHMs
## =============================================================================

# Pit-free thresholds approximate 0, 2, 5, 10, 15, 20 m in feet
pitfree_thresholds <- c(0, 6.6, 16.4, 32.8, 49.2, 65.6)

cat("\n===== Building 2013 CHM =====\n")
chm13_r <- rasterize_canopy(
  las13_clean,
  res = 10.8,
  algorithm = pitfree(
    thresholds = pitfree_thresholds,
    max_edge   = c(0, 2.5)
  )
)

# Apply pulse density mask to 2013 CHM
# Cells below threshold are treated as no-data (gaps in coverage)
density_threshold <- 0  # adjust based on your distribution
density13_r <- resample(density13, chm13_r, method = "near")
mask13 <- ifel(density13_r == 0, NA, 1)
chm13 <- mask(chm13_r, mask13)


writeRaster(
  chm13,
  "D:/MASTERS/THESIS/2026_fresh/CHM_new/Patuxent_2013_chm_10.8ft.tif",
  overwrite = TRUE
)


cat("\n===== Building 2018 CHM =====\n")
chm18_r <- rasterize_canopy(
  las18_clean,
  res = 10.8,
  algorithm = pitfree(
    thresholds = pitfree_thresholds,
    max_edge   = c(0, 2.5)
  )
)

# Apply pulse density mask to 2013 CHM
# Cells below threshold are treated as no-data (gaps in coverage)
density_threshold <- 0  # adjust based on your distribution
density18_r <- resample(density18, chm18_r, method = "near")
mask18 <- ifel(density18_r == 0, NA, 1)
chm18 <- mask(chm18_r, mask18)


writeRaster(
  chm18,
  "D:/MASTERS/THESIS/2026_fresh/CHM_new/Patuxent_2018_chm_10.8ft.tif",
  overwrite = TRUE
)

## =============================================================================
## SECTION 10  Difference CHM
## =============================================================================

# Align grids before differencing
chm18_aligned <- resample(chm18, chm13, method = "bilinear")

# Difference where both have data
valid_mask <- !is.na(chm13) & !is.na(chm18_aligned)
diff_chm <- mask(chm18_aligned - chm13, valid_mask)

plot(
  diff_chm,
  main = "CHM Difference: 2018 - 2013 (ft)",
  col  = colorRampPalette(c("purple", "white", "yellow"))(100)
)

writeRaster(
  diff_chm,
  "D:/MASTERS/THESIS/2026_fresh/CHM_new/Patuxent_diff_2018minus2013.tif",
  overwrite = TRUE
)
cat("Difference raster saved.\n")

# Quick stats
cat("\n===== Difference Raster Summary =====\n")
print(global(diff_chm, c("min", "mean", "max"), na.rm = TRUE))
vals <- values(diff_chm)
cat("% cells with loss (< -2 ft):",
    round(100 * sum(vals < -2, na.rm = TRUE) / sum(!is.na(vals)), 1), "%\n")
cat("% cells with gain (>  2 ft): ",
    round(100 * sum(vals >  2, na.rm = TRUE) / sum(!is.na(vals)), 1), "%\n")

