##### ----- Seneca CHM Creation Workflow ----- #####
#### Author: Leah E. Staub
#### Creation Date: 04/25/2023
#### Update Date: 11/08/26
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

las13 <- readLAS("D:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2013b.laz")
las18a <- readLAS("D:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2018a.laz")
las18b <- readLAS("D:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2018b.laz")
las20a <- readLAS("D:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2020a.laz")
las20b <- readLAS("D:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2020b.laz")


cat("===== 2013 LAS =====\n")
crs(las13)
ext(las13)

cat("===== 2018 LAS =====\n")
crs(las18a)
ext(las18a)

crs(las18b)
ext(las18b)

cat("===== 2018 LAS =====\n")
crs(las20a)
ext(las20a)

crs(las20b)
ext(las20b)

#Looks like 2013 and 2018 don't have crs assigned but 2020 does. 

# Assign CRS - Maryland State Plane NAD83(NSRS2007), US survey feet
# EPSG:3582 confirmed via QGIS/metadata
projection(las13) <- "EPSG:3582"
projection(las18a) <- "EPSG:3582"
projection(las18b) <- "EPSG:3582"

# Reproject 2020 to match the other years (from EPSG:2893 to EPSG:3582)
las20a <- st_transform(las20a, crs = 3582)
las20b <- st_transform(las20b, crs = 3582)

# Verify
crs(las20a)
crs(las20b)

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
print(crs(las18a))
print(ext(las18a))
ps18a <- pointsource_summary(las18a)
cat("PointSourceIDs present (2018):", nrow(ps18a), "strips\n")
print(head(ps18a, 20))

print(crs(las18b))
print(ext(las18b))
ps18b <- pointsource_summary(las18b)
cat("PointSourceIDs present (2018):", nrow(ps18b), "strips\n")
print(head(ps18b, 20))

cat("\n===== Spatial Info: 2020 =====\n")
print(crs(las20a))
print(ext(las20a))
ps20a <- pointsource_summary(las20a)
cat("PointSourceIDs present (2020):", nrow(ps20a), "strips\n")
print(head(ps120a, 20))

print(crs(las20b))
print(ext(las20b))
ps20b <- pointsource_summary(las20b)
cat("PointSourceIDs present (2020):", nrow(ps20b), "strips\n")
print(head(ps20b, 20))

## ==== SECTION 3  Inspect acquisition metadata and combine ================================================
##     - Point density
##     - Number of flight lines
##     - Basic point statistics

cat("\n===== Acquisition Metadata =====\n")
meta13 <- basic_point_stats(las13, "2013")
meta18a <- basic_point_stats(las18a, "2018a")
meta18b <- basic_point_stats(las18b, "2018b")
meta20a <- basic_point_stats(las20a, "2020a")
meta20b <- basic_point_stats(las20b, "2020b")
print(meta13)
print(meta18a)
print(meta18b)
print(meta20a)
print(meta20b)
cat("Number of flight lines (2013):", nrow(ps13), "\n")
cat("Number of flight lines (2018a):", nrow(ps18a), "\n")
cat("Number of flight lines (2018b):", nrow(ps18b), "\n")
cat("Number of flight lines (2020a):", nrow(ps20a), "\n")
cat("Number of flight lines (2020b):", nrow(ps20b), "\n")

# Stitch 2018 and 2020 before doing anything else
las18 <- rbind(las18a, las18b)
las20 <- rbind(las20a, las20b)


## =============================================================================
## SECTION 4  Build year-specific DEMs from raw point clouds
## =============================================================================

#No deduplication of overlapping points, do we need this? Will need to check output! ********

# Then your existing pipeline works unchanged for all three years
dem13 <- rasterize_terrain(las13, res = 10.8, algorithm = tin())
dem18 <- rasterize_terrain(las18, res = 10.8, algorithm = tin())
dem20 <- rasterize_terrain(las20, res = 10.8, algorithm = tin())


writeRaster(dem13, "D:/MASTERS/THESIS/2026_fresh/Data/Seneca_dem13_10.8ft.tif", overwrite = TRUE)
writeRaster(dem18, "D:/MASTERS/THESIS/2026_fresh/Data/Seneca_dem18_10.8ft.tif", overwrite = TRUE)
writeRaster(dem20, "D:/MASTERS/THESIS/2026_fresh/Data/Seneca_dem20_10.8ft.tif", overwrite = TRUE)

cat("DEM extents:\n")
print(ext(dem13))
print(ext(dem18))
print(ext(dem20))
cat("DEM value ranges:\n")
print(global(dem13, c("min", "mean", "max"), na.rm = TRUE))
print(global(dem18, c("min", "mean", "max"), na.rm = TRUE))
print(global(dem20, c("min", "mean", "max"), na.rm = TRUE))

## =============================================================================
## SECTION 6  Normalize point clouds
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

cat("===== Normalizing 2020 =====\n")
las20_n <- normalize_height(las20, dem20)

cat("\n===== Acquisition Metadata =====\n")
meta20 <- basic_point_stats(las20_n, "2018")
print(meta20)

writeLAS(las20_n, "D:/MASTERS/THESIS/2026_fresh/Data/las20_n.las")

# Quick summaries
cat("\n--- 2013 normalized Z summary ---\n")
print(summary(las13_n@data$Z))
print(quantile(las13_n$Z, seq(0, 1, 0.1), na.rm = TRUE))

cat("\n--- 2018 normalized Z summary ---\n")
print(summary(las18_n@data$Z))
print(quantile(las18_n$Z, seq(0, 1, 0.1), na.rm = TRUE))

cat("\n--- 2020 normalized Z summary ---\n")
print(summary(las20_n@data$Z))
print(quantile(las20_n$Z, seq(0, 1, 0.1), na.rm = TRUE))



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

cat("\n===== Cleaning 2020 =====\n")
las20_clean <- clean_las(las20_n)
cat("2020 points remaining:", nrow(las20_clean@data), "\n")
print(summary(las20_clean@data$Z))

writeLAS(las20_clean,"D:/MASTERS/THESIS/2026_fresh/Data/las20_clean.las")

# Quick summaries
cat("\n--- 2013 cleaned Z summary ---\n")
print(summary(las13_clean@data$Z))
print(quantile(las13_clean$Z, seq(0, 1, 0.1), na.rm = TRUE))

cat("\n--- 2018 cleaned Z summary ---\n")
print(summary(las18_clean@data$Z))
print(quantile(las18_clean$Z, seq(0, 1, 0.1), na.rm = TRUE))

cat("\n--- 2020 cleaned Z summary ---\n")
print(summary(las20_clean@data$Z))
print(quantile(las20_clean$Z, seq(0, 1, 0.1), na.rm = TRUE))


## =============================================================================
## Checking removed points
## =============================================================================
#Spatially mapped the removed points after las cleaning
# Get row indices of points remaining after cleaning
# Match on XYZ since there's no persistent point ID in lidR

# kept <- paste(las13_clean@data$X, las13_clean@data$Y, las13_clean@data$Z)
# all  <- paste(las13_n@data$X,     las13_n@data$Y,     las13_n@data$Z)

# removed_idx <- which(!all %in% kept)
# cat("Points removed:", length(removed_idx), "\n")

# df_removed <- as.data.frame(las13_n@data[removed_idx, c("X", "Y", "Z")])

# # Sample if large
# set.seed(42)
# if (nrow(df_removed) > 200000) df_removed <- df_removed[sample(nrow(df_removed), 200000), ]

# ggsave(
#   "D:/MASTERS/THESIS/2026_fresh/Data/removedpts13_hist.png",
#   ggplot(df_removed, aes(x = X, y = Y, color = Z)) +
#     geom_point(size = 0.3, alpha = 0.6) +
#     scale_color_viridis_c(name = "Height (ft)", option = "plasma") +
#     labs(
#       title = "Spatially Removed Points (las13_n → las13_clean)",
#       subtitle = paste0("n = ", nrow(df_removed), " removed points shown"),
#       x = "Easting (ft)", y = "Northing (ft)"
#     ) +
#     coord_equal() +
#     theme_minimal()
# )


# kept <- paste(las18_clean@data$X, las18_clean@data$Y, las18_clean@data$Z)
# all  <- paste(las18_n@data$X,     las18_n@data$Y,     las18_n@data$Z)

# removed_idx <- which(!all %in% kept)
# cat("Points removed:", length(removed_idx), "\n")

# df_removed <- as.data.frame(las18_n@data[removed_idx, c("X", "Y", "Z")])

# # Sample if large
# set.seed(42)
# if (nrow(df_removed) > 200000) df_removed <- df_removed[sample(nrow(df_removed), 200000), ]

# ggsave(
#   "D:/MASTERS/THESIS/2026_fresh/Data/removedpts18_hist.png",
#   ggplot(df_removed, aes(x = X, y = Y, color = Z)) +
#     geom_point(size = 0.3, alpha = 0.6) +
#     scale_color_viridis_c(name = "Height (ft)", option = "plasma") +
#     labs(
#       title = "Spatially Removed Points (las18_n → las18_clean)",
#       subtitle = paste0("n = ", nrow(df_removed), " removed points shown"),
#       x = "Easting (ft)", y = "Northing (ft)"
#     ) +
#     coord_equal() +
#     theme_minimal()
# )

# # Compare raw Z ranges
# summary(las13@data$Z)
# summary(las18@data$Z)


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
density20 <- rasterize_density(las20_clean, res = 10.8)

writeRaster(density13, "D:/MASTERS/THESIS/2026_fresh/Data/Seneca_density13_10.8ft.tif", overwrite = TRUE)
writeRaster(density18, "D:/MASTERS/THESIS/2026_fresh/Data/Seneca_density18_10.8ft.tif", overwrite = TRUE)
writeRaster(density20, "D:/MASTERS/THESIS/2026_fresh/Data/Seneca_density20_10.8ft.tif", overwrite = TRUE)

#understand distribution of densities
# vals13 <- values(density13, na.rm = TRUE)
# summary(vals13)
# ggsave(
#   "D:/MASTERS/THESIS/2026_fresh/Data/density13_hist.png",
#   ggplot(data.frame(density = vals13), aes(x = density)) +
#     geom_histogram(bins = 50, fill = "#377eb8", color = "white") +
#     labs(title = "2013 Pulse Density Distribution", x = "pulses/ft²", y = "Count") +
#     theme_minimal(),
#   width = 8, height = 5, dpi = 150
# )

# vals18 <- values(density18, na.rm = TRUE)
# summary(vals18)
# ggsave(
#   "D:/MASTERS/THESIS/2026_fresh/Data/density18_hist.png",
#   ggplot(data.frame(density = vals18), aes(x = density)) +
#     geom_histogram(bins = 50, fill = "#377eb8", color = "white") +
#     labs(title = "2018 Pulse Density Distribution", x = "pulses/ft²", y = "Count") +
#     theme_minimal(),
#   width = 8, height = 5, dpi = 150
# )

# #understand density distribution based on flight lines
# # --- Summarize mean pulse density per strip ---
# density13_r <- rast(density13)  # ensure it's a SpatRaster

# strip_density_summary <- function(las_clean, density_rast, label) {
#   pts <- as.data.frame(las_clean@data[, c("X", "Y", "PointSourceID")])
#   pts_sf <- st_as_sf(pts, coords = c("X", "Y"), crs = 3582)
  
#   # Extract density value at each point location
#   pts_sf$pixel_density <- terra::extract(density_rast, pts_sf)[, 2]
  
#   pts_sf %>%
#     st_drop_geometry() %>%
#     group_by(PointSourceID) %>%
#     summarise(
#       mean_density = mean(pixel_density, na.rm = TRUE),
#       min_density  = min(pixel_density, na.rm = TRUE),
#       n_pts        = n(),
#       .groups = "drop"
#     ) %>%
#     arrange(mean_density) %>%
#     mutate(year = label)
# }

# strip13 <- strip_density_summary(las13_clean, density13, "2013")
# strip18 <- strip_density_summary(las18_clean, density18,   "2018")

# print(strip13)
# print(strip18)

# #all of the minimums are 0.00857, this must be the equivalent of 0. Let's get a better idea of what the floor looks like
# thresholds <- c(0.00857, 0.02, 0.05, 0.10)
# for (t in thresholds) {
#   cat(sprintf("Below %.5f — 2013: %.1f%%  |  2018: %.1f%%\n",
#       t,
#       100 * mean(vals13 < t, na.rm = TRUE),
#       100 * mean(vals18 < t, na.rm = TRUE)))
# }

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
  "D:/MASTERS/THESIS/2026_fresh/CHM_new/Seneca_2013_chm_10.8ft.tif",
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
  "D:/MASTERS/THESIS/2026_fresh/CHM_new/Seneca_2018_chm_10.8ft.tif",
  overwrite = TRUE
)

cat("\n===== Building 2020 CHM =====\n")
chm20_r <- rasterize_canopy(
  las20_clean,
  res = 10.8,
  algorithm = pitfree(
    thresholds = pitfree_thresholds,
    max_edge   = c(0, 2.5)
  )
)

# Apply pulse density mask to 2013 CHM
# Cells below threshold are treated as no-data (gaps in coverage)
density_threshold <- 0  # adjust based on your distribution
density20_r <- resample(density20, chm20_r, method = "near")
mask20 <- ifel(density20_r == 0, NA, 1)
chm20 <- mask(chm20_r, mask20)


writeRaster(
  chm20,
  "D:/MASTERS/THESIS/2026_fresh/CHM_new/Seneca_2020_chm_10.8ft.tif",
  overwrite = TRUE
)

## =============================================================================
## SECTION 10  Difference CHM
## =============================================================================


### 2020 - 2013 Difference ---------------------------------
# Align grids before differencing
chm20_aligned <- resample(chm20, chm13, method = "bilinear")

# Difference where both have data
valid_mask <- !is.na(chm13) & !is.na(chm20_aligned)
diff_chm <- mask(chm20_aligned - chm13, valid_mask)

plot(
  diff_chm,
  main = "CHM Difference: 2020 - 2013 (ft)",
  col  = colorRampPalette(c("purple", "white", "yellow"))(100)
)

writeRaster(
  diff_chm,
  "D:/MASTERS/THESIS/2026_fresh/CHM_new/Seneca_diff_2020minus2013.tif",
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

### 2018 - 2013 Difference ---------------------------------

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
  "D:/MASTERS/THESIS/2026_fresh/CHM_new/Seneca_diff_2018minus2013.tif",
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
