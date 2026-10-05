##### ----- Patuxent CHM Creation Workflow ----- #####
#### Author: Leah E. Staub
#### Creation Date: 04/25/2023
#### Update Date: 06/26/2026
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

# List present PointSourceIDs (with counts)
pointsource_summary <- function(las) {
  df <- as.data.frame(las@data[, "PointSourceID", drop = FALSE])
  df %>%
    count(PointSourceID, name = "n_pts") %>%
    arrange(desc(n_pts))
}

# Basic acquisition stats

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


# Cleaning workflow: remove negative, extreme outliers, then SOR
clean_las <- function(las_n, max_height = 200, sor_k = 15, sor_m = 3) {
  las_f <- filter_poi(las_n, Z >= 0 & Z <= max_height)
  las_f <- classify_noise(las_f, sor(k = sor_k, m = sor_m))
  las_f <- filter_poi(las_f, Classification != 18)             # remove noise points
  las_f <- filter_poi(las_f, !is.na(Z))                        # remove failed normalization
  las_f
}

# Flightline quality assessment on normalized heights
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

## ==== SECTION 4  Load DTM ==========================================

# Shared DTM for normalizing both years
# CRS: EPSG:3582 (Maryland State Plane NAD83 NSRS2007, US survey feet)
dtm <- rast("D:/MASTERS/THESIS/2026_fresh/DTM/PatuxentTerrain.tif")

cat("===== DTM Info =====\n")
print(crs(dtm))
print(ext(dtm))
print(summary(values(dtm)))

## =============================================================================
## SECTION 5  Prepare DTM
##     - Pad dimensions if necessary
##     - Verify CRS
##     - Verify units
## =============================================================================

# Ensure rows & cols divisible by 3 (for 3× resampling grids like 10.8 ft)
dtm_dims <- dim(dtm)
nrows <- dtm_dims[1]; ncols <- dtm_dims[2]
cat("DTM dimensions — Rows:", nrows, "| Cols:", ncols, "\n")
cat("Rows ÷ 3 remainder:", nrows %% 3, "\n")
cat("Cols ÷ 3 remainder:", ncols %% 3, "\n")

target_rows <- ceiling(nrows / 3) * 3
target_cols <- ceiling(ncols / 3) * 3
rows_to_add <- target_rows - nrows
cols_to_add <- target_cols - ncols
cat("Padding — Rows to add:", rows_to_add, "| Cols to add:", cols_to_add, "\n")

dtm_padded <- extend(dtm, c(cols_to_add, 0, 0, rows_to_add))
cat("Padded DTM dimensions:", dim(dtm_padded), "\n")
cat("Cols ÷ 3:", ncol(dtm_padded) / 3, "| Rows ÷ 3:", nrow(dtm_padded) / 3, "\n")

# CRS & units checks (informational)
stopifnot(grepl("US survey foot", terra::crs(dtm_padded)))
# Double check make sure it's correct
terra::crs(dtm_padded)
cat("DTM CRS verified as EPSG:3582 (US-ft). Z units expected to be feet.\n")

## =============================================================================
## SECTION 6  Normalize point clouds
##     las13_n
##     las18_n
## =============================================================================

cat("===== Normalizing 2013 =====\n")
las13_n <- normalize_height(las13, dtm_padded)

cat("===== Normalizing 2018 =====\n")
las18_n <- normalize_height(las18, dtm_padded)

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
las13_clean <- clean_las(las13_n, max_height = 200)
cat("2013 points remaining:", nrow(las13_clean@data), "\n")
print(summary(las13_clean@data$Z))

cat("\n===== Cleaning 2018 =====\n")
las18_clean <- clean_las(las18_n, max_height = 150)
cat("2018 points remaining:", nrow(las18_clean@data), "\n")
print(summary(las18_clean@data$Z))

## =============================================================================
## SECTION 8  Flightline Quality Assessment
##     - Compute ground-height statistics by PointSourceID
##     - Determine whether flightlines exhibit systematic vertical offsets
##     - If offsets exceed tolerance:
##           compute strip corrections
##           apply corrections
##           regenerate normalized cloud
##           (or apply correction directly to normalized heights)
## =============================================================================

cat("\n===== Flightline QA: 2013 =====\n")
qa13 <- flightline_quality_and_correction(
  las_raw = las13,
  las_norm_clean = las13_clean,
  good_strips = NULL,       # or c(30214, 30215) if known good reference strips
  tol_ft = 0.5,
  ground_window = c(-1, 2),
  apply_to = "normalized",  # default: adjust normalized heights
  dtm = dtm_padded
)
las13_qc <- qa13$corrected_norm
cat("2013 reference near-ground median:", round(qa13$reference, 3), "ft\n")
cat("2013: correction applied?", qa13$applied, "\n")
print(qa13$offsets)

# Re-clean after correction to ensure no negatives/outliers were introduced
las13_qc <- clean_las(las13_qc, max_height = 200)

cat("\n===== Flightline QA: 2018 =====\n")
qa18 <- flightline_quality_and_correction(
  las_raw = las18,
  las_norm_clean = las18_clean,
  good_strips = NULL,
  tol_ft = 0.5,
  ground_window = c(-1, 2),
  apply_to = "normalized",
  dtm = dtm_padded
)
las18_qc <- qa18$corrected_norm
cat("2018 reference near-ground median:", round(qa18$reference, 3), "ft\n")
cat("2018: correction applied?", qa18$applied, "\n")
print(qa18$offsets)

# Re-clean after correction
las18_qc <- clean_las(las18_qc, max_height = 150)

## =============================================================================
## SECTION 9  Build CHMs
## =============================================================================

# Pit-free thresholds approximate 0, 2, 5, 10, 15, 20 m in feet
pitfree_thresholds <- c(0, 6.6, 16.4, 32.8, 49.2, 65.6)

cat("\n===== Building 2013 CHM =====\n")
chm13_r <- rasterize_canopy(
  las13_qc,
  res = 10.8,
  algorithm = pitfree(
    thresholds = pitfree_thresholds,
    max_edge   = c(0, 2.5)
  )
)
chm13_r[chm13_r < 0] <- 0

# Convert to SpatRaster for terra I/O
chm13 <- rast(chm13_r)

writeRaster(
  chm13,
  "D:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_10.8ft.tif",
  overwrite = TRUE
)
cat("2013 CHM saved.\n")

cat("\n===== Building 2018 CHM =====\n")
chm18_r <- rasterize_canopy(
  las18_qc,
  res = 10.8,
  algorithm = pitfree(
    thresholds = pitfree_thresholds,
    max_edge   = c(0, 2.5)
  )
)
chm18_r[chm18_r < 0] <- 0

chm18 <- rast(chm18_r)

writeRaster(
  chm18,
  "D:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2018_chm_10.8ft.tif",
  overwrite = TRUE
)
cat("2018 CHM saved.\n")

# Quick visual check
plot(chm13, main = "2013 CHM (10.8 ft)")
plot(chm18, main = "2018 CHM (10.8 ft)")

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
  "D:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_diff_2018minus2013.tif",
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

## =============================================================================
## SECTION 11  Diagnostic figures
## =============================================================================

# 11a. Flightline vs. CHM-diff comparison (using QA-corrected 2013 cloud)
diag_path <- "D:/MASTERS/THESIS/2026_fresh/CHM/flightline_diagnosis.png"
plot_flightline_diagnosis(
  las = las13_qc,
  diff_raster = diff_chm,
  year_label = "2013 (QC)",
  out_path = diag_path
)

# 11b. Strip offset barplots (saved as PNGs)
offset_plot <- function(offsets_df, year_label, out_path) {
  p <- ggplot(offsets_df, aes(x = factor(PointSourceID), y = offset_ft)) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_col(fill = "#377eb8") +
    labs(
      title = paste("Strip Offsets (near-ground) —", year_label),
      x = "PointSourceID", y = "Offset (ft; ref - median_groundZ)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(out_path, p, width = 12, height = 6, dpi = 200)
  cat("Saved strip offset plot to:", out_path, "\n")
}

if (!is.null(qa13$offsets)) {
  offset_plot(qa13$offsets,
              "2013",
              "D:/MASTERS/THESIS/2026_fresh/CHM/strip_offsets_2013.png")
}
if (!is.null(qa18$offsets)) {
  offset_plot(qa18$offsets,
              "2018",
              "D:/MASTERS/THESIS/2026_fresh/CHM/strip_offsets_2018.png")
}

