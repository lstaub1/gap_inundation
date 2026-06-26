## ---- Flight Line Banding Diagnosis ---- ##

library(pacman)
p_load(lidR, terra, ggplot2, sf, tidyterra, patchwork)

# --- 1) Load data ----
las <- readLAS("D:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2013a.laz")
projection(las) <- "EPSG:3582"

diff_chm <- rast("D:/MASTERS/THESIS/2026_fresh/CHM/dem_Patuxent_diff.tif")

# --- 2) Extract points and convert to sf ----
pts <- as.data.frame(las@data[, c("X", "Y", "PointSourceID")])

# If the file is very large, sample it down for plotting speed
set.seed(42)
if (nrow(pts) > 2e6) {
  pts <- pts[sample(nrow(pts), 2e6), ]
}

pts_sf <- st_as_sf(pts, coords = c("X", "Y"), crs = 3582)

# --- 3) Plot PointSourceID map ----
p_flightlines <- ggplot(pts_sf) +
  geom_sf(aes(color = factor(PointSourceID)), size = 0.1, alpha = 0.4) +
  scale_color_discrete(name = "Flight Line\n(PointSourceID)") +
  labs(title = "2013 Flight Lines by PointSourceID") +
  theme_minimal() +
  theme(legend.position = "right")

# --- 4) Plot difference raster ----
p_diff <- ggplot() +
  geom_spatraster(data = diff_chm) +
  scale_fill_viridis_c(
    name = "Height\nChange (ft)",
    option = "plasma",
    na.value = "transparent",
    limits = c(-20, 20),   # adjust to your data range
    oob = scales::squish
  ) +
  labs(title = "CHM Difference (2018 - 2013)") +
  theme_minimal()

# --- 5) Side by side ----
p_flightlines + p_diff +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Flight Line Banding Diagnosis",
    subtitle = "If bands in the difference raster align with strip boundaries, banding is the cause"
  )

# Save it
ggsave("D:/MASTERS/THESIS/2026_fresh/CHM/flightline_diagnosis.png",
       width = 16, height = 8, dpi = 200)


# --- 6) Optional: check Z bias per flight line ----
# Normalize first if not done, then look at median Z per strip
z_by_strip <- as.data.frame(las@data[, c("PointSourceID", "Z")]) |>
  dplyr::group_by(PointSourceID) |>
  dplyr::summarise(
    median_Z = median(Z, na.rm = TRUE),
    mean_Z   = mean(Z, na.rm = TRUE),
    sd_Z     = sd(Z, na.rm = TRUE),
    n_pts    = dplyr::n()
  )

print(z_by_strip)

# If median_Z varies noticeably between PointSourceIDs (e.g. >0.5 ft),
# that's strong evidence of a strip-to-strip Z bias