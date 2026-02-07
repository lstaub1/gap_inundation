library(terra)

# Load your rasters
rast1 <- rast("F:/MASTERS/THESIS/2026_fresh/DTM/PatuxentTerrain.tif")
rast2 <- rast("F:/MASTERS/THESIS/MD_HECRAS/pax/IQ50/Velocity (Max).PaxTerrain.tif")

# Check resolution
res(rast1)  # Returns c(x_resolution, y_resolution)
res(rast2)

# Compare if resolutions are the same
identical(res(rast1), res(rast2))

# Check if grids align (same resolution AND same origin)
compareGeom(rast1, rast2, stopOnError = FALSE)
# Returns TRUE if they align perfectly, FALSE otherwise