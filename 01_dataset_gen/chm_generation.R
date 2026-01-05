##### ----- CHM Creation----- #####
#### Author: Leah E. Staub
#### Creation Date: 04/25/2023
#### Update Date: 1/4/2026
#### Purpose: This script takes DTM and a raw point cloud and creates a chm.

library(pacman)
p_load(lidR, terra, ggplot2)


## PATUXENT ----
# 1) Read the point cloud and DTM
las <- readLAS("F:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2013a.laz")   # or .las
dtm <- rast("F:/MASTERS/THESIS/2026_fresh/DTM/PatuxentTerrain.tif")             

# Ensure CRS matches
crs(dtm)
#Looks like the crs of the dtm is Maryland State Plane (SPCS) on NAD83(NSRS2007) in US survey feet, a Lambert Conformal Conic (2SP) projection. The matching EPSG is EPSG:3582 (“NAD83 (NSRS2007) / Maryland (ftUS)”). The DTM is in WTK format, which las does not accept.
crs(las)
#Looks like the crs of the las file is Maryland State Plane (SPCS) on NAD83(NSRS2007) using Lambert Conformal Conic (2SP) and US survey feet for the axes. Same as the DTM!

# Tell lidR the projection using EPSG (preferred for LAS headers)
projection(las) <- "EPSG:3582"

# 2) Normalize heights above ground using the DTM
# This subtracts ground elevation from each point’s Z
las_n <- normalize_height(las, dtm)

#Check out height distribution of normalized data
summary(las_n@data$Z)
# Sample up to N points for plotting
set.seed(42)
N <- 5e6
z <- las_n$Z
z <- z[!is.na(z)]
z_sample <- if (length(z) > N) sample(z, N) else z

# Kernel Density plot
ggplot(data.frame(Z = z_sample), aes(Z)) +
  geom_density(fill = "orange", alpha = 0.3) +
  labs(title = "Density of heights", x = "Height above ground (m)", y = "Density") +
  theme_minimal()

#Looks like the majority of sampled points fall within 125. 

#Visualize normalized point cloud.
plot(las_n, color = "Z")

# (Optional) Remove obvious outliers and below-ground points
las_n <- filter_poi(las_n, Z >= 0 & Z <= 200)  # adjust max height for your forest

plot(las_n, color = "Z")

# 3) Build the CHM (choose resolution and algorithm)
# Build CHM that matches DTM's grid by using DTM as the template
chm   <- rasterize_canopy(
  las_n,
  res = dtm,  # <-- use the DTM raster as a template
  algorithm = pitfree(thresholds = c(0,2,5,10,15), max_edge = c(0,30))
)

#evaluate settings
plot(chm)                 # raster view
plot(las_n, color = "Z")  # point cloud colored by height


# Option B: point-to-raster with a height threshold (fast/simple)
chm_p2r <- rasterize_canopy(
  las_n,
  res = 1,
  algorithm = p2r(threshold = 0.2)  # ignore very low returns
)

# 4) Post-process (fill tiny holes, light smoothing)
chm_pf <- terra::focal(chm_pf, w = matrix(1, 3, 3), fun = max, na.policy = "omit")
chm_pf <- terra::subst(chm_pf, NA, 0) # set NAs to 0 if desired

# 5) Save
writeRaster(chm_pf, "outputs/CHM_1m_pitfree.tif", overwrite = TRUE)
writeRaster(chm_p2r, "outputs/CHM_1m_p2r.tif", overwrite = TRUE)





## NO DTM???? generate one!

library(lidR)
library(terra)

las <- readLAS("data/pointcloud.laz")

# (Optional) Classify ground if points aren’t already labeled
las <- classify_ground(las, csf())  # Cloth Simulation Filter

# Build a DTM (choose an algorithm; TIN is good when ground is well classified)
dtm <- grid_terrain(las, res = 1, algorithm = tin())

# Normalize and build CHM as above
las_n <- normalize_height(las, dtm)
chm <- rasterize_canopy(las_n, res = 1, algorithm = pitfree(c(0,2,5,10,15), c(0,30)))
writeRaster(chm, "outputs/CHM_1m.tif", overwrite = TRUE)

