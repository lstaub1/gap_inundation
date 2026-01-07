##### ----- CHM Creation----- #####
#### Author: Leah E. Staub
#### Creation Date: 04/25/2023
#### Update Date: 1/4/2026
#### Purpose: This script takes DTM and a raw point cloud and creates a chm.

library(pacman)
p_load(lidR, terra, ggplot2)


## PATUXENT 2013----
# 1) Read the point cloud and DTM
las <- readLAS("F:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2013a.laz")   # or .las
dtm <- rast("F:/MASTERS/THESIS/2026_fresh/DTM/PatuxentTerrain.tif")             

#Grab las point density
summary(las)
#density: 0.22 points/us-ft²

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

#Visually check for outliers
plot(las_n, color = "Z")

#double check tallest points
#We know that x and y are in ft, let's make sure z is in ft too.
#If 99.9th percentile ≈ 35, your Z is in meters.
#If 99.9th percentile ≈ 125, your Z is in feet.
quantile(las_n$Z, c(0.95, 0.99, 0.999), na.rm = TRUE)

# 3) Build the CHM (choose resolution and algorithm)
# Build CHM that matches DTM's grid by using DTM as the template
#Max_edge has to do with interpolation over areas with no points (gaps). To keep interpolation low, since we are looking for gaps, I've set this to 0.
chm   <- rasterize_canopy(
  las_n,
  res = dtm,  # <-- use the DTM raster as a template
  algorithm = pitfree(thresholds = c(0, 6, 16, 33, 49, 82, 115), max_edge = c(0, 10, 15, 20, 25, 30, 30))
)

# 4) Post-process (fill tiny holes, light smoothing)
chm <- terra::focal(chm, w = matrix(1, 3, 3), fun = max, na.policy = "omit")
chm <- terra::subst(chm, NA, 0) # set NAs to 0 if desired

#evaluate settings
plot(chm)                 # raster view
plot(las, color = "Z")  # point cloud colored by height

# 5) Save
writeRaster(chm, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm.tif", overwrite = TRUE)







## PATUXENT 2018----
# 1) Read the point cloud and DTM
las <- readLAS("F:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Howard2018a.laz")   # or .las
dtm <- rast("F:/MASTERS/THESIS/2026_fresh/DTM/PatuxentTerrain.tif")             

plot(dtm)

#Grab las point density
summary(las)
#density: 0.46 points/us-ft²

# Ensure CRS matches
crs(dtm)
#Looks like the crs of the dtm is Maryland State Plane (SPCS) on NAD83(NSRS2007) in US survey feet, a Lambert Conformal Conic (2SP) projection. The matching EPSG is EPSG:3582 (“NAD83 (NSRS2007) / Maryland (ftUS)”). The DTM is in WTK format, which las does not accept.
crs(las)
#Looks like the las file does not have a crs. After checking the metadata, it looks like it has the same projection as the DTM.

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
  labs(title = "Density of heights", x = "Height above ground (ft)", y = "Density") +
  theme_minimal()

#Looks like the majority of sampled points fall within 125. 

#Visualize normalized point cloud.
plot(las_n, color = "Z")

# (Optional) Remove obvious outliers and below-ground points
las_n <- filter_poi(las_n, Z >= 0 & Z <= 150)  # adjust max height for your forest

#Visually check for outliers
plot(las_n, color = "Z")

#double check tallest points
#We know that x and y are in ft, let's make sure z is in ft too.
#If 99.9th percentile ≈ 35, your Z is in meters.
#If 99.9th percentile ≈ 125, your Z is in feet.
quantile(las_n$Z, c(0.95, 0.99, 0.999), na.rm = TRUE)

# 3) Build the CHM (choose resolution and algorithm)
# Build CHM that matches DTM's grid by using DTM as the template
#thresholds are breakdowns of tree heights. 
#Max_edge has to do with interpolation over areas with no points (gaps). This is dependent on point density. Since this years las file has higher point density, these can be lower. 
chm   <- rasterize_canopy(
  las_n,
  res = dtm,  # <-- use the DTM raster as a template
  algorithm = pitfree(thresholds = c(0, 6, 16, 33, 49, 65, 82, 115, 130), max_edge = c(0, 8, 12, 14, 16, 18, 20, 22, 24, 24))
)

# 4) Post-process (fill tiny holes, light smoothing)
#chm <- terra::focal(chm, w = matrix(1, 3, 3), fun = max, na.policy = "omit")
#chm <- terra::subst(chm, NA, 0) # set NAs to 0 if desired

#evaluate settings
plot(chm)                 # raster view
plot(las, color = "Z")  # point cloud colored by height

# 5) Save
writeRaster(chm, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2018_chm2.tif", overwrite = TRUE)

quantile(values(chm), c(.95, .99, .999), na.rm = TRUE)

## PATUXENT DIFFERENCE

diff_chm <- chm_y2 - chm_y1
diff_chm <- mask(diff_chm, !is.na(chm_y1) & !is.na(chm_y2))



## NO DTM???? generate one! ----

library(lidR)
library(terra)

las <- readLAS("F:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Howard2018a.laz")

# (Optional) Classify ground if points aren’t already labeled
las2 <- classify_ground(las, csf())  # Cloth Simulation Filter

# Build a DTM (choose an algorithm; TIN is good when ground is well classified)
dtm <- grid_terrain(las2, res = 1, algorithm = tin())

# Normalize and build CHM as above
las_n <- normalize_height(las, dtm)
chm <- rasterize_canopy(las_n, res = 1, algorithm = pitfree(c(0,2,5,10,15), c(0,30)))
writeRaster(chm, "outputs/CHM_1m.tif", overwrite = TRUE)

