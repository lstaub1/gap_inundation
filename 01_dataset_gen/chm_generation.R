##### ----- CHM Creation----- #####
#### Author: Leah E. Staub
#### Creation Date: 04/25/2023
#### Update Date: 2/7/2026
#### Purpose: This script takes DTM and a raw point cloud and creates a chm.

library(pacman)
p_load(lidR, terra, ggplot2, sf)


#### PATUXENT 2013----
# 1) Read the point cloud and DTM
las <- readLAS("F:/MASTERS/THESIS/2026_fresh/Clipped_Lidar/Montgomery2013a.laz")   # or .las
dtm <- rast("F:/MASTERS/THESIS/2026_fresh/DTM/PatuxentTerrain.tif")             

# Spatial info for both datasets ----
# Ensure CRS matches
crs(dtm)
#Looks like the crs of the dtm is Maryland State Plane (SPCS) on NAD83(NSRS2007) in US survey feet, a Lambert Conformal Conic (2SP) projection. The matching EPSG is EPSG:3582 (“NAD83 (NSRS2007) / Maryland (ftUS)”). The DTM is in WTK format, which las does not accept.
crs(las)
#Looks like the crs of the las file is Maryland State Plane (SPCS) on NAD83(NSRS2007) using Lambert Conformal Conic (2SP) and US survey feet for the axes. Same as the DTM! I found this out by bringing the pt cloud into QGIS and....guessing! 

# Tell lidR the projection using EPSG (preferred for LAS headers)
projection(las) <- "EPSG:3582"

# DTM data Handling ----
#We know from the crs that our x and y values are in feet, this also was confirmed in QGIS. Let's see if the z value is in feet too.
# If values are ~100-500: likely feet
# If values are ~30-150: likely meters
summary(values(dtm)) #Max is 545, so looks like z value is in feet. 

# We know that 1.8 m grid size is too small for a chm. Let's try to make our CHM 5.4 ft res. 
# To do this, we need to make sure that our source raster (DTM) is able to scale up to 5.4 ft. 
# If it's number of pixels rows and columns are not able to scale up, we will need to add a row of no data to the raster. 
dtm_dims <- dim(dtm)
nrows <- dtm_dims[1]
ncols <- dtm_dims[2]

cat("Rows:", nrows, "\n")
cat("Cols:", ncols, "\n\n")

cat("Rows ÷ 3 =", nrows / 3, "(Remainder:", nrows %% 3, ")\n")
cat("Cols ÷ 3 =", ncols / 3, "(Remainder:", ncols %% 3, ")\n\n")

#If there is a remainder in the output we need to adjust the raster!
#If there is a remainder in rows, we need to pad the rows. Ditto for columns.

#If padding is needed, we need to extend raster
#In this case we just need to pad the columns.

#Check what we need
#target_rows <- ceiling(nrows / 3) * 3
target_cols <- ceiling(ncols / 3) * 3
target_cols #4059
cols_to_add <- target_cols - ncol(dtm)  #4059-4057 = 2 ~ we just need to add 2 columns
cols_to_add

#Specify extent to add on each side (left, right, bottom, top)
# To add 2 columns on the right: c(0, 2, 0, 0)
dtm_padded <- extend(dtm, c(0, cols_to_add, 0, 0))

#Check to make sure padding worked
cat("New dimensions:", dim(dtm_padded), "\n")
cat("Cols ÷ 3 =", ncol(dtm_padded) / 3, "\n")

# Point Cloud Data Handling ----
#Grab las point density
summary(las)
#density: 0.22 points/us-ft²

# Normalize heights above ground using the DTM
# This subtracts ground elevation from each point’s Z
las_n <- normalize_height(las, dtm)

#Check out height distribution of normalized data
summary(las_n@data$Z)
quantile(las_n$Z, seq(0, 1, 0.1), na.rm = TRUE)  # Deciles

# Count points above certain thresholds
sum(las_n$Z > 100, na.rm = TRUE) / length(las_n$Z)  # % above 100
sum(las_n$Z > 150, na.rm = TRUE) / length(las_n$Z)  # % above 150

# Plot histogram instead of density for better intuition
hist(las_n$Z, breaks = 100, main = "Height distribution", 
     xlab = "Height", xlim = c(0, 1035))  #set xlim to max value from summary

#Looks like the majority of sampled points fall below 200. So we need to remove some outliers

#Visualize normalized point cloud.
plot(las_n, color = "Z")

#Filter Method 1: Remove obvious outliers and below-ground points
las_n <- filter_poi(las_n, Z >= 0 & Z <= 200)  # adjust max height for your forest

#Filter Method 2: Statistical outlier removal (SOR)
# Removes isolated points that are far from their neighbors
las_clean <- classify_noise(las_n, sor(k = 15, m = 3))
las_clean <- filter_poi(las_clean, Classification != 18)
#Remove points that couldn't be normalized
las_clean <- filter_poi(las_clean, !is.na(Z))

#Recheck point height distribution

#Visually check for outliers
plot(las_n, color = "Z")

plot(las_clean, color = "Z")

summary(las_clean@data$Z)

#double check tallest points
#We know that x and y are in m, let's see what z is in.

# Check the actual distribution 
summary(las_clean$Z)
quantile(las_clean$Z, seq(0, 1, 0.1), na.rm = TRUE)  # Deciles

# Count points above certain thresholds
sum(las_clean$Z > 100, na.rm = TRUE) / length(las_clean$Z)  # % above 100
sum(las_clean$Z > 150, na.rm = TRUE) / length(las_clean$Z)  # % above 150

# Plot histogram instead of density for better intuition
hist(las_clean$Z, breaks = 100, main = "Height distribution", 
     xlab = "Height", xlim = c(0, 150))


## CHM method investigation ----
# P2R algorithm
chm_p2r <- rasterize_canopy(las_n, 
                            res = res(dtm)[1],  # Match DTM resolution
                            algorithm = p2r())

plot(chm_p2r)   

chm_p2r <- resample(chm_p2r, dtm, method = "near")

writeRaster(chm_p2r, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_p2r.tif", overwrite = TRUE)

# Pitfree algorithm
chm_pit <- rasterize_canopy(las_n, 
                            res = res(dtm)[1],  # Match DTM resolution
                            algorithm = pitfree())

plot(chm_pit_re)  

chm_pit_re <- resample(chm_pit, dtm, method = "bilinear")

writeRaster(chm_pit_re, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_pit.tif", overwrite = TRUE)


# dsmTIN algorithm
chm_tin <- rasterize_canopy(las_n, 
                            res = res(dtm)[1],
                            algorithm = dsmtin())


plot(chm_tin_re)  

chm_tin_re <- resample(chm_tin, dtm, method = "bilinear")

writeRaster(chm_tin_re, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_tin.tif", overwrite = TRUE)

# investigate chm algorithms
difference <- chm_pit_re - chm_tin_re
plot(difference, main = "Difference: Pitfree - TIN")

## Final algorithm choice ---- ###
#Going with pit free for now, modifying some of the settings to account for some stripping and small holes. 

# Thresholds in FEET (for height/Z values):
# 0, 6.6, 16.4, 32.8, 49.2, 65.6 feet
# which correspond to 0, 2, 5, 10, 15, 20 meters

chm <- rasterize_canopy(las_n, 
                        res = res(dtm)[1],
                        algorithm = pitfree(
                          thresholds = c(0, 6.6, 16.4, 32.8, 49.2, 65.6),
                          max_edge = c(0, 2.5)))

#Align grids
chm <- resample(chm, dtm, method = "bilinear")
chm[chm < 0] <- 0

writeRaster(chm, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_pf_t1.tif", overwrite = TRUE)

# Focal smoothing
chm_final <- focal(chm, w = 3, fun = mean, na.rm = TRUE)

writeRaster(chm_final, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_pf_t2.tif", overwrite = TRUE)

# 4) Post-process (fill tiny holes, light smoothing)
chm <- terra::focal(chm, w = matrix(1, 3, 3), fun = max, na.policy = "omit")
chm <- terra::subst(chm, NA, 0) # set NAs to 0 if desired

#evaluate settings
plot(chm)                 # raster view
plot(las, color = "Z")  # point cloud colored by height

# 5) Save
writeRaster(chm, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm.tif", overwrite = TRUE)








## CHM grid size investigation ----
# Create CHM at 5.4 ft resolution (3x the DTM's 1.8 ft)
chm_5.4ft <- rasterize_canopy(las_clean, 
                              res = 5.4,  # 5.4 feet resolution
                              algorithm = pitfree())

writeRaster(chm_5.4ft, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_5.4.tif", overwrite = TRUE)

chm_10.8ft <- rasterize_canopy(las_clean, 
                              res = 10.8,  # 10.8 feet resolution
                              algorithm = pitfree())

writeRaster(chm_10.8ft, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_10.8.tif", overwrite = TRUE)

chm_16.2ft <- rasterize_canopy(las_clean, 
                               res = 16.2,  # 10.8 feet resolution
                               algorithm = pitfree())

writeRaster(chm_16.2ft, "F:/MASTERS/THESIS/2026_fresh/CHM/Patuxent_2013_chm_16.2.tif", overwrite = TRUE)

# 3) Build the CHM (choose resolution and algorithm)
# Build CHM that matches DTM's grid by using DTM as the template
#Max_edge has to do with interpolation over areas with no points (gaps). To keep interpolation low, since we are looking for gaps, I've set this to 0.
chm   <- rasterize_canopy(
  las_n,
  res = dtm,  # <-- use the DTM raster as a template
  algorithm = pitfree(thresholds = c(0, 6, 16, 33, 49, 82, 115), max_edge = c(0, 10, 15, 20, 25, 30, 30))
)





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







#### OLD GUNK ----
# Create CHM, remembering that x and y and z are different units.
chm <- rasterize_canopy(las_clean, 
                        res = res(dtm_padded)[1],
                        algorithm = pitfree(
                          thresholds = c(0, 6.6, 16.4, 32.8, 49.2, 65.6),
                          max_edge = c(0, 2.5)))

plot(chm)