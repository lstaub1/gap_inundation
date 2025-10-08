##### ----- SENECA CREEK CHM ----- #####
#### Author: Leah E. Staub
#### Creation Date: 04/25/2023
#### Update Date: 11/4/2024
#### Purpose: This script grabs the site buffer I created for Seneca Creek, converts from laz to las, creates a DSM, and explores CHM characteristics. 

# Load packages
library(pacman)
p_load(lidR, here, ggplot2, sf, raster, tidyverse, rgdal, ForestTools, ggpubr, ForestGapR, viridis)

sen_buff <- st_read("F:/MASTERS/THESIS/data/buffers/senvel90_buff500_proj.shp")

plot(sen_buff,
     axes=T)

#Montgomery county lidar
mont_lidar <- paste0("F:/MASTERS/THESIS/LIDAR/SEN_MONT2020")

mont_cat<- readLAScatalog(mont_lidar)

st_crs(mont_cat)

plot(mont_cat, mapview=TRUE, map.type="Esri.WorldStreetMap")

#Clip las catalog item to buffer shapefile
opt_output_files(mont_cat) <- paste0(here(),"/LASfiles/seneca_500m")
clipped_las<- clip_roi(mont_cat, sen_buff)

plot(clipped_las)

#Load Clipped Las file
sen_las<- readLAS(paste0(here(), "/LASfiles/seneca_500m.las"))

#Create DSMs (there are 3 different methods)
dsm_pf2<- rasterize_canopy(sen_las, res=1, algorithm=pitfree())

writeRaster(dsm_pf2, paste0(here(),"/DSMs/sen_dsm_pf3.tif"))

dsm_pt <- rasterize_canopy(sen_las, res=1, algorithm=p2r())

writeRaster(dsm_pt, paste0(here(),"/DSMs/p2rv2.tif"))

dsm_tin <- rasterize_canopy(sen_las, res=1, algorithm=dsmtin())

writeRaster(dsm_tin, paste0(here(),"/DSMs/tin1.tif"))

#Load DSM
DSM<- raster("F:/MASTERS/THESIS/Scripts/DSMs/sen_dsm_pf3_clip.tif")
st_crs(DSM)
crs(DSM)<- "EPSG:2893"

#Load DEM
DEM<- raster("E:/MASTERS/THESIS/DEMS/montgomery/2020/Countywide/mont_2020_m.ovr")
st_crs(DEM)
crs(DEM)<- "EPSG:2893"

#Load CHM
chm <- raster("F:/MASTERS/THESIS/CHM/sen500_clip")

###Looking at canopy tops
#Identify the tree tops
lin <- function(x){x * 0.05 + 0.6}
ttops <- vwf(CHM = chm, winFun = lin, minHeight = 2)
plot(chm)

# Create crown map
crowns <- mcws(treetops = ttops, CHM = chm, minHeight = 1.5, verbose = FALSE)
# Plot crowns
plot(chm, main="act-02 CHM")
plot(crowns, col = sample(rainbow(50), length(unique(crowns[])), replace = TRUE), legend = FALSE, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
#Canopy crowns map 
crownsPoly <- mcws(treetops = ttops, CHM = chm, format = "polygons", minHeight = 1.5, verbose = FALSE)
# Plot CHM
plot(chm, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
# Add crown outlines to the plot
plot(crownsPoly, border = "blue", lwd = 0.5, add = TRUE)
crown_stats<- as.data.frame(sp_summarise(crownsPoly, variables = c("crownArea", "height"))) 
crown_stats <- crown_stats %>%
  rownames_to_column('Stat')
crown_stats <-  pivot_wider(crown_stats, names_from = Stat, values_from = Value) 
crown_stats$Site <- i
orig_df<- rbind(orig_df, crown_stats)


#Looking at gaps!
#Let's convert our CHM to a different proj so it's in meters
#Load CHM
chm <- raster("F:/MASTERS/THESIS/CHM/senchm_proj")
crs(chm) <- "ESRI:102039"
st_crs(chm)

hi<- plot(chm, main="seneca")

# set height thresholds (e.g. 10 meters)
threshold <- 150
size <- c(1, 10^4) # m2

# Detecting forest gaps
gaps_duc <- getForestGaps(chm_layer = chm, threshold = threshold)

# Ploting gaps
hi<- plot(gaps_duc, col = "red", add = TRUE, main = "Forest Canopy Gap", legend = FALSE)
hi

data(ALS_CHM_CAU_2012)
# Plotting chm
plot(ALS_CHM_DUC, col = viridis(10), main = "ALS CHM")
grid()

plot(chm, col = viridis(10), main = "seneca")

# set height thresholds (e.g. 10 meters)
threshold <- 200
size <- c(30, 10^4) # m

# Detecting forest gaps
gaps_duc2 <- getForestGaps(chm_layer = chm, threshold = threshold, size = size)

writeRaster(gaps_duc2, "F:/MASTERS/THESIS/GAPS/sen_gaps200.tif")



plot(gaps_duc, col = "red", add = TRUE, main = "Forest Canopy Gap", legend = FALSE)




