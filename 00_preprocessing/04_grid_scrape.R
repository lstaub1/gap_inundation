##### ----- Grid Scraping ----- #####
#### Author: Leah E. Staub
#### Creation Date: 07/28/2025
#### Update Date: 07/28/2025
#### Purpose: This script takes rasterized HECRAS inundation model outputs and scrapes the grid information from the files. This information will be used when we create the DEMs so that both datasets match up smoothly.    

# Set up Environment ----
library(pacman)
p_load(mapview, sf, ggplot2, mapview, lidR, terra, tidyterra, fs, archive, tools, abind, RCSF)

#Load Datasets
#Seneca duration hrcras raster as first example
hecras<- rast("F:/MASTERS/THESIS/MD_HECRAS/sen/IQ50/Duration (hrs).sen_ft.tif")
hecter<- rast("F:/MASTERS/THESIS/TerrainFiles/SenecaTerrainmo/SenecaTerrain/SenecaTerrain.tif")

#Grab grid info
extent<- ext(hecras)
resolution <- res(hecras) #4ft by 4ft!

crs_info <- crs(hecras) #

#Grab grid info
extent<- ext(hecter)
resolution <- res(hecter) #1.8 by 1.8 ft? 

crs_info <- crs(hecter) 


#Load Datasets
#Seneca duration hrcras raster as first example
hecras<- rast("F:/MASTERS/THESIS/MD_HECRAS/sen/IQ50/Duration (hrs).sen_ft.tif")
hecter<- rast("F:/MASTERS/THESIS/2026_fresh/DTM/PatuxentTerrain.tif")

#Grab grid info
extent<- ext(hecras)
resolution <- res(hecras) #4ft by 4ft!

crs_info <- crs(hecras) #

#Grab grid info
extent<- ext(hecter)
resolution <- res(hecter) #1.8 by 1.8 ft? 

crs_info <- crs(hecter) 