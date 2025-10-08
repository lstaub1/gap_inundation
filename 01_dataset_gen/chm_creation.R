##### ----- SENECA CREEK CHM Creation----- #####
#### Author: Leah E. Staub
#### Creation Date: 04/25/2023
#### Update Date: 11/11/2024
#### Purpose: This script takes a DSM and DTM and creates a CHM from them. 

library(pacman)
p_load(lidR, terra)


#Load DSM
DSM<- rast("E:/MASTERS/THESIS/Scripts/DSMs/sen_dsm_pf3_clip.tif")
st_crs(DSM)
crs(DSM)<- "EPSG:2893"
DSM

#Load clipped DTM
DTM<- rast("E:/MASTERS/THESIS/TerrainFiles/SenecaTerrainmo/senterr_clip.vrt")
st_crs(DTM)
crs(DTM)<- "EPSG:2893"
DTM

plot(DSM)
plot(DTM)

#Do these rasters have different pixel sizes? 
#Let's resample to the smaller pixel size to be safe. 
#Hopefully they are on the same grid

######### Kendall's resample code is bath mode, let's modify it to be single.
#Raster to be resampled



datadir <- "F://USGS_FY23//DSWx//WV3//Maxar//Cottonwood_Lakes//ML_output" # Directory containing data to be resampled
setwd(datadir)
dir.create("//ML_resampled") # Create directory to put resampled data
outputdir <- paste0(datadir,"//ML_resampled") # set resampled directory as variable to use later

### Rasters to resample
tiffs <- list.files(pattern = ".tif$", full.names = TRUE, recursive = FALSE) # Extract file names of all tifs in data directory
dates <- substr(tiffs,3,10) # extract dates to use as file names later (assumes yyyymmdd at beginning of file name)

# Example Raster to resample to
resamp_ex <- rast('F://USGS_FY23//DSWx//WV3//Cottonwood_Lakes//HLS.S30.T14TMT.2023160T172859.v2.0.B02.TIF') # Location of example data to resample to
to_resamp_ex <- rast(tiffs[1]) # Select one scene of non-resampled data to extract scene extents
to_resamp_ex <- terra::project(to_resamp_ex, crs(resamp_ex)) # reproject non-resampled to match that of data we want to resample to
extent <- ext(to_resamp_ex) # extract extents of non-resampled data
resamp_ex <- crop(resamp_ex,extent) # crop example data to match extent of non-resampled data

### Loop through and resample data
## Expects format: 0 = Not_water, 1 = Water, 2 = Snow, 3 = Cloud
for (i in 1:length(tiffs)){
  tif <- rast(tiffs[i]) # import data as raster
  un <- unique(tif) # extract unique classification values
  classes <- un$bin_class # convert from dataframe to vector for looping later (may have to change column name from bin_class)
  for (j in 1:length(classes)){ # loop through each classification value
    k = classes[j] # extract classification value for this loop
    classified <- classify(tif, cbind(k,1), others = 0) # reclassify so values for this loop equal 1 and all others equal 0
    resampled <- resample(classified, resamp_ex, method="average") # perform resampling
    
    ###### save resampled tifs as percentage of each classification value (except for Not_Water)
    if (classes[j] == 1){
      setwd(outputdir)
      writeRaster(resampled,paste0(dates[i],"_resampled_pct_water.tif"),overwrite=TRUE)
      setwd(datadir)
    }else if (classes[j] == 2){
      setwd(outputdir)
      writeRaster(resampled,paste0(dates[i],"_resampled_pct_snow.tif"),overwrite=TRUE)
      setwd(datadir)
    }else if (classes[j] == 3){
      setwd(outputdir)
      writeRaster(resampled,paste0(dates[i],"_resampled_pct_cloud.tif"),overwrite=TRUE)
      setwd(datadir)
    }else{
      
    }
    
  }
}



