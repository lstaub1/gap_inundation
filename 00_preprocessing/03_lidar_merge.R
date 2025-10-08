##### ----- Filtered Point Cloud to DEM, DSM, and CHM ----- #####
#### Author: Leah E. Staub
#### Creation Date: 07/16/2025
#### Update Date: 07/16/2025
#### Purpose: This script takes filtered lidar point cloud data and creates DEM, DSM, and CHM rasters. The lidar data extents are first mapped to better understand overlap.    

# Set up Environment ----
library(pacman)
p_load(mapview, sf, ggplot2, mapview, lidR, terra, tidyterra, fs, archive, tools, abind, RCSF)

# Named list of LAScatalogs
las_catalogs <- list(
  Harf2013 = readLAScatalog("F:/MASTERS/THESIS/data/Processed/Harf2013/Harf2013_tile_1466411.8_668030.01.laz"),
  #Harf2020 = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Harf2020_reprojected/HARF2020_ft.laz"),
  Balt2015a = readLAScatalog("F:/MASTERS/THESIS/data/Processed/Balt2015a/Balt2015a_tile_1466411.8_666000.laz")
  # Balt2015b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ31/clipped_chunk_2.laz"),
  # How2011a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ26/clipped_chunk_2.laz"),
  # How2011b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ26/clipped_chunk_3.laz"),
  # How2018a  = readLAScatalog("F:/MASTERS/THESIS/data/Clip/How_2018_BLK_1/clipped_chunk_3.laz"),
  # How2018b  = readLAScatalog("F:/MASTERS/THESIS/data/Clip/How_2018_BLK_2/clipped_chunk_2.laz"),
  # Mont2020a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Montgomery_2020_BLK2/clipped_chunk_4.laz"),
  # Mont2020b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Montgomery_2020_BLK3/clipped_chunk_4.laz"),
  # Mont2018a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Mont_2018_BLK2/clipped_chunk_4.laz"),
  # Mont2018b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Mont_2018_BLK4/clipped_chunk_4.laz"),
  # Mont2013a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ12/clipped_chunk_3.laz"),
  # Mont2013b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ12/clipped_chunk_4.laz")
)
# Loop through list and assign EPSG

for (name in names(las_catalogs)) {
  
  # Assign EPSG:2248 
  projection(las_catalogs[[name]]) <- st_crs(2248)$proj4string
}

#check to see if it worked! 
lapply(las_catalogs, st_crs)



#Visualize footprints for all downloaded laz data to ensure correct files were downloaded.
#Load study areas shapefile
extents<- read_sf("F:/MASTERS/THESIS/data/extents/all_extents.shp")
# Combine the footprints list into one sf object and add a column that indicates which zip folder each group of laz files came from.
#Harford county is a different crs, lets make sure everything is the same crs
target_crs <- st_crs(extents)

all_footprints_named <- lapply(names(all_footprints), function(name) {
  sf_obj <- all_footprints[[name]]
  
  # Reproject if CRS doesn't match target
  if (!is.null(sf_obj) && st_crs(sf_obj) != target_crs) {
    sf_obj <- st_transform(sf_obj, target_crs)
  }
  
  sf_obj$source_zip <- name  # Add zip name as a new column
  return(sf_obj)
})

combined_footprints <- do.call(rbind, all_footprints_named)

#save
st_write(combined_footprints, "F:/MASTERS/THESIS/data/raw_lidar/all_lidar_footprints.shp", append=FALSE)

#read in footprints


# Plot with your study area polygons (`extents`)
mapview(extents, col.regions = "red", alpha.regions = 0.5) + 
  mapview(combined_footprints, 
          color = "lightblue", 
          layer.name = "LAS Catalog",
          zcol = "source_zip",
          alpha.regions = 0.3
  )

mapview(blk31_footprints, col.regions = "red", alpha.regions = 0.5)
















#Function for getting laz footprint
footprint <- function(name, catalog) {
  bbox <- st_as_sfc(st_bbox(catalog), crs = st_crs(catalog))
  sf <- st_sf(geometry = bbox)
  sf$source <- name
  return(sf)
}

# Apply function to each named item
footprints_list <- Map(footprint, names(las_catalogs), las_catalogs)

# Combine all into one sf object
footprints <- do.call(rbind, footprints_list)


# Split your sf object into a list of footprints by `source`
footprints_by_source <- split(footprints, footprints$source)

#Load study area polygons
extents<- read_sf("F:/MASTERS/THESIS/data/extents/all_extents.shp")


# Create a named list of mapview objects (one per source)
map_layers <- lapply(names(footprints_by_source), function(name) {
  mapview(footprints_by_source[[name]],
          col.regions = "blue",
          alpha.regions = 0.4,
          layer.name = name)
})

# Combine all into one map with toggleable layers
final_map <- mapview(extents, 
                     col.regions = "red", 
                     alpha.regions = 0.3, 
                     layer.name = "Study Area") + 
  Reduce(`+`, map_layers)

# View the final map
final_map
