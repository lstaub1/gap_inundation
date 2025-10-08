##### ----- Raw Point Cloud data Processing ----- #####
#### Author: Leah E. Staub
#### Creation Date: 07/06/2025
#### Update Date: 07/16/2025
#### Purpose: This script takes raw lidar point cloud data and classifies & filters it to ground, vegetation, and surface points. It does this by retiling the point cloud data for efficient processing.   

# Set up Environment ----
library(pacman)
p_load(mapview, sf, ggplot2, mapview, lidR, terra, tidyterra, fs, archive, tools, abind, RCSF)

# Named list of LAScatalogs
las_catalogs <- list(
  #Harf2013 = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ34/clipped_chunk_1.laz"),
  #Harf2020 = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Harf2020_reprojected/HARF2020_ft.laz"),
  #Balt2015a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ30/clipped_chunk_1.laz"),
  Balt2015b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ31/clipped_chunk_2.laz"),
  How2011a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ26/clipped_chunk_2.laz"),
  How2011b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ26/clipped_chunk_3.laz"),
  How2018a  = readLAScatalog("F:/MASTERS/THESIS/data/Clip/How_2018_BLK_1/clipped_chunk_3.laz"),
  How2018b  = readLAScatalog("F:/MASTERS/THESIS/data/Clip/How_2018_BLK_2/clipped_chunk_2.laz"),
  Mont2020a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Montgomery_2020_BLK2/clipped_chunk_4.laz"),
  Mont2020b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Montgomery_2020_BLK3/clipped_chunk_4.laz"),
  Mont2018a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Mont_2018_BLK2/clipped_chunk_4.laz"),
  Mont2018b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Mont_2018_BLK4/clipped_chunk_4.laz"),
  Mont2013a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ12/clipped_chunk_3.laz"),
  Mont2013b = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ12/clipped_chunk_4.laz")
)

#For some reason all of the laz files have NA as their crs. Let's assign the correct crs to them all. 
# Loop through list and assign EPSG
for (name in names(las_catalogs)) {
  
  # Assign EPSG:2248 
  projection(las_catalogs[[name]]) <- st_crs(2248)$proj4string
}

#check to see if it worked! 
lapply(las_catalogs, st_crs)




#### Potential New workflow -----
filter_tiles <- function(chunk, output_dir, buffer, base_name, remove_duplicates = TRUE) {
  las <- readLAS(chunk)
  
  if (is.null(las) || npoints(las) == 0) {
    message("⚠️ Empty tile found, skipping")
    return(NULL)
  }
  
  message("📦 Processing tile with ", npoints(las), " points")
  
  # Optionally remove exact duplicate points
  if (remove_duplicates) {
    las@data <- distinct(las@data)
  }
  
  # Classify ground
  las <- classify_ground(las, algorithm = csf())
  
  # Keep ground and surface (last return) points
  las <- las[las@data$Classification == 2 | las@data$ReturnNumber == las@data$NumberOfReturns, ]
  
  # Trim buffer — keep only the core area
  core <- ext(las)  # Get extent using terra
  core[1] <- core[1] + buffer  # xmin
  core[2] <- core[2] - buffer  # xmax
  core[3] <- core[3] + buffer  # ymin
  core[4] <- core[4] - buffer  # ymax
  
  # Clip the LAS object to the new extent
  las <- clip_rectangle(las, core[1], core[3], core[2], core[4])
  
  # Output filename based on chunk's core origin
  #tile_name <- paste0("tile_", las@header$X[1], "_", las@header$Y[1], ".laz")
  #out_file <- file.path(output_dir, tile_name)
  
  # Output filename based on the provided base name and tile coordinates
  tile_name <- paste0(base_name, "_tile_", round(core[1], 3), "_", round(core[3], 3), ".laz")
  out_file <- file.path(output_dir, tile_name)
  
  writeLAS(las, out_file)
  
  #writeLAS(las)
  message("✅ Processed tile written to: ", out_file)  
  return(out_file)
}

retile <- function(ctg, processed_dir, tile_size, buffer, crs_target, base_name) {
  if (!inherits(ctg, "LAScatalog")) stop("Input must be a LAScatalog object.")
  
  if (!is.null(crs_target)) {
    projection(ctg) <- crs_target
  }
  
  # Configure tiling
  opt_chunk_size(ctg) <- tile_size
  opt_chunk_buffer(ctg) <- buffer
  opt_output_files(ctg) <- ""  # Don’t write by default, we control writing
  opt_progress(ctg) <- TRUE
  opt_filter(ctg) <- "-drop_withheld"
  #opt_output_files(ctg) <- file.path(processed_dir, paste0("{filename}_{XLEFT}_{YBOTTOM}"))
  
  dir_create(processed_dir)
  
  catalog_apply(ctg, function(chunk, ...) {
    filter_tiles(chunk, processed_dir, buffer, base_name)
  })
  
  message("🎉 All tiles processed and saved in: ", processed_dir)
}

for (name in names(las_catalogs)) {
  ctg <- las_catalogs[[name]]  # Get the LAScatalog object
  processed_dir <- file.path("F:/MASTERS/THESIS/data/Processed", name)  # Create output directory
  
  # Call the retile function
  retile(
    ctg = ctg,
    processed_dir = processed_dir,
    tile_size = 1000,
    buffer = 30,
    crs_target = "EPSG:2248",
    base_name = name  # Use the name of the catalog as the base name
  )
}



##Need to rerun HARF2020....its extent is showing up in west virginia???

##Next steps: I need to combine the processed tiles based on study area. Then I need to remove duplicates.I also need to think about which points to keep in the areas that counties overlap. I should make sure I apply county/year to each point

#Update- I think I should just make CHMs with each of the county lidar datasets and stitch together later...?

