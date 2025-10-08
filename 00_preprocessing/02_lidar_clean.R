##### ----- Raw Point Cloud data Processing ----- #####
#### Author: Leah E. Staub
#### Creation Date: 07/06/2025
#### Update Date: 08/09/2025
#### Purpose: This script takes raw lidar point cloud data that is clipped to study area extents and ensures it is in the proper projection. 

library(pacman)
p_load(mapview, sf, ggplot2, mapview, lidR, terra, tidyterra, fS)

# Named list of LAScatalogs
las_catalogs <- list(
  Harf2013 = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ34/clipped_chunk_1.laz"),
  Harf2020 = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Harford_2020_BLK1/clipped_chunk_1.laz"),
  #Harf2020 = readLAScatalog("F:/MASTERS/THESIS/data/Clip/Harf2020_reprojected/HARF2020_ft.laz"),
  Balt2015a = readLAScatalog("F:/MASTERS/THESIS/data/Clip/LAZ30/clipped_chunk_1.laz"),
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


