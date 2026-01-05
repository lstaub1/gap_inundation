##### ----- Raw Point Cloud data Processing ----- #####
#### Author: Leah E. Staub
#### Creation Date: 07/06/2025
#### Update Date: 10/16/2025
#### Purpose: This script takes raw lidar point cloud data that is clipped to study area extents and ensures it is in the proper projection. 

library(pacman)
p_load(lidR)

#Need to fix HARF2020 projection before loading this list

# Named list of LAScatalogs
las_catalogs <- list(
  Harf2013 = readLAScatalog("F:/LEAH/Clipped_Lidar/Harford2013.laz"),
  Harf2020 = readLAScatalog("F:/LEAH/Clipped_Lidar/Harford2020.laz"),
  Balt2015a = readLAScatalog("F:/LEAH/Clipped_Lidar/Baltimore2015a.laz"),
  Balt2015b = readLAScatalog("F:/LEAH/Clipped_Lidar/Baltimore2015b.laz"),
  How2011a = readLAScatalog("F:/LEAH/Clipped_Lidar/Howard2011a.laz"),
  How2011b = readLAScatalog("F:/LEAH/Clipped_Lidar/Howard2011b.laz"),
  How2018a  = readLAScatalog("F:/LEAH/Clipped_Lidar/Howard2018a.laz"),
  How2018b  = readLAScatalog("F:/LEAH/Clipped_Lidar/Howard2018b.laz"),
  Mont2020a = readLAScatalog("F:/LEAH/Clipped_Lidar/Montgomery2020a.laz"),
  Mont2020b = readLAScatalog("F:/LEAH/Clipped_Lidar/Montgomery2020b.laz"),
  Mont2018a = readLAScatalog("F:/LEAH/Clipped_Lidar/Montgomery2018a.laz"),
  Mont2018b = readLAScatalog("F:/LEAH/Clipped_Lidar/Montgomery2018b.laz"),
  Mont2013a = readLAScatalog("F:/LEAH/Clipped_Lidar/Montgomery2013a.laz"),
  Mont2013b = readLAScatalog("F:/LEAH/Clipped_Lidar/Montgomery2013b.laz")
)

#Load counties for 1 study site
las_catalogs <- list(
  
)


#Test with just one file
file<- "Montgomery2020a"

# Path to your LAZ file
Mont2020a = readLAS("F:/LEAH/Clipped_Lidar/Montgomery2020a.laz")
#Took approx. 3 min to load

# Path to your LAZ file
Mont2018b = readLAS("F:/LEAH/Clipped_Lidar/Montgomery2018b.laz")


# Create new attribute column in data itself
test<- add_attribute(Mont2020a, file, "Source")

#Update header
test<- add_lasattribute(Mont2020a, file, "Source", "County and Year")
print(header(las)) 


# Write back out with extra byte fields (LAZ 1.4 supports this)
writeLAS(test, "F:/LEAH/Data/test.laz")

test2 = readLAS("F:/LEAH/Data/test.laz")

