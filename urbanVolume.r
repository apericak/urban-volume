# Urban Density from Lidar, DRIVE E:/
# Andrew Pericak, aap35@duke.edu

# This script uses FUSION (http://forsys.cfr.washington.edu/fusion/fusionlatest.html) to give a rough estimate of 
# "urban volume" for an input city. The script clips a list of raw lidar files to a polygon boundary (Step 1); removes
# obvious error in the raw data (Step 2); generates rough statistics using LAStools (Step 3); determines actual "ground"
# points, since sometimes pre-classified ground points can't differentiate building roofs from the ground (Step 4);
# creates a ground surface, aka a digital elevation model (Step 5); creates an "urban canopy" surface, aka a digital
# surface model (Step 6); determines volume between the DEM and DSM (Step 7); re-clips the data to the polygon extent
# (Step 8); and deletes intermediate data (Step 9). Note, in this current iteration, trees are not filtered out.

# The output products are three volume rasters, as described in the Fusion manual. "Potential volume" represents the
# maximum possible volume per raster cell; "surface volume" represents a best-guess true volume per raster cell; and 
# "surface volume ratio" is the ratio between the two. Note that Fusion produces other volume products; consult the 
# manual and modify the script if you want those volume products.

# NOTE: if you are using .laz files, you must copy the LASzip.dll from LASzip or LAStools to the Fusion install folder.
# Consult the FUSION manual for additional details about this.



### LIBRARIES
library(rgdal)
library(raster)


### INPUT VARIABLES

# Fusion, lastools, and data folders
fusionInstall = "E:/Fusion"  # Path to the Fusion folder
lastoolsInstall = "E:/LAStools/bin" # Path to bin folder of lastools
rawDataFolder = "E:/ENV857_project/data/RawLidar/BelowBP/AboveAvgDensity/HibbingMN" # Path to raw data folder (.las files)
processedDataFolder = "E:/ENV857_project/data/ProcessedLidar/BelowBP/AboveAvgDensity/HibbingMN" # Path to processed data folder

# Shapefile of urban area
cityShapeFolder = "E:/ENV857_project/data/TargetUACs/BelowBP/AboveAvgDensity" # Folder containing shapefile of city boundary
cityShape = "HibbingMN_proj.shp" # Shapefile polygon of city boundary
proj4 = "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" # PROJ-4 string for input shapefile; see spatialreference.org to get these

# Parameters used for Fusion tools; check lidar metadata and see descriptions below. Also refer to Fusion manual.
cellsize = 1 # Output cell size in same units as input data
xyunits = "M" # Use M for meters or F for feet
zunits = "M" # Use M for meters or F for feet
coordsys = 1 # Use 0 for unknown, 1 for UTM, 2 for state plane
zone = 15 # Coordinate system zone, if known; otherwise use 0
horizdatum = 2 # Use 0 for unknown, 1 for NAD27, 2 for NAD83
vertdatum = 2 # Use 0 for unknown, 1 for NGVD29, 2 for NAVD88, 3 for GRS80


### PROCESSING

# STEP 1: CLIP LAS FILES TO CITY BOUNDARY (note: shapefile MUST be in same projection as lidar data)
system2(paste(fusionInstall,"PolyClipData", sep="/"), 
        args = c(
          paste(cityShapeFolder, cityShape, sep="/"), 
          paste(rawDataFolder, "clipped.las", sep="/"), 
          paste(rawDataFolder, "fileIndex.txt", sep="/"))) # fileIndex.txt is a text list of raw lidar data files, in .laz or .las
# Note: if a lot of .laz/.las files and you need a .txt list:
# CMD: dir /s/b *.laz > fileIndex.txt

# STEP 2: REMOVE POINT ERRORS / OUTLYING DATA
system2(paste(fusionInstall,"FilterData", sep="/"),
        args = c("outlier", "5.0", "100", 
                 paste(rawDataFolder, "cleaned.las", sep="/"), 
                 paste(rawDataFolder, "clipped.las", sep="/")))

# STEP 3: GENERATE STATS FOR DATA; use these stats to determine an appropriate output cellsize (see above)
system2(paste(lastoolsInstall, "lasinfo.exe", sep="/"),
        args = c("-v", "-i", 
                 paste(rawDataFolder, "cleaned.las", sep="/"), "-o", 
                 paste(rawDataFolder, "cleaned_stats.csv", sep="/"), "-cd"))

# STEP 4: FILTER TRUE GROUND POINTS
system2("C:/Program Files/Fusion/GroundFilter",   # note: for some reason, doesn't run with external drive version...
        args = c(
          paste(rawDataFolder, "groundpts.las", sep="/"), "200",
          paste(rawDataFolder, "cleaned.las", sep="/")))

# STEP 5: CREATE GROUND SURFACE / DEM
system2(paste(fusionInstall, "TINSurfaceCreate", sep="/"),
        args = c(
          paste(processedDataFolder, "/", "DEM_", toString(cellsize), "m.dtm", sep=""),
          cellsize, xyunits, zunits, coordsys, zone, horizdatum, vertdatum,
          paste(rawDataFolder, "groundpts.las", sep="/")))
system2(paste(fusionInstall, "DTM2ASCII", sep="/"),
        args = c("/raster", 
                 paste(processedDataFolder, "/", "DEM_", toString(cellsize), "m.dtm", sep="")))

# STEP 6: CREATE 'URBAN CANOPY' SURFACE / DSM
# Non-normalized DSM
system2(paste(fusionInstall, "CanopyModel", sep="/"),
        args = c("/ascii",
                 paste("/align:", processedDataFolder, "/", "DEM_", toString(cellsize), "m.dtm", sep=""),
                 paste(processedDataFolder, "/", "DSM_", toString(cellsize), "m.dtm", sep=""),
                 cellsize, xyunits, zunits, coordsys, zone, horizdatum, vertdatum,
                 paste(rawDataFolder, "cleaned.las", sep="/")))

# STEP 7: DETERMINE VOLUME
system2(paste(fusionInstall, "GridSurfaceStats", sep="/"),
        args = c(
          paste("/ground:", processedDataFolder, "/", "DEM_", toString(cellsize), "m.dtm", sep=""), "/ascii",
          paste("/align:", processedDataFolder, "/", "DEM_", toString(cellsize), "m.dtm", sep=""),
          paste(processedDataFolder, "/", "DSM_", toString(cellsize), "m.dtm", sep=""),
          paste(processedDataFolder, "/", "volume", toString(cellsize), "m.dtm", sep=""), "1"))

# STEP 8: MASK THE VOLUME SURFACES TO URBAN EXTENT
# Read shapefile (again) and read three desired volume rasters
uac <- readOGR(dsn = cityShapeFolder, layer = substr(cityShape, 1, nchar(cityShape)-4), p4s = proj4)
pv <- raster(paste(processedDataFolder, "/", "volume", toString(cellsize), "m_potential_volume.asc", sep=""))
sv <- raster(paste(processedDataFolder, "/", "volume", toString(cellsize), "m_surface_volume.asc", sep=""))
svr <- raster(paste(processedDataFolder, "/", "volume", toString(cellsize), "m_surface_volume_ratio.asc", sep=""))
# Peform initial clip and then mask; save as tif
# Potential volume
pv.sub <- crop(pv, extent(uac))
pv.sub <- mask(pv.sub, uac)
writeRaster(pv.sub, paste(processedDataFolder, "clipped_potentialVolume.tif", sep="/"), overwrite = TRUE)
# Surface volume
sv.sub <- crop(sv, extent(uac))
sv.sub <- mask(sv.sub, uac)
writeRaster(sv.sub, paste(processedDataFolder, "clipped_surfaceVolume.tif", sep="/"), overwrite = TRUE)
# Surface volume ratio
svr.sub <- crop(svr, extent(uac))
svr.sub <- mask(svr.sub, uac)
writeRaster(svr.sub, paste(processedDataFolder, "clipped_surfaceVolumeRatio.tif", sep="/"), overwrite = TRUE)

# STEP 9: DELETE INTERMEDIATE FILES
file.remove(c(
  paste(rawDataFolder, "/", "clipped.las", sep=""),
  paste(processedDataFolder, "/", "volume", toString(cellsize), "m_max_height.asc", sep=""),
  paste(processedDataFolder, "/", "volume", toString(cellsize), "m_mean_height.asc", sep=""),
  paste(processedDataFolder, "/", "volume", toString(cellsize), "m_potential_volume.asc", sep=""),
  paste(processedDataFolder, "/", "volume", toString(cellsize), "m_stddev_height.asc", sep=""),
  paste(processedDataFolder, "/", "volume", toString(cellsize), "m_surface_area_ratio.asc", sep=""),
  paste(processedDataFolder, "/", "volume", toString(cellsize), "m_surface_volume.asc", sep=""),
  paste(processedDataFolder, "/", "volume", toString(cellsize), "m_surface_volume_ratio.asc", sep="")
))
