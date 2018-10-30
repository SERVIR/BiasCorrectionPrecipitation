# Script: config.R
# Purpose: Provide basic input parameters (input & output file locations, spatial reference system, geographic extent, time extent)
# Author: Francisco Delgado (francisco.delgadoolivares@nasa.gov), SERVIR Science Coordination Office
# Date: September 7, 2018, revised October 30, 2018
# Based on code prepared by Begum Rushi (begumrabeya.rushi@nasa.gov), Regional Associate, HKH Region, SERVIR, as of August 4, 2018

# ------------ Geospatial Extent and CRS ------------

# Spatial Coordinate Reference System
crs_definition = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
# Ouptut file extent
raster_output_extent <- c(34.75, 36.0, -0.5, 0.1)

# ------------ Time Period for Analysis ------------
start_year <- 2015
end_year <- 2016
start_month = 1
end_month = 12
start_day = 1
end_day = 31 # Has to be consistent with end_month

# ----------------- Base directory -----------------
# NOTE: - assign the full path to your base data directory
#       - use / instead of \ to separate directories (including Windows machines)
#       - finish with /
#       - If data directories are not under a common path, use full paths for the variables below
base_directory = "c:/precipdata/biascorrection/"
#
# For the sample dataset, the directory structure is as follows:
# base_directory/
#    - CHIRPS/
#    - CHIRPS_corrected/
#    - observed/
#    - chirps_obs_daily_comparison/
#    - persiann/
#    - persiann_corrected_linear/
#    - persiann_corrected_quantile/

# ------------ For Adjusted_CHIRPS.R ---------------------------

# Location of the Observed rainfall files (INPUT)
# Yearly rainfall files follow the pattern: STATION_NAME/YYYY/YYYY_daily.csv
observed_rainfall_dir <- paste(base_directory, "observed/", sep = "")

# Location for daily comparisson files (OUTPUT)
# Files contain the following columns (with column titles): chirps_prcp, obs_prcp
# Files named following the pattern: YYYY/YYYY.mm.dd_chirps_obs_prcp.csv
comparison_dir <- paste(base_directory, "chirps_obs_daily_comparison/", sep = "")

# Location of the ground stations file (INPUT)
# File contains the following columns (with column titles): station_name, Lat, Lon
ground_stations_file <- paste(base_directory, "Stations.csv", sep = "")

# ---------- For linear_method.R & quantile_method.R  ---------------

# Location of the CHIRPS files (INPUT)
# Files should be organized in subdirectories by year, following the pattern: YYYY/YYYY.mm.dd.tif
chirps_dir <- paste(base_directory, "chirps/", sep = "")

# Location for corrected CHIRPS files (OUTPUT for adjusted_CHIRPS.R, INPUT for linear_method.R and quantile_method.R if use_corrected_chirps set to TRUE)
# Files will be written to subdirectories by year, following the pattern: YYYY/YYYY.mm.dd.tif
corrected_chirps_dir <- paste(base_directory,"chirps_corrected/", sep = "")

# Indicate whether to use Ground Corrected CHIRPS (TRUE) or basic CHIRPS (FALSE)
use_corrected_chirps <- FALSE

# Location for Satellite Precipitation Product to be bias-corrected (INPUT)
# Files should be organized in subdirectories by year, following the pattern: YYYY/YYYY.mm.dd.tif
spp_dir <- paste(base_directory,"persiann/", sep = "")

# Location for corrected Satellite Precipitation Product files (OUTPUT)
# Files will be written to subdirectories by year, following the pattern: YYYY/YYYY.mm.dd.tif 
corrected_spp_dir <- paste(base_directory, "persiann_corrected_linear/", sep = "")
