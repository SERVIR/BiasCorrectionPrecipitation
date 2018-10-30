# Script: linear_method.R
# Purpose: Linear adjustment for Satellite Precipitation Products using CHIRPS or Adjusted CHIRPS
# Author: Francisco Delgado (francisco.delgadoolivares@nasa.gov), SERVIR Science Coordination Office
# Date: September 7, 2018, revised October 30, 2018
# Based on code prepared by Begum Rushi (begumrabeya.rushi@nasa.gov), Regional Associate, HKH Region, SERVIR, as of August 4, 2018

# Clear up in memory variables
rm(list = ls())

# load Required Libraries
library(rgdal)
library(raster)
library(sp)
require(raster)

# Read parameters from configuration file
source("./config.R")
print("Initialization Complete - Check config values")

# Using config value to select CHIRPS or corrected CHIRPS
if (use_corrected_chirps) {
  chirps_dataset_dir <- corrected_chirps_dir
} else {
  chirps_dataset_dir <- chirps_dir
}

# Make mean monthly matrix for chirps and SPP

data_names <- c(chirps_dataset_dir, spp_dir)

for (d in 1:length(data_names)) {
  for (mon in start_month:end_month) {
    myList <- list()
    for (yyyy in start_year:end_year) {
      inputfolder <- paste(data_names[d], yyyy, "/", sep = "")
      setwd(inputfolder)
      files <- list.files(path = ".", pattern = paste(yyyy, ".", sprintf("%02d", mon), ".*", sep = ""))
      for (file in files) {
        # directory <- paste(inputfolder, file, sep = "")
        myList[[length(myList)+1]] <- as.matrix(raster(file))
      }
    }
    MonthlyMean <- Reduce("+", myList) / length(myList)
    mean_monthly <- as.matrix(MonthlyMean)
    rb <- raster(mean_monthly)
    # Use class(rb) to verify that rb is a "RasterLayer"
    crs(rb) <- CRS(crs_definition)
    
    # Extent coordinates provided in the configuration file
    extent(rb) <- raster_output_extent

    assign(paste(mon, data_names[d], "monthly_mean.tif", sep = "_"), rb)
    monthly_folder <- data_names[d]
    setwd(monthly_folder)
    monthly_tif <- paste(monthly_folder, mon, "_monthly_mean.tif", sep = "")
    writeRaster(rb, filename=monthly_tif, format="GTiff", overwrite=TRUE)
  }
}
print("Monthly Matrices generated")

# Calculate factors and write resulting rasters

for (yyyy in start_year:end_year) {
  for (mon in start_month:end_month) {
    
    # Calculate Monhtly Bias Factor
    monthly_chirps_factor <- raster(paste(chirps_dataset_dir, mon, "_monthly_mean.tif", sep = ""))
    monthly_spp_factor <- raster(paste(spp_dir, mon, "_monthly_mean.tif", sep = ""))
    monthly_bias_factor <- overlay(monthly_chirps_factor, monthly_spp_factor, fun=function(x,y){ x/y })
    
    inputfolder_spp <- paste(spp_dir, yyyy, "/", sep = "")
    setwd(inputfolder_spp)
    files_spp <- list.files(path=".", pattern = paste(yyyy, ".", sprintf("%02d", mon), ".*", sep = ""))
    for (file in files_spp) {
      spp_daily <- raster(paste(inputfolder_spp, file, sep = ""))
      # Calculate corrected SPP values
      corrected_spp <- overlay(spp_daily, monthly_bias_factor, fun=function(x,y){x*y})
      correct_folder <- paste(corrected_spp_dir, yyyy, "/", sep = "")

      if (!file.exists(correct_folder)) {
        dir.create(correct_folder)
      }
      setwd(correct_folder)
      corrected_file <- paste(correct_folder, file, sep = "")
      
      crs(corrected_spp) <- CRS(crs_definition)
      writeRaster(corrected_spp, filename=corrected_file, format="GTiff", overwrite=TRUE)
    }
  }
}

print("Finished running: Linear Method SPP Bias Correction b")
