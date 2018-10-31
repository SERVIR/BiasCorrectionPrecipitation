# Script: adjusted_CHIRPS.R
# Purpose: Linear adjustment for CHIRPS using Rain Gauge Data
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

# Loading Ground Observation Point file for Area of Interest
groundObs_points <- read.csv(ground_stations_file,header = T) # read the station names
groundObs_coords <-data.frame(groundObs_points$Lon,groundObs_points$Lat) # Only Lat and lon have been kept. 
colnames(groundObs_coords) <- c("Lon","Lat")
groundObs_co <- data.matrix(groundObs_coords,rownames.force = T)
ground_location_count <-nrow(groundObs_co)
groundObs_points <- as.matrix(groundObs_points) 

print("Ground Observation Points Loaded")


# Extract cell information for CHIRPS data on observed locations
# For every year in the Time Extent
for (yyyy in start_year:end_year) {

  date_begin <- paste(yyyy,start_month,start_day,sep = "-") # Start of the time period 
  date_end <- paste(yyyy,end_month,end_day,sep = "-") # End of the time period

  # Create day sequence
  date1 <- seq(as.Date(date_begin), as.Date(date_end), "days")

  # Replace "-" with "."
  date_seq<- gsub("-",".",date1)

  # Create output sub-folder for the year
  chirps_obs_outfolder <- paste(comparison_dir, yyyy, "/", sep="")
  if(!file.exists(chirps_obs_outfolder)) {
	dir.create(chirps_obs_outfolder)
  }
  
  # For every day of the sequence
  for (d in 1:length(date_seq)) {
    # Create empty matrix
    st_chirps <- data.frame(matrix(NA, nrow = ground_location_count, ncol = 5))
    # Assign column names
    colnames(st_chirps) <- c("StationName","Lon","Lat","CHIRPSPrcp","ObsPrcp")
    # Re-read Station Names and coordinates
    st_chirps$StationName <- groundObs_points[,1]
    st_chirps[,2:3]<-data.frame(coordinates(groundObs_co))
    
    # Input file for the selected day
    inputfile<- paste(chirps_dir, yyyy, "/", date_seq[d],".tif",sep = "") 
    
    chirps_global_crop <- raster(inputfile)
    # Extract CHIRPS value at ground observation locations
    st_chirps[,4] <- data.frame(extract(chirps_global_crop,groundObs_co))
    
    # For each station
    for(p in 1:ground_location_count) {
      # print( paste("Station ", groundObs_points[p,1], sep="" ))
      obsinfolder <- paste(observed_rainfall_dir, groundObs_points[p,1], "/", yyyy,sep="")
      inputfile <- paste(obsinfolder, "/", yyyy, "_daily.csv", sep="")
      if ( file.exists(inputfile) ) {
        obs_data <- read.csv( paste( inputfile, sep = ""), header = T)
        # Checking if a particular date is found on the inputfile
        obs_data_index <- which( obs_data[,1]==gsub("[.]","", date_seq[d]) )
        if ( length(obs_data_index) > 0 ) {
          st_chirps[p,5] <- obs_data[obs_data_index,4]
          # print( paste("Date, station, st_chirps: ", date_seq[d], ", ", groundObs_points[p,1], ", ", st_chirps[p,5], sep="" ) )
        } else {
            print( paste("Observed data for ", groundObs_points[p,1], " for ", date_seq[d], " is not available", sep="") )
        }
        # Replace NA values with 0
        st_chirps[is.na(st_chirps)] <- 0
        
        chirps_obs_prcp <- as.data.frame(st_chirps)
        chirps_obs_prcp_final <- data.frame(chirps_obs_prcp$CHIRPSPrcp, chirps_obs_prcp$ObsPrcp)
        colnames(chirps_obs_prcp_final) <- c("chirps_prcp","obs_prcp")
        # Note: Why make all values lower than 1 equal to 0????
        chirps_obs_prcp_final$chirps_prcp <- ifelse(chirps_obs_prcp_final$chirps_prcp<1, 0, chirps_obs_prcp_final$chirps_prcp)
        chirps_obs_prcp_final$obs_prcp <- ifelse(chirps_obs_prcp_final$obs_prcp<1, 0, chirps_obs_prcp_final$obs_prcp)
        
      }
      chirps_obs_outfile <- paste(chirps_obs_outfolder, date_seq[d], "_chirps_obs_prcp.csv", sep="")
      write.table(chirps_obs_prcp_final, chirps_obs_outfile, row.names=FALSE, na="", col.names=T, sep=",")
    }
  }
}

month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Calculate Monthyly Average Component
for (mon in start_month:end_month) {
  myList <- list()
  # For each year in the specified time period
  for (yyyy in start_year:end_year) {
    inputfolder <- paste(comparison_dir, yyyy, "/", sep="")
    file_list <- list.files(path=inputfolder, pattern=paste(yyyy,".",sprintf("%02d", mon), ".*", sep=""))
    for(comparison_file in file_list) {
      myList[[length(myList)+1]] <- as.matrix(read.csv(comparison_file))
    }
  }
  m_mean_temp <- Reduce("+", myList)/length(myList)
  m_mean <- as.data.frame( m_mean_temp )
  m_mean$Ratio <- m_mean$obs_prcp/m_mean$chirps_prcp
  m_mean_garbage_removal <- m_mean[!is.infinite(m_mean$Ratio)]
  m_mean_final <- m_mean_garbage_removal[!is.na(m_mean_garbage_removal$Ratio)]
  monthly_bias_factor <- as.numeric(mean(m_mean_final$Ratio))
  assign(paste(mon, "monthly_obs_bias_factor", sep="_"), monthly_bias_factor)
}

for (yyyy in start_year:end_year) {
  # Create YYYY folder
  correct_folder <- paste(corrected_chirps_dir, yyyy, "/", sep="")
  if(!file.exists(correct_folder)) {
    dir.create(correct_folder)
  }
  
  inputfolder_chirps <- paste(chirps_dir, yyyy, "/",sep = "")
  for (mon in start_month:end_month) {
    chirps_files <- list.files(path=inputfolder_chirps, pattern = paste(yyyy,".",sprintf("%02d", mon),".*",sep="")) #path=inputfolder
    
    for (chirps_file in chirps_files) {
      chirps_tif_file <- paste(inputfolder_chirps, chirps_file, sep = "")
      chirps_daily <- as.matrix(raster(chirps_tif_file))
      monthly_factor <- as.numeric(lapply(paste(mon,"monthly_obs_bias_factor", sep = "_"),get))
      corrected_chirps <- data.matrix(chirps_daily*monthly_factor,rownames.force = T)
      
      rb <- raster(corrected_chirps)
      # Use class(rb) to verify that rb is a "RasterLayer"
      crs(rb) <- CRS(crs_definition)
      
      # Extent coordinates provided in the configuration file
      extent(rb) <- raster_output_extent
      
      correct_chirps_file <- paste(correct_folder, chirps_file, sep="")
      writeRaster(rb, filename=correct_chirps_file, format="GTiff", overwrite=TRUE)
    }
  }
}
  
print("Finished running Linear adjustment for CHIRPS using Rain Gauge Data")
