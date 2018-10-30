# Script: quantile_method.R
# Purpose: Quantile adjustment for Satellite Precipitation Products using CHIRPS or Adjusted CHIRPS
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

# Load Statistical Libraries 

# Package: abind
# Combine multidimensional arrays into a single array. Generalization of 'cbind' and 'rbind'.
# Works with vectors, matrices, and higher-dimensional arrays.
require(abind)
# Package: MASS
# Functions and datasets to support Venables and Ripley
library(MASS)
# Package: pscl
# Political Science Computational Laboratory. Bayesian analysis of item-response theory (IRT) models, 
# roll call analysis, computing highest density regions, maximum likelihood estimation of zero-inflated
# and hurdle models for count data, goodness-of-fit measures for GLMs
library(pscl) 
# Package: EDISON
# Estimation of Directed Interaction from Sequences of Non-homogeneous gene expressions. MCMC simulation 
# to recontruct networks from time series data, using a non-homogeneous, time-varying dynamic Bayesian network
library(EDISON)
# Package: MCMCpack
# Markov Chain Monte Carlo (MCMC) package. Perform Bayesian inference using posterios simulation for a number
# of statistical models. Returns coda mcmc objects that can be summarized with the coda package
library(MCMCpack)
# Package: invgamma
# Inverse Gamma Distribution
library(invgamma)

# Memory Cleanup/Initialization
rm(list = ls())

# Read parameters from configuration file
source("./config.R")

month_names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Determine wether to use original CHIRPS or corrected CHIRPS (if available). 
# Set "use_corrected_chirps"" in config.R
if (use_corrected_chirps) {
  chirps_dataset_dir <- corrected_chirps_dir
} else {
  chirps_dataset_dir <- chirps_dir
}

print("Quantile Method. Initialization Complete - Check config values")

# Main process

for (mon in start_month:end_month){
  # Clear values from the previous month
  myList <- list()

  chirpsCounter <- 0
  sppCounter <- 0

  for (yyyy in start_year:end_year) {

    print(paste("Processing month ", mon, ", Year ", yyyy, sep=""))
    # Loading CHIRPS and SPP data for all years in the current month (mon)
    chirpsInputfolder <- paste(chirps_dataset_dir,yyyy,"/",sep = "")
    sppInputfolder<- paste(spp_dir,yyyy,"/",sep = "")
    chirpsFiles <- list.files(path=chirpsInputfolder,pattern = paste(yyyy,".",sprintf("%02d", mon),".*",sep=""))
    sppFiles <- list.files(path=sppInputfolder,pattern = paste(yyyy,".",sprintf("%02d", mon),".*",sep=""))
    
    # Create CHIRPS array
    for(aFile in chirpsFiles) { 
      chirpsCounter <- chirpsCounter + 1
      chirpsFile <- paste(chirpsInputfolder,aFile,sep = "")
      if (chirpsCounter==1){
    	  chirpsPrecip <- as.matrix(raster(chirpsFile))
      } else {
    	  chirpsPrecip <- cbind(chirpsPrecip,as.matrix(raster(chirpsFile)))
      }
      #myList[[length(myList)+1]]<-as.matrix(raster(i)) 
    }
  
  	# Create SPP array
  	for(aFile in sppFiles) { 
  	  sppCounter <- sppCounter + 1
  	  sppFile <- paste(sppInputfolder,aFile,sep = "")
  	  if (sppCounter==1){
  			sppPrecip<-as.matrix(raster(sppFile))
  			sppNames <- as.list(aFile)
  	  } else {
  			sppPrecip<-cbind(sppPrecip,as.matrix(raster(sppFile)))
  			sppNames<- cbind(sppNames,as.list(aFile))
  	  }
  	}
  }
  # print(chirpsCounter)
  # print(sppCounter)

  # Drizzle value, Less than 1 mm rain is considered drizzle
  # NOTE: Check that the values for SPP are in mm (for example, IMERG tiff files are originally produced in Tenths of a milimiter)
  drizzle <- 1
  # Make 3-D Matrix for Each Month Considering All Years and discard Drizzle values
  # Precipitation arrays for CHIRPS (chirpsPrecip) and SPP (sppPrecip) are converted into 3-dimensional arrays. 
  # The dimensions are Rows, Columns, Days
  
  dim(chirpsPrecip) <- c(dim(raster(chirpsFile))[1], dim(raster(chirpsFile))[2], chirpsCounter)
  chirpsPrecip[which(chirpsPrecip<drizzle)] <- 0
  
  dim(sppPrecip) <- c(dim(raster(sppFile))[1], dim(raster(sppFile))[2], sppCounter)
  sppPrecip[which(sppPrecip<drizzle)] <- 0
  
  GammaCDF_chirps <- array(0, dim=c(dim(raster(chirpsFile))[1], dim(raster(chirpsFile))[2], chirpsCounter))
  GammaCDF_SPP    <- array(0, dim=c(dim(raster(sppFile))[1],    dim(raster(sppFile))[2],    sppCounter))
  
  CHIRPSParmsLambda <-matrix(0, dim(raster(chirpsFile))[1], dim(raster(chirpsFile))[2])
  CHIRPSParmsTheta  <-matrix(0, dim(raster(chirpsFile))[1], dim(raster(chirpsFile))[2])
   
  GammaParmsLambda  <-matrix(0, dim(raster(sppFile))[1], dim(raster(sppFile))[2])
  GammaParmsTheta   <-matrix(0, dim(raster(sppFile))[1], dim(raster(sppFile))[2])
  
  for (m in 1:dim(sppPrecip)[1]) {
    for (n in 1:dim(sppPrecip)[2]) {
      
      # Extract non-Zero values from both arrays
      # From CHIRPS
      IndexNonZeroCHIRPS <- which(chirpsPrecip[m,n,]>0)
      NonZeroCHIRPS <- chirpsPrecip[m,n,][which(chirpsPrecip[m,n,]>0)]
      # From SPP
      IndexNonZeroSPP <- which(sppPrecip[m,n,]>0)
      NonZeroSPP <- sppPrecip[m,n,][which(sppPrecip[m,n,]>0)]
      
	  # Only apply correction if more than 5 non-zero unique values are present in both datasets for the current pixel
      if (length(IndexNonZeroCHIRPS)>5 & length(IndexNonZeroSPP)>5  & length(unique(NonZeroCHIRPS))>5 & length(unique(NonZeroSPP))>5) {
        CHIRPSParmsLambda[m,n] <- fitdistr(NonZeroCHIRPS, "gamma")$estimate[1] # lambda OR SHAPE
        CHIRPSParmsTheta[m,n]<-fitdistr(NonZeroCHIRPS, "gamma")$estimate[2] # theta or rate
        
        GammaParmsLambda[m,n]<-fitdistr(NonZeroSPP, "gamma")$estimate[1] # lambda
        GammaParmsTheta[m,n]<-fitdistr(NonZeroSPP, "gamma")$estimate[2] # theta
        
        NonZeroGammaCDF <- pgamma(NonZeroSPP, GammaParmsLambda[m,n], rate=GammaParmsTheta[m,n], log=FALSE)
        GammaCDF_SPP[m,n,IndexNonZeroSPP] <- qgamma(NonZeroGammaCDF, CHIRPSParmsLambda[m,n], CHIRPSParmsTheta[m,n]) # Inverse
      } else {
        print(NonZeroSPP)
        GammaCDF_SPP[m,n,IndexNonZeroSPP] <- NonZeroSPP # No bias correction is done if only 2 points are available
      }
    }
  }
  for (kk in 1:sppCounter) {
    corrected_spp_daily <- as.matrix(GammaCDF_SPP[,,kk])
    rb <- raster(corrected_spp_daily)
    crs(rb) <- CRS(crs_definition)
    class(rb)
    
    # replace with correct coordinates
    extent(rb) <- raster_output_extent
    if (!file.exists(corrected_spp_dir)) {
      dir.create(corrected_spp_dir)
    }
    corrected_spp_file <- paste(corrected_spp_dir, sppNames[kk], sep="")
    writeRaster(rb, filename = corrected_spp_file, format="GTiff", overwrite=TRUE)
  }
}

print("Finished reading data into CHIRPS and SPP arrays")  


