# Bias Correction of Satellite Precipitation 

The scripts in this repository are used to bias-correct satellite-observed precipitation using [CHIRPS](http://chg.geog.ucsb.edu/data/chirps/).  These methods are derived from algorithms written in MATLAB, developed by the SWAAT research group at the University of Arizona, lead by Professor Juan Valdes (Roy et al. 2016).

The code offers two different techniques for bias correction:
a. Linear Scaling,
b. Quantile Mapping

The process is divided into several scripts, as follows:
  - `config.r`: Contains basic configuration information (e.g., base path to file locations, time extent). This file is not executable on itself, it is called from all the other scripts to gather this parameters.
  - `adjusted_CHIRPS.r`:  This step is optional, use it when reliable ground observations are available for the area of interest, during the time period to be bias corrected. Allows to correct CHIRPS data using observed values (rain gauge observations). The output is ground-corrected CHIRPS  data, which can be used in place of base CHIRPS for both the Linear and the Quantile bias-correction methods. 
  - `linear_method.r`: This script calculates monthly averages across all years for both CHIRPS and the Satellite Precipitation Product (SPP) to correct. These averages are then used to calculate the bias of the SPP compared to CHIRPS and applied to each daily SPP file.
  - `quantile_method.r`:  This script uses a probability density function (PDF) to calculate the two parameter of the distribution: gamma (λ), and theta (θ). The process ignores values lower than 1 mm (drizzle), and creates a cumulative density function (CDF) for each matrix. Next, it map CHIRPS PDF to the SPP being bias-corrected. The process is applied only when there are more than 5 non-zero unique values in each month (for a given grid cell).

The scripts were translated into Open Source code (in R) by Begum Rabeya Rushi, Faith Mitheu and Stella Masese.

For more information on the **SERVIR** Program, click [here](https://servirglobal.net) 

# Requirements
- Important: Edit the `config.r` file to match the locations of the files to process before running any of the scripts.
- An R interpreter. Please install R from the [CRAN website](https://cran.r-project.org/). As an alternative, you can use [Microsoft Open R](https://mran.microsoft.com/open). R 3.3.1 or later is needed.
- A text editor to make changes to the configuration file. An IDE supporting R is preferred, such as [RStudio](https://www.rstudio.com/).
- R packages: Raster, sp, rgdal, abind, MASS, Pscl, EDISON, MCMCpack, invgamma.
- CHIRPS data should be organized in subfolders per year. Each subfolder must be named "YYYY", using the four digit year value. Individual files should be named "YYYY.mm.dd.tif". For example: <<base path to CHIRPS data>>/2015/2015.01.01.tif
- Observed (Rain gauge) data, if available: A CSV file with the Station IDs and their location is required. The observation data should be presented in CSV format, with a header row for columns Date, Lat, Lon and PRCP (precipitation). Date should be written in YYYYmmdd format, Lat and Lon in decimal degrees and PRCP in mm. The files should be organized by folders named after the station identifier and the observations should be broken into a file per year. For example: <<path to observed data>>/StationID1/YYYY_daily.csv. 
- Input SPP files should be organized in a similar manner as the CHIRPS files. For example: <<base path to SPP data>>/2015/2015.01.01.tif

# SERVIR Credits

> Rushi, B.R., Science Coordination Office Regional Science Associate for Hindu-Kush Himalaya, UAH
> Mithue, F., Thematic Lead, Water Resources & Disasters, RCMRD
> Francisco Delgado, Science Coordination Office Geospatial IT Lead, USRA

# References

If you would like to use this resource, here is the reference:

> Rushi, B.R., Adams, E.C., Roy, T. Valdes, R.M., Valdes,J.B., Ellenburg,W.L., Anderson, E., Markert, K.N., Florescordova, A., Limaye, A. (In Preparation). Bias Correction of Satellite Precipitation using Open Source Integrated Development Environment . The Earth Observer. 
