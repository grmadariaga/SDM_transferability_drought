if(!require(sf)){
  install.packages("sf")
  library(sf)
}

if(!require(sp)){
  install.packages("sp")
  library(sp)
}


if(!require(raster)){
  install.packages("raster")
  library(raster)
}

if(!require(dismo)){
  install.packages("dismo")
  library(dismo)
}

if(!require(lubridate)){
  install.packages("lubridate")
  library(lubridate)
}

if(!require(foreach)){
  install.packages("foreach")
  library(foreach)
}

if(!require(parallel)){
  install.packages("parallel")
  library(parallel)
}

if(!require(doParallel)){
  install.packages("doParallel")
  library(doparallel)
}

year <- as.character(2005:2021)
month <- c("01_Jan", "02_Feb", "03_Mar", "04_Apr", "05_May", "06_Jun", "07_Jul", "08_Aug", "09_Sep", "10_Oct", "11_Nov", "12_Dec")
submonth <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

###First download rasters from DWD for Temperature minima, Temperature Maxima, and Precipitation for each month in the period considered

###Tmax

for (m in month) {
  for (y in year) {
    subm <- substr(m, start = 1, stop = 2)  # Extract the submonth from the month
    
    if (subm %in% submonth) {  # Check if the submonth is valid
      #tryCatch({
      download.file(paste0("https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_max/",
                           m, "/grids_germany_monthly_air_temp_max_", y, subm, ".asc.gz", sep = ""),
                    destfile = paste0("Your/Folder/Tmax/", y,"_",m,"_", subm, ".asc.gz"), ###Change to your own directory
      )
      #}, error = function(e) {})
    }
  }
}

##Tmin

for (m in month) {
  for (y in year) {
    subm <- substr(m, start = 1, stop = 2)  # Extract the submonth from the month
    
    if (subm %in% submonth) {  # Check if the submonth is valid
      #tryCatch({
      download.file(paste0("https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_min/",
                           m, "/grids_germany_monthly_air_temp_min_", y, subm, ".asc.gz", sep = ""),
                    destfile = paste0("Your/Folder/Tmin/", y,"_",m,"_", subm, ".asc.gz"), ###your directory
      )
      #}, error = function(e) {})
    }
  }
}

###Precipitation

for (m in month) {
  for (y in year) {
    subm <- substr(m, start = 1, stop = 2)  # Extract the submonth from the month
    
    if (subm %in% submonth) {  # Check if the submonth is valid
      #tryCatch({
      download.file(paste0("https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/precipitation/",
                           m, "/grids_germany_monthly_precipitation_", y, subm, ".asc.gz", sep = ""),
                    destfile = paste0("Your/Folder/Precipitation/", y,"_",m,"_", subm, ".asc.gz"), ###Your directory
      )
      #}, error = function(e) {})
    }
  }
}


###DWD come in different projections compared to Kinzig, we can convert all rasters with the code below


# Set the source and destination coordinate reference systems
src_crs <- "+init=EPSG:31467"
dst_crs <- "+init=EPSG:25832"

# Set the input folder where the downloaded ASC files are located
input_tmax <- "Your/Folder/Tmax"

# Set the output folder for the reprojected rasters
output_tmax <- "Your/output/Folder/Tmax/25832"

# List the ASC files in the input folder
asc_tmax <- list.files(input_tmax, pattern = "//.asc$", full.names = TRUE)

# Create the output folder if it doesn't exist
dir.create(output_tmax, showWarnings = FALSE)

# Loop through each ASC file
for (i in asc_tmax) {
  # Read the raster from the ASC file
  r <- raster(i)
  crs(r) <- src_crs
  
  # Reproject the raster to EPSG:25832
  r_proj <- projectRaster(r, crs = dst_crs)
  
  # Create the output file path with .tif extension
  output_tmax_i <- file.path(output_tmax, paste0(tools::file_path_sans_ext(basename(i)), ".tif"))
  
  # Save the reprojected raster as a GeoTIFF
  writeRaster(r_proj, filename = output_tmax_i, format = "GTiff", overwrite = TRUE)
}


# Set the input folder where the downloaded ASC files are located
input_tmin <- "Your/Folder/Tmin"

# Set the output folder for the reprojected rasters
output_tmin <- "Your/output/Folder/Tmin/25832"

# List the ASC files in the input folder
asc_tmin <- list.files(input_tmin, pattern = "//.asc$", full.names = TRUE)

# Create the output folder if it doesn't exist
dir.create(output_tmin, showWarnings = FALSE)

# Loop through each ASC file
for (i in asc_tmin) {
  # Read the raster from the ASC file
  r <- raster(i)
  crs(r) <- src_crs
  
  # Reproject the raster to EPSG:25832
  r_proj <- projectRaster(r, crs = dst_crs)
  
  # Create the output file path with .tif extension
  output_tmin_i <- file.path(output_tmin, paste0(tools::file_path_sans_ext(basename(i)), ".tif"))
  
  # Save the reprojected raster as a GeoTIFF
  writeRaster(r_proj, filename = output_tmin_i, format = "GTiff", overwrite = TRUE)
}

# Set the input folder where the downloaded ASC files are located
input_prec <- "Your/Folder/Precipitation"

# Set the output folder for the reprojected rasters
output_prec <- "Your/output/Folder/Precipitation/25832"

# List the ASC files in the input folder
asc_prec <- list.files(input_prec, pattern = "//.asc$", full.names = TRUE)

# Create the output folder if it doesn't exist
dir.create(output_prec, showWarnings = FALSE)

# Loop through each ASC file
for (i in asc_prec) {
  # Read the raster from the ASC file
  r <- raster(i)
  crs(r) <- src_crs
  
  # Reproject the raster to EPSG:25832
  r_proj <- projectRaster(r, crs = dst_crs)
  
  # Create the output file path with .tif extension
  output_prec_i <- file.path(output_prec, paste0(tools::file_path_sans_ext(basename(i)), ".tif"))
  
  # Save the reprojected raster as a GeoTIFF
  writeRaster(r_proj, filename = output_prec_i, format = "GTiff", overwrite = TRUE)
}



# Set the site coordinates dataframe
sites <- st_read("G:/RESIST/MedinaMadariaga/M2_Modeltransferability/Codes_data/Kinzig_IHApoints.gpkg") ###Your sites' coordinates
sites$x <- st_coordinates(sites)[,1]
sites$y <- st_coordinates(sites)[,2]
sites<- sites[,c("site_code","x", "y", "year", "sampling")]



###Set list of tmax

tmax_tif<- list.files(output_tmax, pattern = "TADXMM_\\d{2}_\\d{4}_01.tif", full.names = TRUE)


###Create
tmax_list_site <- list()

for (site_index in 1:nrow(sites)) {
  
  s <- sites$site_code[site_index]
  sampling<- sites$sampling[site_index]
  x <- sites[sites$site_code== s,"x"]
  y <- sites[sites$site_code== s,"y"]
  
  last_12_months <- list()
  # Get the last 12 months' month and year
  for (i in 0:11) {
    month_number <- sprintf("%02d", month(sampling %m-% months(i)))
    yearx <- year(sampling %m-% months(i))
    last_12_months[[12 - i]] <- c(month_number, as.character(yearx))
  }
  
  year_files <- c()
  
  for (m in last_12_months) {
    
      # Get files for the current year
    month_file <- grep(paste0(m[1],"_", m[2], "_01.tif$"), tmax_tif, value = TRUE)
    year_files <- c(year_files, month_file) 
  }
    # Stack the rasters for the current year
    raster_stack <- stack(year_files)
    
    # Extract the values for each month of the current year at the current site
    raster_values <- extract(raster_stack, cbind(x, y))
    
    print(paste("Processing", site_index, sampling))
    # Save the values for the current site and year in the list
    samp_key <- as.character(sampling)
    tmax_list_site[[s]][[samp_key]] <- raster_values
  }


 



###Set list of tmin

tmin_tif<- list.files(output_tmin, pattern = "TADNMM_\\d{2}_\\d{4}_01.tif", full.names = TRUE)



###Create
tmin_list_site <- list()

for (site_index in 1:nrow(sites)) {
  
  s <- sites$site_code[site_index]
  sampling<- sites$sampling[site_index]
  x <- sites[sites$site_code== s,"x"]
  y <- sites[sites$site_code== s,"y"]
  
  last_12_months<- list()
  # Get the last 12 months' month and year
  for (i in 0:11) {
    month_number <- sprintf("%02d", month(sampling %m-% months(i)))
    yearx <- year(sampling %m-% months(i))
    last_12_months[[12 - i]] <- c(month_number, as.character(yearx))
  }
  
  year_files <- c()
  
  for (m in last_12_months) {
    
    # Get files for the current year
    month_file <- grep(paste0(m[1],"_", m[2], "_01.tif$"), tmin_tif, value = TRUE)
    year_files <- c(year_files, month_file) 
  }
  # Stack the rasters for the current year
  raster_stack <- stack(year_files)
  print(paste("Processing", site_index, sampling))
  # Extract the values for each month of the current year at the current site
  raster_values <- extract(raster_stack, cbind(x, y))
  
 
  # Save the values for the current site and year in the list
  samp_key <- as.character(sampling)
  tmin_list_site[[s]][[samp_key]] <- raster_values
}



###Set list of prec

prec_tif<- list.files(output_prec, pattern = "RSMS_\\d{2}_\\d{4}_01.tif", full.names = TRUE)



###Create
prec_list_site <- list()

for (site_index in 1:nrow(sites)) {
  
  s <- sites$site_code[site_index]
  sampling<- sites$sampling[site_index]
  x <- sites[sites$site_code== s,"x"]
  y <- sites[sites$site_code== s,"y"]
  
  last_12_months<- list()
  # Get the last 12 months' month and year
  for (i in 0:11) {
    month_number <- sprintf("%02d", month(sampling %m-% months(i)))
    yearx <- year(sampling %m-% months(i))
    last_12_months[[12 - i]] <- c(month_number, as.character(yearx))
  }
  
  year_files <- c()
  
  for (m in last_12_months) {
    
    # Get files for the current year
    month_file <- grep(paste0(m[1],"_", m[2], "_01.tif$"), prec_tif, value = TRUE)
    year_files <- c(year_files, month_file) 
  }
  # Stack the rasters for the current year
  raster_stack <- stack(year_files)
  print(paste("Processing", site_index, sampling))
  # Extract the values for each month of the current year at the current site
  raster_values <- extract(raster_stack, cbind(x, y))
  
  
  # Save the values for the current site and year in the list
  samp_key <- as.character(sampling)
  prec_list_site[[s]][[samp_key]] <- raster_values
}
####Create a dataframe with all bioclim per year per site

biovar_df <- data.frame(
  site_code = character(),
  sample = Date(),
  stringsAsFactors = FALSE
)


for (site_index in 1:nrow(sites)) {
  s <- sites$site_code[site_index]
  sample_date <- sites$sampling[site_index]
  
  # Loop through each year
  # Get the lists of values for the current site and year
    Tmin_vec <- tmin_list_site[[s]][[as.character(sample_date)]]
    Tmax_vec <- tmax_list_site[[s]][[as.character(sample_date)]]
    prec_vec <- prec_list_site[[s]][[as.character(sample_date)]]
    
    # Calculate the biovars for the current site and year
    biovars_values <- biovars(prec_vec, Tmin_vec, Tmax_vec)
    
    # Create a data frame with the results for the current site and year
    result_row <- data.frame(
      site_code = s,
      sampling = sample_date,
      biovars_values,
      stringsAsFactors = FALSE
    )
    
    # Add the row to the results dataframe
    biovar_df <- rbind(biovar_df, result_row)
  }



####Because the values extracted from the rasters of DWD are 1/10
#### we need to divide all results

biovar_cols <- colnames(biovar_df)[startsWith(colnames(biovar_df), "bio")]

biovar_df <- biovar_df %>%
  mutate(across(all_of(biovar_cols), ~ ./10))

biovar_df <- biovar_df[!duplicated(biovar_df[1:2]),]

biovar_df$Year <- year(biovar_df$sampling)

###write csv

write.csv(biovar_df, "G:/RESIST/MedinaMadariaga/M2_Modeltransferability/Codes_data/biovarfullK.csv")

Cut_biovar <- biovar_df%>%
  mutate(time_frame = cut(Year, ###We create intervals of 3 years to find real absences
                          breaks = c(2007, 2009, 2012, 2015, 2018, 2021), 
                          labels = c('2007-2009', '2010-2012', '2013-2015', '2016-2018', '2019-2021'),
                          include.lowest = TRUE))

Cut_biovar <- Cut_biovar %>%
  group_by(time_frame, site_code) %>%
  mutate(across(bio1:bio19, ~ mean(.))) 

Cut_biovar <- Cut_biovar[!duplicated(Cut_biovar[c("site_code","time_frame")]),]

write.csv(Cut_biovar, "G:/RESIST/MedinaMadariaga/M2_Modeltransferability/Codes_data/biovarcutK.csv")
