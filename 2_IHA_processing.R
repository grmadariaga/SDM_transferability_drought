library("devtools")
library("EflowStats")
library("readxl")
library("tidyverse")
library("reshape2")
library(data.table)
library(dplyr)
library(purrr)
library(sf)
library(leaflet)
library(magrittr)
library(stringr)
library(readr)
library(doParallel)
library(foreach)
library(SSN)

Myfiles<- "Your/folder/"

setwd(Myfiles)

# Function to calculate all IHA metrics
run_IHA <- function (sfd, wyMonth, drainArea){
  
  # Import data and cut to water year (starting from the months of the "wyMonth" variable)
  sfd=as.data.frame(sfd)
  sfd$Date=as.Date(sfd$Date)  
  sfd <- sfd[, c('Date', 'flo_out')]
  head(sfd)
  sfd_cut <- EflowStats::cut_dataToWaterYear(sfd, yearType= "water", wyMonth = wyMonth)
  
  # Get the annual maximum peak values of the defined years
  peakValues <- sfd_cut
  #print(peakValues)
  peakValues$wat_yr <- EflowStats::get_waterYear(sfd_cut$Date, wyMonth = wyMonth)
  head(peakValues)
  peakValues <- aggregate(peakValues$flo_out, by = list(peakValues$wat_yr), max)
  peakValues <- setNames(peakValues, c("peak_dt","peak_va"))
  peakValues$peak_dt <- as.Date(paste(peakValues$peak_dt,wyMonth,"01", sep = "-"))
  
  # Calculate the peak threshold if more than 1 year of data
  if (nrow(peakValues)>1) {
    sfd_for_threshold <- sfd_cut
    sfd_for_threshold$flo_out <- ifelse(sfd_for_threshold$flo_out == 0, 0.001, sfd_for_threshold$flo_out)
    peakThreshold <- EflowStats::get_peakThreshold(sfd_for_threshold, peakValues, perc = 0.6, yearType = "water", wyMonth = wyMonth)
  }else{
    peakThreshold <- NULL
  }
  
  ## Calculate 171 metrics
  allHIT <- EflowStats::calc_allHIT(sfd_cut, yearType="water", wyMonth=wyMonth, stats="all", drainArea=drainArea, floodThreshold=peakThreshold)
  return(allHIT)
}

########################
# Import SWAT data that contains channel parameters
sd <- fread("channel_sd_day.txt", skip = 3) 
# Extract column names 
col_name_sd <- read_lines("channel_sd_day.txt", n_max = 2)[2] %>%
  str_split_1(pattern = " ") %>%
  magrittr::extract(. != "")
# Add column names
names(sd) <- col_name_sd
#Subset columns of interest
sd_sel <- sd   %>%
  dplyr::select(yr, gis_id, mon, day, unit, flo_out)%>%
  mutate(Date= as.Date(with(., paste(yr, mon, day,sep="-")), "%Y-%m-%d"))%>%
  filter(yr>2005 & yr <2022)


rm(sd)
gc()
#Load SSN for Drainage (H2O Area)
Kpoints <- importSSN(file.path("KinzigIHA.ssn"))

Kpoints <- getSSNdata.frame(Kpoints, "Obs")

#### Subset to sites of interest
channels <- unique(Kpoints$ChaUnit)
unique_id <- unique(Kpoints$id)


###Merge to Strm order from river gpkg for sampling month
Rivs <- st_read("rivs1.gpkg")%>%st_drop_geometry()

#### If strmOrder 1 or 2, sampling in April. If 3 or 4 sampling in June
Kpoints<- left_join(Kpoints, Rivs[,c("ChaUnit","strmOrder")], by="ChaUnit")

#### Create list with dataframes and inputs per site for EFlowStats

dflist <- list()

for (i in unique_id){
  g <- Kpoints[Kpoints$id==i, "ChaUnit"]
    df <- sd_sel%>%
    filter(,gis_id== g)%>%
    dplyr::select("Date", "flo_out")
  
  WaterMonth <- ifelse(Kpoints[Kpoints$ChaUnit == g, "strmOrder"] %>% unique() %in% c(1, 2), 4, 6)
  
  drainArea <- Kpoints[Kpoints$id== i, "H2OAreaA" ]%>%unique()%>%as.numeric()%>%max()
  
  dflist[[i]] <-list(df = df, Wmonth = WaterMonth, drainArea = drainArea)
}


dflist_WY <- list()
for (i in unique_id){
  # Calculate the water year and cut data
  print(i)
  if (Kpoints[Kpoints$id==i, "ChaUnit"] != 48)
  dfWY <- EflowStats::cut_dataToWaterYear(dflist[[i]][["df"]], yearType = "water", wyMonth = dflist[[i]][["Wmonth"]]) %>%as.data.frame()
  dfWY$wat_year <- EflowStats::get_waterYear(dfWY$Date, wyMonth = dflist[[i]][["Wmonth"]])
  
  dflist_WY[[i]] <- dfWY
}


# Calculate the IHA for each year in data
# Set the number of cores to use
num_cores <- 6  # Adjust the number of cores as needed

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)
#parallel::mcaffinity(1:20)

iha_all <- foreach(i = unique_id, .combine = 'c', .errorhandling= "remove") %do% {
  wydf_subset <- dflist_WY[[i]]
  wyMonth <- dflist[[i]][["Wmonth"]]
  drainArea <- dflist[[i]][["drainArea"]]
  years <- unique(wydf_subset$wat_year)
  iha_list <- vector("list",length = 15)
  x <- 1
  foreach(year = years, .combine = 'c') %do% {
    print(paste("Executing year: ", year, "point", i))
    dfwy_year <- wydf_subset[wydf_subset$wat_year == year,]
    dfwy_year <- dfwy_year[, c('Date', 'flo_out')]
    iha <- run_IHA(dfwy_year, wyMonth=wyMonth, drainArea=drainArea)
    iha$id <- i
    iha$sampling <- Kpoints[Kpoints$id==i, "sampling"]
    iha <- setNames(iha, c("IHA", year, "id", "sampling"))
    iha_list[[x]] <- iha
    x <- x + 1
    
  }
  iha_list
}

# Stop the parallel backend
stopCluster(cl)

dfarr <- foreach(i = 1:length(iha_all), .combine = 'rbind') %do% {
  df1 <- as.data.frame(do.call(cbind, iha_all[[i]]))
  df2 <- df1[, !duplicated(names(df1))]
  u<- unique(df1$id)
  sam <-unique(df1$sampling)
  df3 <- df2 %>%
    dplyr::select(-c(id, sampling))%>%
    pivot_longer(cols = -IHA, names_to = "Year", values_to = "Value") %>%
    spread(key = IHA, value = Value)
  df3$id <- u
  df3$sampling <- sam
  as.data.frame(df3)
}


write.csv(dfarr, "IHA_Kinzig.csv")




