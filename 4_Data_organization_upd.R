library(dplyr)
library(sf)
library(tidyr)
library(random)
library(stringr)



###Load full dataset
Myfiles<- "Your/folder/"
setwd("Myfiles")

Full_data <- read.csv("Full_dataset.csv")[-1]

Knames <- colnames(Full_data[,8:49])

saveRDS(Knames, "Knames.rds")

write.csv(Full_data,"Kinzig_occurrences.csv" )

data_train <- Full_data%>%filter(year<=2016)

write.csv(data_train,"data_train.csv")

data_eval <- Full_data%>%filter(year>2016)


write.csv(data_eval,"data_eval.csv")

Kinzig_train_locs <- st_as_sf(data_train, 
                              coords= c("x","y"), crs= 4326 )

Kinzig_eval_locs <- st_as_sf(data_eval, 
                              coords= c("x","y"), crs= 4326 )

st_write(Kinzig_train_locs[,c(1:6)], "Kinzig_PA_train.gpkg", overwrite=TRUE)

st_write(Kinzig_train_locs[,c(1:6)], "Kinzig_PA_evaluation.gpkg", overwrite=TRUE)



