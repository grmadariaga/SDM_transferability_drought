library(tidymodels)
library(SSN)
library(sirad)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(MuMIn)




###Merge SSN object from STARS to dataframe on environmental variables

###Load SSN object Kinzig
###Create a variable with the name of the folderwhere the files were downloaded

Myfiles<- "Your/folder/"

setwd(Myfiles)



#createDistMat(K, predpts="Kinzig_PA_evaluation",  o.write = TRUE)
#distPred1km <- getStreamDistMat(K, Name = "Kinzig_PA_evaluation")

#dmats <- getStreamDistMat(K)


K<-additive.function(K, "H2OArea", "computed.afv")

names(K)

###Get observation dataframe

obs_df <- getSSNdata.frame(K,"Obs")


###Load Kinzig occurrence

Kinzig_pa <- read.csv("data_train.csv")[,-1]
filtered_cols <- Kinzig_pa %>%
  dplyr::select(8:49)
Kinzig_df <- Kinzig_pa %>%
  dplyr::select(1:7, all_of(names(filtered_cols)))

Kinzig_df <- left_join(obs_df, Kinzig_df[, c(1,8:49)], by= c("id"))

#### Load Bioclim variables
Biovar <- read.csv("biovarfullK.csv")[,-1]


Kinzig_bio<- left_join(Kinzig_df, Biovar[, c(1,3:21)], by= "id")%>%
  mutate(wat_year= year-1)

####Load IHA

IHA_cut <- read.csv("IHA_Kinzig.csv")%>%
  mutate(wat_year= Year)%>%
  dplyr::select(-Year)%>%
  select(-sampling)

####

Kinzig_all <- left_join(
  Kinzig_bio, 
  IHA_cut[,c(2:174)],
  by= c("id", "wat_year"))%>% 
  #filter(Channel != 48)%>% ### Channel 48 is empty we remove it and we return it after variable selection
  dplyr::select_if(~ !any(is.na(.)))%>%  ####Remove columns with NA
  #dplyr::select_if(~ !all(. == .[[1]]))%>%###Remove columns with the same value for all sites
  #dplyr::select(c(colnames(filtered_cols), where(~ n_distinct(.) >= 4)))%>%###Remove columns with less than 4 different values (except for species)
  dplyr::select(where(~ !any(is.infinite(.))))%>%### Remove columns with inf
  dplyr::select(-wat_year)


new_columns <- readRDS("new_cols_red.rds")


Knames <- readRDS("Knames.rds")

obs_df <- Kinzig_all[,new_columns]



####Extract names to an object to avoid for rerunning this process

#saveRDS(Knames, "SDM_Kinzig_updated/Knames.rds")


#cv_ids <- list()

#for (k in Knames){
#  folds <- ssn_list[[k]][["folds"]]
#  cv_ids[[k]] <- folds
#}


####Extract cv_ids to an object

#saveRDS(cv_ids,"/SDM_Kinzig_updated/cv_ids.rds")



### Create object list for training in a dataset and evaluating in another
#Knames_ev <- paste0(Knames, ".ev")

#obs_ev <- obs_df

#obs_ev <- obs_ev%>%
#  rename_with(~ paste0(., ".ev"), any_of(Knames))

#for (i in Knames_ev){
#  obs_ev[,i]<- ifelse(obs_ev$dataset=="Senckenberg", NA,obs_ev[,i])
#}


ssn_objects <- list()

# Loop through each name in Knames
for (name in Knames) {
  # Subset the obs_df dataframe
 subset <- obs_df[, c(names(obs_df[,c(1:25, 68:125)]), name)]
  
  
  # Create a new SSN object with the subset
  ssn_obj <- putSSNdata.frame(subset, K)  
  
  # Assign the SSN object to the list with the corresponding name
  ssn_objects[[name]] <- ssn_obj
}


#saveRDS(ssn_objects, "SDM_Kinzig_updated/ssn_objects.rds")

preds <- getSSNdata.frame(K, "Kinzig_PA_evaluation")


Kinzig_pa <- read.csv("Kinzig_occurrences.csv")[,-1]
filtered_cols <- Kinzig_pa %>%
  dplyr::select(8:49)
Kinzig_preds <- Kinzig_pa %>%
  dplyr::select(1:7, all_of(names(filtered_cols)))

Kinzig_preds <- left_join(preds, Kinzig_preds[, c(1,8:49)], by= c("id"))

#### Load Bioclim variables
Biovar <- read.csv("biovarfullK.csv")[,-1]


Kinzig_bio_preds<- left_join(Kinzig_preds, Biovar[, c(1,3:21)], by= "id")%>%
  mutate(wat_year= year-1)

####Load IHA

IHA_cut <- read.csv("IHA_full.csv")%>%
  mutate(wat_year= Year)%>%
  dplyr::select(-Year)

####

Kinzig_all_preds <- left_join(
  Kinzig_bio_preds, 
  IHA_cut[,c(2:174)],
  by= c("id", "wat_year"))



Kinzig_all_preds <- Kinzig_all_preds[,new_columns]



ssn.obj <- list()

for (name in Knames) {
  # Subset the obs_df dataframe
  
  ssn_objects[[name]]@obspoints@SSNPoints[[1]]@point.data<- obs_df[, c(names(obs_df[,c(1:25, 68:125)]), name)]
  
  
  # Create a new SSN object with the subset
  ssn_objects[[name]]@predpoints@SSNPoints[[1]]@point.data <-Kinzig_all_preds[, c(names(Kinzig_all_preds[,c(1:25, 68:125)]), name)]
  
  
  # Assign the SSN object to the list with the corresponding name
  ssn.obj[[name]] <- ssn_objects[[name]]
}


saveRDS(ssn.obj, "ssn_object_eval_red.rds")


Kinzig_full_data <- rbind(obs_df,Kinzig_all_preds)

write.csv(Kinzig_full_data, "Kinzig_full_data_red.csv")
