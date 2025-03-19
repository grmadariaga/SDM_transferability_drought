library(dplyr)


###Create a variable with the name of the folderwhere the files were downloaded

Myfiles<- "Your/folder/"

setwd(Myfiles)


glm <- readRDS("glm_evaluation_results.rds")
ssn<- readRDS("ssn_evaluation_results.rds")
rf <- readRDS("rf_models_results.rds")
maxent <- readRDS("maxent_models.rds")

Knames <- readRDS("Knames.rds")

###Compile all results in a dataframe

#Start indicidually to check results

##MAXENT FIRST
# Create an empty list to store the results
maxent_list <- list()

# Loop through each name in maxent_models
for (name in Knames) {
      for (r in 1:5) {
        tryCatch({
    # Extract relevant information
    TSS_tr <- maxent[[name]][[r]]$TSS_tr
    AUC_tr <- maxent[[name]][[r]]$AUC_tr
    TSS_ev <- maxent[[name]][[r]]$TSS_ev
    AUC_ev <- maxent[[name]][[r]]$AUC_ev
    
    # Create a dataframe for the current iteration
    maxent_df <- data.frame(Species = name,
                            TSS_tr = TSS_tr,
                            AUC_tr = AUC_tr,
                            TSS_ev = TSS_ev,
                            AUC_ev = AUC_ev,
                            Fold=r,
                            model = "Maxent")
    
    # Append the dataframe to the result list
    maxent_list <- c(maxent_list, list(maxent_df))
        }, error=function(e){})
      }
  }

# Combine all dataframes into a single dataframe
maxent_dataframe <- do.call(rbind, maxent_list)


#### RANDOM FOREST

# Create an empty list to store the results
rf_list <- list()

# Loop through each name in maxent_models
for (name in Knames) {
  for (r in 1:5) {
    tryCatch({
      # Extract relevant information
      TSS_tr <- rf[[name]][[r]]$TSS_tr
      AUC_tr <- rf[[name]][[r]]$AUC_tr
      TSS_ev <- rf[[name]][[r]]$TSS_ev
      AUC_ev <- rf[[name]][[r]]$AUC_ev
      
      # Create a dataframe for the current iteration
      rf_df <- data.frame(Species = name,
                              TSS_tr = TSS_tr,
                              AUC_tr = AUC_tr,
                              TSS_ev = TSS_ev,
                              AUC_ev = AUC_ev,
                              Fold=r,
                              model = "RF")## Change model name
      
      # Append the dataframe to the result list
      rf_list <- c(rf_list, list(rf_df)) ##change arguments
    }, error=function(e){})
  }
}

# Combine all dataframes into a single dataframe
rf_dataframe <- do.call(rbind, rf_list)

####GLM

# Create an empty list to store the results
glm_list <- list()

# Loop through each name in maxent_models
for (name in Knames) {
  for (r in 1:5) {
    tryCatch({
      # Extract relevant information
      TSS_tr <- glm[[name]][[r]]$TSS_tr
      AUC_tr <- glm[[name]][[r]]$AUC_tr
      TSS_ev <- glm[[name]][[r]]$TSS_ev
      AUC_ev <- glm[[name]][[r]]$AUC_ev
      
      # Create a dataframe for the current iteration
      glm_df <- data.frame(Species = name, ###change variable name
                              TSS_tr = TSS_tr,
                              AUC_tr = AUC_tr,
                              TSS_ev = TSS_ev,
                              AUC_ev = AUC_ev,
                              Fold=r,
                              model = "glm") ###change model name
      
      # Append the dataframe to the result list
      glm_list <- c(glm_list, list(glm_df)) ###change arguments
    }, error=function(e){})
  }
}

# Combine all dataframes into a single dataframe
glm_dataframe <- do.call(rbind, glm_list)

#### SSN

# Create an empty list to store the results
ssn_list <- list()

# Loop through each name in maxent_models
for (name in Knames) {
  for (r in 1:5) {
    tryCatch({
      # Extract relevant information
      TSS_tr <- ssn[[name]][[r]]$TSS_tr
      AUC_tr <- ssn[[name]][[r]]$AUC_tr
      TSS_ev <- ssn[[name]][[r]]$TSS_ev
      AUC_ev <- ssn[[name]][[r]]$AUC_ev
      
      # Create a dataframe for the current iteration
      ssn_df <- data.frame(Species = name, ###change variable name
                           TSS_tr = TSS_tr,
                           AUC_tr = AUC_tr,
                           TSS_ev = TSS_ev,
                           AUC_ev = AUC_ev,
                           Fold=r,
                           model = "ssn") ###change model name
      
      # Append the dataframe to the result list
      ssn_list <- c(ssn_list, list(ssn_df)) ###change arguments
    }, error=function(e){})
  }
}

# Combine all dataframes into a single dataframe
ssn_dataframe <- do.call(rbind, ssn_list)



###Create a dataframe with all models 

Models_dataframe <- rbind(maxent_dataframe, 
                          ssn_dataframe, 
                          glm_dataframe,
                          rf_dataframe)



saveRDS(Models_dataframe, "df_models_all.rds")
