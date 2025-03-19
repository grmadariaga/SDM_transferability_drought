library(tidymodels)
library(randomForest)
library(dplyr)
library(sirad)
library(parallel)
library(doParallel)
library(foreach)

###Create a variable with the name of the folderwhere the files were downloaded

Myfiles<- "Your/folder/"

setwd(Myfiles)

###Read species names

Knames <- readRDS("Knames.rds")

###Read df for each species

ssn_list <-readRDS("ssn_list_train_or.rds")

###Read formulas selected

formula_list <- readRDS("formulas_or.rds")


### Load ssn_objects 



ssn_object <- readRDS("ssn_object_eval_or.rds")
###Load functions

source("Functions.R")





rf_models <- list()

rf_dfs <- list()
for (name in Knames) {
  time_sp <- system.time({
    tryCatch({
      
      ssn_p_i <- ssn_object[[name]]
      
      print(paste0("Modelling ", name))
      #1. Import data set with presence - absences and predictors of one species 
      dataset <- getSSNdata.frame(ssn_p_i, "Obs")
      predict_sites <- getSSNdata.frame(ssn_p_i, "Kinzig_PA_evaluation")
      # output_path3 <- paste0("./FBAC/data/SDM/output/ClimateScenarios/",files_names[sp])
      # ifelse(!dir.exists(output_path3), dir.create(output_path3), FALSE)
      #Matrix to save results
      thres <- matrix(nrow = 1, ncol = 5)
      perform <- matrix(nrow = 2, ncol = 5)
      PredProb <- matrix(nrow = nrow(predict_sites), ncol = 6)
      PredProb[,6] <- predict_sites$id 
      colnames(PredProb) <- c( paste0("prob_split", 1:5),"id")
      
      predPresAbs <- matrix(nrow = nrow(predict_sites), ncol = 6)
      predPresAbs[,6] <- predict_sites$id
      colnames(predPresAbs) <- c( paste0("pres_abs_fold", 1:5),"id")
      
      print("prep_ready")
      
      # # Import raster with distances from occurrence points
      # distance_raster <- raster(paste0("./sdm/distance_point/", files_names[sp], ".tif"))
      # # Extract distance values
      # distance <- terra::extract(distance_raster, points_prediction) 
      # # Re-scale distance
      # distance_rescal <- rescale(distance, to = c(1, 0)) 
      # # Add a column with subcatchment ids to rescaled distance
      # distance_rescal <- data.frame('id' = points_prediction$stream, distance_rescal)
      
      set.seed(1234)
      
      
      #2. Run 10 replicates of the model
      for (r in 1:5) {
        
        print(paste("Processing", name,  "Fold", 1))
        
        #Observations used for calibration
        CalIds <- ssn_list[[name]][["folds"]][[1]][[r]][["in_id"]]
        CalibData <- dataset[dataset$pid %in% CalIds,]
        #Observations selected for evaluation 
        
        EvaData <- dataset[!dataset$pid %in% CalIds,]
        
        # Subset of distances. Will be used to multiply probabilities obtained 
        # from each calibration split
        # distance_subset <- distance_rescal |>
        #   filter(id %in% EvaIds)
        #
        mod_rf <- NULL
        # convert the response to factor for RF model to return probabilities
        CalibData[, name] <- as.factor(CalibData[,name])
        #CalibData <- CalibData[,c(2:11, 14)]
        
        # with covariate bias correction
        #CalibData <- CalibData[,c(2:14)]
        
        print("Calculating class weights")
        # calculating the class weights and sample size
        prNum <- as.numeric(table(CalibData[,name])["1"]) # number of presences
        bgNum <- as.numeric(table(CalibData[,name])["0"]) # number of backgrounds
        # cwt <- c("1" = 1, "0" = prNum / bgNum)
        samsize <- c("0" = bgNum, "1" = prNum)
        
        frml <- formula_list[[name]]%>%as.character()%>%sub("as\\.factor\\((.*)\\).*", "\\1", .)
        frml <- as.formula(paste(frml[2], frml[1], frml[3], collapse = " "))
        
        
        print(frml)
        # start modeling! We use the "try" notation so if a species fails to fit, the loop will continue.
        if(inherits(try(
          mod_rf <- randomForest(formula= frml, 
                                 data = CalibData,
                                 ntree = 1000,
                                 sampsize = samsize,
                                 replace = TRUE,
                                 importance = TRUE)
          
          
        ), "try-error")){
          print(paste("Error for species", s, "from", r))
        }
        
        print("Model running")
        
        #Predict on evaluation data set
        #PredProbEval <- as.numeric(predict(mod_rf, EvaData[,2:11], type = "prob")[,"1"])
        
        # correction bias covariate to a constant value across all sub-catchments
        EvaData_bias_cor <- EvaData 
        
        PredProbEval <- as.numeric(predict(mod_rf, EvaData_bias_cor, type = "prob")[,"1"])
        print("Metric evaluation")
        #Evaluation metrics
        idsPresabsProb <- data.frame(EvaData[,c("pid",name)],PredProbEval)
        
        TSS_tr <- eval_mod(idsPresabsProb)$performance[5,2]
        AUC_tr <- eval_mod(idsPresabsProb)$performance[6,2]
        th_tr <- eval_mod(idsPresabsProb)$threshold[1,2]
        thres[,r] <- th_tr
        perform[,r] <- c(TSS_tr, AUC_tr)
        
        print(paste("TSS", TSS_tr, "AUC_tr", AUC_tr))
        if (AUC_tr >= 0.7 & TSS_tr >= 0.5){
          
          # variable importance 
          # rf_var_importance[,r] <- importance(mod_rf, type = 1)[,1]
          
          #Predict on complete data set
          #pred <- as.numeric(predict(mod_rf, predict_sites[,1:10], type = "prob")[,"1"])
          
          pred_var_bias_cor <- predict_sites 
          
          pred <- as.numeric(predict(mod_rf, pred_var_bias_cor, type = "prob")[,"1"])
          
          # Multiply by distance
          #pred <- pred * distance_rescal$distance_rescal
          PredProb[,r] <- pred
          # pred_ACCESS1 <- as.numeric(predict(mod_rf, ACCESS1[,2:15], type = "prob")[,"1"])
          # pred_CESM1 <- as.numeric(predict(mod_rf, CESM1[,2:15], type = "prob")[,"1"])
          # pred_CMCC <- as.numeric(predict(mod_rf, CMCC[,2:15], type = "prob")[,"1"])
          # pred_MIROC5 <- as.numeric(predict(mod_rf, MIROC5[,2:15], type = "prob")[,"1"])
          # pred_MPI <- as.numeric(predict(mod_rf, MPI[,2:15], type = "prob")[,"1"])
          # pred_scenarios <- data.frame(pred_ACCESS1, pred_CESM1, pred_CMCC, pred_MIROC5, pred_MPI)
          print("Binary transformation")
          #Binary transformation 
          presabs <- ifelse(pred >= th_tr,1,0)
          predPresAbs[,r] <- presabs 
          
          obsPresAbsEval <- predict_sites %>%
            dplyr::select(pid,name)
          
          print("Presence and absences binary processing")
          
          idsPresabsProbEval <- data.frame(obsPresAbsEval,pred)
          TSS_ev <- eval_mod(idsPresabsProbEval)$performance[5,2]
          AUC_ev <- eval_mod(idsPresabsProbEval)$performance[6,2]
          th_ev <- eval_mod(idsPresabsProbEval)$threshold[1,2]
          
          print(paste("TSS_EV", TSS_ev))
          print(paste("AUC_EV", AUC_ev))
          
          
          rf_models[[name]][[r]]<- list("Ids"=CalIds,
                                        "AUC_tr"=AUC_tr,
                                        "TSS_tr"=TSS_tr,
                                        "th_tr"=th_tr,
                                        "Prob"=PredProb,
                                        "PredOcc"=presabs,
                                        "AUC_ev"=AUC_ev,
                                        "TSS_ev"=TSS_ev,
                                        "th_ev"=th_ev)
          
          # presabs_ACCESS1 <- ifelse(pred_ACCESS1 >= th,1,0)
          # presabs_CESM1 <- ifelse(pred_CESM1 >= th,1,0)
          # presabs_CMCC <- ifelse(pred_CMCC >= th,1,0)
          # presabs_MIROC5 <- ifelse(pred_MIROC5 >= th,1,0)
          # presabs_MPI <- ifelse(pred_MPI >= th,1,0)
          # presabs_scenarios <- data.frame(presabs_ACCESS1, presabs_CESM1, presabs_CMCC, presabs_MIROC5, presabs_MPI)
          
          # write.csv(pred_scenarios, row.names = F, paste0(output_path3,"/pred_rf_",r,".csv"))
          # write.csv(presabs_scenarios, row.names = F, paste0(output_path3,"/presabs_rf_",r,".csv"))
        }
        
        
      }
      #Export results
      # print(paste0("Export results for ", files_names[sp]))
      # outputpath2 <- paste0("./FBAC/data/SDM/output/rf/",files_names[sp])
      # ifelse(!dir.exists(outputpath2), dir.create(outputpath2), FALSE)
      # write.csv(thres, row.names = F, paste0("./FBAC/data/SDM/output/rf/",files_names[sp], "/thres.csv"))
      # write.csv(perform, row.names = F, paste0("./FBAC/data/SDM/output/rf/",files_names[sp], "/perform.csv"))
      # write.csv(PredProb, row.names = F, paste0("./FBAC/data/SDM/output/rf/",files_names[sp], "/pred.csv"))
      
      #mean predictor importance
      #rf_var_importance[,11] <- apply( rf_var_importance[,1:10], 1, mean, na.rm=T)
      # write.csv(rf_var_importance, paste0(output_path,"/importance_rf_.csv"))
      
      #  Mean probability  
      
      mean_pro = apply(PredProb[,1:5], 1, mean, na.rm=T)
      
      # CV
      cv <- function (x) (sqrt(var(x, na.rm = T) / length(x)))/ mean(x, na.rm = T)
      
      cv_mean_pro = apply(PredProb[,1:5], 1, cv)
      
      # mean threshold
      th_mean <- mean(thres) 
      
      # presence/absence of the mean probability
      mean_pres_abs <- ifelse(mean_pro >= th_mean,1,0) 
      
      # probabilities for each split, mean and cv of probabilities
      # across splits and presence-absence in a data frame
      
      
      
      prediction_df <- data.frame(PredProb,  mean_pro,
                                  cv_mean_pro, predPresAbs[,1:5],
                                  mean_pres_abs)
      
      rf_dfs[[name]] <- prediction_df
      
      print("done")
      
      
    }, error=function(e){})
    
  }) 
}

saveRDS(rf_models, "rf_models_results_or.rds")
saveRDS(rf_dfs, "presabs_rf_or.rds")






