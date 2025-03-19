library(SSN)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)


###Create a variable with the name of the folderwhere the files were downloaded

Myfiles<- "Your/folder/"

setwd(Myfiles)


#### Load functions

source("Functions.R")


options(na.action = "na.omit")

####Load Knames

Knames <- readRDS("Knames.rds")


####Load list with Folds created for Training cv

ssn_list <- readRDS("ssn_list_train_or.rds")


####Load object with preds 

ssn_object <- readRDS("ssn_object_eval_or.rds")

####Load list with formulas created in Object_creation.R

null_AIC_list <-readRDS("formulas_or.rds")


###Load selected models
CorModels <- readRDS("ssn_selection_or.rds")


###START MODELLING



glm_null<- list()

for (name in Knames) {
  for (i in 1:5) {
    tryCatch({
      Fold <- paste0("Fold.", i)  
      
      ssn_p_i <- ssn_object[[name]]
      #Data frame with observation points
      ssn_test_dataDF <- getSSNdata.frame(ssn_p_i, "Obs")
      
      #Ids of observations selected for evaluation
      EvaIds <- ssn_list[[name]][["folds"]][[1]][[i]][["in_id"]]
      
      #Insert NAs in the column of the response variable (pres_abs) in the observations
      #selected for evaluation
      
      print("EvaIds correct")
      ssn_test_dataDF[, paste0(name,Fold)]<- ssn_test_dataDF[,name]
      ssn_test_dataDF[!ssn_test_dataDF$pid %in% EvaIds,paste0(name, Fold)] <- NA
      
      thres <- matrix(nrow = 1, ncol = 5)
      perform <- matrix(nrow = 2, ncol = 5)
      predict_sites <- getSSNdata.frame(ssn_p_i, Name = "Kinzig_PA_evaluation")
      PredProb <- matrix(nrow = nrow(predict_sites), ncol = 6)
      PredProb[,6] <- predict_sites$pid
      
      #Put the data frame with NAs in the SSN object.
      ssn_calib_data <- putSSNdata.frame(ssn_test_dataDF, ssn_p_i)
      
      
      print(paste("Processing", name, "Fold", i))
      
      
      frml <- null_AIC_list[[name]]%>%as.character()%>%sub("as\\.factor\\((.*)\\).*", "\\1", .)
      frml <- as.formula(paste(paste0(frml[2],Fold), frml[1], frml[3], collapse = " "))
      
      #With this, the regression omits observations with NAs and the model is fitted
      # only on the calibration data set
      options(na.action = "na.omit")
      
      #Fit the model with the calibration data set and parameters
      CalibSsn <- try(glmssn(formula =  frml,
                             family = "binomial",
                             ssn_calib_data,
                             CorModels = NULL,
                             addfunccol = "computed.afv",
                             control = list(trunc.pseudo=100)))
      
      
      print("Model running")
      
      # # correction bias covariate to a constant value across all sub-catchments
      #Insert 1 in the columns of the covariates used to correct observer bias
      
      eval_temp <- getSSNdata.frame(CalibSsn$ssn.object, "_MissingObs_") 
      
      CalibSsn_cor_eval <- putSSNdata.frame(eval_temp, CalibSsn, "_MissingObs_")
      
      #Predict points used for evaluation with the model trained with the
      #calibration data set
      print("Predicting...")
      
      PredEval <- predict.glmssn(CalibSsn_cor_eval, "_MissingObs_")
      PredProbEval <- SSNProb(PredEval)
      
      #Extract ids and observed presence-absences of points used for evaluation
      obsPresAbs <- getSSNdata.frame(ssn_p_i) %>%
        filter(!pid %in% EvaIds) %>%
        dplyr::select(pid,name)
      
      #Evaluation metrics
      idsPresabsProb <- data.frame(obsPresAbs,PredProbEval)
      TSS <- eval_mod(idsPresabsProb)$performance[5,2]
      AUC <- eval_mod(idsPresabsProb)$performance[6,2]
      th <- eval_mod(idsPresabsProb)$threshold[1,2]
      thres[,i] <- th
      perform[,i] <- c(TSS, AUC)
      
      print(paste("Fold training finished","AUC", AUC, "TSS", TSS))
      
      predPresAbs <- matrix(nrow = nrow(predict_sites), ncol = 6)
      predPresAbs[,6] <- predict_sites$id
      colnames(predPresAbs) <- c( paste0("Fold", 1:5),"pid")
      
      if (AUC >= 0.7 & TSS >= 0.5){
        
        # # correction bias covariate to a constant value across all sub-catchments
        #Insert 1 in the columns of the covariates used to correct observer bias
        
        temp_pred_df <- getSSNdata.frame(CalibSsn, "Kinzig_PA_evaluation")
        
        print("Evaluation starting")
        
        
        CalibSsn_cor_eval <- putSSNdata.frame(temp_pred_df, CalibSsn, "Kinzig_PA_evaluation")
        
        # predict on complete data set
        print("Predicting on Evaluation")
        
        PredAll <- predict.glmssn(CalibSsn_cor_eval, "Kinzig_PA_evaluation")
        
        # Extract probabilities
        PredProb[,i] <- SSNProb(PredAll)
        
        # pred_ACCESS1 <- predict.glmssn(CalibSsn, "ACCESS1")%>%
        #   SSNProb()
        # pred_CESM1 <- predict.glmssn(CalibSsn, "CESM1")%>%
        #   SSNProb()
        # pred_CMCC <- predict.glmssn(CalibSsn, "CMCC")%>%
        #   SSNProb()
        # pred_MIROC5 <- predict.glmssn(CalibSsn, "MIROC5")%>%
        #   SSNProb()
        # pred_MPI <- predict.glmssn(CalibSsn, "MPI")%>%
        #   SSNProb()
        # pred_scenarios <- data.frame(pred_ACCESS1, pred_CESM1, pred_CMCC, pred_MIROC5, pred_MPI)
        
        #Binary transformation
        presabs <- ifelse(PredProb[,i] >= th,1,0)
        predPresAbs[,i] <- presabs
        
        
        obsPresAbsEval <- getSSNdata.frame(ssn_p_i, "Kinzig_PA_evaluation") %>%
          dplyr::select(id,name)
        
        print("Presence and absences binary processing")
        
        idsPresabsProbEval <- data.frame(obsPresAbsEval,PredProb[,i])
        TSS_ev <- eval_mod(idsPresabsProbEval)$performance[5,2]
        AUC_ev <- eval_mod(idsPresabsProbEval)$performance[6,2]
        th_ev <- eval_mod(idsPresabsProbEval)$threshold[1,2]
        
        
        # presabs_ACCESS1 <- ifelse(pred_ACCESS1 >= th,1,0)
        # presabs_CESM1 <- ifelse(pred_CESM1 >= th,1,0)
        # presabs_CMCC <- ifelse(pred_CMCC >= th,1,0)
        # presabs_MIROC5 <- ifelse(pred_MIROC5 >= th,1,0)
        # presabs_MPI <- ifelse(pred_MPI >= th,1,0)
        # presabs_scenarios <- data.frame(presabs_ACCESS1, presabs_CESM1, presabs_CMCC, presabs_MIROC5, presabs_MPI)
        
        #Save  AUC, TSS, thresholds and predictions of models with AUC and TSS > 0.7
        
        glm_null[[name]][[Fold]]<- list("Ids"=EvaIds,
                                        "AUC_tr"=AUC,
                                        "TSS_tr"=TSS,
                                        "th_tr"=th,
                                        "Prob"=PredProb[,i],
                                        "Orig"=presabs,
                                        "AUC_ev"=AUC_ev,
                                        "TSS_ev"=TSS_ev,
                                        "th_ev"=th_ev
        )
        # write.csv(pred_scenarios, row.names = F, paste0(output_path3,"/pred_ssn_",r,".csv"))
        
        # write.csv(presabs_scenarios, row.names = F, paste0(output_path3,"/presabs_ssn_",r,".csv"))
        
        
        
      }
      
    } , error=function(e){})
  }
}


saveRDS(glm_null, "glm_evaluation_results_or.rds")




##### Train and test in independent datasets
### Load ssn_list for evaluation

ssn_models<- list()

for (name in Knames) {
  for (i in 1:5) {
    tryCatch({
      Fold <- paste0("Fold.", i)  
      
      ssn_p_i <- ssn_object[[name]]
      #Data frame with observation points
      ssn_test_dataDF <- getSSNdata.frame(ssn_p_i, "Obs")
      
      #Ids of observations selected for evaluation
      EvaIds <- ssn_list[[name]][["folds"]][[1]][[i]][["in_id"]]
      
      #Insert NAs in the column of the response variable (pres_abs) in the observations
      #selected for evaluation
      
      print("EvaIds correct")
      ssn_test_dataDF[, paste0(name,Fold)]<- ssn_test_dataDF[,name]
      ssn_test_dataDF[!ssn_test_dataDF$pid %in% EvaIds,paste0(name, Fold)] <- NA
      
      thres <- matrix(nrow = 1, ncol = 5)
      perform <- matrix(nrow = 2, ncol = 5)
      predict_sites <- getSSNdata.frame(ssn_p_i, Name = "Kinzig_PA_evaluation")
      PredProb <- matrix(nrow = nrow(predict_sites), ncol = 6)
      PredProb[,6] <- predict_sites$pid
      
      #Put the data frame with NAs in the SSN object.
      ssn_calib_data <- putSSNdata.frame(ssn_test_dataDF, ssn_p_i)
      
      
      print(paste("Processing", name, "Fold", i))
      
      
      
      frml <- null_AIC_list[[name]]%>%as.character()%>%sub("as\\.factor\\((.*)\\).*", "\\1", .)
      frml <- as.formula(paste(paste0(frml[2],Fold), frml[1], frml[3], collapse = " "))
      crmdl <- CorModels[[name]]%>%as.character() 
      #With this, the regression omits observations with NAs and the model is fitted
      # only on the calibration data set
      options(na.action = "na.omit")
      
      #Fit the model with the calibration data set and parameters
      CalibSsn <- try(glmssn(formula =  frml,
                             family = "binomial",
                             ssn_calib_data,
                             CorModels = crmdl,
                             addfunccol = "computed.afv",
                             control = list(trunc.pseudo=100)))
      
      
      print("Model running")
      
      # # correction bias covariate to a constant value across all sub-catchments
      #Insert 1 in the columns of the covariates used to correct observer bias
      
      eval_temp <- getSSNdata.frame(CalibSsn$ssn.object, "_MissingObs_") 
      
      CalibSsn_cor_eval <- putSSNdata.frame(eval_temp, CalibSsn, "_MissingObs_")
      
      #Predict points used for evaluation with the model trained with the
      #calibration data set
      print("Predicting...")
      
      PredEval <- predict.glmssn(CalibSsn_cor_eval, "_MissingObs_")
      PredProbEval <- SSNProb(PredEval)
      
      #Extract ids and observed presence-absences of points used for evaluation
      obsPresAbs <- getSSNdata.frame(ssn_p_i) %>%
        filter(!pid %in% EvaIds) %>%
        dplyr::select(pid,name)
      
      #Evaluation metrics
      idsPresabsProb <- data.frame(obsPresAbs,PredProbEval)
      TSS <- eval_mod(idsPresabsProb)$performance[5,2]
      AUC <- eval_mod(idsPresabsProb)$performance[6,2]
      th <- eval_mod(idsPresabsProb)$threshold[1,2]
      thres[,i] <- th
      perform[,i] <- c(TSS, AUC)
      
      print(paste("Fold training finished","AUC", AUC, "TSS", TSS))
      
      predPresAbs <- matrix(nrow = nrow(predict_sites), ncol = 6)
      predPresAbs[,6] <- predict_sites$id
      colnames(predPresAbs) <- c( paste0("Fold", 1:5),"pid")
      
      if (AUC >= 0.7 & TSS >= 0.5){
        
        # # correction bias covariate to a constant value across all sub-catchments
        #Insert 1 in the columns of the covariates used to correct observer bias
        
        temp_pred_df <- getSSNdata.frame(CalibSsn, "Kinzig_PA_evaluation")
        
        print("Evaluation starting")
        
        
        CalibSsn_cor_eval <- putSSNdata.frame(temp_pred_df, CalibSsn, "Kinzig_PA_evaluation")
        
        # predict on complete data set
        print("Predicting on Evaluation")
        
        PredAll <- predict.glmssn(CalibSsn_cor_eval, "Kinzig_PA_evaluation")
        
        # Extract probabilities
        PredProb[,i] <- SSNProb(PredAll)
        
        # pred_ACCESS1 <- predict.glmssn(CalibSsn, "ACCESS1")%>%
        #   SSNProb()
        # pred_CESM1 <- predict.glmssn(CalibSsn, "CESM1")%>%
        #   SSNProb()
        # pred_CMCC <- predict.glmssn(CalibSsn, "CMCC")%>%
        #   SSNProb()
        # pred_MIROC5 <- predict.glmssn(CalibSsn, "MIROC5")%>%
        #   SSNProb()
        # pred_MPI <- predict.glmssn(CalibSsn, "MPI")%>%
        #   SSNProb()
        # pred_scenarios <- data.frame(pred_ACCESS1, pred_CESM1, pred_CMCC, pred_MIROC5, pred_MPI)
        
        #Binary transformation
        presabs <- ifelse(PredProb[,i] >= th,1,0)
        predPresAbs[,i] <- presabs
        
        
        obsPresAbsEval <- getSSNdata.frame(ssn_p_i, "Kinzig_PA_evaluation") %>%
          dplyr::select(pid, name)
        
        print("Presence and absences binary processing")
        
        idsPresabsProbEval <- data.frame(obsPresAbsEval,PredProb[,i])
        TSS_ev <- eval_mod(idsPresabsProbEval)$performance[5,2]
        AUC_ev <- eval_mod(idsPresabsProbEval)$performance[6,2]
        th_ev <- eval_mod(idsPresabsProbEval)$threshold[1,2]
        
        
        # presabs_ACCESS1 <- ifelse(pred_ACCESS1 >= th,1,0)
        # presabs_CESM1 <- ifelse(pred_CESM1 >= th,1,0)
        # presabs_CMCC <- ifelse(pred_CMCC >= th,1,0)
        # presabs_MIROC5 <- ifelse(pred_MIROC5 >= th,1,0)
        # presabs_MPI <- ifelse(pred_MPI >= th,1,0)
        # presabs_scenarios <- data.frame(presabs_ACCESS1, presabs_CESM1, presabs_CMCC, presabs_MIROC5, presabs_MPI)
        
        #Save  AUC, TSS, thresholds and predictions of models with AUC and TSS > 0.7
        
        ssn_models[[name]][[Fold]]<- list("Ids"=EvaIds,
                                          "AUC_tr"=AUC,
                                          "TSS_tr"=TSS,
                                          "th_tr"=th,
                                          "Prob"=PredProb[,i],
                                          "Pred_occ"=presabs,
                                          "AUC_ev"=AUC_ev,
                                          "TSS_ev"=TSS_ev,
                                          "th_ev"=th_ev
        )
        # write.csv(pred_scenarios, row.names = F, paste0(output_path3,"/pred_ssn_",r,".csv"))
        
        # write.csv(presabs_scenarios, row.names = F, paste0(output_path3,"/presabs_ssn_",r,".csv"))
        
        
        
      }
      
    } , error=function(e){})
  }
}


saveRDS(ssn_models, "ssn_evaluation_results_or.rds")
