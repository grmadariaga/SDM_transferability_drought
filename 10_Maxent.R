if (!require(dplyr)) install.packages("dplyr")

if (!require(dismo)) install.packages("dismo")

if (!require(biomod2)) install.packages("biomod2")

if (!require(sf)) install.packages("sf")

if (!require(SSN)) install.packages("SSN")

Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8')

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


###Function to select maxent param in https://rvalavi.github.io/Presence-Only-SDM/
maxent_param <- function(data, occ, pred, k = 5){
  
  
  if (!require(dismo)) install.packages("dismo")
  if (!require(caret)) install.packages("caret")
  if (!require(precrec)) install.packages("precrec")
  if (!require(rJava)) install.packages("rJava")
  
  # generate balanced CV folds
  folds <- caret::createFolds(y = as.factor(data[,occ]), k = k)
  
  # regularisation multipliers
  ms <- c(0.5, 1, 2, 3, 4)
  grid <- expand.grid(
    regmult = paste0("betamultiplier=", ms),
    features = list(
      c("noautofeature", "nothreshold"), # LQHP
      c("noautofeature", "nothreshold", "noproduct"), # LQH
      c("noautofeature", "nothreshold", "nohinge", "noproduct"), # LQ
      c("noautofeature", "nothreshold", "nolinear", "noquadratic", "noproduct"), # H
      c("noautofeature", "nothreshold", "noquadratic", "nohinge", "noproduct")), # L
    stringsAsFactors = FALSE
  )
  AUCs <- c()
  for(n in seq_along(grid[,1])){
    full_pred <- data.frame()
    for(i in seq_len(length(folds))){
      trainSet <- unlist(folds[-i])
      testSet <- unlist(folds[i])
      if(inherits(try(
        maxmod <- dismo::maxent(x = data[trainSet, pred],
                                p = data[trainSet,occ],
                                removeDuplicates = FALSE,
                                args = as.character(unlist(grid[n, ]))
        )
      ), "try-error")){
        next
      }
      modpred <- predict(maxmod, data[testSet, pred], args = "outputformat=cloglog")
      pred_df <- data.frame(score = modpred, label = data[testSet,occ])
      full_pred <- rbind(full_pred, pred_df)
    }
    AUCs[n] <- precrec::auc(precrec::evalmod(scores = full_pred$score, 
                                             labels = full_pred$label))[1,4]
  }
  best_param <- as.character(unlist(grid[which.max(AUCs), ]))
  return(best_param)
}

set.seed(123)

maxent_models <- list()
maxent_dfs <- list()

for (name in Knames) {
  time_sp <- system.time({
    tryCatch({
      
      ssn_p_i <- ssn_object[[name]]
      
      print(paste("Modelling ", name))
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
      
      print("prep ready")
      
      
      # Import data frames with sub-catchment ids of calibration and evaluation
      # data sets
      CalDf <- dataset
      EvalDf <- predict_sites
      
      
      for (r in 1:5) {
        tryCatch({
          
          print(paste("Fold", r))
          #Observations used for calibration
          CalIds <- ssn_list[[name]][["folds"]][[1]][[r]][["in_id"]]
          CalibData <- dataset[dataset$pid %in% CalIds,]
          
          #Observations selected for evaluation 
          EvaData <- dataset[!dataset$pid %in% CalIds,]
          
          #Tun Maxent model parameters
          frml <- formula_list[[name]]%>%as.character()
          frml <- unlist(strsplit(gsub("\\s*\\+\\s*1", "", frml[[3]]), "\\s*\\+\\s*"))
          
          
          pred <- frml 
          
          print(pred)
          #response var name
          params <- maxent_param (data = CalibData,
                                  pred = pred,
                                  occ = name,
                                  k = 5)
          
          
          
          #Run model with tuned parameters
          print("Running Maxent")
          
          
          CalibMaxent <- maxent(x = CalibData[, pred], p = CalibData[, name],
                                args = params)
          
          # # correction bias covariate to a constant value across all sub-catchments
          
          
          EvaData_bias_cor <- EvaData
          #Predict on evaluation dataset
          PredProbEval <- predict(CalibMaxent, EvaData_bias_cor[, pred], args = "outputformat=cloglog")
          #Evaluation metrics
          idsPresabsProb <- data.frame(EvaData[,c("id",name)],PredProbEval)
          
          
          TSS_tr <- eval_mod(idsPresabsProb)$performance[5,2]
          AUC_tr <- eval_mod(idsPresabsProb)$performance[6,2]
          th_tr <- eval_mod(idsPresabsProb)$threshold[1,2]
          thres[,r] <- th_tr
          perform[,r] <- c(TSS_tr, AUC_tr)
          
          print(paste("TSS_tr", TSS_tr, "AUC_tr", AUC_tr))
          if (AUC_tr >=0.7 & TSS_tr >= 0.5){
            
            print("Evaluating")
            # Correct for observer bias
            pred_var_bias_cor <- predict_sites  
            
            #Predict on complete dataset
            predicted <- predict(CalibMaxent, pred_var_bias_cor[,pred], args = "outputformat=cloglog")
            PredProb[,r] <- predicted
            
            # pred_ACCESS1 <- predict(CalibMaxent, ACCESS1[,pred], args = "outputformat=cloglog")
            # pred_CESM1 <- predict(CalibMaxent, CESM1[,pred], args = "outputformat=cloglog")
            # pred_CMCC <- predict(CalibMaxent, CMCC[,pred], args = "outputformat=cloglog")
            # pred_MIROC5 <- predict(CalibMaxent, MIROC5[,pred], args = "outputformat=cloglog")
            # pred_MPI <- predict(CalibMaxent, MPI[,pred], args = "outputformat=cloglog")
            # pred_scenarios <- data.frame(pred_ACCESS1, pred_CESM1, pred_CMCC, pred_MIROC5, pred_MPI)
            print("Binary transformation")
            #Binary transformation 
            presabs <- ifelse(predicted >= th_tr,1,0)
            predPresAbs[,r] <- presabs 
            
            obsPresAbsEval <- predict_sites %>%
              dplyr::select(pid,name)
            
            print("Presence and absences binary processing")
            
            idsPresabsProbEval <- data.frame(obsPresAbsEval,predicted)
            TSS_ev <- eval_mod(idsPresabsProbEval)$performance[5,2]
            AUC_ev <- eval_mod(idsPresabsProbEval)$performance[6,2]
            th_ev <- eval_mod(idsPresabsProbEval)$threshold[1,2]
            
            print(paste("TSS_EV", TSS_ev))
            print(paste("AUC_EV", AUC_ev))
            
            
            maxent_models[[name]][[r]]<- list("Ids"=CalIds,
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
            
            #Save  AUC, TSS, thresholds and predictions of models with AUC and TSS > 0.7
            
            #write.csv(EvaIds, row.names = F, paste0(output_path,"/eva_max_",r,".csv"))
            #write.csv(AUC, row.names = F, paste0(output_path,"/auc_max_",r,".csv"))
            #write.csv(TSS, row.names = F, paste0(output_path,"/tss_max_",r,".csv"))
            #write.csv(th, row.names = F, paste0(output_path,"/th_max_",r,".csv"))
            #write.csv(pred, row.names = F, paste0(output_path,"/pred_max_",r,".csv"))
            #write.csv(presabs, row.names = F, paste0(output_path,"/presabs_max_",r,".csv"))
            #write.csv(pred_scenarios, row.names = F, paste0(output_path3,"/pred_max_",r,".csv"))
            #write.csv(presabs_scenarios, row.names = F, paste0(output_path3,"/presabs_max_",r,".csv"))
            
          }
          
          
          
        }, error=function(e){})  
      }
      #Export results
      # print(paste0("Export results for ", files_names[sp]))
      # outputpath2 <- paste0("./FBAC/data/SDM/output/maxent/",files_names[sp])
      # ifelse(!dir.exists(outputpath2), dir.create(outputpath2), FALSE)
      # write.csv(thres, row.names = F, paste0("./FBAC/data/SDM/output/maxent/",files_names[sp], "/thres.csv"))
      # write.csv(perform, row.names = F, paste0("./FBAC/data/SDM/output/maxent/",files_names[sp], "/perform.csv"))
      # write.csv(PredProb, row.names = F, paste0("./FBAC/data/SDM/output/maxent/",files_names[sp], "/pred.csv"))
      
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
      
      maxent_dfs[[name]] <- prediction_df
      
      print("done")
      
    }, error=function(e){})
  })
}

saveRDS(maxent_models, "maxent_models_or.rds")
saveRDS(maxent_dfs, "maxent_dfs_or.rds")
