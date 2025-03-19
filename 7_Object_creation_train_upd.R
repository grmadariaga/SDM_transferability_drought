library(tidymodels)
library(SSN)
#library(sirad)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(MuMIn)
library(fuzzySim)
library(mecofun)




###Merge SSN object from STARS to dataframe on environmental variables

###Load SSN object Kinzig
###Create a variable with the name of the folderwhere the files were downloaded

Myfiles<- "Your/folder/"

setwd(Myfiles)

K <- importSSN( "train.ssn", predpts = "Kinzig_loc")


#createDistMat(K, predpts="Kinzig_loc",  o.write = TRUE)
#distPred1km <- getStreamDistMat(K, Name = "Kinzig_loc")

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
  dplyr::select(-Year)

####

Kinzig_all <- left_join(
  Kinzig_bio, 
  IHA_cut[,c(2:175)],
  by= c("id", "wat_year"))%>% 
  dplyr::select_if(~ !any(is.na(.)))%>%  ####Remove columns with NA
  dplyr::select_if(~ !all(. == .[[1]]))%>%###Remove columns with the same value for all sites
  dplyr::select(where(~ !any(is.infinite(.))))%>%### Remove columns with inf
  dplyr::select(-wat_year)


Kinzig_all$ma11<- ifelse(Kinzig_all$ma11<0, 0, Kinzig_all$ma11)
Kinzig_all$ma10<- ifelse(Kinzig_all$ma10<0, 0, Kinzig_all$ma10)
###Check variable colinearity within each variable categories


####Create subsets depending on categories
duration <- Kinzig_all%>%
  dplyr::select(starts_with("d") & matches("[[:digit:]]"))

magnitude <- Kinzig_all%>%
  dplyr::select(starts_with("m") & matches("[[:digit:]]"))

timing <- Kinzig_all%>%
  dplyr::select(starts_with("t") & matches("[[:digit:]]"))

rate <- Kinzig_all%>%
  dplyr::select(starts_with("r") & matches("[[:digit:]]"))

frequency <- Kinzig_all%>%
  dplyr::select(starts_with("f") & matches("[[:digit:]]"))

bioclim <-Kinzig_all%>%
  dplyr::select(starts_with("bio") & matches("[[:digit:]]"))

luse <-Kinzig_all%>%
  dplyr::select(starts_with("LU") & matches("[[:digit:]]"))



###Remove correlated variables duration

dur_red <- corSelect(duration,coeff=TRUE,cor.thresh = 0.7, var.cols = names(duration))

dur_red_ex <- dur_red$excluded.vars


Kinzig_prep <- Kinzig_all[, !(colnames(Kinzig_all) %in% dur_red_ex)]


###Remove correlated variables magnitude


mag_red<- corSelect(magnitude,coeff=TRUE,cor.thresh = 0.7, var.cols = names(magnitude))

mag_red_ex <- mag_red$excluded.vars


Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% mag_red_ex)]


##### Remove correlated variables frequency
freq_red<- corSelect(frequency,coeff=TRUE,cor.thresh = 0.7, var.cols = names(frequency))

freq_red_ex <- freq_red$excluded.vars


Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% freq_red_ex)]


####Remove correlated variables rate
rate_red<- corSelect(rate,coeff=TRUE,cor.thresh = 0.7, var.cols = names(rate))

rate_red_ex <- rate_red$excluded.vars


Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% rate_red_ex)]

#### timing only has 3 variables
timing_red<- corSelect(timing,coeff=TRUE,cor.thresh = 0.7, var.cols = names(timing))

timing_red_ex <- timing_red$excluded.vars


Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% timing_red_ex)]


###Timing is not autocorrelated, we take all


#### Remove correlated variables bioclim

bioclim_red<- corSelect(bioclim,coeff=TRUE,cor.thresh = 0.7, var.cols = names(bioclim))

bioclim_red_ex <- bioclim_red$excluded.vars


Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% bioclim_red_ex)]

####Remove correlate variables Land use

luse_red<- corSelect(luse,coeff=TRUE,cor.thresh = 0.7, var.cols = names(luse))

luse_red_ex <- luse_red$excluded.vars


Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% luse_red_ex)]




###We remove land use 5 as it's others
LU5_vec <- "LU5"

Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% LU5_vec)]

Samp <- "sampling"
Kinzig_prep <- Kinzig_prep[,!(colnames(Kinzig_prep)%in% Samp)]

####Select variables from Kakouei et al. 2018 available

Kvars <- Kinzig_all[,c("id", "dh4")]

Kinzig_prep<- left_join(Kinzig_prep,Kvars, by="id")


##### Use select 07 to select variables

Knames <- names(filtered_cols)

#saveRDS(Knames, "Knames.rds")

selected_vars <- list()


for (i in Knames) {
  vars<- Kinzig_prep[, c(20:24, 68:126)]
  occ<- Kinzig_prep[[i]]
  
  covar_sel <- select07(X= vars, y=occ, threshold= 0.7, univar= "glm1")
  
  selected_vars[[i]]<- covar_sel
}

final_sel <- list()

for (i in Knames){
  
  s <- sum(Kinzig_prep[[i]])
  
  z <- ceiling(s/10)
  
  variables <- selected_vars[[i]][["pred_sel"]]
  
  final_sel[[i]]<- variables[1:z]
  
  
}


Formula_list <- list()


for (i in Knames){
  
  sp <- i
  
  preds <- final_sel[[i]]
  
  formula <- as.formula(paste(paste0("as.factor(",i,") ~"), paste(preds, collapse = "+"), "+ 1"))
  
  Formula_list[[i]]<- formula
  
}


#readRDS("formulas.rds")

saveRDS(Formula_list, "formulas.rds")

####

obs_df <- Kinzig_prep

####Check for correlation in all variables now

#All_corr <- corSelect(Kinzig_prep,coeff=TRUE,cor.thresh = 0.6, var.cols = colnames(Kinzig_prep[,c(20:23,65:108)])) ###Retain elev
  
#All_corr_ex<- All_corr$excluded.vars
  

#Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% All_corr_ex)]

####dredge accepts a maximum of 30 variables, so we have to reduce it a little bit more.
###we remove further the 4 remaining variables with highest VIF

#extra_filter <- All_corr$remaining.multicollinearity%>%filter(VIF>20)%>%row.names()

#Kinzig_prep <- Kinzig_prep[, !(colnames(Kinzig_prep) %in% extra_filter)]






###Return Channel48 to dataframe

#ch48cols <- colnames(Kinzig_prep)

#Channel48<- Channel48[,ch48cols]

#Kinzig_prep <- rbind(Kinzig_prep,Channel48)%>%
#  arrange(.,pid)

#Kinzig_prep <- Kinzig_prep[, c(40:84, 1:39)]


###Extract species names to create separate dataframes (evaluation sets)




new_columns <- colnames(Kinzig_prep)

saveRDS(new_columns, "new_cols.rds")

#saveRDS()


####Extract names to an object to avoid for rerunning this process

saveRDS(Knames, "Knames.rds")


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
  subset <- obs_df[, c(names(obs_df[,c(1:25, 68:126)]), name)]
  
  
  # Create a new SSN object with the subset
  ssn_obj <- putSSNdata.frame(subset, K)  # Replace K with the appropriate value
  
  # Assign the SSN object to the list with the corresponding name
  ssn_objects[[name]] <- ssn_obj
}


saveRDS(ssn_objects, "ssn_objects.rds")

##### Now we select the best variables for each species
####Best model selection for null models

#var_names <- colnames(obs_df[,c(20:23,64:89)])

#null_AIC_list <- list()

#nCores <- detectCores()-1
#clus <- makeCluster(nCores, type= "PSOCK", outfile=""
#                       )
#parallel::mcaffinity(1:15)
#doParallel::registerDoParallel(clus)

#options(na.action = "na.fail")
# Loop through each response variable
#for (i in Knames) {
  
 # ssn_i <- ssn_objects[[i]]
  #df_ssn_i<- getSSNdata.frame(ssn_i,"Obs")
  #df_ssn_i <- df_ssn_i
  
  # Initialize variables to keep track of the best model and its AIC
  #best_model_formula <- NULL
  #lowest_aic <- Inf
  
  # Create the formula for the current combination
  #formula <- as.formula(paste("as.factor(",i,")", "~", paste(var_names, collapse = " + ")))
  
  # Fit the model
  #model_1 <- try(glm(formula = formula, 
                     #data = df_ssn_i,
                    # family = binomial,
                     #na.action = "na.fail")) 
  
  #print("here")
  ###create cluster for dredge
  
 #clusterExport(clus, "df_ssn_i")
  
  #ListModelNsp <- MuMIn::dredge(model_1, m.lim = c(0,5), cluster = clus, trace=2)
  
 # BestModelNsp <- eval(attributes(ListModelNsp)$model.calls[[1]])
  
 # forml <- BestModelNsp$formula
  
  
  
  # Save the best model formula with the lowest AIC for the current response variable
 # null_AIC_list[[i]] <- list("formula" = forml)
#}



#stopCluster(clus)


#saveRDS(null_AIC_list, "new_dredge_formulas.rds")


#null_AIC_list <- readRDS("new_dredge_formulas.rds")


# Create an empty list to store the SSN objects
ssn_list <- list()

set.seed(432)

# Create 5 balanced folds for each species according to specific prevalence, this folds are going to be used to create five submodels to be tested on test df
for (name in Knames) {
  
  subset <- obs_df[, c(names(obs_df[,c(1:25, 68:126)]), name)]
  subset_folded <- vfold_cv(subset,v=5, repeats = 1,
                            strata= name)
  
  
  ssn_obj <- putSSNdata.frame(subset, K)  
  
  # Assign the SSN object to the list with the corresponding name
  ssn_list[[name]][["obj"]] <- ssn_obj
  ssn_list[[name]][["folds"]]<-subset_folded
}

saveRDS(ssn_list, "ssn_list_train.rds")

####load all possible combinations of CorMdls

load("CorMdls.RData")

CorMdls <-CorMdls[1:179] ###Remove NULL

######Select spatial model for each species

AIC_ssn_select <- list()


###remove previous cluster to avoid conflicts
#rm(cluster)

#### First run glmssn without spatial model included, in the next steps we will use the residuals of those models to find a spatial model
glm_null_models <- list()
# Loop through each name in Knames
for (name in Knames) {
  
  tryCatch(
    {
      
      # Get the ssn object from the ssn_list
      ssn_obj <- ssn_objects[[name]]
      ssn_df_sel <- getSSNdata.frame(ssn_obj)
      
      frml <- Formula_list[[name]]
      # Run the glmssn model without spatial model
      model_2 <- try(glm(formula = frml, 
                         data = ssn_df_sel,
                         family = binomial,
                         na.action = "na.exclude")) 
      
      # Save the model in the glm_null_models list with the name
      glm_null_models[[name]] <- model_2},  
    error =function(e) {
      # Store the error message in the list with the name from Knames
      glm_null_models[[name]] <<- paste("Failed model:", name, conditionMessage(e))
    }
  )
}


#saveRDS(glm_null_models, file = "SDM_Kinzig_updated/glm_null_models.rds")
#glm_null_models <- readRDS("Kinzig_upd/Results/glm_null_models2011.rds")
####Create evaluation results for each species with null_models


#### Function to unregister_dopar() in each loop

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

unregister_dopar()


######Select spatial model for each species

AIC_ssn_select <- list()


for (k in Knames) {
  tryCatch(
    {
      
      ###Select null model for the species
      start<-Sys.time()
      print(paste("Start time", start))
      i<- k
      print(i)
      mod <- glm_null_models[[i]]
      print("mod ok")
      SSNobj <- ssn_objects[[i]]
      print("SSNobj ok")
      ###add residuals from null model to ssn object
      
      SSNobj@obspoints@SSNPoints[[1]]@point.data$RES <- resid(mod)
      print("Selecting a spatial autocorrelation model...")
      
      unregister_dopar()
      #parallel::mcaffinity(1:25)
      cl <- makeCluster(5)
      doParallel::registerDoParallel(cl)
      getDoParWorkers()
      print("ok2")
      
      ssn_cor_test <- foreach (m = 1:length(CorMdls), .packages = "SSN", .errorhandling="remove") %dopar%  {
        print(paste(eval(CorMdls[[m]])))
        
        start1<-Sys.time()
        print <- start1
        glmssn(RES ~ 1, SSNobj, CorModels= eval(CorMdls[[m]]),
               addfunccol = "computed.afv")
      }
      parallel::stopCluster(cl)
      unregister_dopar()
      
      gc()
      print("Performing model comparison")
      # #Models AIC, selects model with lower AIC describing the residuals
      cor_modl_AIC <- InfoCritCompare(keep(ssn_cor_test, is.list))
      model <- eval(CorMdls[[which.min(cor_modl_AIC$AIC)]])
      AIC_ssn_select[[i]] <- model
      End_time<- Sys.time()
      print(paste("End time", End_time))
      print(paste("Duration", End_time-start))
      #saveRDS(AIC_ssn_select, "ssn_selection.rds")
    },
    error = function(err) {
      message(paste("Error occurred for:", i))
      AIC_ssn_select[[i]] <<- paste("Failed:", name, conditionMessage(err))  
    }
  )
  saveRDS(AIC_ssn_select, file = "ssn_selection.rds")
}



gc()

