# Macroinvertebrate SDM Transferability to Drought Conditions

## Overview
This repository contains scripts and data processing functions for the manuscript "Transferability of stream benthic macroinvertebrate distribution models to drought-related conditions". The provided R scripts facilitate data preparation, model training, evaluation, and analysis. The external files can be found in (https://fred.igb-berlin.de/data/package/977). The codes must be executed in the correct order indicated by the number on the file name.

## Repository Structure

### **1. Data Processing & Preparation**
- **`1_STARSiha.R`** – Computes STAR statistics needed to calculate IHAs
- **`2_IHA_processing.R`** – Processes Indices of Hydrologic Alteration (IHA) from flow data.
- **`3_BioclimExtraction.R`** – Extracts bioclimatic variables for their use in SDM.
- **`4_Data_organization_upd.R`** – Organizes and structures the full dataset for analysis.
- **`5_SSN_train.R`** – Creates .ssn files used for the fitting of SSN SDM models
- **`6_SSN_eval.R`** – Creates .ssn files used for the fexternal evaluation of SSN SDM models
- **`7_Object_creation_train_upd.R`** – Creates training dataset objects.
- **`8_Object_creation_eval_upd.R`** – Creates evaluation dataset objects.

### **2. Model Training & Evaluation**
- **`9_SSN_model.R`** – Implements the SSN model framework for all macroinvertebrate species.
- **`10_Maxent.R`** – Runs species distribution models for all species using MaxEnt.
- **`11_RF.R`** – Implements a Random Forest species distribution model for each species.


### **3. Model Analysis & Interpretation**
- **`12_Dataframe_models.R`** – Structures the model outputs into dataframes for further analysis.
- **`13_Model_analysis.R`** – Analyzes and compares model results.

