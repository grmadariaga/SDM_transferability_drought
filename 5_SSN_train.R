###SSN_Kinzig

## Loading required package: tmap
library(tmap)
library(tmaptools)
library(raster)
library(terra)


### Load openSTARS
library(openSTARS)

Sys.setenv(PATH=paste0("C:/Program Files/QGIS 3.28.7/bin",";", Sys.getenv("PATH")))

grass_program_path <- "C:/Program Files/QGIS 3.28.7/apps/grass/grass82/"

working_dir <- file.path("C:/") ##Your grass path file
grass_db_path <- file.path(working_dir)
dir.create(working_dir)
setwd(working_dir)

###Create a variable with the name of the folderwhere the files were downloaded

Myfiles<- "Your/folder/"

dem_path <- paste0(Myfiles,"Kinzig_filled.tif")



setup_grass_environment(dem = dem_path,
                        gisBase = grass_program_path,
                        gisDbase = grass_db_path,
                        location = "SDM_EE_train", ###decide your location
                        remove_GISRC = TRUE,
                        override=TRUE)

Sys.setenv(GRASS_PYTHON="C:/Program Files/QGIS 3.28.7/bin/python.exe")
Sys.setenv(PYTHONHOME="C:/Program Files/QGIS 3.28.7/apps/Python39")
Sys.setenv(GDAL_DATA="C:/Program Files/QGIS 3.28.7/apps/share/gdal")





gmeta()
#### Paths for predictor rasters
Elevation<- paste0(Myfiles,"Kinzig_filled.tif")

#Paths for sampling points
sites_path <- paste0(Myfiles,"Kinzig_PA_train.gpkg")

#Paths for prediction points
preds_path <- paste0(Myfiles,"Kinzig_loc.gpkg")

#Path stream network
streams_path <- paste0(Myfiles,"rivs1.gpkg")

LU <- paste0(Myfiles,"NNZART_2016.tif")

import_data(dem = dem_path, sites = sites_path, streams = streams_path,
            predictor_raster = c(Elevation,LU), predictor_r_names = c("Elevation","LU"),
            pred_sites = preds_path )

dem <- readRAST("dem", ignore.stderr = TRUE)
sites_o <- readVECT("sites_o", ignore.stderr = TRUE)
streams_o <- readVECT("streams_o", ignore.stderr = TRUE)
Elevation <- readRAST("Elevation", ignore.stderr = TRUE)
#biom <- readRAST("bioms", ignore.stderr = TRUE)
preds_o <- readVECT("Kinzig_loc_o", ignore.stderr = TRUE)
LU_r <- readRAST("LU", ignore.stderr = TRUE)

sites_o$col <- 1
streams_o$col <- 1


tmap_mode("plot")
tm_shape(LU_r) +
  tm_raster("LU", palette = terrain.colors(n = 5), n = 5,
            title = "digital elevation \nmodel [m asl]", style = "cont") +
  tm_shape(streams_o) +
  tm_lines(col = "col", lwd = 1.5, legend.col.show = TRUE, palette = "blue4",
           title.col = "", legend.format = list(fun=function(x) "input streams")) +
  tm_shape(sites_o) +
  tm_symbols(col = "col", palette = "red", size = 0.4,
             legend.col.show = TRUE, title.col = "",
             legend.format = list(fun=function(x) "input sampling sites")) +
  tm_view(symbol.size.fixed = TRUE) + # helpful for tmap_mode('view')
  tm_layout(scale = 1, legend.bg.color = T, legend.position = c("left", "bottom"))

#tm_shape(pw) +
#  tm_polygons(col = "tmx_mean", palette = terrain.colors(n=15),
#              title = "Max temp", n = 8) +
#  tm_layout(scale = 1, legend.bg.color = T, legend.position = c("left", "bottom"))


derive_streams(accum_threshold = 2400, condition = TRUE, clean = TRUE, burn = 10)


streams <- readVECT("streams_v", ignore.stderr = TRUE)
streams$col <- 1

#tm_shape(dem) +
#  tm_raster("dem", palette = terrain.colors(n = 9), n = 9,
#            title = "digital elevation \nmodel [m asl]", style = "cont") +
#  tm_shape(streams_o) +
#  tm_lines(col = "col", lwd = 1, legend.col.show = TRUE, palette = "blue4",
#           title.col = "", legend.format = list(fun=function(x) "input streams")) +
#  tm_shape(streams) +
#  tm_lines(col = "col", palette = "dodgerblue", legend.col.show = TRUE, lwd = 2,
#           title.col = "", legend.format = list(fun=function(x) "derived streams")) +
#  tm_layout(scale = 1, legend.bg.color = T, legend.position = c("left", "bottom"))


cp <- check_compl_confluences()
if (cp) correct_compl_confluences()

streams <- readVECT('streams_v', ignore.stderr = TRUE)
#streams_orig <- readVECT('streams_v_o3', ignore.stderr = TRUE)
streams$col <- 1
#streams_orig$col <- 1

#tmap_mode("view")
#bb <- cbind(c(339723.7, 5720582), c(364103.7, 5708772))

#tm_shape(streams_orig, bbox = bb, unit = "m") +
#  tm_lines(col = "col", lwd = 2, legend.col.show = TRUE, palette = "dodgerblue",
#           title.col = "", legend.format = list(fun=function(x) "original streams")) +
#  tm_shape(streams) +
#  tm_lines(col = "col", lwd = 3, lty = 2, palette = "darkblue", legend.col.show = TRUE,
#           title.col = "", legend.format = list(fun=function(x) "corrected streams")) +
#  tm_scale_bar()

#streams@data$changed.str <- ifelse(streams@data$changed == 1, "changed", "unchanged")

#tmap_mode("view")
#tm_shape(streams) +
#  tm_lines(col = "changed.str", lwd = 2, legend.col.show = TRUE,
#           palette = c("red", "grey"),
#           title.col = "changed segments") +
#  tm_shape(streams_orig) +
#  tm_lines(col = "col", lwd = 1, lty = 2, legend.col.show = TRUE, palette = "dodgerblue",
#           title.col = "", legend.format = list(fun=function(x) "before correction")) +
#  tm_basemap(NULL) +
#  tm_layout(scale = 1, legend.bg.color = T) +
#  tm_scale_bar()

calc_edges()

edges <- readVECT("edges", ignore.stderr = TRUE)
head(edges@data)

#calc_sites()


#tmap_mode("plot")
## tmap mode set to plotting
#tm_shape(edges) +
#  tm_lines(col = "col", palette = "darkblue", legend.col.show = TRUE,
#           title.col = "", legend.format = list(fun=function(x) "edges")) +
#  tm_shape(sites_o) +
#  tm_symbols(col = "col", palette = "grey35", size = 0.4, legend.col.show = TRUE,
#             legend.format = list(fun=function(x) "input sampling sites"),
#             title.col = "") +
#  tm_shape(sites) +
#  tm_symbols(col = "col", palette = "red", size = 0.2, legend.col.show = TRUE,
#             title.col = "", legend.format = list(fun=function(x) "snapped sites")) +
#  tm_layout(scale = 1, legend.bg.color = T, legend.position = c("left", "bottom"))

calc_sites(predictions = "Kinzig_loc_o")
pred_sites <- readVECT("Kinzig_loc", ignore.stderr = TRUE)
sites <- readVECT("sites", ignore.stderr = TRUE)
head(sites@data, n = 4)
sites$col <- 1
edges$col <- 1
head(pred_sites@data, n = 4)

#tmap_mode("plot")
## tmap mode set to plotting
#tm_shape(edges) +
#  tm_lines(col = "col", palette = "darkblue", legend.col.show = TRUE,
#           title.col = "", legend.format = list(fun=function(x) "edges")) +
#  tm_shape(sites_o) +
#  tm_symbols(col = "col", palette = "grey35", size = 0.4, legend.col.show = TRUE,
#             legend.format = list(fun=function(x) "input sampling sites"),
#             title.col = "") +
#  tm_shape(sites) +
#  tm_symbols(col = "col", palette = "red", size = 0.2, legend.col.show = TRUE,
#             title.col = "", legend.format = list(fun=function(x) "snapped sites")) +
#  tm_layout(scale = 1, legend.bg.color = T, legend.position = c("left", "bottom"))

#tmap_mode("view")
#tm_shape(edges) +
#  tm_lines(col = "col", palette = "darkblue", legend.col.show = TRUE,
#           title.col = "", legend.format = list(fun=function(x) "edges")) +
#  tm_shape(sites_o) +
#  tm_symbols(col = "col", palette = "grey35", size = 0.08, legend.col.show = TRUE,
#             title.col = "", legend.format = list(fun=function(x) "input sampling sites")) +
#  tm_shape(sites) +
#  tm_symbols(col = "col", palette = "red", size = 0.04, legend.col.show = TRUE,
#             legend.format = list(fun=function(x) "snapped sites"),
#             title.col = "") +
#  tm_basemap(NULL) +
#  tm_layout(scale = 1, legend.bg.color = T) +
#  tm_scale_bar() +
#  tm_view(symbol.size.fixed = T)



calc_attributes_edges(input_raster =c("LU","Elevation"),
                      stat_rast = c("percent", "mean"),
                      attr_name_rast= c("LU","Elev"), round_dig = 2)
calc_attributes_sites_approx(sites_map = "sites",
                             input_attr_name = c('LUp_1_1_015686', 'LUp_1_988235_2_003922', 'LUp_2_992157_3_007843', 'LUp_3_996078_4_011765', 'LUp_4_984314_5', "Elev"),
                             output_attr_name = c(paste0(rep("LU",5),1:5), "Elev"),
                             stat = c(rep("percent",5), "mean"), overwrite = TRUE)
calc_attributes_sites_approx(sites_map = "Kinzig_loc",
                             input_attr_name = c('LUp_1_1_015686', 'LUp_1_988235_2_003922', 'LUp_2_992157_3_007843', 'LUp_3_996078_4_011765', 'LUp_4_984314_5', "Elev"),
                             output_attr_name = c(paste0(rep("LU",5),1:5), "Elev"),
                             stat = c(rep("percent",5), "mean"), overwrite = TRUE)
sites <- readVECT("sites", ignore.stderr = TRUE)
pred_sites<-readVECT("Kinzig_loc", ignore.stderr = TRUE)
head(sites@data, n = 4)
head(pred_sites@data, n = 4)

calc_edges()
edges <- readVECT("edges", ignore.stderr = TRUE)
head(edges@data)

dirs <- readRAST("dirs", ignore.stderr = TRUE)

#calc_attributes_sites_ap(sites_map = "sites",
#                           input_raster = "Elevation",
#                          stat_rast = "mean",
#                         attr_name_rast = "Elv",
#                        round_dig = 4,
#                       overwrite = TRUE)

#### Do not RUN takes too long better to use approximates
#####calc_attributes_sites_exact(sites_map = "Kinzig_loc",
#                           input_raster = "Elevation",
#                          stat_rast = "mean",
#                         attr_name_rast = "Elv",
#                        round_dig = 4,
#                       overwrite = TRUE)


#sites <- readVECT("sites", ignore.stderr = TRUE)
#pred_sites <- readVECT("Kinzig_loc", ignore.stderr = TRUE)

#head(sites@data, n = 4)
#head(pred_sites@data)


execGRASS("g.copy",
          parameters = list(
            vector = "sites,sites_s"
          ))

ssn_dir_kinzig <- file.path(Myfiles, "train.ssn")
export_ssn(ssn_dir_kinzig, predictions = "Kinzig_loc", delete_directory = TRUE)
