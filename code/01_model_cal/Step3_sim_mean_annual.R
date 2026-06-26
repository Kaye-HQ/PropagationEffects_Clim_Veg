
rm(list=ls())
library(data.table)
library(lubridate)
library(doParallel)

source("01_model_cal/functions.R")

exps <- c('SimulationRaw', 'State_LAI', 'State_P', 'State_T', 'State_P_T_LAI')

pathdata0 <- "../Data/sample_data/"
LAIpath <- "../Data/sample_data/LAI/"

infos <- fread("../Data/info_2252.csv")

folderPath <-  "../Result/HBV_PML/02_cal_log/"
fileList <- list.files(folderPath, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)
Name_all_str <- gsub("\\.RData$", "", basename(fileList))
csv_names_str <- infos$Name
Name_all_com <- intersect(Name_all_str, csv_names_str)


no_cores <- 1 # 2#detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

for (iexp in seq_along(exps)) {
  exp <- exps[iexp]
  
  output_main_path <- "../Result/Result_Exps_detrend/Exp_data/"
  sim_path <- sprintf("%s/01_Simulate/%02d_%s/", output_main_path,iexp, exp)
  dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
  
  output_path <- paste0(sim_path, "output/"); dir.create(output_path, showWarnings = FALSE)
  input_path <- paste0(sim_path, "input/"); dir.create(input_path, showWarnings = FALSE)
  
  Done_files <- list.files(input_path, pattern = "*.csv")
  Done_Name <- gsub(".csv", "", Done_files)
  
  if (length(Done_Name) >= 3) {
      names_to_run <- setdiff(Name_all_com, Done_Name[3:length(Done_Name)])
  } else {
      names_to_run <- Name_all_com
  }
  
  fileList_in <- fileList[basename(fileList) %in% paste0(names_to_run, ".RData")]
  

  foreach(i = seq_along(fileList_in), 
          .packages = c("data.table", "lubridate"),
          .export = c("parread_exps_mean_annual", "parread", "parrun_HBV_PML", 
                      "TraToPar_HBV_PML", "prelocate", "snow_routine", "Run_hbv_PML", 
                      "mytrimf", "get_Route", "userf", "methods_evaluate", "detrend_Data")) %dopar% {
                        
   filePath <- fileList_in[i]
    basin_name <- gsub("\\.RData$", "", basename(filePath))
    FileName <- paste0(basin_name, ".csv")
    
    load(filePath)
    
    ifRoute <- 1; ifPML <- 1
    
    basin_info_n <- infos[Name %in% basin_name,]
    res_read <- parread_exps_mean_annual(pathdata0, basin_name, basin_info_n, LAIpath, exp)
    INPUTS <- res_read$INPUTS
    
    missingLAI <- is.na(INPUTS$in_df$LAI_GLASS)
    INPUTS$in_df$LAI_GLASS[missingLAI] <- mean(INPUTS$in_df$LAI_GLASS[!missingLAI], na.rm=TRUE)
    
    missingEmiss <- is.na(INPUTS$in_df$Emiss)
    if("Emiss" %in% names(INPUTS$in_df)) INPUTS$in_df$Emiss[missingEmiss] <- mean(INPUTS$in_df$Emiss[!missingEmiss], na.rm=TRUE)
    
    missingAlbedo <- is.na(INPUTS$in_df$Albedo)
    if("Albedo" %in% names(INPUTS$in_df)) INPUTS$in_df$Albedo[missingAlbedo] <- mean(INPUTS$in_df$Albedo[!missingAlbedo], na.rm=TRUE)
    
    if (length((INPUTS$Date)) > 0) {
      v_res <- parrun_HBV_PML(INPUTS, Best1_par, 'v', ifRoute, ifPML, list())
      output <- v_res$output
      output$Name <- basin_name
      output$exp <- exp
      
      fwrite(output, paste0(output_path, FileName))
      fwrite(INPUTS$in_df, paste0(input_path, FileName))
      print(basin_name)
    }
  }
}

stopCluster(cl)