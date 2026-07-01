rm(list=ls())
library(data.table)
library(GA)
library(Rcpp)

source("./01_model_cal/functions.R")
Rcpp::sourceCpp("./01_model_cal/HBV_PML.cpp") 

basin_all  <- fread("../Data/info_2252.csv")
basin_info <- basin_all[1,]# here we only run the example basin

col_mapping <- list(
  lai_essential = c("Name", "Date", "LAI_GLASS", "Emiss", "Albedo"),
  met_essential = c("Name", "Date", "Flow", "P_MSWEP", "Temp", "Tmax", "Tmin", 
                    "SWd", "LWd", "Pres", "Wind", "CO2", "SpecHum")
)

# Units of inputs:
# P_MSWEP           mm/d
# Temp, Tmax, Tmin  degree celecous,   2-m air mean/max/min temperature
# Pres              Pa,                surface pressure
# CO2               umol/mol
# Wind              m/s,               wind speed at 10m height
# SWd LWd           W/m^2,             Downward shortwave & longwave radiation
# SpecHum           g/g,               2-m specific humidity
# PET_FAO           mm/d               potential evapotranspiration 
# Albedo            -
# Emiss             -
# LAI_GLASS         m2/m2

pathdata0 <- "../Data/sample_data/"
LAIpath <- "../Data/sample_data/LAI/"

lb <- c(0, 0, 0, 0, 0.5010, 0.1000, 0.6800, 0, 0, 0.6570, 3.5000, 0.0100, -5, 0, 0, 0, 0, 0.3, 1, 0, 0.05, 0.01, 0.02, 0, 1)
ub <- c(0.0480, 0.0690, 25.7550, 46.3750, 2.0000, 1.0000, 0.8990, 0.1740, 0.1580, 1.4920, 6.5000, 10.0000, 5, 20, 1, 0.8, 7, 1, 2000, 100, 2, 1, 0.1, 300, 6)
numberOfVariables <- 25; ifRoute <- 1; ifPML <- 1

schemes <- list(
  scheme1 = list(
    folder_name = "01_cal_daily",
    target_f = list(
      list("Daily", "KGE", 1.0)
    )
  ),
  scheme2 = list(
    folder_name = "02_cal_log",
    target_f = list(
      list("Daily", "KGE", 0.5),
      list("Q_90", "KGE", 0.1),
      list("Q_log_10", "Abias", 0.3),
      list("LFD_Q10", "Abias", 0.1)
    )
  )
)


for (s_name in names(schemes)) {
  current_scheme <- schemes[[s_name]]
  target_f <- current_scheme$target_f
 
  destination <- paste0("../Result/HBV_PML/", current_scheme$folder_name, "/")
  dir.create(destination, showWarnings = FALSE, recursive = TRUE)
  sapply(c("mat/", "output/", "Evaluation/"), 
         function(x) dir.create(paste0(destination, x), showWarnings = FALSE))
  
  print(paste("================ Starting Scheme:", current_scheme$folder_name, "================"))
  
  for (ibasin in 1:nrow(basin_info)) {
    basin_name <- basin_info$Name[ibasin]
    basin_info_n <- basin_info[ibasin, ]
    
    res_read <- parread(pathdata0, basin_name, basin_info_n, LAIpath, "SimulationRaw",col_map = col_mapping)
    INPUTS <- res_read$INPUTS
    
    if (length(INPUTS$Date) > 0) {
      print(paste0("[", current_scheme$folder_name, "] ", basin_name, "_start"))
      
      fitness_func <- function(par) {
        score <- parrun_HBV_PML_Rcpp(INPUTS, par, 'c', ifRoute, ifPML, target_f)
        if(is.na(score) || is.infinite(score)) return(-9999)
        return(-score) 
      }
      
      ptm <- proc.time()
      
      ga_res <- ga(type = "real-valued", 
                   fitness = fitness_func, 
                   lower = lb, upper = ub, 
                   popSize = 100, maxiter = 50,
                   run = 50,
                   monitor = TRUE)
      
      pop <- ga_res@population
      fits <- ga_res@fitness
      top1_idx <- order(fits, decreasing = TRUE)[1]
      Best_par <- pop[top1_idx, , drop = FALSE]
      
      v_res <- parrun_HBV_PML_Rcpp(INPUTS, as.numeric(Best_par), 'v', ifRoute, ifPML, target_f)
      
      output <- v_res$output
      Evalue <- v_res$Evalue
      
      save(Best_par, INPUTS, v_res, file = paste0(destination, "mat/", basin_name, ".RData"))
      fwrite(output, paste0(destination, "output/", basin_name, ".csv"))
      fwrite(Evalue, paste0(destination, "Evaluation/", basin_name, ".csv"))
      
      print(paste0(basin_name, " done"))
      print(paste("Elapsed time in running ga:", (proc.time() - ptm)["elapsed"]))
    }
  }
}