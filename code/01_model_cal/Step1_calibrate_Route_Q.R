rm(list=ls())
library(data.table)
library(GA)

source("01_model_cal/functions.R")

destination <- "../Result/HBV_PML/01_cal_10/"
dir.create(destination, showWarnings = FALSE, recursive = TRUE)
sapply(c("mat/", "output/", "Evaluation/"), 
       function(x) dir.create(paste0(destination, x), showWarnings = FALSE))


basin_all  <- fread("../Data/info_2252.csv")
basin_info <- basin_all

pathdata0 <- "../Data/sample_data/"
LAIpath <- "../Data/sample_data/LAI/"

lb <- c(0, 0, 0, 0, 0.5010, 0.1000, 0.6800, 0, 0, 0.6570, 3.5000, 0.0100, -5, 0, 0, 0, 0, 0.3, 1, 0, 0.05, 0.01, 0.02, 0, 1)
ub <- c(0.0480, 0.0690, 25.7550, 46.3750, 2.0000, 1.0000, 0.8990, 0.1740, 0.1580, 1.4920, 6.5000, 10.0000, 5, 20, 1, 0.8, 7, 1, 2000, 100, 2, 1, 0.1, 300, 6)

target_f <- list(c("Mean", "KGE"))
numberOfVariables <- 25; ifRoute <- 1; ifPML <- 1

for (ibasin in 1:nrow(basin_info)) {
  basin_name <- basin_info$Name[ibasin]
  basin_info_n <- basin_info[ibasin, ]
  
  res_read <- parread(pathdata0, basin_name, basin_info_n, LAIpath, "SimulationRaw")
  INPUTS <- res_read$INPUTS
  
  if (length(INPUTS$Date) > 0) {
    print(paste0(basin_name, "_start"))
    

    fitness_func <- function(par) {
        score <- parrun_HBV_PML(INPUTS, matrix(par, nrow=1), 'c', ifRoute, ifPML, target_f)
        if(is.na(score)) return(-9999)
        return(-score) 
    }
    
    ptm <- proc.time()
    ga_res <- ga(type = "real-valued", 
                 fitness = fitness_func, 
                 lower = lb, upper = ub, 
                 popSize = 500, maxiter = 100, 
                 run = 10,
                 monitor = TRUE)
    
    pop <- ga_res@population
    fits <- ga_res@fitness
    top1_idx <- order(fits, decreasing = TRUE)[1]
    Best1_par <- pop[top1_idx, , drop = FALSE]
    Best_par <- Best1_par
    
    v_res <- parrun_HBV_PML(INPUTS, Best1_par, 'v', ifRoute, ifPML, target_f)
    
    output <- v_res$output
    output_sim_median <- v_res$output_sim_median
    output_sim_sd <- v_res$output_sim_sd
    wb_v <- "WB_valied"
    
    save(Best1_par, Best_par,INPUTS, output, output_sim_median, output_sim_sd, wb_v, 
         file = paste0(destination, "mat/", basin_name, ".RData"))
    fwrite(output, paste0(destination, "output/", basin_name, ".csv"))
  
    print(paste0(basin_name, " done"))

    print(paste("Elapsed time in running ga:", (proc.time() - ptm)["elapsed"]))
    
  }
}