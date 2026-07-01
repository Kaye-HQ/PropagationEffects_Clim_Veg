library(data.table)
Rcpp::sourceCpp("./01_model_cal/utils_fast.cpp")

userf <- function(data_table, target_f) {
  dt <- as.data.table(data_table)
  sim_cols <- grep("Simulated_", names(dt), value = TRUE)
  
  if (is.list(target_f) && !is.data.frame(target_f)) {
    target_f <- do.call(rbind, lapply(target_f, function(x) {
       ux <- unlist(x)
      c(ux[1], ux[2], ifelse(length(ux) >= 3, ux[3], "1.0"))
    }))
    target_f <- as.data.frame(target_f, stringsAsFactors = FALSE)
  }
  
  if(nrow(target_f) > 0){
    
    fit_all <- matrix(NA, nrow = nrow(target_f), ncol = length(sim_cols))
    Evalue_all <- matrix(NA, nrow = nrow(target_f), ncol = length(sim_cols))
    
    for (n in 1:nrow(target_f)) {
      Sigsn <- as.character(unlist(target_f[n, 1]))
      methodn <- as.character(unlist(target_f[n, 2]))
      
      if (Sigsn == 'Daily') {
        res <- dt
        temp_num <- which(!is.na(res$Obs) & !is.infinite(res$Obs) & complete.cases(res[, ..sim_cols]))
        XX <- res$Obs
        YY <- as.matrix(res[, ..sim_cols])
      } else {
        
        counts <- dt[, .N, by = Year]
        valid_years <- counts[N >= 346]$Year
        dt_valid <- dt[Year %in% valid_years]
        
        if (Sigsn == 'Q_10') {
          res <- dt_valid[, lapply(.SD, function(x) quantile(x, 0.10, na.rm = TRUE)), by = Year, .SDcols = c("Obs", sim_cols)]
          XX <- res$Obs; YY <- as.matrix(res[, ..sim_cols])
        } else if (Sigsn == 'Q_log_10') {
          res <- dt_valid[, lapply(.SD, function(x) quantile(x, 0.10, na.rm = TRUE)), by = Year, .SDcols = c("Obs", sim_cols)]
          XX <- log(res$Obs); YY <- log(as.matrix(res[, ..sim_cols]))
        } else if (Sigsn == 'Q_90') {
          res <- dt_valid[, lapply(.SD, function(x) quantile(x, 0.90, na.rm = TRUE)), by = Year, .SDcols = c("Obs", sim_cols)]
          XX <- res$Obs; YY <- as.matrix(res[, ..sim_cols])
        } else if (Sigsn == 'Mean') {
          res <- dt_valid[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = Year, .SDcols = c("Obs", sim_cols)]
          XX <- res$Obs; YY <- as.matrix(res[, ..sim_cols])
        } else if (Sigsn == 'LFD') {
          TLFD <- 0.2 * mean(dt_valid$Obs, na.rm = TRUE)
          res <- dt_valid[, lapply(.SD, function(x) sum(x < TLFD, na.rm = TRUE)), by = Year, .SDcols = c("Obs", sim_cols)]
          XX <- res$Obs; YY <- as.matrix(res[, ..sim_cols])
        } else if (Sigsn == 'LFD_Q10') {
          TLFD <- quantile(dt_valid$Obs, 0.10, na.rm = TRUE)
          res <- dt_valid[, lapply(.SD, function(x) sum(x < TLFD, na.rm = TRUE)), by = Year, .SDcols = c("Obs", sim_cols)]
          XX <- res$Obs; YY <- as.matrix(res[, ..sim_cols])
        }
        
        
        temp_num <- which(!is.na(XX) & !is.infinite(XX) & complete.cases(YY) & !is.infinite(rowSums(YY)))
      }
      
      observed <- XX[temp_num]
      simulated <- YY[temp_num, , drop = FALSE]
      
      eval_res <- methods_evaluate_cpp(methodn, observed, simulated)
      Evalue_all[n, ] <- eval_res$value
      fit_all[n, ] <- eval_res$fit
    }
    
    
    combinedData <- cbind(target_f[, 1:2, drop = FALSE], Evalue_all)
    col_names <- c("Sig", "Method", paste0("Simulated_", 1:length(sim_cols)))
    Evalue_dt <- as.data.frame(combinedData)
    names(Evalue_dt) <- col_names
    
    return(list(fit_all = fit_all, Evalue_dt = Evalue_dt))
  }else{
    return(list(fit_all = NULL, Evalue_dt = NULL))
  }
}

parread <- function(pathdata0, basin_name, basin_info_n, LAIpath, usage, col_map ) {
  print(paste0('read ', basin_name))
  
  in_dt <- fread(paste0(pathdata0, basin_name, '.csv'))
  LAI <- fread(paste0(LAIpath, basin_name, '.csv'))
  
  in_dt[, Date := as.Date(Date)]
  LAI[, Date := as.Date(Date)]
  
  existing_lai_cols <- intersect(names(LAI), col_map$lai_essential)
  LAI_sub <- LAI[, ..existing_lai_cols]
  
  in_dt <- merge(in_dt, LAI_sub, all = TRUE)
  
  INPUTS <- list()
  if(is.na(basin_info_n$Area)) {
    INPUTS$area <- basin_info_n$Area_shp
  } else {
    INPUTS$area <- basin_info_n$Area
  }
  
  
  if (usage == 'SimulationRaw') {
    desiredLength <- nrow(in_dt[Year >= 1981])
    LAI[, Year := year(Date)]
    
    lai_1981 <- LAI[Year == 1981]$LAI_GLASS
    RepeatedLAI <- rep(lai_1981, length.out = desiredLength)
    
    LAI_GLASSraw <- in_dt$LAI_GLASS
    in_dt[Year >= 1981, LAI_GLASS := RepeatedLAI]
  }
  
  
  in_dt <- add_Rn(in_dt)
  in_dt <- add_VPD(in_dt)
  in_dt <- add_fvalSoil(in_dt)
  
  INPUTS$in_df <- in_dt   
  if (usage == 'SimulationRaw') {
    INPUTS$in_df$LAI_GLASSraw <- LAI_GLASSraw 
  }
  
  warmupyears <- 2
  
  
  in_dt[is.na(Emiss), Emiss := 0.95]
  temp <- which(!is.na(in_dt$Flow) & !is.na(in_dt$LAI_GLASS) & 
                  !is.na(in_dt$Emiss) & !is.na(in_dt$Albedo) & 
                  !is.na(year(in_dt$Date)) & !is.na(in_dt$P_MSWEP))
  
  INPUTS$Date <- in_dt$Date[temp]
  
  if (length(temp) > 0) {
    syear <- max(1979, year(in_dt$Date[temp[1]]) - warmupyears)
    s1 <- which(year(in_dt$Date) == syear & !is.na(in_dt$P_MSWEP))
    Nstart <- s1[1]
    if (Nstart < 50) Nstart <- Nstart + 50
    
    
    idx_range <- (Nstart - 50):nrow(in_dt)
    tempLAI <- in_dt$LAI_GLASS[idx_range]
    if (any(is.na(tempLAI))) {
      tempLAI[is.na(tempLAI)] <- mean(in_dt$LAI_GLASS, na.rm = TRUE)
      in_dt$LAI_GLASS[idx_range] <- tempLAI
    }
    
    
    cal_index <- temp[1:floor(length(temp) / 3 * 2)]
    val_index <- setdiff(temp, cal_index)
    
    INPUTS$yearn <- year(in_dt$Date)
    INPUTS$Nstart <- Nstart
    INPUTS$Nend <- cal_index[length(cal_index)]
    INPUTS$Nstart_c <- temp[1]
    INPUTS$Nend_c <- cal_index[length(cal_index)]
    INPUTS$Nstart_v <- val_index[1]
    INPUTS$Nend_v <- temp[length(temp)]
    
    
    INPUTS$KeyDates <- in_dt$Date[c(INPUTS$Nstart, INPUTS$Nstart_c, INPUTS$Nend, INPUTS$Nstart_v, INPUTS$Nend_v)]
    
    Nstart <- INPUTS$Nstart      
    Nend <- INPUTS$Nend
  } else {
    Nstart <- NA
    Nend <- NA
  }
  
  return(list(Nstart = Nstart, Nend = Nend, INPUTS = INPUTS))
}


add_Rn <- function(in_dt) {
  
  Stefan <- 4.903e-9 #  [MJ K-4 m-2 day-1]
  
  Rns <- (1 - in_dt$Albedo) * in_dt$SWd
  RLout <- in_dt$Emiss * Stefan * ((in_dt$Tmax + in_dt$Tmin)/2 + 273.15)^4
  RLout <- RLout / 0.0864 #  [W/m2]
  Rnl <- in_dt$LWd - RLout 
  
  
  in_dt$Rn <- pmax(Rns + Rnl, 0) 
  return(in_dt)
}

add_VPD <- function(in_dt) {
  
  vapor_pressure <- function(x) { 0.6108 * exp(17.27 * x / (x + 237.3)) } # kPa
  
  p <- in_dt$Pres / 1000 #  [kPa]
  q <- in_dt$SpecHum     # [g/g or g/kg]
  ea <- q * p / (0.622 + 0.378 * q) #  [kPa]
  
  es_tmax <- vapor_pressure(in_dt$Tmax)
  es_tmin <- vapor_pressure(in_dt$Tmin) # ℃
  es <- (es_tmax + es_tmin) / 2 
  
  in_dt$VPD <- pmax(es - ea, 0.001)
  return(in_dt)
}

add_fvalSoil <- function(in_dt) {
  Cp <- 4.2 * 0.242 
  lamada <- 2500 - 2.2 * in_dt$Temp     
  
  frame <- 33  
  
  gama <- Cp * (in_dt$Pres / 1000) / (0.622 * lamada) # [kPa/℃]
  slop <- 4098 * 0.6108 * exp((17.27 * in_dt$Temp) / (in_dt$Temp + 237.3)) / (in_dt$Temp + 237.3)^2 # [kPa/℃]
  epsilon <- slop / gama
  
  Eeq <- epsilon / (epsilon + 1) * in_dt$Rn / lamada * 86400 * 10^-3
  Eeq <- pmax(0.0001, Eeq, na.rm = TRUE)  
  
  
  kA <- 0.9 
  Tou <- exp(-kA * in_dt$LAI_GLASS)
  Es_eq <- Eeq * Tou   
  
  mov_P <- roll_mean_right_partial_cpp(in_dt$P_MSWEP, frame)
  mov_E <- roll_mean_right_partial_cpp(Es_eq, frame)
  
  fval_soil <- mov_P / mov_E
  fval_soil <- pmin(fval_soil, 1) 
  fval_soil <- pmax(fval_soil, 0) 
  
  in_dt$Eeq <- Eeq       
  in_dt$fval_soil <- fval_soil
  
  return(in_dt)
}


parrun_HBV_PML_Rcpp <- function(INPUTS, PAR, period, ifRoute, ifPML, target_f) {
  
  par_df <- data.frame(
    Alpha=PAR[1], Thelta=PAR[2], m=PAR[3], Am_25=PAR[4], D0=PAR[5],
    kQ=PAR[6], kA=PAR[7], S_sls=PAR[8], fER0=PAR[9], VPDmin=PAR[10],
    VPDmax=PAR[11], hc=PAR[12], Ts=PAR[13], CFMAX=PAR[14], CFR=PAR[15],
    CWH=PAR[16], BETA=PAR[17], LP=PAR[18], FC=PAR[19], PERC=PAR[20],
    K0=PAR[21], K1=PAR[22], K2=PAR[23], UZL=PAR[24], MAXBAS=PAR[25]
  )
  
  Nstart <- INPUTS$Nstart
  Nend <- if(period == 'c') INPUTS$Nend else INPUTS$Nend_v
  input_data <- INPUTS$in_df[Nstart:Nend, ]
  
  sim <- HBV_PML_Rcpp(pars = par_df, input = input_data, ifPML = ifPML, area = INPUTS$area, ifRoute = ifRoute)
  
  if (period == 'c') {
    start_idx <- INPUTS$Nstart_c - Nstart + 1
    end_idx   <- INPUTS$Nend_c - Nstart + 1
    
    cal_obs <- INPUTS$in_df$Flow[INPUTS$Nstart_c:INPUTS$Nend_c]
    cal_sim <- sim$QMOD[start_idx:end_idx]
    
    cal_dt <- data.frame(Year = INPUTS$in_df$Year[INPUTS$Nstart_c:INPUTS$Nend_c], 
                         Obs = cal_obs, Simulated_1 = cal_sim)
    
    res_cal <- userf(cal_dt, target_f)
    fit_all <- res_cal$fit_all
    fit_all[is.na(fit_all)] <- 0
    
    weights <- sapply(target_f, function(x) as.numeric(x[3]))
    weights[is.na(weights)] <- 0
    
    fit <- as.numeric(weights %*% fit_all)
    return(fit)
    
  } else {
    
    n_len <- nrow(sim)
    P_eff <- sim$P_after_s[2:n_len]; ET_sim <- sim$ET_sim[2:n_len]; QT_mm <- sim$QT_mm[2:n_len]
    dSM <- sim$SM[1:(n_len-1)] - sim$SM[2:n_len]
    dUZ <- sim$UZ[1:(n_len-1)] - sim$UZ[2:n_len]
    dLZ <- sim$LZ[1:(n_len-1)] - sim$LZ[2:n_len]
    
    R_check <- P_eff - ET_sim - QT_mm + dSM + dUZ + dLZ
    lag_check <- sim$QT_mm - (sim$QMOD / INPUTS$area * 86.4)
    
    wb_v <- ifelse(abs(sum(lag_check, na.rm = TRUE)) < 5, "WB_valied", "WB_not_valied")
    
    
    # --- Calibration Evaluation Block ---
    cal_start <- INPUTS$Nstart_c - Nstart + 1
    cal_end   <- INPUTS$Nend_c - Nstart + 1
    cal_dt <- data.frame(Year = input_data$Year[cal_start:cal_end], 
                         Obs = input_data$Flow[cal_start:cal_end], 
                         Simulated_1 = sim$QMOD[cal_start:cal_end])
    res_cal <- userf(cal_dt, target_f)
    Evalue_cal <- res_cal$Evalue_dt
    Evalue_cal$Type <- 'cal'
    
    # Add exact time periods for the calibration phase
    Evalue_cal$StartDate <- input_data$Date[cal_start]
    Evalue_cal$EndDate   <- input_data$Date[cal_end]
    
    # --- Validation Evaluation Block ---
    val_start <- INPUTS$Nstart_v - Nstart + 1
    val_end   <- INPUTS$Nend_v - Nstart + 1
    val_dt <- data.frame(Year = input_data$Year[val_start:val_end], 
                         Obs = input_data$Flow[val_start:val_end], 
                         Simulated_1 = sim$QMOD[val_start:val_end])
    res_val <- userf(val_dt, target_f)
    Evalue_val <- res_val$Evalue_dt
    Evalue_val$Type <- 'val'
    
    # Add exact time periods for the validation phase
    Evalue_val$StartDate <- input_data$Date[val_start]
    Evalue_val$EndDate   <- input_data$Date[val_end]
    
    # Combine the dataframes
    Evalue <- rbind(Evalue_cal, Evalue_val)
    
    output <- data.frame(Date = input_data$Date, Qobs = input_data$Flow, QMOD = sim$QMOD, 
                         QT_m3 = sim$QT_m3, QT_mm = sim$QT_mm, ET_sim = sim$ET_sim)
    
    return(list(output = output, Best_par = PAR, R_check = R_check, wb_v = wb_v, Evalue = Evalue))
  }
}

detrend_Data <- function(dailydata, dates, var_name) {
  
  yrs <- year(dates)
  df <- data.table(val = dailydata, Year = yrs)
  annualStats <- df[, .(annual_mean = mean(val, na.rm = TRUE)), by = Year]
  setorder(annualStats, Year) 
  
  uniqueYears <- annualStats$Year
  annualdata <- annualStats$annual_mean
  yearsNumeric <- uniqueYears - uniqueYears[1] 
  
  fit <- lm(annualdata ~ yearsNumeric)
  p_intercept <- coef(fit)[["(Intercept)"]]
  p_slope <- coef(fit)[["yearsNumeric"]]
  
  if(is.na(p_slope)) p_slope <- 0
  
  predictedAnnualdata1 <- p_intercept + (p_slope * yearsNumeric)
  predictedAnnualdata2 <- rep(p_intercept, length(yearsNumeric)) # [0, p(2)]
  
  # Detrending
  detrendedAnnual <- annualdata + (predictedAnnualdata2 - predictedAnnualdata1)
  
  if (var_name == 'LAI') {
    detrendedAnnual[detrendedAnnual < 0] <- 0
  }
  
  if (var_name == 'T') {
    scalingFactor <- detrendedAnnual - annualdata
    # Linear interpolation/extrapolation mapping to daily dates
    daily_scaling <- approx(x = uniqueYears, y = scalingFactor, xout = yrs, rule = 2)$y
    detrendeddata <- dailydata + daily_scaling
  } else {
    scalingFactor <- detrendedAnnual / annualdata
    scalingFactor[is.nan(scalingFactor) | is.infinite(scalingFactor)] <- 1
    daily_scaling <- approx(x = uniqueYears, y = scalingFactor, xout = yrs, rule = 2)$y
    detrendeddata <- dailydata * daily_scaling
  }
  
  return(detrendeddata)
}

parread_exps <- function(pathdata0, basin_name, basin_info_n, LAIpath, exp_type, col_map,
                         Year_base =1981 ,
                         Year_end = 2018  ) {
  print(paste0('read ', basin_name))
  
  in_dt <- fread(paste0(pathdata0, basin_name, '.csv'))
  LAI <- fread(paste0(LAIpath, basin_name, '.csv'))
  
  in_dt[, Date := as.Date(Date)]
  LAI[, Date := as.Date(Date)]
  
  existing_lai_cols <- intersect(names(LAI), col_map$lai_essential)
  LAI_sub <- LAI[, ..existing_lai_cols]
  
  in_dt <- merge(in_dt, LAI_sub, all = TRUE)
  
  INPUTS <- list()
  if(is.na(basin_info_n$Area)) {
    INPUTS$area <- basin_info_n$Area_shp
  } else {
    INPUTS$area <- basin_info_n$Area
  }
  
 
  in_dt[, `:=`(
    LAI_GLASSraw = LAI_GLASS,
    Praw = P_MSWEP,
    Traw = Temp,
    Year = year(Date)
  )]
  
  idx_period <- which(in_dt$Year >= Year_base & in_dt$Year <= Year_end)
  
  if (exp_type == 'State_LAI') {
    detrended <- detrend_Data(in_dt$LAI_GLASS[idx_period], in_dt$Date[idx_period], 'LAI')
    in_dt$LAI_GLASS[idx_period] <- detrended
    
  } else if (exp_type == 'State_P') {
    detrended <- detrend_Data(in_dt$P_MSWEP[idx_period], in_dt$Date[idx_period], 'P')
    in_dt$P_MSWEP[idx_period] <- detrended
    
  } else if (exp_type == 'State_T') {
    detrended <- detrend_Data(in_dt$Temp[idx_period], in_dt$Date[idx_period], 'T')
    in_dt$Temp[idx_period] <- detrended
    in_dt$Tmax[idx_period] <- in_dt$Tmax[idx_period] + in_dt$Temp[idx_period] - in_dt$Traw[idx_period]
    in_dt$Tmin[idx_period] <- in_dt$Tmin[idx_period] + in_dt$Temp[idx_period] - in_dt$Traw[idx_period]
    
  } else if (exp_type == 'State_P_T_LAI') {
    # 1. LAI
    detrended_LAI <- detrend_Data(in_dt$LAI_GLASS[idx_period], in_dt$Date[idx_period], 'LAI')
    in_dt$LAI_GLASS[idx_period] <- detrended_LAI
    
    # 2. P
    detrended_P <- detrend_Data(in_dt$P_MSWEP[idx_period], in_dt$Date[idx_period], 'P')
    in_dt$P_MSWEP[idx_period] <- detrended_P
    
    # 3. T
    detrended_T <- detrend_Data(in_dt$Temp[idx_period], in_dt$Date[idx_period], 'T')
    in_dt$Temp[idx_period] <- detrended_T
    
    in_dt$Tmax[idx_period] <- in_dt$Tmax[idx_period] + in_dt$Temp[idx_period] - in_dt$Traw[idx_period]
    in_dt$Tmin[idx_period] <- in_dt$Tmin[idx_period] + in_dt$Temp[idx_period] - in_dt$Traw[idx_period]
  }
  
  in_dt <- add_Rn(in_dt)
  in_dt <- add_VPD(in_dt)
  in_dt <- add_fvalSoil(in_dt)
  
  INPUTS$in_df <- in_dt   
  
  warmupyears <- 2
  
  temp <- which(!is.na(in_dt$Flow) & !is.na(in_dt$LAI_GLASS) & 
                  !is.na(in_dt$Albedo) & !is.na(in_dt$Year) & !is.na(in_dt$P_MSWEP))
  
  idx_emiss <- intersect(which(is.na(in_dt$Emiss)), temp)
  if(length(idx_emiss) > 0) in_dt$Emiss[idx_emiss] <- 0.95
  
  INPUTS$Date <- in_dt$Date[temp]
  
  if (length(temp) > 0) {
    syear <- max(1979, in_dt$Year[temp[1]] - warmupyears)
    s1 <- which(in_dt$Year == syear & !is.na(in_dt$P_MSWEP))
    Nstart <- s1[1]
    if (Nstart < 50) Nstart <- Nstart + 50
    
    idx_range <- (Nstart - 50):nrow(in_dt)
    tempLAI <- in_dt$LAI_GLASS[idx_range]
    if (any(is.na(tempLAI))) {
      tempLAI[is.na(tempLAI)] <- mean(in_dt$LAI_GLASS, na.rm = TRUE)
      in_dt$LAI_GLASS[idx_range] <- tempLAI
    }
    
    cal_index <- temp[1:floor(length(temp) / 3 * 2)]
    val_index <- setdiff(temp, cal_index)
    
    INPUTS$yearn <- in_dt$Year
    INPUTS$Nstart <- Nstart
    INPUTS$Nend <- cal_index[length(cal_index)]
    INPUTS$Nstart_c <- temp[1]
    INPUTS$Nend_c <- cal_index[length(cal_index)]
    INPUTS$Nstart_v <- val_index[1]
    INPUTS$Nend_v <- temp[length(temp)]
    
    INPUTS$KeyDates <- in_dt$Date[c(INPUTS$Nstart, INPUTS$Nstart_c, INPUTS$Nend, INPUTS$Nstart_v, INPUTS$Nend_v)]
    
    Nstart <- INPUTS$Nstart      
    Nend <- INPUTS$Nend
  } else {
    Nstart <- NA
    Nend <- NA
  }
  
  return(list(Nstart = Nstart, Nend = Nend, INPUTS = INPUTS))
}

parread_exps_mean_annual <- function(pathdata0, basin_name, basin_info_n, LAIpath, exp_type, col_map,
                                     Year_base = 1981 ) {
  print(paste0('read ', basin_name))
  
  in_dt <- fread(paste0(pathdata0, basin_name, '.csv'))
  LAI <- fread(paste0(LAIpath, basin_name, '.csv'))
  
  in_dt[, Date := as.Date(Date)]
  LAI[, Date := as.Date(Date)]
  
  existing_lai_cols <- intersect(names(LAI), col_map$lai_essential)
  LAI_sub <- LAI[, ..existing_lai_cols]
  
  in_dt <- merge(in_dt, LAI_sub,  all = TRUE)
  
  INPUTS <- list()
  if(is.na(basin_info_n$Area)) {
    INPUTS$area <- basin_info_n$Area_shp
  } else {
    INPUTS$area <- basin_info_n$Area
  }
  

  in_dt$LAI_GLASSraw <- in_dt$LAI_GLASS
  in_dt$Praw <- in_dt$P_MSWEP
  in_dt$Traw <- in_dt$Temp
  
  if (exp_type == 'State_LAI') {
    desiredLength <- nrow(in_dt[Year >= Year_base])
    LAI[, Year := year(Date)]
    state_data <- LAI[Year %in% c(Year_base, Year_base+1, Year_base+2)]$LAI_GLASS
    RepeatedLAI <- rep(state_data, length.out = desiredLength)
    
    in_dt[Year >= Year_base, LAI_GLASS := RepeatedLAI]
    
  } else if (exp_type == 'State_P') {
    desiredLength <- nrow(in_dt[Year >= Year_base])
    state_data <- in_dt[Year %in% c(Year_base, Year_base+1, Year_base+2)]$P_MSWEP
    Repeated_data <- rep(state_data, length.out = desiredLength)
    
    in_dt[Year >= Year_base, P_MSWEP := Repeated_data]
    
  } else if (exp_type == 'State_T') {
    desiredLength <- nrow(in_dt[Year >= Year_base])
    state_data <- in_dt[Year %in% c(Year_base, Year_base+1, Year_base+2)]$Temp
    Repeated_data <- rep(state_data, length.out = desiredLength)
    
    
    in_dt[Year >= Year_base, Temp := Repeated_data]
    in_dt$Tmax =  in_dt$Tmax +  in_dt$Temp - in_dt$Traw
    in_dt$Tmin =  in_dt$Tmin +  in_dt$Temp - in_dt$Traw
    
    
  } else if (exp_type == 'State_P_T_LAI') {
    desiredLength <- nrow(in_dt[Year >= Year_base])
    LAI[, Year := year(Date)]
    
    state_data_T <- in_dt[Year %in% c(Year_base, Year_base+1, Year_base+2)]$Temp
    
    in_dt[Year >= Year_base, Temp := rep(state_data_T, length.out = desiredLength)]
    in_dt$Tmax =  in_dt$Tmax +  in_dt$Temp - in_dt$Traw
    in_dt$Tmin =  in_dt$Tmin +  in_dt$Temp - in_dt$Traw
    
    
    state_data_P <- in_dt[Year %in% c(Year_base, Year_base+1, Year_base+2)]$P_MSWEP
    
    in_dt[Year >= Year_base, P_MSWEP := rep(state_data_P, length.out = desiredLength)]
    
    
    state_data_LAI <- LAI[Year %in% c(Year_base, Year_base+1, Year_base+2)]$LAI_GLASS
    in_dt[Year >= Year_base, LAI_GLASS := rep(state_data_LAI, length.out = desiredLength)]
    
  }
  
  in_dt <- add_Rn(in_dt)
  in_dt <- add_VPD(in_dt)
  in_dt <- add_fvalSoil(in_dt)
  
  INPUTS$in_df <- in_dt   
  
  
  warmupyears <- 2
  
  temp <- which(!is.na(in_dt$Flow) & !is.na(in_dt$LAI_GLASS) & 
                  !is.na(in_dt$Albedo) & !is.na(year(in_dt$Date)) & !is.na(in_dt$P_MSWEP))
  
  idx_emiss <- intersect(which(is.na(in_dt$Emiss)), temp)
  if(length(idx_emiss) > 0) in_dt$Emiss[idx_emiss] <- 0.95
  
  INPUTS$Date <- in_dt$Date[temp]
  
  if (length(temp) > 0) {
    syear <- max(1979, year(in_dt$Date[temp[1]]) - warmupyears)
    s1 <- which(year(in_dt$Date) == syear & !is.na(in_dt$P_MSWEP))
    Nstart <- s1[1]
    if (Nstart < 50) Nstart <- Nstart + 50
    
    idx_range <- (Nstart - 50):nrow(in_dt)
    tempLAI <- in_dt$LAI_GLASS[idx_range]
    if (any(is.na(tempLAI))) {
      tempLAI[is.na(tempLAI)] <- mean(in_dt$LAI_GLASS, na.rm = TRUE)
      in_dt$LAI_GLASS[idx_range] <- tempLAI
    }
    
    cal_index <- temp[1:floor(length(temp) / 3 * 2)]
    val_index <- setdiff(temp, cal_index)
    
    INPUTS$yearn <- year(in_dt$Date)
    INPUTS$Nstart <- Nstart
    INPUTS$Nend <- cal_index[length(cal_index)]
    INPUTS$Nstart_c <- temp[1]
    INPUTS$Nend_c <- cal_index[length(cal_index)]
    INPUTS$Nstart_v <- val_index[1]
    INPUTS$Nend_v <- temp[length(temp)]
    
    INPUTS$KeyDates <- in_dt$Date[c(INPUTS$Nstart, INPUTS$Nstart_c, INPUTS$Nend, INPUTS$Nstart_v, INPUTS$Nend_v)]
    
    Nstart <- INPUTS$Nstart      
    Nend <- INPUTS$Nend
  } else {
    Nstart <- NA
    Nend <- NA
  }
  
  return(list(Nstart = Nstart, Nend = Nend, INPUTS = INPUTS))
}
