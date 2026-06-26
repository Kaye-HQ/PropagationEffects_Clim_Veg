library(zoo)
library(data.table)
library(lubridate)


TraToPar_HBV_PML <- function(X) {
  if(is.vector(X)) X <- matrix(X, nrow=1)
  par <- list(
    Alpha   = X[, 1],  Thelta  = X[, 2],  m       = X[, 3],
    Am_25   = X[, 4],  D0      = X[, 5],  kQ      = X[, 6],
    kA      = X[, 7],  S_sls   = X[, 8],  fER0    = X[, 9],
    VPDmin  = X[, 10], VPDmax  = X[, 11], hc      = X[, 12],
    Ts      = X[, 13], CFMAX   = X[, 14], CFR     = X[, 15],
    CWH     = X[, 16], BETA    = X[, 17], LP      = X[, 18],
    FC      = X[, 19], PERC    = X[, 20], K0      = X[, 21],
    K1      = X[, 22], K2      = X[, 23], UZL     = X[, 24],
    MAXBAS  = X[, 25]
  )
  return(par)
}

prelocate <- function(nsize, PAR, INPUTS) {
  state <- list(
    v   = rep(0, nsize), vl  = rep(0, nsize),
    SM  = rep(0, nsize), UZ  = rep(0, nsize),
    LZ  = rep(0, nsize), fip = rep(0, nsize),
    fop = rep(0, nsize)
  )
  par <- TraToPar_HBV_PML(PAR)
  
  target_vars <- c("Temp", "P_MSWEP", "LAI_GLASS", "Pres", "Wind", 
                   "SWd", "CO2", "Rn", "VPD", "Eeq", "fval_soil", "PET_FAO")
  
  
  valid_vars <- intersect(target_vars, names(INPUTS$in_df))
  

  Table0 <- INPUTS$in_df[, valid_vars,with = F]
  In_N <- array(0, dim = c(nrow(Table0), ncol(Table0), nsize))
  for (isize in 1:nsize) {
    In_N[,,isize] <- as.matrix(Table0)
  }
  return(list(state = state, par = par, In_N = In_N, valid_vars = valid_vars))
}

snow_routine <- function(in_data, par, state) {
  Ts    <- par$Ts; CFMAX <- par$CFMAX; CFR <- par$CFR; CWH <- par$CWH
  temp  <- in_data$Temp; prec <- in_data$P_MSWEP
  v     <- state$v; vl <- state$vl
  N     <- length(prec)
  P <- rep(0, N)
  
  rain <- prec; rain[temp < Ts] <- 0
  snow <- prec; snow[temp >= Ts] <- 0
  Ta <- temp - Ts; Ta[temp < Ts] <- 0
  Tn <- Ts - temp; Tn[temp >= Ts] <- 0
  
  m   <- pmin(CFMAX * Ta, v)
  rfz <- pmin(CFR * CFMAX * Tn, vl)
  
  v   <- v - m + snow + rfz
  vl  <- vl + m + rain - rfz
  
  index  <- which(vl > (CWH * v))
  index2 <- which(vl <= (CWH * v))
  
  if(length(index) > 0) {
    P[index]  <- vl[index] - (CWH[index] * v[index])
    vl[index] <- CWH[index] * v[index]
  }
  if(length(index2) > 0) { P[index2] <- 0 }
  
  state$v <- v; state$vl <- vl
  out <- list(rfz = rfz, m = m, snow = snow, rain = rain)
  return(list(P = P, state = state, out = out))
}

Run_hbv_PML <- function(P, in_data, par, state, Case, area, ifPML) {
  if (ifPML == 1) {
    Prcp <- P
    LAI <- in_data$LAI_GLASS; tmean <- in_data$Temp; Pa <- in_data$Pres / 1000
    u2 <- in_data$Wind; Rs <- in_data$SWd; PAR_rad <- 0.45 * Rs
    PAR_mol <- PAR_rad * 4.57; Ca <- in_data$CO2; Rn <- in_data$Rn
    VPD <- in_data$VPD; Eeq <- in_data$Eeq
    
    hc <- par$hc; Alpha <- par$Alpha; Thelta <- par$Thelta; m <- par$m
    Am_25 <- par$Am_25; D0 <- par$D0; kQ <- par$kQ; kA <- par$kA
    S_sls <- par$S_sls; fER0 <- par$fER0; Dmin <- par$VPDmin; Dmax <- par$VPDmax
    LAIref <- 5; kmar <- 0.40; Zob <- 15
    rou_a <- 3.846 * 10^3 * Pa / (tmean + 273.15)
    
    f_VPD_gc <- 1 / (1 + VPD / D0)
    f_VPD <- ifelse(VPD > Dmax, 0, ifelse(VPD >= Dmin & VPD <= Dmax, (Dmax - VPD) / (Dmax - Dmin), 1))
    fT2 <- pmin(exp(0.031 * (tmean - 25)) / (1 + exp(0.115 * (tmean - 41))), 1)
    
    Am <- Am_25 * fT2
    P1 <- Am * Alpha * Thelta * PAR_mol; P2 <- Am * Alpha * PAR_mol
    P3 <- Am * Thelta * Ca; P4 <- Alpha * PAR_mol * Thelta * Ca
    
    Ags <- Ca * P1 / (P2 * kQ + P4 * kQ) * (kQ * LAI + log((P2 + P3 + P4) / (P2 + P3 * exp(kQ * LAI) + P4)))
    Ag <- Ags * f_VPD
    Gc <- m / Ca * Ag * 1.6 * f_VPD_gc
    Gc <- pmax(Gc * 1e-2 / (0.446 * (273 / (273 + tmean)) * (Pa / 101.3)), 10^-6)
    
    Tou <- exp(-kA * LAI)
    d <- 0.667 * hc; zom <- 0.125 * hc; zoh <- 0.135 * zom
    Cp <- 4.2 * 0.242; lamada <- 2500 - 2.4 * tmean
    gama <- Cp * Pa / (0.622 * lamada)
    slop <- 4098 * 0.6108 * exp((17.27 * tmean) / (tmean + 237.3)) / (tmean + 237.3)^2
    epsilon <- slop / gama
    
    uz <- log(67.8 * Zob - 5.42) / 4.87 * u2
    Ga <- uz * kmar^2 / (log((Zob - d) / zom) * log((Zob - d) / zoh))
    
    LEcr <- epsilon * Rn * (1 - Tou) / (epsilon + 1 + Ga / Gc)
    LEca <- (rou_a * Cp * Ga * VPD / gama) / (epsilon + 1 + Ga / Gc)
    Ecr <- LEcr / lamada * 86400 * 10^-3
    Eca <- LEca / lamada * 86400 * 10^-3
    
    fveg <- 1 - exp(-LAI / LAIref); Sveg <- S_sls * LAI; fER <- fER0 * fveg
    Pwet <- -log(1 - fER0) / fER0 * Sveg / fveg
    Ei <- ifelse(Prcp < Pwet, fveg * Prcp, fveg * Pwet + fER * (Prcp - Pwet))
    Ei[is.na(Ei)] <- 0
    
    Es_eq <- Eeq * Tou; 
    frame <- 32
    fval_soil <- zoo::rollapply(Prcp - Ei, width=frame, FUN=mean, na.rm=TRUE, align="right", fill=NA) / 
                 zoo::rollapply(Es_eq, width=frame, FUN=mean, na.rm=TRUE, align="right", fill=NA)
    fval_soil[is.na(fval_soil)] <- 0
    fval_soil <- pmax(pmin(fval_soil, 1), 0)
    Es <- fval_soil * Es_eq; Est <- Es + Eca + Ecr
    
    E_st <- rep(NA, length(P)); Pei <- rep(NA, length(P))
    Pei[(P - Ei) >= 0] <- P[(P - Ei) >= 0] - Ei[(P - Ei) >= 0]
    Pei[(P - Ei) < 0] <- 0
    Ei[(P - Ei) < 0] <- P[(P - Ei) < 0]
    
    UZ <- state$UZ; LZ <- state$LZ
    cond1 <- Est <= Pei + UZ + LZ; cond2 <- Est > Pei + UZ + LZ
    E_st[cond1] <- Est[cond1]
    E_st[cond2] <- Pei[cond2] + LZ[cond2] + UZ[cond2]
    ept <- E_st + Ei
  } else {
    Prcp <- P; LAI <- in_data$LAI_GLASS; ept <- in_data$PET_FAO
  }
  
  BETA <- par$BETA; LP <- par$LP; FC <- pmax(1e-10, par$FC)
  PERC <- par$PERC; K0 <- par$K0; K1 <- par$K1; K2 <- par$K2; UZL <- par$UZL
  SM <- state$SM; UZ <- state$UZ; LZ <- state$LZ
  
  R <- P * (SM / FC)^BETA
  SM_dummy <- pmax(pmin(SM + P - R, FC), 0)
  R <- R + pmax(SM + P - R - FC, 0) + pmin(SM + P - R, 0)
  
  EA <- ept * pmin(SM_dummy / (FC * LP), 1)
  SM <- pmax(pmin(SM_dummy - EA, FC), 0)
  EA <- EA + pmax(SM_dummy - EA - FC, 0) + pmin(SM_dummy - EA, 0)
  
  if (Case == 1) {
    Q0 <- pmax(pmin(K1 * UZ + K0 * pmax(UZ - UZL, 0), UZ), 0)
    RL <- pmax(pmin(UZ - Q0, PERC), 0)
  } else {
    RL <- pmax(pmin(PERC, UZ), 0)
    Q0 <- pmax(pmin(K1 * UZ + K0 * pmax(UZ - UZL, 0), UZ - RL), 0)
  }
  
  UZ <- UZ + R - Q0 - RL
  Q1 <- pmax(pmin(K2 * LZ, LZ), 0)
  LZ <- LZ + RL - Q1
  Q <- Q0 + Q1
  
  state$SM <- SM; state$UZ <- UZ; state$LZ <- LZ
  out <- list(EA=EA, ept=ept, R=R, RL=RL, Q0=Q0, Q1=Q1, Qsim_mm=Q, Qsim_m3s=Q*area/86.4)
  
  if (ifPML == 1) {
    out$Ei <- Ei; out$Est <- Est; out$Es <- Es; out$Et <- Eca + Ecr
    out$Ei_act <- Ei / (Es + Eca + Ecr + Ei) * EA
    out$Es_act <- Es / (Es + Eca + Ecr + Ei) * EA
    out$Et_act <- (Eca + Ecr) / (Es + Eca + Ecr + Ei) * EA
  } else { out$Ei <- NA; out$Est <- NA }
  
  return(list(state = state, out = out))
}

mytrimf <- function(x, param) {
  f <- rep(0, length(x))
  idx <- (x > param[1]) & (x <= param[2])
  f[idx] <- (x[idx] - param[1]) / (param[2] - param[1])
  idx <- (x > param[2]) & (x <= param[3])
  f[idx] <- (param[3] - x[idx]) / (param[3] - param[2])
  return(f)
}

get_Route <- function(cn, MAXBAS, c_list, N, QT_m3) {
  Qsim_out <- rep(NA, N)
  maxb <- round(MAXBAS[cn])
  if(maxb < 1) maxb <- 1
  for (t in maxb:N) {
    seq_idx <- (t - maxb + 1):t
    Qsim_out[t] <- sum(c_list[[cn]] * QT_m3[cn, seq_idx])
  }
  return(Qsim_out)
}

methods_evaluate <- function(method, observed, simulated) {
  simulated <- as.matrix(simulated)
  
  if (method == 'KGE') {
    mu_obs <- mean(observed, na.rm=TRUE)
    sigma_obs <- sd(observed, na.rm=TRUE)
    
    sim_means <- colMeans(simulated, na.rm=TRUE)
    sim_sds <- apply(simulated, 2, sd, na.rm=TRUE)
    
    obs_c <- observed - mu_obs
    sim_c <- sweep(simulated, 2, sim_means, "-")
    
    r <- colSums(obs_c * sim_c, na.rm=TRUE) / 
         (sqrt(sum(obs_c^2, na.rm=TRUE)) * sqrt(colSums(sim_c^2, na.rm=TRUE)))
    
    alpha <- sim_sds / sigma_obs
    beta <- sim_means / mu_obs
    
    value <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
    fit <- 1 - value
  } else if (method == 'NSE') {
    SS_total <- sum((observed - mean(observed, na.rm=TRUE))^2, na.rm=TRUE)
    SS_residual <- colSums(sweep(simulated, 1, observed, "-")^2, na.rm=TRUE)
    value <- 1 - (SS_residual / SS_total)
    fit <- 1 - value
  } else if (method == 'Rsq') {
    SS_total <- sum((observed - mean(observed, na.rm=TRUE))^2, na.rm=TRUE)
    SS_residual <- colSums(sweep(simulated, 1, observed, "-")^2, na.rm=TRUE)
    value <- 1 - (SS_residual / SS_total)
    fit <- 1 - value
  } else if (method == 'Abias') {
    diffs <- sweep(simulated, 1, observed, "-")
    value <- abs(colMeans(diffs, na.rm=TRUE) / mean(observed, na.rm=TRUE))
    fit <- value
  }
  
  value[is.infinite(value)] <- NA
  fit[is.infinite(fit)] <- NA
  
  return(list(value=value, fit=fit))
}

userf <- function(data_table, target_f) {
  fit_all <- matrix(NA, nrow=length(target_f), ncol=ncol(data_table)-2)
  Evalue_all <- matrix(NA, nrow=length(target_f), ncol=ncol(data_table)-2)
  
  dt <- as.data.table(data_table)
  sim_cols <- grep("Simulated_", names(dt), value=TRUE)
  
  for (n in seq_along(target_f)) {
    Sigsn <- target_f[[n]][1]
    methodn <- target_f[[n]][2]
    
    if (Sigsn == 'Daily') {
      res <- dt
      temp_num <- which(!is.na(res$Obs) & !is.infinite(res$Obs))
      XX <- res$Obs
      YY <- as.matrix(res[, ..sim_cols])
    } else {
      if (Sigsn == 'Q_10' || Sigsn == 'Q_log_10') {
        res <- dt[, lapply(.SD, function(x) quantile(x, 0.10, na.rm=TRUE)), by=Year, .SDcols=c("Obs", sim_cols)]
      } else if (Sigsn == 'Q_90') {
        res <- dt[, lapply(.SD, function(x) quantile(x, 0.90, na.rm=TRUE)), by=Year, .SDcols=c("Obs", sim_cols)]
      } else if (Sigsn == 'Mean') {
        res <- dt[, lapply(.SD, function(x) mean(x, na.rm=TRUE)), by=Year, .SDcols=c("Obs", sim_cols)]
      } else if (Sigsn == 'LFD') {
        TLFD <- 0.2 * mean(dt$Obs, na.rm=TRUE)
        res <- dt[, lapply(.SD, function(x) sum(x < TLFD, na.rm=TRUE)), by=Year, .SDcols=c("Obs", sim_cols)]
      } else if (Sigsn == 'LFD_Q10') {
        TLFD <- quantile(dt$Obs, 0.1, na.rm=TRUE)
        res <- dt[, lapply(.SD, function(x) sum(x < TLFD, na.rm=TRUE)), by=Year, .SDcols=c("Obs", sim_cols)]
      }
      
      
      counts <- dt[, .N, by=Year]
      valid_years <- counts[N >= 346]$Year
      res <- res[Year %in% valid_years]
      
      if (Sigsn == 'Q_log_10') {
        XX <- log(res$Obs)
        YY <- log(as.matrix(res[, ..sim_cols]))
      } else {
        XX <- res$Obs
        YY <- as.matrix(res[, ..sim_cols])
      }
      temp_num <- which(!is.na(XX) & !is.infinite(XX))
    }
    
    observed <- XX[temp_num]
    simulated <- YY[temp_num, , drop=FALSE]
    
    eval_res <- methods_evaluate(methodn, observed, simulated)
    Evalue_all[n, ] <- eval_res$value
    fit_all[n, ] <- eval_res$fit
  }
  
  return(list(fit_all = fit_all, Evalue_all = Evalue_all))
}


detrend_Data <- function(dailydata, dates, var) {
  dt <- data.table(Date = as.Date(dates), Data = dailydata)
  dt[, Year := year(Date)]
  
  annualStats <- dt[, .(mean_data = mean(Data, na.rm=TRUE)), by = Year]
  annualStats[, yearsNumeric := Year - min(Year)]
  
  
  model <- lm(mean_data ~ yearsNumeric, data = annualStats)
  p <- coef(model) # p[1] intercept, p[2] slope
  
  predicted1 <- p[1] + p[2] * annualStats$yearsNumeric
  predicted2 <- 0 + p[2] * annualStats$yearsNumeric
  
  annualStats[, detrendedAnnual := mean_data + (predicted2 - predicted1)]
  if(var == 'LAI') annualStats[detrendedAnnual < 0, detrendedAnnual := 0]
  
  dt <- merge(dt, annualStats[, .(Year, mean_data, detrendedAnnual)], by = "Year")
  
  if(var == 'T') {
    dt[, ScalingFactor := detrendedAnnual - mean_data]
    dt[, detrendeddata := Data + ScalingFactor]
  } else {
    dt[, ScalingFactor := detrendedAnnual / mean_data]
    dt[, detrendeddata := Data * ScalingFactor]
  }
  return(dt$detrendeddata)
}



parrun_HBV_PML <- function(INPUTS, PAR, period, ifRoute, ifPML, target_f) {
  
  area <- INPUTS$area
  PAR <- matrix(PAR, ncol=25) 
  nsize <- nrow(PAR)
  
    if (length(target_f) > 0) {
    weights <- sapply(target_f, function(x) {
      if (length(x) >= 3) {
        return(as.numeric(x[[3]]))
      } else {
        return(1.0)
      }
    })
  } else {
    weights <- numeric(0)
  }
  
  res_pre <- prelocate(nsize, PAR, INPUTS)
  state <- res_pre$state
  par <- res_pre$par
  In_N <- res_pre$In_N
  valid_vars <- res_pre$valid_vars 
  
  Nstart <- INPUTS$Nstart
  Nend <- if(period == 'c') INPUTS$Nend else INPUTS$Nend_v
  
  QT_m3 <- matrix(NA, nrow=nsize, ncol=Nend-Nstart+1)
  QT_mm <- matrix(NA, nrow=nsize, ncol=Nend-Nstart+1)
  ET_sim <- matrix(NA, nrow=nsize, ncol=Nend-Nstart+1)
  
 
  for (dloop in Nstart:Nend) {
    jj <- dloop - Nstart + 1
    inloop <- as.list(In_N[dloop, , 1])
    
      names(inloop) <- valid_vars 
    
    sn_res <- snow_routine(inloop, par, state)
    P <- sn_res$P; state <- sn_res$state
    
    md_res <- Run_hbv_PML(P, inloop, par, state, Case=1, area=area, ifPML=ifPML)
    state <- md_res$state
    
    QT_m3[, jj] <- md_res$out$Qsim_m3s
    QT_mm[, jj] <- md_res$out$Qsim_mm
    ET_sim[, jj] <- md_res$out$EA
  }
  

  MAXBAS <- round(par$MAXBAS)
  c_list <- lapply(1:nsize, function(n) mytrimf(1:MAXBAS[n], c(0, (MAXBAS[n]+1)/2, MAXBAS[n]+1)))
  
  QMOD <- matrix(NA, nrow=nsize, ncol=ncol(QT_m3))
  for(i in 1:nsize) {
    if(ifRoute == 1) {
      QMOD[i,] <- get_Route(i, MAXBAS, c_list, ncol(QT_m3), QT_m3)
    } else {
      QMOD[i,] <- QT_m3[i,]
    }
  }

  if (period == 'c') {
    N_start_c <- INPUTS$Nstart_c
    N_end_c <- INPUTS$Nend_c
    
    cal_obs <- INPUTS$in_df$Flow[N_start_c:N_end_c]
    cal_sim <- t(QMOD[, (N_start_c - Nstart + 1):(N_end_c - Nstart + 1), drop=FALSE])
    
    cal_dt <- data.frame(Year = INPUTS$in_df$Year[N_start_c:N_end_c], Obs = cal_obs)
    sim_df <- as.data.frame(cal_sim)
    colnames(sim_df) <- paste0("Simulated_", 1:ncol(sim_df))
    cal_dt <- cbind(cal_dt, sim_df)
    
    res_cal <- userf(cal_dt, target_f)
    fit_all <- res_cal$fit_all
    

    fit_all[is.na(fit_all)] <- 0
    weights[is.na(weights)] <- 0
    
    if (length(target_f) == 1) {
      fit <- fit_all[1, ]
    } else {
      fit <- as.numeric(weights %*% fit_all)
    }
    return(fit) 
    
  } else {
    N_start_v <- INPUTS$Nstart_v
    N_end_v <- INPUTS$Nend_v
    
    val_obs <- INPUTS$in_df$Flow[N_start_v:N_end_v]
    val_sim <- t(QMOD[, (N_start_v - Nstart + 1):(N_end_v - Nstart + 1), drop=FALSE])
    
    output <- data.frame(
      Date = INPUTS$in_df$Date[Nstart:Nend],
      Qobs = INPUTS$in_df$Flow[Nstart:Nend],
      QMOD = QMOD[1, ],
      QT_m3 = QT_m3[1, ],
      QT_mm = QT_mm[1, ],
      ET_sim = ET_sim[1, ]
    )
    

    output_sim_median <- data.frame(
      Date = INPUTS$in_df$Date[Nstart:Nend],
      QMOD = apply(QMOD, 2, median, na.rm=TRUE)
    )
    output_sim_sd <- data.frame(
      Date = INPUTS$in_df$Date[Nstart:Nend],
      QMOD = apply(QMOD, 2, sd, na.rm=TRUE)
    )
    
    return(list(
      output = output, 
      output_sim_median = output_sim_median, 
      output_sim_sd = output_sim_sd,
      Best10_par = PAR
    ))
  }
}


library(data.table)
library(lubridate)
library(zoo)


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
  
  mov_P <- rollapply(in_dt$P_MSWEP, width = frame, FUN = function(x) mean(x, na.rm=TRUE), align = "right", fill = NA, partial = TRUE)
  mov_E <- rollapply(Es_eq, width = frame, FUN = function(x) mean(x, na.rm=TRUE), align = "right", fill = NA, partial = TRUE)
  
  fval_soil <- mov_P / mov_E
  fval_soil <- pmin(fval_soil, 1) 
  fval_soil <- pmax(fval_soil, 0) 
  
  in_dt$Eeq <- Eeq       
  in_dt$fval_soil <- fval_soil
  
  return(in_dt)
}


parread <- function(pathdata0, basin_name, basin_info_n, LAIpath, usage) {
  print(paste0('read ', basin_name))
  
   in_dt <- fread(paste0(pathdata0, basin_name, '.csv'))
  LAI <- fread(paste0(LAIpath, basin_name, '.csv'))
  
 
  lai_cols <- names(LAI)[c( 2, 4, 7, 8,9)] 
  LAI_sub <- LAI[, lai_cols,with = F]
  

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

parread_exps_mean_annual <- function(pathdata0, basin_name, basin_info_n, LAIpath, exp_type) {
  print(paste0('read ', basin_name))
  
  in_dt <- fread(paste0(pathdata0, basin_name, '.csv'))
  LAI <- fread(paste0(LAIpath, basin_name, '.csv'))
  
  
  lai_cols <- names(LAI)[c( 2, 4, 7, 8,9)] 
  LAI_sub <- LAI[, lai_cols,with = F]
  
  in_dt <- merge(in_dt, LAI_sub,  all = TRUE)
  
  INPUTS <- list()
  if(is.na(basin_info_n$Area)) {
    INPUTS$area <- basin_info_n$Area_shp
  } else {
    INPUTS$area <- basin_info_n$Area
  }
  
  Year_base <- 1981 
 
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
