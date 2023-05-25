model{
  # note: cannot redefine parameters in jags.model
  
  #### these are the proxy-derived data
  d13Cc.data ~ dnorm(d13Cc, 1 / 0.03 ^ 2)
  d18Oc.data ~ dnorm(d18Oc, 1 / 0.03 ^ 2)
  d13Co.data ~ dnorm(d13Co, 1 / 0.2 ^ 2)
  Dp17c.data ~ dnorm(Dp17c, 1 / 5 ^ 2)
  
  # isotope ratio constants
  R13_VPDB <- 0.011237
  R18_VSMOW <- 0.0020052
  R17_VSMOW <- 0.0003799
  R18_VPDB <- 0.0020672
  R17_VPDB <- 0.0003860
  alpha18_diff <- 1.028489
  alpha17_diff <- 1.014672
  theta_v_l_eq <- 0.529
  
  
  PPCQ <- MAP * PfPCQ # precipitation at pedogenic carbonate quarter
  TmPCQ <- MAT + TmPCQ_min_a # temperature at pedogenic carbonate quarter
  
  ### Depth to carbonate formation
  # top of Bk equation based on Retallack 2005 data - consider incorporate the uncertainty of the coefficients?
  z_min <- MAP * 0.0925 + 13.4
  # thickness of Bk, CMP as proxy for seasonality
  z_thick <- abs(PPCQ - MAP / 4) * 0.74 + 17.4
  # find middle of Bk in cm
  z_mean <- z_min + z_thick / 2
  #gamma scale parameter, using 22cm as variance from Retallack fit
  z_theta <- 22 ^ 2 / z_mean
  #gamma shape parameter
  z_shp <- z_mean / z_theta
  z ~ dgamma(z_shp, z_theta)
  
  # z <- 20 # for fixed value
  
  # depth in meters
  z_m <- z / 100
  
  ### Soil temperatures at depth z
  # Assume thermal conductivity = 0.0007 cal / cm2  s  *C, volumetric heat capacity of 0.3 cal / cm2 *C, Quade 2013
  d <- sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  # T is highest (avg summer temps) at t = 0.29, mean at t = 0 (avg spring and fall temps), low at t = 0.79 (avg winter temps)
  t <- 0.29 # assume pedogenic carbonate grew during warmest season
  Tsoil <- MAT + (TmPCQ_min_a * sin(2 * 3.1415 * t - z / d) / exp(z / d)) 
  Tsoil_K <- Tsoil + 273.15
  
  ### Relative humidity
  ## based on quarterly precipitation
  # h_m <- 0.25 + 0.7 * (PPCQ / 900) # Simmons (2010)
  # # Large variance
  # h_var <- 0.1 ^ 2
  # size <- h_m * (1 - h_m) / h_var - 1
  # h_alpha <- h_m * size
  # h_beta <- (1 - h_m) * size
  # h ~ dbeta(h_alpha, h_beta)
  RH <- h_pre * 100 # prescribed humidity
  
  ### Potential Evapotranspiration - Hargreaves and Samani (1982)
  # Solar radiation estimate
  Rs <- Ra * 0.16 * sqrt(12)
  # PET (mm/day) - Turc (1961)
  # Annual
  PET_A_D_m1 <- ifelse (RH < 50, 0.013 * (MAT / (MAT + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (MAT / (MAT + 15)) * (23.885 * Rs + 50))
  # PET_A_D_m <- max(PET_A_D_m1, 0.01)
  # PET_theta <- 0.26 ^ 2 / PET_A_D_m  
  # PET_shp <- PET_A_D_m / PET_theta
  # PET_A_D ~ dgamma(PET_shp, PET_theta)
  PET_A_A <- PET_A_D_m1 * 365 
  # PCQ
  PET_PCQ_D_m1 <- ifelse (RH < 50, 0.013 * (TmPCQ / (TmPCQ + 15)) * (23.8856 * Rs + 50)* (1 + ((50 - RH) / 70)), 0.0133 * (TmPCQ / (TmPCQ + 15)) * (23.885 * Rs + 50))
  PET_PCQ_D_m <- max(PET_PCQ_D_m1, 0.01)
  PET_PCQ_theta <- 0.26 ^ 2 / PET_PCQ_D_m  
  PET_PCQ_shp <- PET_PCQ_D_m / PET_PCQ_theta
  PET_PCQ_D ~ dgamma(PET_PCQ_shp, PET_PCQ_theta)
  PET_PCQ <- PET_PCQ_D * 90
  ## AET - actual evapotranspiration
  AET_var ~ dnorm(1, 1 / 0.2 ^ 2) # noise parameter - Gentine (2012)
  # AET in mm/quarter from Budyko curve - Pike (1964)
  AET_PCQ <- PPCQ * (1 / (sqrt(1 + (1 / ((PET_PCQ / (PPCQ)) * AET_var)) ^ 2)))
  
  # Free air porosity
  #Scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
  FAP1 <- min((pores - ((PPCQ - AET_PCQ) / (L * 10 * pores))), pores - 0.05)
  #At least 1% free air porosity
  FAP <- max(FAP1, 0.01)
  
  ### Respiration
  ## calculate CO2 production depth - Yang (2016)
  # assume equal to average rooting depth (cm) - Quade (2007)
  # AI <- PET_A_A / MAP
  # L <- ifelse(AI > 1.4, 60, -200 * AI ^ 2 + 250 * AI + 100) # yield too high S(z)
  # characteristic production depth (cm) - Quade (2007)
  k <- L / 2 / log(2)
  ## calculate soil respiration rate (gC/m2/d) - Raich (2002)
  #R_PCQ_D_m1 <- 1.25 * exp(0.05452 * TmPCQ) * PPCQ / (127.77 + PPCQ)  # Raich (2002)
  R_PCQ_D_m1 <- 1.25 * exp(0.07987 * MAT) * MAP / (29.86 + MAP) # regional model based on CLP data
  R_PCQ_D_m <- R_PCQ_D_m1 * f_R # optimized model - Fischer-Femal & Bowen (2020)
  R_PCQ_D <- R_PCQ_D_m
  # R_theta <- (R_PCQ_D_m * 0.5) ^ 2 / R_PCQ_D_m
  # R_shp <- R_PCQ_D_m / R_theta
  # R_PCQ_D ~ dgamma(R_shp, R_theta)
  # convert to molC/cm3/s
  R_PCQ_D1 <- R_PCQ_D / (12.01 * 100^2)  # from gC/m2/d to molC/cm2/d
  R_PCQ_S <- R_PCQ_D1 / (24 * 3600)  # molC/ cm2 / s
  R_PCQ_S_0 <- R_PCQ_S / (L * pores) # Quade et al. (2007)
  # soil CO2 diffusion coefficient
  Dair <- 0.1369 * (Tsoil_K / 273.15) ^ 1.958
  DIFC <- FAP * tort * Dair
  # S(z) 
  S_z_mol <- k ^ 2 * R_PCQ_S_0 / DIFC * (1 - exp(-z / k)) # (mol/cm3)
  S_z <- S_z_mol * (0.08206 * Tsoil_K * 10^9) # ppmv
  
  # fractionation of soil organic matter - PBUQ (Breecker, 2013)
  d13Co <- d13Cr + 1 + SOM.frac
  
  ### d13C of pedogenic carbonate
  d13Cs <- (pCO2 * (-6.5) + S_z * (1.0044 * d13Cr + 4.4))/(S_z + pCO2) # assumed a uniform d13Ca of -6.5‰
  # T-dependent fractionation equation - Romanek et al., 1992
  d13Cc <- ((1 + (11.98 - 0.12 * Tsoil) / 1000) * (d13Cs + 1000)) - 1000
  
  
  #### d18O of soil water
  # d18O - MAP based on modern observations across the CLP - Da (2023)
  # d18Op_m <- -0.0082 * MAP - 3.639
  d18p_m <- -7 # assign a fixed value
  d18p ~ dnorm(d18p_m, 1 / 3 ^ 2)
  R18p <- (d18p / 1000 + 1) * R18_VSMOW
  
  ### soil water evaporation
  # evaporation is 6+/-4% of AET - Good (2015)
  e_mean <- 0.06
  e_var <- 0.04^2
  e_size <- e_mean*(1-e_mean)/e_var - 1
  e_alpha <- e_mean * e_size
  e_beta <- (1-e_mean) * e_size
  ETR ~ dbeta(e_alpha, e_beta)
  E1 <- ETR * AET_PCQ
  # minimum of 1 mm
  E <- max(E1, 1)
  # Soil water diffusion evaporation balance
  # evaporation in m/sec
  E_s <- E / (1000 * 90 * 24 * 3600)
  
  # # saturated water vapor content 
  # P_sat <- 10^(8.07131 - 1730.63/(Tsoil + 233.426)) # Antoine Equation (mmHg)
  # n <- (133.322368 * P_sat * 10^-6)/(8.31 * Tsoil_K) # moles of water in the air (PV/RT: mol/cm3)
  # N_sat <- 18*n # (g/cm3)
  
  # self diffusivity coefficient of water - Holz (2000)
  DIFW <- 1.637e-8 * (Tsoil_K / 216.25 - 1) ^ 2.074
  # Diffusion, scaled to temperature, soil water content and tortuosity - Barnes & Allison (1983)
  DIFO <- DIFW * (pores - 0.05) * tort # assuming minimum water content of 0.05
  # mean penetration depth of evap, in m
  # density_water <- 1 # g/cm3
  # z_i <- DIFO * N_sat / (E_s * density_water)
  z_i <- DIFO / E_s
  # Diffusion ratio factor - Barnes & Allison (1983)
  frac_t <- 0.1 # the relative amount of turbulence above the saturated soil column
  DRF <- frac_t + (1 - frac_t) * (1 / 0.9723)
  
  ## calculating the triple oxygen isotopes - Passey & Levin (2021)
  
  # equilibrium fractionation factor between water vapor and liquid
  alpha18_v_l_eq <- exp(((1.137e6 / (Tsoil_K ^ 2) - 0.4156e3/Tsoil_K - 2.0667) /1000))
  alpha17_v_l_eq <- alpha18_v_l_eq ^ theta_v_l_eq
  
  # the steady-state isotopic composition that the water may reach at low fraction of remaining soil water
  dp18p <- log(d18p/1000 + 1) * 1000
  dp17p <- 0.5267 * dp18p + 0.013 # Aron (2021)
  R17p <- exp(dp17p / 1000) * R17_VSMOW
  
  ### water vapor pressure and relative humidity
  AI <- MAP / PET_A_A
  Tdew <- ifelse( # dew point temperature - Qiu et al. (2021)
    AI < 0.05, MAT - TmPCQ_min_a - 4,
    ifelse(
      AI < 0.2, MAT - TmPCQ_min_a - 2,
      ifelse(
        AI < 0.65, MAT - TmPCQ_min_a - 1,
        ifelse(
          AI < 1, MAT - TmPCQ_min_a, MAT - 2
        )
      )
    )
  )
  e_a <- 0.6108 * exp(17.27 * Tdew / (Tdew + 237.3)) * 10^3 # actual water vapor pressure (kPa)
  Tair_PCQ <- MAT + TmPCQ_min_a * sin(2 * 3.1415 * t)
  e_s <- 0.6112 * exp(17.67 * Tair_PCQ / (Tair_PCQ + 243.5)) * 10^3 # saturated water vapor pressure (kPa)
  h <- e_a / e_s
  
  # evaporation from isolated bodies of water
  R18v_ini <- R18p / alpha18_v_l_eq
  R17v_ini <- R17p / alpha17_v_l_eq
  h_z <- h + z_m / z_i # relative humidity at depth z
  R18wss <- (alpha18_v_l_eq * h * R18v_ini) / (1 - alpha18_v_l_eq * alpha18_diff * (1 - h))
  R17wss <- (alpha17_v_l_eq * h * R17v_ini) / (1 - alpha17_v_l_eq * alpha17_diff * (1 - h))
  
  u18 <- (1 - alpha18_v_l_eq * alpha18_diff * (1 - h)) / (alpha18_v_l_eq * alpha18_diff * (1 - h))
  u17 <- (1 - alpha17_v_l_eq * alpha17_diff * (1 - h)) / (alpha17_v_l_eq * alpha17_diff * (1 - h))
  
  R18surf <- f_soil ^ u18 * (R18p - R18wss) + R18wss
  R17surf <- f_soil ^ u17 * (R17p - R17wss) + R17wss
  
  # R18Ov <- R18Op / alpha18_v_l_eq
  # d18Ov_m <- (R18Ov / R18_VSMOW - 1) * 1000
  # d18Ov ~ dnorm(d18Ov_m, 1)
  # Surface water isotopes - Barnes & Allison (1983)
  # R18Osurf <- ((1 - h) * DRF * R18Op + h * R18Ov) * alpha18_v_l_eq # assuming infinite soil water reservoir
  # # Soil water isotopes at depth
  R18soil <- (R18surf - R18p) * exp(-z_m / z_i) + R18p
  R17soil <- (R17surf - R17p) * exp(-z_m / z_i) + R17p
  
  ### soil carbonate
  # Soil carbonate O isotope fractionation - Kim & O'neal (1997)
  # A_H20_Carb <- exp((1.803e4 / Tsoil_K - 32.42) / 1000)
  # Soil carbonate O isotope fractionation - Wostbrock (2020)
  alpha18_c_w_eq <- exp((1.61e4 / Tsoil_K - 24.6) / 1000)
  theta_c_w_eq <- 0.5305 - 1.39/Tsoil_K
  alpha17_c_w_eq <- alpha18_c_w_eq ^ theta_c_w_eq
  R18c <- R18soil * alpha18_c_w_eq
  R17c <- R17soil * alpha17_c_w_eq
  d18Oc <-(R18c/R18_VPDB - 1) * 1000
  dp18c <- log(R18c / R18_VSMOW) * 1000
  dp17c <- log(R17c / R17_VSMOW) * 1000
  Dp17c <- (dp17c - 0.528 * dp18c) * 1000
  
  #### set priors
  # z ~ dunif(20, 80) # soil carbonate formation depth (cm)
  Ra ~ dunif(28, 31) # radiation at the top of the atmosphere - based on modern CLP latitude
  MAP ~ dunif(150, 700) # mean annual precipitation (mm) - based on modern CLP data
  PfPCQ ~ dunif(0.5, 0.8) # fraction of precipitation at the pedogenic carbonate quarter
  MAT ~ dunif(3, 12) # mean annual temperature (celcius) - based on modern CLP data
  TmPCQ_min_a ~ dunif(10, 17) # difference between PCQ temperature and MAT - based on modern CLP data
  # Tsoil ~ dunif(10, 30) # carbonate formation temperature
  pCO2 ~ dunif(150, 600) # atmospheric CO2 level (ppm) - based on boron estimates between 8-2 Ma
  # d13Ca ~ dunif(-6.4, -6.6) # d13C of atmospheric CO2
  L ~ dunif(30, 60) # mean rooting depth
  # S_z ~ dunif(100, 2000) # soil-respired CO2 level at depth z
  d13Cr ~ dunif(-30, -20)
  SOM.frac ~ dunif(-0.5, 0.5)
  pores ~ dunif(0.35, 0.5) # porosity
  tort ~ dunif(0.6, 0.8) # tortuosity
  f_soil ~ dunif(0.1, 0.9) # the fraction of remaining soil water from evaporation
  f_R ~ dunif(0.1, 0.2) # the proportion of respired CO2 during PCQ
  h_pre ~ dunif(0.3, 0.7) # air humidity
}