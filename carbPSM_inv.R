model{
  # note: cannot redefine parameters in jags.model
  
  #### these are the proxy-derived data
  for (i in 1:length(d13Cc.data)){
    d13Cc.data[i] ~ dnorm(d13Cc[i], d13Cc.p[i])
    d13Cc.p[i] = 1/d13Ccsd.data[i]^2 # Gaussian precision for d13Cc measurements from se of replicate analyses 
    d18Oc.data[i] ~ dnorm(d18Oc[i], d18Oc.p[i])
    d18Oc.p[i] <- 1/d18Ocsd.data[i]^2
    d13Co.data[i] ~ dnorm(d13Co[i], d13Co.p[i])
    d13Co.p[i] <- 1/d13Cosd.data[i]^2
  }

  # isotope ratio constants
  R13_VPDB <- 0.011237
  R18_VSMOW <- 0.0020052
  R17_VSMOW <- 0.0003799
  R18_VPDB <- 0.0020672
  R17_VPDB <- 0.0003860
  alpha18_diff <- 1.028489
  alpha17_diff <- 1.014672
  theta_v_l_eq <- 0.529
  # Damping term d for calculating soil temperature 
  # Assume thermal conductivity = 0.0007 cal / cm2  s  *C, volumetric heat capacity of 0.3 cal / cm2 *C, Quade 2013
  d <- sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  
  for (i in 1:length(d13Cc.data)) {
    PPCQ[i] <- MAP[i] * PfPCQ[i] # precipitation at pedogenic carbonate quarter
    TmPCQ[i] <- MAT[i] + TmPCQ_min_a[i] # temperature at pedogenic carbonate quarter
    
    ### Depth to carbonate formation
    # top of Bk equation based on Retallack (2005) data 
    # consider incorporate the uncertainty of the coefficients?
    # z_min[i] <- MAP[i] * 0.0925 + 13.4
    # # thickness of Bk, CMP as proxy for seasonality
    # z_thick[i] <- abs(PPCQ[i] - MAP[i] / 4) * 0.74 + 17.4
    # # find middle of Bk in cm
    # z_mean[i] <- z_min[i] + z_thick[i] / 2
    # #gamma scale parameter, using 22cm as variance from Retallack fit
    # z_theta[i] <- 22 ^ 2 / z_mean[i]
    # #gamma shape parameter
    # z_shp[i] <- z_mean[i] / z_theta[i]
    # z[i] ~ dgamma(z_shp[i], z_theta[i])
    
    z[i] <- 20 # for fixed value
    
    # depth in meters
    z_m[i] <- z[i] / 100
    
    ### Soil temperatures at depth z
    # T is highest (avg summer temps) at t = 0.29, mean at t = 0 (avg spring and fall temps), low at t = 0.79 (avg winter temps)
    # t[i] <- 0.29 # assume pedogenic carbonate grew during warmest season
    Tsoil[i] <- MAT[i] + (TmPCQ_min_a[i] * sin(2 * 3.1415 * t[i] - z[i] / d) / exp(z[i] / d)) 
    Tsoil_K[i] <- Tsoil[i] + 273.15
    
    ### Relative humidity
    ## based on quarterly precipitation
    # h_m[i] <- 0.25 + 0.7 * (PPCQ[i] / 900) # Simmons (2010)
    # # Large variance
    # h_var[i] <- 0.1 ^ 2
    # size[i] <- h_m[i] * (1 - h_m[i]) / h_var[i] - 1
    # h_alpha[i] <- h_m[i] * size[i]
    # h_beta[i] <- (1 - h_m[i]) * size[i]
    # h[i] ~ dbeta(h_alpha[i], h_beta[i])
    
    RH[i] <- h_pre[i] * 100 # prescribed humidity
    
    ### Potential Evapotranspiration - Hargreaves and Samani (1982)
    # Solar radiation estimate
    Rs[i] <- Ra[i] * 0.16 * sqrt(12)
    
    # PET (mm/day) - Turc (1961)
    PET_A_D_m[i] <- ifelse (RH[i] < 50, 0.013 * (MAT[i] / (MAT[i] + 15)) * (23.8856 * Rs[i] + 50)* (1 + ((50 - RH[i]) / 70)),
                            0.0133 * (MAT[i] / (MAT[i] + 15)) * (23.885 * Rs[i] + 50))
    # PET_A_D_m <- max(PET_A_D_m1, 0.01)
    # PET_theta <- 0.26 ^ 2 / PET_A_D_m  
    # PET_shp <- PET_A_D_m / PET_theta
    # PET_A_D ~ dgamma(PET_shp, PET_theta)
    PET_A_A[i] <- PET_A_D_m[i] * 365 # Annual
    # PCQ
    PET_PCQ_D_m1[i] <- ifelse (RH[i] < 50, 0.013 * (TmPCQ[i] / (TmPCQ[i] + 15)) * (23.8856 * Rs[i] + 50)* (1 + ((50 - RH[i]) / 70)),
                               0.0133 * (TmPCQ[i] / (TmPCQ[i] + 15)) * (23.885 * Rs[i] + 50))
    PET_PCQ_D_m[i] <- max(PET_PCQ_D_m1[i], 0.01)
    PET_PCQ_theta[i] <- 0.26 ^ 2 / PET_PCQ_D_m[i] 
    PET_PCQ_shp[i] <- PET_PCQ_D_m[i] / PET_PCQ_theta[i]
    PET_PCQ_D[i] ~ dgamma(PET_PCQ_shp[i], PET_PCQ_theta[i])
    PET_PCQ[i] <- PET_PCQ_D[i] * 90 # Seasonal
    
    ## AET - actual evapotranspiration
    AET_var[i] ~ dnorm(1, 1 / 0.2 ^ 2) # noise parameter - Gentine (2012)
    # AET in mm/quarter from Budyko curve - Pike (1964)
    AET_PCQ[i] <- PPCQ[i] * (1 / (sqrt(1 + (1 / ((PET_PCQ[i] / (PPCQ[i])) * AET_var[i])) ^ 2)))
    
    # Free air porosity
    #Scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
    FAP1[i] <- min((pores[i] - ((PPCQ[i] - AET_PCQ[i]) / (L[i] * 10 * pores[i]))), pores[i] - 0.05)
    #At least 1% free air porosity
    FAP[i] <- max(FAP1[i], 0.01)
    
    ### Respiration
    ## calculate CO2 production depth - Yang (2016)
    # assume equal to average rooting depth (cm) - Quade (2007)
    # AI <- PET_A_A / MAP
    # L <- ifelse(AI > 1.4, 60, -200 * AI ^ 2 + 250 * AI + 100) # yield too high S(z)
    # characteristic production depth (cm) - Quade (2007)
    k[i] <- L[i] / 2 / log(2)
    
    ## calculate soil respiration rate (gC/m2/d)
    #R_PCQ_D_m1 <- 1.25 * exp(0.05452 * TmPCQ) * PPCQ / (127.77 + PPCQ)  # Raich (2002)
    R_PCQ_D_m1[i] <- 1.25 * exp(0.07987 * MAT[i]) * MAP[i] / (29.86 + MAP[i]) # regional model based on CLP data
    R_PCQ_D_m[i] <- R_PCQ_D_m1[i] * f_R[i] # optimized model - Fischer-Femal & Bowen (2020)
    R_PCQ_D[i] <- R_PCQ_D_m[i]
    # R_theta <- (R_PCQ_D_m * 0.5) ^ 2 / R_PCQ_D_m
    # R_shp <- R_PCQ_D_m / R_theta
    # R_PCQ_D ~ dgamma(R_shp, R_theta)
    
    # convert to molC/cm3/s
    R_PCQ_D1[i] <- R_PCQ_D[i] / (12.01 * 100^2)  # from gC/m2/d to molC/cm2/d
    R_PCQ_S[i] <- R_PCQ_D1[i] / (24 * 3600)  # molC/ cm2 / s
    R_PCQ_S_0[i] <- R_PCQ_S[i] / (L[i] * pores[i]) # Quade et al. (2007)
    
    # soil CO2 diffusion coefficient
    Dair[i] <- 0.1369 * (Tsoil_K[i] / 273.15) ^ 1.958
    DIFC[i] <- FAP[i] * tort[i] * Dair[i]
    
    # S(z) 
    S_z_mol[i] <- k[i] ^ 2 * R_PCQ_S_0[i] / DIFC[i] * (1 - exp(-z[i] / k[i])) # (mol/cm3)
    S_z[i] <- S_z_mol[i] * (0.08206 * Tsoil_K[i] * 10^9) # ppmv
    
    ## estimate the d13Cr of soil-respired CO2
    # d13Ca[i] <- -6.5 # assumed a fixed value
    # DD13_water[i] <- 25.09 - 1.2 * (MAP[i] + 975) / (27.2 + 0.04 * (MAP[i] + 975))
    # D13C_plant[i] <- (28.26 * 0.22 * (pCO2[i] + 23.9)) / (28.26 + 0.22 * (pCO2[i] + 23.9)) - DD13_water # schubert & Jahren (2015)
    # d13Cr[i] <- d13Ca[i] - D13C_plant[i]
    # fractionation of soil organic matter - PBUQ (Breecker, 2013)
    d13Co[i] <- d13Cr[i] + 1 + SOM.frac[i]
    
    ### d13C of pedogenic carbonate
    d13Cs[i] <- (pCO2[i] * d13Ca[i] + S_z[i] * (1.0044 * d13Cr[i] + 4.4))/(S_z[i] + pCO2[i])
    # T-dependent fractionation equation - Romanek et al., 1992
    d13Cc[i] <- ((1 + (11.98 - 0.12 * Tsoil[i]) / 1000) * (d13Cs[i] + 1000)) - 1000
    
    
    #### d18O of soil water
    # d18Op_m[i] <- -0.0082 * MAP[i] - 3.639  #CLP observations - Da (2023)
    d18p_m[i] <- -7 # assign a fixed value
    d18p[i] ~ dnorm(d18p_m[i], 1 / 3 ^ 2)
    R18p[i] <- (d18p[i] / 1000 + 1) * R18_VSMOW
    
    ### soil water evaporation
    # evaporation is 6+/-4% of AET - Good (2015)
    e_mean[i] <- 0.06
    e_var[i] <- 0.04^2
    e_size[i] <- e_mean[i]*(1-e_mean[i])/e_var[i] - 1
    e_alpha[i] <- e_mean[i] * e_size[i]
    e_beta[i] <- (1-e_mean[i]) * e_size[i]
    ETR[i] ~ dbeta(e_alpha[i], e_beta[i])
    E1[i] <- ETR[i] * AET_PCQ[i]
    E[i] <- max(E1[i], 1) #minimum of 1 mm
    E_s[i] <- E[i] / (1000 * 90 * 24 * 3600) #evaporation in m/sec
    
    # # saturated water vapor content 
    # P_sat <- 10^(8.07131 - 1730.63/(Tsoil + 233.426)) # Antoine Equation (mmHg)
    # n <- (133.322368 * P_sat * 10^-6)/(8.31 * Tsoil_K) # moles of water in the air (PV/RT: mol/cm3)
    # N_sat <- 18*n # (g/cm3)
    
    # self diffusivity coefficient of water - Holz (2000)
    DIFW[i] <- 1.637e-8 * (Tsoil_K[i] / 216.25 - 1) ^ 2.074
    # Diffusion, scaled to temperature, soil water content and tortuosity - Barnes & Allison (1983)
    DIFO[i] <- DIFW[i] * 0.05 * tort[i] # assuming minimum water content of 0.05 during PCQ
    # mean penetration depth of evaporation
    # density_water <- 1 # g/cm3
    # z_i <- DIFO * N_sat / (E_s * density_water)
    z_i[i] <- DIFO[i] / E_s[i] # m
    
    ### calculating the triple oxygen isotopes
    
    # equilibrium fractionation factor between water vapor and liquid
    alpha18_v_l_eq[i] <- exp(((1.137e6 / (Tsoil_K[i] ^ 2) - 0.4156e3/Tsoil_K[i] - 2.0667) /1000))
    # alpha17_v_l_eq[i] <- alpha18_v_l_eq[i] ^ theta_v_l_eq
    
    # the steady-state isotopic composition that the water may reach at low fraction of remaining soil water
    dp18p[i] <- log(d18p[i]/1000 + 1) * 1000
    # dp17p[i] <- 0.5267 * dp18p[i] + 0.013 # Aron (2021)
    # R17p[i] <- exp(dp17p[i] / 1000) * R17_VSMOW
    
    ## water vapor pressure and relative humidity - Qiu et al. (2021)
    # determine Aridity Index first
    AI[i] <- MAP[i] / PET_A_A[i]
    # determine dew point temperature
    Tdew[i] <- ifelse(AI[i] < 0.05, MAT[i] - TmPCQ_min_a[i] - 4,
                      ifelse(AI[i] < 0.2, MAT[i] - TmPCQ_min_a[i] - 2,
                             ifelse(AI[i] < 0.65, MAT[i] - TmPCQ_min_a[i] - 1,
                                    ifelse(AI[i] < 1, MAT[i] - TmPCQ_min_a[i], MAT[i] - 2))))
    # actual water vapor pressure (kPa)
    e_a[i] <- 0.6108 * exp(17.27 * Tdew[i] / (Tdew[i] + 237.3)) * 10^3
    # saturated water vapor pressure (kPa)
    Tair_PCQ[i] <- MAT[i] + TmPCQ_min_a[i] * sin(2 * 3.1415 * t[i])
    e_s[i] <- 0.6112 * exp(17.67 * Tair_PCQ[i] / (Tair_PCQ[i] + 243.5)) * 10^3
    # air humidity
    h[i] <- e_a[i] / e_s[i]
    
    ## calculating triple oxygen isotopes at the evaporation front
    # assume soil water as isolated water body - Passey & Levin (2021)
    R18v_ini[i] <- R18p[i] / alpha18_v_l_eq[i]
    # R17v_ini[i] <- R17p[i] / alpha17_v_l_eq[i]
    h_z[i] <- h[i] + z_m[i] / z_i[i] # relative humidity at depth z
    R18wss[i] <- (alpha18_v_l_eq[i] * h[i] * R18v_ini[i]) / (1 - alpha18_v_l_eq[i] * alpha18_diff * (1 - h[i]))
    # R17wss[i] <- (alpha17_v_l_eq[i] * h[i] * R17v_ini[i]) / (1 - alpha17_v_l_eq[i] * alpha17_diff * (1 - h[i]))
    
    u18[i] <- (1 - alpha18_v_l_eq[i] * alpha18_diff * (1 - h[i])) / (alpha18_v_l_eq[i] * alpha18_diff * (1 - h[i]))
    # u17[i] <- (1 - alpha17_v_l_eq[i] * alpha17_diff * (1 - h[i])) / (alpha17_v_l_eq[i] * alpha17_diff * (1 - h[i]))
    
    R18surf[i] <- f_soil[i] ^ u18[i] * (R18p[i] - R18wss[i]) + R18wss[i]
    # R17surf[i] <- f_soil[i] ^ u17[i] * (R17p[i] - R17wss[i]) + R17wss[i]
    
    # R18Ov[i] <- R18Op[i] / alpha18_v_l_eq[i]
    # d18Ov_m[i] <- (R18Ov[i] / R18_VSMOW - 1) * 1000
    # d18Ov[i] ~ dnorm(d18Ov_m[i], 1)
    
    ## Soil water isotopes at depth z
    R18soil[i] <- (R18surf[i] - R18p[i]) * exp(-z_m[i] / z_i[i]) + R18p[i]
    # R17soil[i] <- (R17surf[i] - R17p[i]) * exp(-z_m[i] / z_i[i]) + R17p[i]
    
    ### soil carbonate isotopes
    # Soil carbonate O isotope fractionation - Kim & O'neal (1997)
    # alpha18_c_w_eq[i] <- exp((1.803e4 / Tsoil_K[i] - 32.42) / 1000)
    # Soil carbonate O isotope fractionation - Wostbrock (2020)
    alpha18_c_w_eq[i] <- exp((1.61e4 / Tsoil_K[i] - 24.6) / 1000)
    # theta_c_w_eq[i] <- 0.5305 - 1.39 / Tsoil_K[i]
    # alpha17_c_w_eq[i] <- alpha18_c_w_eq[i] ^ theta_c_w_eq[i]
    R18c[i] <- R18soil[i] * alpha18_c_w_eq[i]
    # R17c[i] <- R17soil[i] * alpha17_c_w_eq[i]
    d18Oc[i] <-(R18c[i] / R18_VPDB - 1) * 1000
    dp18c[i] <- log(R18c[i] / R18_VSMOW) * 1000
    # dp17c[i] <- log(R17c[i] / R17_VSMOW) * 1000
    # Dp17c[i] <- (dp17c[i] - 0.528 * dp18c[i]) * 1000
    
    #### set priors
    # z[i] ~ dunif(20, 80) # soil carbonate formation depth (cm)
    t[i] ~ dunif(0.25, 0.3) # timing of soil carbonate formation
    Ra[i] ~ dunif(28, 31) # radiation at the top of the atmosphere - based on modern CLP latitude
    MAP[i] ~ dunif(150, 700) # mean annual precipitation (mm) - based on modern CLP data
    PfPCQ[i] ~ dunif(0.5, 0.8) # fraction of precipitation at the pedogenic carbonate quarter
    MAT[i] ~ dunif(3, 12) # mean annual temperature (celcius) - based on modern CLP data
    TmPCQ_min_a[i] ~ dunif(10, 17) # difference between PCQ temperature and MAT - based on modern CLP data
    # Tsoil[i] ~ dunif(10, 30) # carbonate formation temperature
    pCO2[i] ~ dunif(150, 600) # atmospheric CO2 level (ppm) - based on boron estimates between 8-2 Ma
    L[i] ~ dunif(30, 60) # mean rooting depth
    # S_z[i] ~ dunif(100, 2000) # soil-respired CO2 level at depth z
    d13Cr[i] ~ dunif(-30, -20)
    SOM.frac[i] ~ dunif(-0.5, 0.5)
    pores[i] ~ dunif(0.35, 0.5) # porosity
    tort[i] ~ dunif(0.6, 0.8) # tortuosity
    f_soil[i] ~ dunif(0.1, 0.9) # the fraction of remaining soil water from evaporation
    f_R[i] ~ dunif(0.1, 0.2) # the proportion of respired CO2 during PCQ
    h_pre[i] ~ dunif(0.3, 0.7) # air humidity
  }
}