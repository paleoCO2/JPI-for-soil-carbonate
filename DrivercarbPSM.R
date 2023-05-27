setwd("/Users/jiawei/Documents/GitHub/JPI-for-soil-carbonate/")

library(tidyverse)
library(R2jags)

############################################################################################
#    INPUT TO PASS TO JAGS
############################################################################################  

### Read in proxy time series data
set.seed(42)
prox.in <- read.csv("data/Jiaxian.csv")
prox.in <- prox.in[, c(1:9)]
names(prox.in) <- c("age", "d13Cc", "d13Ccsd", "d18Oc", "d18Ocsd", "d13Co", "d13Cosd", "dp17Oc", "dp17Ocsd")
prox.in <- prox.in %>%
  drop_na(age, d13Cc, d18Oc, d13Co)
prox.in$age <- prox.in$age / 1000 # ka to Ma

# read in the d13C of atmospheric CO2 from Tipple et al. (2010)
d13Ca <- read.csv("data/tipple_et_al_2010.csv") # age in Ma
d13Ca_interp <- approx(d13Ca$age, d13Ca$d13Ca, xout = prox.in$age, method = "linear")
prox.in$d13Ca <- d13Ca_interp$y

### other inputs
# Damping term d for calculating soil temperature at depth
# Assume thermal conductivity = 0.0007 cal / cm2  s  *C, volumetric heat capacity of 0.3 cal / cm2 *C, Quade 2013
d <- sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))

# These are the parameters you want to record in the output
params <- c("pCO2", "h", "MAP", "S_z", "PfPCQ", "MAT", "f_soil", "z", "L", "f_R", "R_PCQ_D", "PET_A_A")

# Data inputs to jags 
data <- list("d13Cc.data" = prox.in$d13Cc,
             "d13Ccsd.data" = prox.in$d13Ccsd,
             "d18Oc.data" = prox.in$d18Oc,
             "d18Ocsd.data" = prox.in$d18Ocsd,
             "d13Co.data" = prox.in$d13Co,
             "d13Cosd.data" = prox.in$d13Cosd,
             "d13Ca" = prox.in$d13Ca
)

# Run the inversion
jout <- jags(model.file = "carbPSM_inv.R", parameters.to.save = params,
                     data = data, inits = NULL, n.chains = 3, n.iter = 1e4,
                     n.burnin = 1e3, n.thin = 10)


# Look at the results, first check visually for burn-in effect and convergence
traceplot(jout, col = c("blue", "red", "green"))

# Look at the summary stats for the posterior
# looking here at the Rhat and sample size, which we'd like to be 
# < about 1.03 and > several hundred, ideally. Suggests we should
# run more iterations before publishing!
View(jout$BUGSoutput$summary)

MCMC_mean <- do.call(cbind, jout$BUGSoutput$mean)
MCMC_sd <- do.call(cbind, jout$BUGSoutput$sd)
write.csv(MCMC_mean, "data/MCMC_mean.csv")
write.csv(MCMC_sd, "data/MCMC_sd.csv")
