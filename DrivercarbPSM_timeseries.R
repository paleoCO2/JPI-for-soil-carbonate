setwd("/Users/jiawei/Documents/GitHub/JPI-for-soil-carbonate/data/")

library(tidyverse)
library(R2jags)

# Read in proxy time series data
set.seed(42)
prox.in <- read.csv("Jiaxian.csv")
prox.in <- prox.in[, c(1:9)]
names(prox.in) <- c("age", "d13Cc", "d13Ccsd", "d18Oc", "d18Ocsd", "d13Co", "d13Cosd", "dp17Oc", "dp17Ocsd")
prox.in$age <- prox.in$age / 1000 # ka to Ma

# d13Cc <- measurements$d13Cc
# d18Oc <- measurements$d18Oc
# d13Co <- measurements$d13Co
# Dp17Oc <- measurements$Dp17O

# read in the d13C of atmospheric CO2 from Tipple et al. (2010)
d13Ca <- read.csv("tipple_et_al_2010.csv") # age in Ma
d13Ca_interp <- approx(d13Ca$age, d13Ca$d13Ca, xout = prox.in$age, method = "linear")
prox.in$d13Ca <- d13Ca_interp$y

# These are the parameters you want recorded in the output
params <- c("pCO2", "h", "MAP", "S_z", "PfPCQ", "MAT", "f_soil", "z", "L", "f_R", "R_PCQ_D", "PET_A_A")


# run the inversion
# Run the inversion
jout <- list()
for (i in 1:length(d13Cc)) {
  data <- list("d13Cc.data" = d13Cc[i], "d18Oc.data" = d18Oc[i], "d13Co.data" = d13Co[i], "Dp17c.data" = Dp17Oc[i])
  jout[[i]] <- jags(model.file = "carbPSM_inv.R", parameters.to.save = params,
                    data = data, inits = NULL, n.chains = 3, n.iter = 1e4,
                    n.burnin = 1e3, n.thin = 10)
}

# Dp17O
for (i in 1:length(d13Cc)) {
  data <- list("d13Cc.data" = d13Cc[i], "d18Oc.data" = d18Oc[i], "d13Co.data" = d13Co[i], "Dp17c.data" = Dp17Oc[i])
  jout[[i]] <- jags(model.file = "d13CcarbPSM_inv_17O.R", parameters.to.save = params,
                    data = data, inits = NULL, n.chains = 4, n.iter = 1e4,
                    n.burnin = 1e3, n.thin = 10)
}

# Look at the results, first check visually for burn-in effect and 
# convergence
traceplot(jout[[1]], col = c("blue", "red", "green"))

# Look at the summary stats for the posterior
# looking here at the Rhat and sample size, which we'd like to be 
# < about 1.03 and > several hundred, ideally. Suggests we should
# run more iterations before publishing!
View(jout[[1]]$BUGSoutput$summary)

s1 <- data.frame(matrix(nrow = length(d13Cc), ncol = 7))
names(s1) <- c("CO2", "sd", "Sz", "MAP", "PET_A_A", "MAT", "h")
for (i in 1:length(d13Cc)) {
  s1$CO2[i] <- jout[[i]]$BUGSoutput$summary[12,1]
  s1$sd[i] <- jout[[i]]$BUGSoutput$summary[12,2]
  s1$Sz[i] <- jout[[i]]$BUGSoutput$summary[7,1]
  s1$MAP[i] <- jout[[i]]$BUGSoutput$summary[2,1]
  s1$PET_A_A[i] <- jout[[i]]$BUGSoutput$summary[4,1]
  s1$MAT[i] <- jout[[i]]$BUGSoutput$summary[3,1]
  s1$h[i] <- jout[[i]]$BUGSoutput$summary[11,1]
}
write.csv(s1, "posterior.csv")

# plot
plot_u = function(l, u, d){
  dens.pri = density(runif(1e6, l, u))
  dens.post = density(d)
  xrange = range(c(dens.pri$x, dens.post$x))
  yrange = range(c(dens.pri$y, dens.post$y))
  plot(c(l, l, u, u), c(0, rep(1 / (u-l), 2), 0), type = "l", 
       xlim = xrange, ylim = yrange, lty = 3)
  lines(dens.post)
}

for (i in 1:10) {
  s2 <- jout[[i]]$BUGSoutput$sims.list
  plot_u(150, 650, s2$pCO2)
}
