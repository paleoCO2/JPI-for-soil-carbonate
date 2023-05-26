# setwd("/Users/jiawei/Desktop/R_code/PSM/PSM for soil carbonates/")

library(tidyverse)
library(R2jags)
set.seed(42)

#  These are the observed data
#data <- list("d13Cc.data" = d13Cc)

# These are the parameters you want recorded in the output
params <- c("pCO2", "S_z", "MAP", "f_soil", "h")

# run the inversion
# Run the inversion
data <- list("d13Cc.data" = -1.91, "d18Oc.data" = -9.12, "d13Co.data" = -22.3, "Dp17c.data" = -102)
jout <- jags(model.file = "carbPSM_inv.R", parameters.to.save = params,
                  data = data, inits = NULL, n.chains = 4, n.iter = 1e4,
                  n.burnin = 1e3, n.thin = 10)

traceplot(jout, col = c("blue", "red", "green", "purple"))

View(jout$BUGSoutput$summary)

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

s1 <- jout$BUGSoutput$sims.list
plot_u(100, 500, s1$pCO2)
plot_u(100, 1000, s1$S_z)
plot_u(0.1, 0.9, s1$f_soil)
plot_u(20, 50, s1$k)
