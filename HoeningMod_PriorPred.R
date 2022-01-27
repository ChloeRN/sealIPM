########################################################
#### SEAL NATURAL MORTALITY PRIORS FOR IPM ANALYSIS ####
########################################################

library(coda)
library(ggplot2)
library(nimble)

mySeed <- 88
set.seed(mySeed)

#**********#
# 1) SETUP #
#**********#

## Function to calculate prior from max age (Tom Porteus' Hoening Model)
mOpred_hoening = function(maxAge, param.a, param.b, sigma.z, N, mu.t.max){
  
  # Set the number of posterior samples
  nsamples <- length(param.a)
  
  # Set the number of simulations per posterior sample
  nsim <- N
  
  # Prepare a data frame to store results
  data <- data.frame(sampleNo = NA, simNo = NA, mO_est = NA)
  
  # Loop over posterior samples and simulate mO
  for(i in 1:nsamples){
    
    log_mu.mO <- param.a[i] + param.b[i] * (log(maxAge) - log(mu.t.max))
    
    mO <- rlnorm(nsim, meanlog = log_mu.mO, sdlog = sigma.z[i])
    data <- rbind(data, data.frame(sampleNo = rep(i, nsim), 
                                    simNo = 1:nsim, 
                                    mO_est = mO))
  } # i
  
  # Remove first row (NA) and return data
  data <- data[-1,]
  return(data)
}

## Load posterior samples for Hoening model parameters
#  (obtained from Tom Porteus 01/12/2020)
PostSamples <- read.table('Hoenig_Posteriors_NoTitle.txt', header = T)
mu.t.max <- 22.61062


#*******************************#
# 3) SIMULATE PRIOR FOR RED FOX #
#*******************************#

## Set maximum age for Svalbard ringed seal females (3 periods)
maxAge_RS <- c(45, 29, 33)

## Prepare vectors to save parameter values
log_mean <- log_sd <- rep(NA, 3)

## Simulate natural mortality prior for period 1
mO_RS_P1 <- mOpred_hoening(maxAge = maxAge_RS[1], 
                           param.a = PostSamples[,'a_1'], 
                           param.b = PostSamples[,'b_1'], 
                           sigma.z = PostSamples[,'sigma.z_1'], 
                           N = 30, mu.t.max = mu.t.max)

## Extract log-mean and log-sd for simulation (period 1)
log_mean[1] <- mean(log(mO_RS_P1$mO_est))
log_sd[1] <- sd(log(mO_RS_P1$mO_est))
rm(mO_RS_P1)

## Simulate natural mortality prior for period 2
mO_RS_P2 <- mOpred_hoening(maxAge = maxAge_RS[2], 
                           param.a = PostSamples[,'a_1'], 
                           param.b = PostSamples[,'b_1'], 
                           sigma.z = PostSamples[,'sigma.z_1'], 
                           N = 30, mu.t.max = mu.t.max)

## Extract log-mean and log-sd for simulation (period 1)
log_mean[2] <- mean(log(mO_RS_P2$mO_est))
log_sd[2] <- sd(log(mO_RS_P2$mO_est))
rm(mO_RS_P2)

## Simulate natural mortality prior for period 3
mO_RS_P3 <- mOpred_hoening(maxAge = maxAge_RS[3], 
                           param.a = PostSamples[,'a_1'], 
                           param.b = PostSamples[,'b_1'], 
                           sigma.z = PostSamples[,'sigma.z_1'], 
                           N = 30, mu.t.max = mu.t.max)

## Extract log-mean and log-sd for simulation (period 1)
log_mean[3] <- mean(log(mO_RS_P3$mO_est))
log_sd[3] <- sd(log(mO_RS_P3$mO_est))
rm(mO_RS_P3)

## Save distribution parameter estimates
mO.params <- list(log_mean = log_mean, log_sd = log_sd, 
                  Period = c('1981-1982', '2002-2004', '2012-2020'),
                  MaxAge = maxAge_RS)
saveRDS(mO.params, file = '220119_mO_HoeningParameters_Seal.rds')

## Simulate prior values using MCMC
Hoening.Sim <- nimbleCode({
  for(p in 1:3){
    mO[p] ~ dlnorm(meanlog = mO.logmean[p], sdlog = mO.logsd[p])
  }
})

init.sim <- function(){list(mO = runif(3, 0.01, 1))}
test.inits <- list(init.sim(), init.sim(), init.sim())

mO.sim <- nimbleMCMC(code = Hoening.Sim, 
                     constants = list(), 
                     data = list(mO.logmean = log_mean, mO.logsd = log_sd), 
                     inits = test.inits, 
                     monitors = c('mO'), 
                     niter = 10000, nburnin = 2000, nchains = 3, thin = 1,
                     samplesAsCodaMCMC = TRUE)

## Save MCMC samples
saveRDS(mO.sim, file = '220119_mO_HoeningSim.rds')
