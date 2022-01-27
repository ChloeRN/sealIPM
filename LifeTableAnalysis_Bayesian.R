library(nimble)

mySeed <- 88
set.seed(mySeed)


## Load data
LT.data.all <- readRDS('220118_LifeTables.rds')

## Make a list to store MCMC samples for all different model runs
LTA.MCMC <- list(
  FM.All = NULL,
  FM.P1 = NULL,
  FM.P2 = NULL,
  FM.P3 = NULL,
  F.All = NULL,
  F.P1 = NULL,
  F.P2 = NULL,
  F.P3 = NULL
)

## Nimble code for Bayesian life table analysis
nimble.LTA <- nimbleCode({
  
  # Data likelihood #
  #-----------------#
  
  for(a in 1:Amax_LT){
    dx[a] ~ dbin(qx[a], nx[a])
  }

  # Priors and constraints #
  #------------------------#
  
  ## Age-specific finite mortality rates (qx)
  qx[1:SA_Amax] <- 1-S_SA[1:SA_Amax] # Subadult age classes
  qx[(SA_Amax+1):Amax_LT] <- 1-S_MA # Adult age classes
  
  # Priors
  for(a in 1:SA_Amax){
    S_SA[a] ~ dunif(0, 1)
  }
  S_MA ~ dunif(0, 1)
  
})


## Function for sampling initial values
init.fun <- function(){list(
  S_SA = runif(5, 0, 1),
  S_MA = runif(1, 0, 1))}

## Sample initial values
#test.inits <- list(init.fun())
test.inits <- list(init.fun(), init.fun(), init.fun())

## Parameters to monitor
params <- c('S_SA', 'S_MA')

## MCMC setup
niter <- 10000
nburnin <- 2000
nthin <- 1
nchains <- 3

for(i in 1:8){
  ## Select dataset to use
  LT.data <- LT.data.all[[i]]
  
  ## Remove first and last row
  #  First row is removed because we expect sampling og pups of the year to be non-representative
  #  Last row is removed because it contains no information
  LT.data <- LT.data[2:(dim(LT.data)[1]-1),]
  
  ## Prepare data and constants for nimble
  LTA.data <- list(dx = LT.data[,'dx'], nx = LT.data[,'nx'])
  
  LTA.constants <- list(SA_Amax = 5, Amax_LT = dim(LT.data)[1])
  
  ## Test run 
  testRun <- nimbleMCMC(code = nimble.LTA, 
                        constants = LTA.constants, 
                        data = LTA.data, 
                        inits = test.inits, 
                        monitors = params, 
                        niter = niter, nburnin = nburnin, nchains = nchains, thin = nthin,
                        samplesAsCodaMCMC = TRUE, setSeed = mySeed)
  
  ## Store MCMC samples
  LTA.MCMC[[i]] <- testRun
  
}

## Save MCMC data
saveRDS(LTA.MCMC, file = '220118_LifeTableAnalysis_MCMC.rds')
#pdf('Traceplots.pdf', height = 10, width = 8)
#plot(testRun)
#dev.off()
