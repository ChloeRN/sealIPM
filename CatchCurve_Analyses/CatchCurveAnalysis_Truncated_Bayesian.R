library(nimble)

mySeed <- 88
set.seed(mySeed)


## Load data
AaH.data.all <- readRDS('CatchCurve_Analyses/AaHData.rds')

## Make a list to store MCMC samples for all different model runs
CCA.MCMC <- list(
  FM.All = NULL,
  FM.P1 = NULL,
  FM.P2 = NULL,
  FM.P3 = NULL,
  F.All = NULL,
  F.P1 = NULL,
  F.P2 = NULL,
  F.P3 = NULL
)

## Nimble code for Bayesian catch curve analysis
nimble.CCA <- nimbleCode({
  
  # Data likelihood #
  #-----------------#
  
  for(a in 1:Amax_AaH){
    obsH[a] ~ dpois(expH[a])
  }

  # Priors and constraints #
  #------------------------#
  
  ## Log-linear relationship of harvest numbers and age
  log(expH[1:Amax_AaH]) <- log(Mu) + beta.AaH*Age[1:Amax_AaH]
  
  ## Relationship of age slope and survival probability
  beta.AaH <- log(S_MA)
  
  # Priors
  Mu ~ dunif(0, 250)
  S_MA ~ dunif(0, 1)
  
})


## Function for sampling initial values
init.fun <- function(){list(
  Mu = runif(1, 40, 70),
  S_MA = runif(1, 0, 1))}

## Sample initial values
#test.inits <- list(init.fun())
test.inits <- list(init.fun(), init.fun(), init.fun())

## Parameters to monitor
params <- c('Mu', 'S_MA')

## MCMC setup
niter <- 10000
nburnin <- 2000
nthin <- 1
nchains <- 3

for(i in 1:8){
  ## Select dataset to use
  AaH.data <- AaH.data.all[[i]]
  
  ## Truncate dataset
  #  Start age = 10, because that is after the "saddle" of the data
  AaH.data <- subset(AaH.data, Age > 9)
  
  ## Prepare data and constants for nimble
  CCA.data <- list(obsH = AaH.data$Number, Age = AaH.data$Age)
  
  CCA.constants <- list(Amax_AaH = dim(AaH.data)[1])
  
  ## Test run 
  testRun <- nimbleMCMC(code = nimble.CCA, 
                        constants = CCA.constants, 
                        data = CCA.data, 
                        inits = test.inits, 
                        monitors = params, 
                        niter = niter, nburnin = nburnin, nchains = nchains, thin = nthin,
                        samplesAsCodaMCMC = TRUE, setSeed = mySeed)
  
  ## Store MCMC samples
  CCA.MCMC[[i]] <- testRun
  
}

## Save MCMC data
saveRDS(CCA.MCMC, file = 'CatchCurve_Analyses/CatchCurveAnalysis_Truncated_MCMC.rds')
#pdf('Traceplots.pdf', height = 10, width = 8)
#plot(testRun)
#dev.off()
