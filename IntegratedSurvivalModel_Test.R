library(nimble)

mySeed <- 88
set.seed(mySeed)

## WHAT'S NEW:
# - In this model, I am properly integrating the information from the Hoening

####################
# DATA PREPARATION #
####################

## Load different data sources
LT.data.all <- readRDS('220118_LifeTables.rds')
AaH.data.all <- readRDS('220118_AaHData.rds')
mO.param.data <- readRDS('220119_mO_HoeningParameters_Seal.rds')

## Select relevant data subset (females only, all time)
LT.data <- LT.data.all[[5]]
AaH.data <- AaH.data.all[[5]]
mO.data <- c(mO.param.data$log_mean[1], mO.param.data$log_sd[1])

## Truncate data as necessary
LT.data <- LT.data[2:(dim(LT.data)[1]-1),]
AaH.data <- subset(AaH.data, Age > 9)

## Prepare data and constants for nimble
ISM.data <- list(
  dx = LT.data[,'dx'], nx = LT.data[,'nx'],
  obsH = AaH.data$Number, Age = AaH.data$Age,
  mO.logmean = mO.data[1], mO.logsd = mO.data[2])

ISM.constants <- list(
  SA_Amax = 5, Amax_LT = dim(LT.data)[1],
  Amax_AaH = dim(AaH.data)[1])


##############
# MODEL CODE #
##############

## Nimble code for Bayesian life table analysis
nimble.ISM <- nimbleCode({
  
  # Data likelihoods #
  #------------------#
  
  ## Life Table Data
  for(a in 1:Amax_LT){
    dx[a] ~ dbin(qx[a], nx[a])
  }
  
  ## AaH Data (Catch curve analysis)
  for(a in 1:Amax_AaH){
    obsH[a] ~ dpois(expH[a])
  }
  
  # Priors and constraints #
  #------------------------#
  
  ## LTA: Age-specific finite mortality rates (qx)
  qx[1:SA_Amax] <- 1-S_SA[1:SA_Amax] # Subadult age classes
  qx[(SA_Amax+1):Amax_LT] <- 1-S_MA # Adult age classes
  
  ## CCA: Log-linear relationship of harvest numbers and age
  log(expH[1:Amax_AaH]) <- log(Mu) + beta.AaH*Age[1:Amax_AaH]
  
  ## CCA : Relationship of age slope and survival probability
  beta.AaH <- log(S_MA)
  
  ## mO Sim: Relationship between harvest mortality, natural mortality, and survival
  S_MA <- exp(-(mO_MA + mH_MA))
  
  ## Priors
  for(a in 1:SA_Amax){
    S_SA[a] ~ dunif(0, 1)
  }
  
  Mu ~ dunif(0, 250)
  
  mH_MA ~ dunif(0, 5)
  mO_MA ~ dlnorm(meanlog = mO.logmean, sdlog = mO.logsd)

})


##############
# MCMC SETUP #
##############

## Function for sampling initial values
init.fun <- function(){list(
  S_SA = runif(5, 0, 1),
  Mu = runif(1, 40, 70),
  mO_MA = runif(1, 0.01, 1),
  mH_MA = runif(1, 0.01, 0.5))}

## Sample initial values
#test.inits <- list(init.fun())
test.inits <- list(init.fun(), init.fun(), init.fun())

## Parameters to monitor
params <- c('S_SA', 'S_MA', 'Mu', 'mO_MA', 'mH_MA')

## MCMC setup
niter <- 20000
nburnin <- 4000
nthin <- 2
nchains <- 3


############
# TEST RUN #
############

## Run nimble model
testRun <- nimbleMCMC(code = nimble.ISM, 
                      constants = ISM.constants, 
                      data = ISM.data, 
                      inits = test.inits, 
                      monitors = params, 
                      niter = niter, nburnin = nburnin, nchains = nchains, thin = nthin,
                      samplesAsCodaMCMC = TRUE, setSeed = mySeed)

## Plot traces
plot(testRun, ask = T)

## Save MCMC samples
saveRDS(testRun, file = '220119_ISM_Exp2.rds')


