library(nimble)

mySeed <- 0
set.seed(mySeed)

#################
# DATA ASSEMBLY #
#################

## Set paths
#DataPath <- '/data/P-Prosjekter/41201625_sustainable_harvesting_of_seals_in_svalbard/Code/'
DataPath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/Code/'

#CodePath <- '/data/P-Prosjekter/41201625_sustainable_harvesting_of_seals_in_svalbard/SealIPM/'
CodePath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/sealIPM/'

## Load data
IceData <- readRDS(paste0(DataPath, '220509_SealIPM_IceData.rds'))

## Set model switches

# Ice model type switch
# FALSE = period
# TRUE = Trend
TrendIceModel <- TRUE
#TrendIceModel <- FALSE

# Ice/lair habitat covariate switch
# TRUE = pup survival is a function of sea ice availability
# FALSE = pup survival is a function of lair habitat availability
IceCov <- TRUE
# IceCov <- FALSE

# Environmental change scenario switch
# TRUE = Simulated future years will have same env. conditions on average as last year with data (2020)
# FALSE = Simulated future years will have env. conditions that continue deteriorating
StableFuture <- TRUE
#StableFuture <- FALSE


## Set parameters that are fixed a priori (but may be varied)

# Number of years to simulate model beyond data
sim_extraYears <- 10

# Start and end years for population model
sim_Tmin <- 2002 - 1981 + 1 
sim_Tmax <- length(1981:2020) + sim_extraYears

# Pup survival
Mu.S_pup.ideal <- 0.80 # pup survival under ideal conditions
sdlog.m_pup <- 0.20 # standard deviation of uncertainty in pup mortality hazard rate (on the log scale)

## Set env. covariate incl. ideal conditions
if(IceCov){
  env.cov <- c(IceData$ice.GL.mean, rep(NA, sim_extraYears-1))
}else{
  env.cov <- c(IceData$lairH.GL.mean, rep(NA, sim_extraYears-1))
}

env.period <- c(IceData$ice.period, rep(3, sim_extraYears))

if(sim_extraYears < 1){
  env.cov <- env.cov[1:sim_Tmax]
  env.period <- env.period[1:sim_Tmax]
}

env.ideal <- quantile(env.cov[which(IceData$ice.years < 2006)], probs = 0.75, na.rm = TRUE) # threshold env. conditions which are considered ideal
sdlink.env_ideal <- 0.1 # standard deviation of uncertainty in ideal env. conditions on the link scale (= log)


## Make dummy year covariate
if(StableFuture){
  year.idx <- c(1:length(1981:2020), rep(length(1981:2020), sim_extraYears))
}else{
  year.idx <- c(1:sim_Tmax)
}


## Assemble data and constants
seal.data <- list(
  env = env.cov,
  year.idx = year.idx
)

seal.constants <- list(
   
  Tmax = length(1981:2020), # Total length of the study period
  sim_Tmin = sim_Tmin, # Start year for population model
  sim_Tmax = sim_Tmax, # End year for population model
  
  Pmax = 3, 
  
  Mu.S_pup.ideal = Mu.S_pup.ideal,
  sdlog.m_pup = sdlog.m_pup,
  
  env.period = env.period,
  Mu.env.ideal = env.ideal,
  sdlink.env_ideal = sdlink.env_ideal
)


####################
# NIMBLE FUNCTIONS #
####################


## Nimble function for calculating pup survival
calc.S_pup <- nimbleFunction(
  
  run = function(S_pup.ideal = double(0),
                 env.ideal = double(0),
                 env = double(0)) {
    
    # Calculate environment factor
    if(env <= env.ideal){
      P_env <- env/env.ideal
    }else{
      P_env <- 1
    }
    
    # Calculate survival probability and multiply by environment factor
    S_pup <- S_pup.ideal*P_env
    
    # Return matrix
    return(S_pup)
    returnType(double(0))
  }
)


##############
# MODEL CODE #
##############

ice.Mod <- nimbleCode({

  
  #*************#
  # CONSTRAINTS #
  #*************#
  
  # Pup survival #
  #--------------#
  
  # Year-specific pup survival (dependent on sea ice)
  for(t in 1:sim_Tmax){
    S_pup[t] <- calc.S_pup(S_pup.ideal = S_pup.ideal, env.ideal = env.ideal, env = env[t])
  }
  
  
  #*********#
  # PRIORS  #                                                          
  #*********#
  
  ## Vital rate averages
  
  # Pup survival under ideal conditions
  m_pup.ideal ~ dlnorm(meanlog = log(-log(Mu.S_pup.ideal)), sdlog = sdlog.m_pup)
  S_pup.ideal <- exp(-m_pup.ideal)
  
  
  ## Fixed effects
  env.ideal ~ dlnorm(meanlog = log(Mu.env.ideal), sdlog = sdlink.env_ideal)


  
  #********************************#
  # ENVIRONMENTAL COVARIATE MODEL  #                                                          
  #********************************#
  
  ## Data likelihood (imputation of missing values)
  for(t in 1:sim_Tmax){
    env[t] ~ dlnorm(meanlog = log(env.pred[t]), sdlog = sigmaY.env)
  }
  
  if(TrendIceModel){
    
    ## Ice model - Trend
    log(env.pred[1:sim_Tmax]) <- log(Mu.env) + beta.env*(year.idx[1:sim_Tmax])
    
    Mu.env ~ dunif(0, Mu.env.ideal*2)
    beta.env ~ dunif(-5, 0)
    
  }else{

    ## Ice model - Period
    for(t in 1:sim_Tmax){
      env.pred[t] <- Mu.env[env.period[t]]
    }
    
    for(p in 1:3){
      if(IceCov){
        Mu.env[p] ~ dunif(0, 100)
      }else{
        Mu.env[p] ~ dunif(0, Mu.env.ideal*2)
      }
    }
  }
  
  sigmaY.env ~ dunif(0, 10)
  
})


###########################
# INITAL VALUE SIMULATION #
###########################

## Function for simulating initial values for icemodel
initIceSim <- function(data, constants){
  
  # 0) Set parameters for ice model #
  #---------------------------------#
  
  # 0.1) Missing covariate values
  env <- rep(NA, length(data$env))
  env.cov <- data$env
  for(t in 1:length(env)){
    if(is.na(data$env[t])){
      env.upper <- ifelse(IceCov, 1, max(data$env, na.rm = TRUE))
      env[t] <- env.cov[t] <- EnvStats::rlnormTrunc(1, mean = mean(log(data$env), na.rm = T), sd = sd(log(data$env), na.rm = T), max = env.upper)
    }
  }
  
  
  # 0.2) Ice model parameters
  if(TrendIceModel){
    Mu.env <- runif(1, max(data$env, na.rm = T)*0.9, max(data$env, na.rm = T)*1.1)
  }else{
    Mu.env <- runif(3, mean(data$env, na.rm = T)*0.9, mean(data$env, na.rm = T)*1.1)
  }
  
  beta.env <- runif(1, -0.1, 0)
  sigmaY.env <- runif(1, 0, 2)
  
  #-----------------------------------------------------------------------
  
  # 1) Simulate vital rates #
  #-------------------------#
  
  # 1.1) Vital rate averages
  
  ## Pup survival
  S_pup.ideal <- constants$Mu.S_pup.ideal
  
  
  # 1.3) Time-dependent vital rates
  
  ## Pup survival
  S_pup <- rep(NA, constants$sim_Tmax)
  env.ideal <- constants$Mu.env.ideal
  
  for(t in 1:constants$sim_Tmax){
    S_pup[t] <- calc.S_pup(S_pup.ideal = S_pup.ideal, env.ideal = env.ideal, env = env.cov[t])
  }
  
  
  #-----------------------------------------------------------------------
  
  # 7) Assemble and return initial values #
  #---------------------------------------#
  
  InitVals <- list(
    S_pup.ideal = S_pup.ideal, m_pup.ideal = -log(S_pup.ideal),

    S_pup = S_pup, env.ideal = unname(env.ideal),
    
    env = env,
    Mu.env = Mu.env, beta.env = beta.env, sigmaY.env = sigmaY.env
  )
  
  return(InitVals)
}
  
  
## Sample initial values for ice model
#Inits <- list(initValSim(data = seal.data, constants = seal.constants))
Inits <- list(initIceSim(data = seal.data, constants = seal.constants),
              initIceSim(data = seal.data, constants = seal.constants),
              initIceSim(data = seal.data, constants = seal.constants))


############
# TEST RUN #
############

## Set parameters to monitor
params <- c('S_pup.ideal', 'S_pup', 'env.ideal')

## Add parameters for sea ice model
if(TrendIceModel){
  params <- c(params, c('env', 'Mu.env', 'beta.env', 'sigmaY.env'))
}else{
  params <- c(params, c('env', 'Mu.env', 'sigmaY.env'))
}

## MCMC specs
#niter <- 10
#nburnin <- 0
#nthin <- 1
#nchains <- 3
#nchains <- 1

niter <- 10000
nburnin <- 3000
nthin <- 1
nchains <- 3

## Testrun
testRun <- nimbleMCMC(code = ice.Mod, 
                      constants = seal.constants, 
                      data = seal.data, 
                      inits = Inits, 
                      monitors = params, 
                      niter = niter, nburnin = nburnin, nchains = nchains, thin = nthin,
                      samplesAsCodaMCMC = TRUE, setSeed = mySeed)


plot(testRun, ask = TRUE)

