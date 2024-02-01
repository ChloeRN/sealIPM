library(nimble)

mySeed <- 0
set.seed(mySeed)

#################
# DATA ASSEMBLY #
#################

## Set paths
#DataPath <- '/data/P-Prosjekter/41201625_sustainable_harvesting_of_seals_in_svalbard/Code/'
DataPath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/Code/'
IceData <- readRDS(paste0(DataPath, '220428_SealIPM_IceData.rds'))

## Set parameters that are fixed a priori (but may be varied)
# Pup survival ~ sea ice
Mu.S_pup.ideal <- 0.60 # pup survival under ideal conditions
sdlog.m_pup <- 0.20 # standard deviation of uncertainty in pup mortality hazard rate (on the log scale)
ice.ideal <- median(IceData$fastice.mean[which(IceData$ice.period == 1)])# threshold ice extent above which conditions are considered ideal

# First-year survival
S_YOY.fix <- 0.75

# Ice model type switch (FALSE = period, TRUE = Trend)
#TrendIceModel <- FALSE
TrendIceModel <- TRUE

## Make ice covariate
ice.cov <- IceData$fastice.mean
ice.cov[which(IceData$ice.dataKeep==0)] <- NA

## Assemble data and constants
ice.data <- list(
  ice = ice.cov
)

ice.constants <- list(
  sim_Tmax = length(1981:2021), # End year for population model
  Mu.S_pup.ideal = Mu.S_pup.ideal,
  sdlog.m_pup = sdlog.m_pup,
  ice.period = IceData$ice.period,
  ice.ideal = ice.ideal
)

####################
# NIMBLE FUNCTIONS #
####################

## Nimble function for calculating pup survival
calc.S_pup <- nimbleFunction(
  
  run = function(S_pup.ideal = double(0),
                 ice.ideal = double(0),
                 ice = double(0)) {
    
    # Calculate sea ice factor
    if(ice <= ice.ideal){
      P_ice <- ice/ice.ideal
    }else{
      P_ice <- 1
    }
    
    # Calculate survival probability and multiply by ice factor
    S_pup <- S_pup.ideal*P_ice
  
    # Return matrix
    return(S_pup)
    returnType(double(0))
  }
)

#############
# FIT MODEL #
#############

## Model code
ice.model <- nimbleCode({
  
  #********************************#
  # PUP SURVIVAL PREDICTION MODEL  #                                                          
  #********************************#
  
  # Pup survival under ideal conditions
  m_pup.ideal ~ dlnorm(meanlog = log(-log(Mu.S_pup.ideal)), sdlog = sdlog.m_pup)
  S_pup.ideal <- exp(-m_pup.ideal)
  
  # Year-specific pup survival
  for(t in 1:sim_Tmax){
    S_pup[t] <- calc.S_pup(S_pup.ideal = S_pup.ideal, ice.ideal = ice.ideal, ice = ice[t])
  }

  #**************************#
  # SEA ICE COVARIATE MODEL  #                                                          
  #**************************#
  
  ## Data likelihood (imputation of missing values)
  for(t in 1:sim_Tmax){
    ice[t] ~ dlnorm(meanlog = log(ice.pred[t]), sdlog = sigmaY.ice)
  }
  
  if(TrendIceModel){
    ## Ice model - Trend
    log(ice.pred[1:sim_Tmax]) <- log(Mu.ice) + beta.ice*(1:sim_Tmax)
    Mu.ice ~ dunif(0, ice.ideal*2)
    beta.ice ~ dunif(-5, 0)
  
  } else {
    ## Ice model - Period
    for(t in 1:sim_Tmax){
      log(ice.pred[t]) <- log(Mu.ice[ice.period[t]])
    }
    
    for(p in 1:3){
      Mu.ice[p] ~ dunif(0, ice.ideal*2)
    }
  }
  
  sigmaY.ice ~ dunif(0, 10)
})

## Initial values
ice.inits <- function(ice){
  ice.inits <- rep(NA, length(ice))
  for(t in 1:length(ice)){
    if(is.na(ice[t])){
      ice.inits[t] <- rnorm(1, mean = mean(ice, na.rm = T), sd = sd(ice, na.rm = t))
    }
  }
  return(ice.inits)
}

initVal.sim <- function(ice.data, TrendIceModel){
  
  if(TrendIceModel){
    Inits <- list(ice = ice.inits(ice.data$ice),
                  Mu.ice = runif(1, max(ice.data$ice, na.rm = T)*0.9, max(ice.data$ice, na.rm = T)*1.1),
                  beta.ice = runif(1, -0.1, 0),
                  sigmaY.ice = runif(1, 0, 2),
                  S_pup.ideal = ice.constants$Mu.S_pup.ideal,
                  m_pup.ideal = -log(ice.constants$Mu.S_pup.ideal))
  }else{
    Inits <- list(ice = ice.inits(ice.data$ice),
                  Mu.ice = runif(3, mean(ice.data$ice, na.rm = T)*0.9, mean(ice.data$ice, na.rm = T)*1.1),
                  sigmaY.ice = runif(1, 0, 2),
                  S_pup.ideal = ice.constants$Mu.S_pup.ideal,
                  m_pup.ideal = -log(ice.constants$Mu.S_pup.ideal))
  }
  
  
  return(Inits)
}

Inits <- list(initVal.sim(ice.data, TrendIceModel), initVal.sim(ice.data, TrendIceModel), initVal.sim(ice.data, TrendIceModel))

# Parameters
if(TrendIceModel){
  params <- c('ice', 'Mu.ice', 'beta.ice', 'sigmaY.ice', 'S_pup', 'S_pup.ideal')
}else{
  params <- c('ice', 'Mu.ice', 'sigmaY.ice', 'S_pup', 'S_pup.ideal')
}

## Test run
testRun <- nimbleMCMC(code = ice.model, 
                      constants = ice.constants, 
                      data = ice.data, 
                      inits = Inits, 
                      monitors = params,
                      niter = 5000, nburnin = 1000, nchains = 3, thin = 1,
                      samplesAsCodaMCMC = TRUE, setSeed = mySeed)
plot(testRun)

out.mat <- as.matrix(testRun)

####################
# MAKE PREDICTIONS #
####################

## Predicting sea ice change over time
ice.pred <- data.frame(Year = 1981:2021, Ice_lCI = NA, Ice_Median = NA, Ice_uCI = NA)

for(t in 1:nrow(ice.pred)){
  
  if(TrendIceModel){
    sam.pred <- exp(log(out.mat[,'Mu.ice']) + out.mat[,'beta.ice']*t)
  }else{
    sam.pred <- out.mat[,paste0('Mu.ice[', ice.constants$ice.period[t], ']')]
  }
  
  ice.pred[t, 2:4] <- quantile(sam.pred, probs = c(0.025, 0.5, 0.975))
}


## Predicting pup survival as a function of sea ice
test.ice <- seq(0, max(ice.data$ice, na.rm = T), length.out = 100)
S_pup.pred <- data.frame(Ice = test.ice, S_pup_lCI = NA, S_pup_Median = NA, S_pup_uCI = NA)

for(i in 1:length(test.ice)){
  
  sam.pred <- calc.S_pup(S_pup.ideal = out.mat[,'S_pup.ideal'], ice.ideal = ice.constants$ice.ideal, ice = S_pup.pred$Ice[i])
  S_pup.pred[i, 2:4] <- quantile(sam.pred, probs = c(0.025, 0.5, 0.975))
}

## Predicting sea ice change over time
S_pup_t.pred <- data.frame(Year = 1981:2021, S_pup_lCI = NA, S_pup_Median = NA, S_pup_uCI = NA)

for(t in 1:nrow(S_pup_t.pred)){
  sam.pred <- out.mat[,paste0('S_pup[', t, ']')]
  S_pup_t.pred[t, 2:4] <- quantile(sam.pred, probs = c(0.025, 0.5, 0.975))
}

## Rename predictions
if(TrendIceModel){
  ice.pred$Model <- 'Trend'
  S_pup.pred$Model <- 'Trend'
  S_pup_t.pred$Model <- 'Trend'
  trend.ice.pred <- ice.pred
  trend.S_pup.pred <- S_pup.pred
  trend.S_pup_t.pred <- S_pup_t.pred
}else{
  ice.pred$Model <- 'Period'
  S_pup.pred$Model <- 'Period'
  S_pup_t.pred$Model <- 'Period'
  period.ice.pred <- ice.pred
  period.S_pup.pred <- S_pup.pred
  period.S_pup_t.pred <- S_pup_t.pred
}

##########################################
# VISUALIZE PREDICTIONS FROM BOTH MODELS #
##########################################
library(ggplot2)

## Combine data
ice.pred <- rbind(trend.ice.pred, period.ice.pred)
S_pup.pred <- rbind(trend.S_pup.pred, period.S_pup.pred)
S_pup_t.pred <- rbind(trend.S_pup_t.pred, period.S_pup_t.pred)

## Plot: Ice predictions
obs.data <- data.frame(Year = year, Ice = ice)
pdf('220428_IceModel_IcePred.pdf', width = 8, height = 4)
ggplot(ice.pred) + geom_line(aes(x = Year, y = Ice_Median, color = Model)) +
  geom_ribbon(aes(x = Year, ymin = Ice_lCI, ymax = Ice_uCI, fill = Model), alpha = 0.5) +
  geom_point(data = obs.data, aes(x = Year, y = Ice), color = 'cornflowerblue') + 
  geom_line(data = obs.data, aes(x = Year, y = Ice), color = 'cornflowerblue', linetype = 'dashed') + 
  ylab('Ice extent (average)') +
  scale_color_manual(values = c('hotpink', 'orange')) + 
  scale_fill_manual(values = c('hotpink', 'orange')) + 
  theme_bw()
dev.off()

## Plot: Pup survival ~ ice predictions
pdf('220428_IceModel_S_pupPred.pdf', width = 6, height = 4)
ggplot(S_pup.pred) + geom_line(aes(x = Ice, y = S_pup_Median, color = Model)) +
  geom_ribbon(aes(x = Ice, ymin = S_pup_lCI, ymax = S_pup_uCI, fill = Model), alpha = 0.5) +
  ylab('Pup survival') +
  scale_color_manual(values = c('hotpink', 'orange')) + 
  scale_fill_manual(values = c('hotpink', 'orange')) + 
  theme_bw()
dev.off()

## Plot: Pup survival ~ year predictions
pdf('220428_IceModel_S_pupPred_t.pdf', width = 8, height = 4)
ggplot(S_pup_t.pred) + geom_line(aes(x = Year, y = S_pup_Median, color = Model)) +
  geom_ribbon(aes(x = Year, ymin = S_pup_lCI, ymax = S_pup_uCI, fill = Model), alpha = 0.5) +
  ylab('Pup survival') +
  scale_color_manual(values = c('hotpink', 'orange')) + 
  scale_fill_manual(values = c('hotpink', 'orange')) + 
  theme_bw()
dev.off()