library(nimble)
library(gridExtra)

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
load(paste0(DataPath, '220114_SealIPM_Data.RData'))
IceData <- readRDS(paste0(DataPath, '220509_SealIPM_IceData.rds'))

# Pup survival
Mu.S_pup.ideal <- 0.80 # pup survival under ideal conditions
sdlog.m_pup <- 0.20 # standard deviation of uncertainty in pup mortality hazard rate (on the log scale)

# Ice/lair habitat covariate switch
# TRUE = pup survival is a function of sea ice availability
# FALSE = pup survival is a function of lair habitat availability
IceCov <- TRUE
# IceCov <- FALSE

## Set env. covariate incl. ideal conditions
if(IceCov){
  env.cov <- c(IceData$ice.GL.mean)
}else{
  env.cov <- c(IceData$lairH.GL.mean)
}

Mu.env.ideal <- quantile(env.cov[which(IceData$ice.years < 2006)], probs = 0.75, na.rm = TRUE) # threshold env. conditions which are considered ideal
sdlink.env_ideal <- 0.1 # standard deviation of uncertainty in ideal env. conditions on the link scale (= log)


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

##########################
# S_PUP ~ ENV (a priori) #
##########################

## Simulate Pup survival under ideal conditions
m_pup.ideal <- rlnorm(100000, meanlog = log(-log(Mu.S_pup.ideal)), sdlog = sdlog.m_pup)
S_pup.ideal <- exp(-m_pup.ideal)

## Simulate ideal conditions
env.ideal.sim <- rlnorm(100000, meanlog = log(Mu.env.ideal), sdlog = sdlink.env_ideal) 

## Set range over which to calculate S_pup
env.sim <- seq(0, 100, length.out = 200)

## Simulate S_pup over range
Sim.S_pup <- data.frame(SeaIce = env.sim, Median = NA, lCI = NA, uCI = NA)

for(x in 1:length(env.sim)){
  
  S_pup.sample <- rep(NA, length(S_pup.ideal))
  
  for(i in 1:length(S_pup.sample)){
    S_pup.sample[i] <- calc.S_pup(S_pup.ideal = S_pup.ideal[i], 
                                  env.ideal = env.ideal.sim[i], 
                                  env = env.sim[x])
  }
  Sim.S_pup[x, 2:4] <- quantile(S_pup.sample, probs = c(0.5, 0.025, 0.975))
}


########
# PLOT #
########

## Plot sea ice time-series
ice.data <- data.frame(SeaIce = IceData$ice.GL.mean, Year = IceData$ice.years)
ice.data <- subset(ice.data, !is.na(SeaIce))

p.ice <- ggplot(ice.data, aes(x = Year, y = SeaIce)) + 
  geom_line() + 
  geom_hline(aes(yintercept = Mu.env.ideal), linetype = 'dashed', color = 'slateblue2') + 
  geom_hline(aes(yintercept = exp(log(Mu.env.ideal) + sdlink.env_ideal)), linetype = 'dotted', color = 'turquoise3') + 
  geom_hline(aes(yintercept = exp(log(Mu.env.ideal) - sdlink.env_ideal)), linetype = 'dotted', color = 'turquoise3') + 
  scale_x_continuous(breaks = seq(1990, 2020, by = 5)) + 
  ylab('Sea ice extent') + 
  ggtitle('A) Sea ice extent (Mar-May average) over time') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())
p.ice

## Plot S_pup simulation
p.S_pup <- ggplot(Sim.S_pup) + 
  geom_line(aes(x = SeaIce, y = Median)) + 
  geom_ribbon(aes(x = SeaIce, ymin = lCI, ymax = uCI), alpha = 0.5) + 
  geom_vline(aes(xintercept = Mu.env.ideal), linetype = 'dashed', color = 'slateblue2') + 
  geom_vline(aes(xintercept = exp(log(Mu.env.ideal) + sdlink.env_ideal)), linetype = 'dotted', color = 'turquoise3') + 
  geom_vline(aes(xintercept = exp(log(Mu.env.ideal) - sdlink.env_ideal)), linetype = 'dotted', color = 'turquoise3') + 
  scale_x_continuous(breaks = seq(0, 100, by = 20)) + 
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.2)) + 
  xlab('Sea ice extent (Mar-May average)') + 
  ylab('Pup survival') + 
  ggtitle('B) Simulated relationship of pup survival and sea ice extent') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())
p.S_pup

## Print to pdf
pdf('S_pup&SeaIce.pdf', width = 6, height = 6)
grid.arrange(p.ice, p.S_pup, ncol = 1, heights = c(0.6, 1))
dev.off()
