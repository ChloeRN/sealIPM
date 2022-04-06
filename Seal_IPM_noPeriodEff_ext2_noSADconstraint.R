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
load(paste0(DataPath, '220114_SealIPM_Data.RData'))
mN.param.data <- readRDS(paste0(CodePath, '220119_mO_HoeningParameters_Seal.rds'))
HarvestData <- readRDS(paste0(CodePath, '220405_HarvestCountData.rds'))
  
## Assemble data and constants
seal.data <- list(
  Mature = SealDataF_mat$Maturity, # Individual maturity status
  
  noOvl = c(112, 77, 63), # Number of ovluating females (per time period)
  noMat = c(123, 90, 67), # Number of mature females (per time period)
  noPrg.P3 = 37, # Number of pregnant females (Aug-Oct of period 3)
  noMat.P3 = 52, # Number of mature females (Aug-Oct of period 3)
  
  S_YOY.fix = 0.75,
  #S_SA.fix = c(0.82, 0.84, 0.86, 0.88, 0.90),
  S_SA.fix = c(0.80, 0.82, 0.84, 0.86, 0.88),
  S_MA.fix = 0.92,
  # NOTE: The survival probabilities are fixed to similar numbers as used by
  # Reimer et al. 2019, Ecological Applications 29(3), e01855
  
  estN.mean = 3242, # Mean of estimated number of seals (males and females) from 2002 aerial survey
  estN.sd = 266, # Standard deviation of estimated number of seals (males and females) from 2002 aerial survey (NOTE: reported as se, but seems to actually be sd)

  mN.logmean = mN.param.data$log_mean[1], mN.logsd = mN.param.data$log_sd[1],
  
  no.H = HarvestData$no.H,
  r = HarvestData$r,
  
  ACaH_k = HarvestData$ACaH_k,
  uMatA = HarvestData$uMatA,
  no.ACaH = HarvestData$No.ACaH.yr,
  
  count.Isfj = HarvestData$count.Isfj,
  count.all = HarvestData$count.all
)

seal.constants <- list(
  Amax = 8, # Total number of age classes
  SA_Amax = 5, # Number of subadult age classes
  
  Tmax = length(1981:2020), # Total length of the study period
  sim_Tmin = 2002 - 1981 + 1, # Start year for population model
  sim_Tmax = length(1981:2020), # End year for population model
  
  Pmax = 3, 
  
  X = length(seal.data$Mature),
  min.matA = min(unique(SealDataF_mat$AgeClass)),
  max.matA = max(unique(SealDataF_mat$AgeClass)),
  Age =  SealDataF_mat$AgeClass,
  #SamYear = as.numeric(as.factor(SealDataF_mat$Year)),
  SamYear = SealDataF_mat$Year - min(SealDataF_mat$Year) + 1,

  estN.yr.idx = 2002 - min(SealDataF_mat$Year) + 1, # Time index for the year of aerial survey (2002)

  S_pup.lLimit = 0.4, # Lower limit for baseline pup survival
  S_pup.uLimit = 0.6, # Upper limit for baseline pup survival
  
  SAD.alpha = rep(1, 8)
)


####################
# NIMBLE FUNCTIONS #
####################

## Nimble function for assembling seal projection matrix
make.sealPM <- nimbleFunction(
  
  run = function(S_YOY = double(0),
                 S_SA = double(1),
                 S_MA = double(0),
                 pMat = double(1),
                 pOvl = double(0),
                 pPrg = double(0),
                 S_pup = double(0)) {
    
    
    # Specify matrix of correct dimensions
    A <- matrix(0, nrow = 8, ncol = 8)
    
    # Add reproduction components
    A[1,7] <- pPrg*S_MA*0.5*S_pup # Offspring production by newly matured females
    A[1,8] <- pOvl*pPrg*S_MA*0.5*S_pup # Offspring production by previously matured females
    
    # Add survival-maturation components
    A[2,1] <- S_YOY # Young-of-the year surviving to next year
    A[3,2] <- S_SA[1] # Age 1 immature subadults surviving to next year
    
    for(a in 3:5){ # Age 2-4 immature subadults...
      A[a+1,a] <- S_SA[a-1]*(1-pMat[a]) # ... surviving to next year and remaining immature
      A[7,a] <- S_SA[a-1]*pMat[a] # ... surviving to next year and maturing
    }
    
    A[7,6] <- S_SA[5] # Age 5 immature subadults surviving to next year and maturing
    A[8,7:8] <- S_MA # (Newly) mature individuals surviving to next year
    
    # Return matrix
    return(A)
    returnType(double(2))
  }
)


##############
# MODEL CODE #
##############

seal.IPM <- nimbleCode({

  # (ST)AGE CLASSES:
  # YOY[t] = young-of-the-year in year t(age 0)
  # SubA[a,t] = age a immature subadult in year t (ages 1-8)
  # nMatA[t] = newly mature adult / first-time reproducer in year t (ages 3-9+)
  # MatA[t] = mature adult in year t(ages 4-9+)
  # Surv_Sub[a,t+1] = Immature individuals that survived from age a-1 to a, and
  #                   over the interval t to t+1.
  # Mat_SubA[a,t+1] = Immature survivors that matured just prior to t+1 census
  
  # VITAL RATE PARAMETERS: 
  # S_YOY[t] = first-year survival (interval t to t+1)
  # S_SA[a,t] = age a subadult survival (interval t to t+1)
  # S_MA[t] = adult survival (interval to t+1)
  # pMat[a,t] = maturation rate of age class a individual in year t
  # pOvl[t] = ovulation rate of mature females in year t
  # pPrg[t] = pregnancy rate of mature females in year t
  # S_pup[t] = survival of pup from conception to weaning in year t
  
  # OTHER PARAMETERS:
  # R_nMat[t] = Expected reproductive output of newly matured adults in year t
  # R_Mat[t] = Expected reproductive output of previously mature adults in year t
  # Ntot[t] = total population at census (early June) in year t
  # lambda[t] = realized population growth rate from year t to t+1 (census)
  
  #******************#
  # POPULATION MODEL #                                                          
  #******************#
  
  # Initialization of population model #
  #------------------------------------#
  
  # This will have to divide Ntot[estN.yr.idx], which has a likelihood below,
  # into the different age classes. 
  # There is no information on age structure in the population survey. 
  # This means we have to use either age structure from the harvest sample 
  # (known to underestimate proportion in youngest age classes) or stable age 
  # distribution, which can be extracted from the matrix (standardized eigenvector)
  
  ## Approximation of initial numbers of females per age class (deterministic)
  #YOY[estN.yr.idx] <- round(SAD[1] * estN.2002 * 0.5)
  #SubA[1:SA_Amax, estN.yr.idx] <- round(SAD[2:(SA_Amax+1)] * estN.2002 * 0.5)
  #nMatA[estN.yr.idx] <- round(SAD[(Amax-1)] * estN.2002 * 0.5)
  #MatA[estN.yr.idx] <- round(SAD[Amax] * estN.2002 * 0.5)
  
  ## Approximation of initial numbers of females per age class (stochastic)
  YOY[estN.yr.idx] ~ dpois(SAD[1]*estN.2002*0.5)
  for(a in 1:SA_Amax){
    SubA[a, estN.yr.idx] ~ dpois(SAD[a+1]*estN.2002*0.5)
  }
  nMatA[estN.yr.idx] ~ dpois(SAD[(Amax-1)]*estN.2002*0.5)
  MatA[estN.yr.idx] ~ dpois(SAD[Amax]*estN.2002*0.5)
  

  for(t in sim_Tmin:(sim_Tmax-1)){
    
    # Survival to the next year #
    #---------------------------#
    
    ## Young-of-the-year --> age 1 subadults
    SubA[1, t+1] ~ dbin(S_YOY, YOY[t])
    
    ## Age 1-5 subadults
    for(a in 1:SA_Amax){
      Surv_SubA[a+1, t+1] ~ dbin(S_SA[a], SubA[a, t])
    }
    
    ## (Newly) mature adults --> mature adults
    Surv_nMatA[t+1] ~ dbin(S_MA, nMatA[t])
    Surv_MatA[t+1] ~ dbin(S_MA, MatA[t])
    
    MatA[t+1] <- Surv_nMatA[t+1] + Surv_MatA[t+1]

    
    # Deaths due to harvest #
    #-----------------------#
    
    ## Number that died in each age class
    #  Age classes re-defined to match AaH data: 
    #  1 = YOY, 2-6 = Sub[1]-Sub[5], 7 = nMatA+MatA
    D[1,t] <- YOY[t] - SubA[1, t+1]
    
    for(a in 1:SA_Amax){
      D[1+a,t] <- SubA[a, t] - Surv_SubA[a+1, t+1]
    }
    
    D[Amax-1,t] <- nMatA[t] - Surv_nMatA[t+1]
    D[Amax,t] <- MatA[t] -  Surv_MatA[t+1]
    
    ## Number that died due to harvest
    for(a in 1:Amax){
      H[a,t] ~ dbin(alpha[a], D[a,t])
    }

    
    # Maturation to the next year #
    #-----------------------------#
  
    # NOTES: 
    # - With the census placed between ovluation and mating,
    #   maturation is conditional on survival (first survive, then mature).
    # - Maturation rates need to be set to 0 for ages 1 and 2, and to 1 for age 9
  
    
    ## Surviving age 1-5 (2-6) subadults maturing
    Mat_SubA[1, t+1] <- 0
    
    for(a in 1:SA_Amax){
      Mat_SubA[a+1, t+1] ~ dbin(pMat[a+1, t+1], Surv_SubA[a+1, t+1])
    }
    
    nMatA[t+1] <- sum(Mat_SubA[1:(SA_Amax+1), t+1])
    # NOTE: Double-check what happens with the subadults in the last age class that HAVE TO mature
    
    ## Surviving age 1-5 (2-6) subadults remaining immature
    SubA[2:SA_Amax, t+1] <- Surv_SubA[2:SA_Amax, t+1] - Mat_SubA[2:SA_Amax, t+1]
  
  
    
    # Reproduction prior to next year's census #
    #------------------------------------------#
  
    ## Expected reproduction from surviving newly mature adults
    R_nMat[t+1] <- Surv_nMatA[t+1] * pPrg * 0.5 * S_pup[t+1]

    ## Expected reproduction from surviving mature adults
    R_Mat[t+1] <- Surv_MatA[t+1] * pOvl * pPrg * 0.5 * S_pup[t+1]
    
    ## Realized reproduction
    YOY[t+1] ~ dpois(R_nMat[t+1] + R_Mat[t+1])
    
  }
  
  # Calculation of population size and growth #
  #-------------------------------------------#
  
  ## Total population size at census
  for(t in sim_Tmin:sim_Tmax){
    Ntot[t] <- YOY[t] + sum(SubA[1:SA_Amax, t]) + nMatA[t] + MatA[t]
  }
  
  ## Realized population growth rate
  lambda_real[sim_Tmin:(sim_Tmax-1)] <- Ntot[(sim_Tmin+1):sim_Tmax] / Ntot[sim_Tmin:(sim_Tmax-1)]

  
  # Book-keeping #
  #--------------#
  
  ## Fixing numbers prior to simulation start to 0
  YOY[1:(sim_Tmin-1)] <- 0
  SubA[1:SA_Amax, 1:(sim_Tmin-1)] <- 0
  nMatA[1:(sim_Tmin-1)] <- 0
  MatA[1:(sim_Tmin-1)] <- 0
  
  SurvSubA[1:(SA_Amax+1), 1:sim_Tmin] <- 0
  MatSubA[1:(SA_Amax+1), 1:sim_Tmin] <- 0
  Surv_nMatA[1:sim_Tmin] <- 0
  Surv_MatA[1:sim_Tmin] <- 0
  R_Mat[1:sim_Tmin] <- 0
  R_nMat[1:sim_Tmin] <- 0
  Ntot[1:(sim_Tmin-1)] <- 0
  lambda_real[1:(sim_Tmin-1)] <- 0
  
  D[1:Amax,1:(sim_Tmin-1)] <- 0
  H[1:Amax,1:(sim_Tmin-1)] <- 0
  
  #******************#
  # DATA LIKELIHOODS #                                                          
  #******************#
  
  # Likelihood for population count data #
  #--------------------------------------#
  
  ## Estimated population size during 2002 aerial survey
  estN.2002 ~ dnorm(mean = estN.mean, sd = estN.sd)
  
  
  # Likelihood for harvest count data #
  #-----------------------------------#
  
  for(t in (sim_Tmin+1):(sim_Tmax-1)){
    no.H[t] ~ dpois((2*sum(H[1:(Amax-1),t])*r[t])/prop.Isfj)
  }
    
  
  # Likelihood for harvest age structure data #
  #-------------------------------------------#
  
  ## Likelihood for age assignment of unknown age class ACaH data (uMatA)
  u_nMat_H[1:3,1:(sim_Tmin-1)] <- 0
  extra_ACaH[1:8,1:(sim_Tmin-1)] <- 0
  
  for(t in sim_Tmin:sim_Tmax){
    
    # Age 4-6 mature females that are newly matured
    u_nMat_H[1,t] ~ dbin((1-pMat[3,t-1])*pMat[4,t], uMatA[1,t]) # Age 4
    u_nMat_H[2,t] ~ dbin((1-pMat[3,t-2])*(1-pMat[4,t-1])*pMat[5,t], uMatA[2,t]) # Age 5
    u_nMat_H[3,t] ~ dbin((1-pMat[3,t-3])*(1-pMat[4,t-4])*(1-pMat[5,t-1]), uMatA[3,t]) # Age 6

    # Age class distribution of additional females to add to ACaH data
    extra_ACaH[1:6,t] <- 0
    extra_ACaH[7,t] <- sum(u_nMat_H[1:3,t])
    extra_ACaH[8,t] <- sum(uMatA[1:3,t]) - extra_ACaH[7,t]
  }
  
  ## Likelihood for known age class ACaH data
  for(t in sim_Tmin:(sim_Tmax-1)){
    for(a in 1:Amax){
      ACaH_k[a,t] ~ dpois((H[a,t]*pACaH[t])-extra_ACaH[a,t]) 
    }
  }
  
  # Likelihood for maturation data #
  #--------------------------------#
  
  for(x in 1:X){
    Mature[x] ~ dbern(pMat[Age[x], SamYear[x]])
  }
  
  
  # Likelihood for ovulation data #
  #-------------------------------#
  
  for(p in 1:Pmax){
    noOvl[p] ~ dbin(pOvl, noMat[p])
  }
  
  # Likelihood for pregnancy data #
  #-------------------------------#
  
  noPrg.P3 ~ dbin(pPrg, noMat.P3)
  
  
  #*************#
  # CONSTRAINTS #
  #*************#
  
  # Maturation rates #
  #------------------#

  ## Age classes 1 & 2 (constrained immature)
  pMat[1:(min.matA-1), 1:Tmax] <- 0
    
  ## Age classes 3-5 (estimated)
  for(t in 1:Tmax){
    for(a in min.matA:max.matA){
      cloglog(pMat[a, t]) <- cloglog(Mu.pMat[a]) + epsilonY.pMat[t]
    }
  }

  ## Age class 6+ (constrained mature)
  pMat[max.matA+1, 1:Tmax] <- 1
  
  ## Random year effects
  for(t in 1:Tmax){
    epsilonY.pMat[t] ~ dnorm(0, sd = sigmaY.pMat)
  }
  
  
  # Pup survival #
  #--------------#
  
  #S_pup[1:Tmax] <- Mu.S_pup
  for(t in 1:Tmax){
    
    S_pup[t] <- exp(-(m_pup[t]))
    log(m_pup[t]) <- log(-log(Mu.S_pup)) + epsilonY.m_pup[t]
    
    epsilonY.m_pup[t] ~ dnorm(0, sd = sigmaY.m_pup)
  }
  
  # Mortality (hazard) rates #
  #--------------------------#
  
  ## Proportion dying due to harvest
  alpha[1] <- mH_YOY/(mH_YOY+mN_YOY) 
  alpha[2:(max.matA+1)] <- mH_SA[1:max.matA]/(mH_SA[1:max.matA]+mN_SA[1:max.matA]) 
  alpha[(Amax-1):Amax] <- mH_MA/(mH_MA+mN_MA)
  
  ## Harvest mortality
  
  # First-year
  mH_YOY <- -log(S_YOY) - mN_YOY
  
  # Subadults
  #mH_SA[1:max.matA] <- -log(S_SA[1:max.matA]) - mN_SA[1:max.matA]
  mH_SA[1:max.matA] <- mH_MA
  
  # Adults
  mH_MA <- -log(S_MA) - mN_MA
  
  
  # Proportion harvested in Isfjorden area #
  #----------------------------------------#
  
  prop.Isfj <- sum(count.Isfj[1:3])/sum(count.all[1:3])
  
  # NOTE: This could alternatively be formulated as a stochastic relationship:
  #       e.g. count.Isfj[x] ~ dbin(prop.Isfj, count.all[x])
  
  
  # Proportion of the harvest represented in ACaH data #
  #----------------------------------------------------#
  
  pACaH[1:(sim_Tmin-1)] <- 0
  
  for(t in sim_Tmin:(sim_Tmax-1)){
    pACaH[t] <- no.ACaH[t]/sum(H[1:(Amax-1),t])
  }
  
  #*********#
  # PRIORS  #                                                          
  #*********#
  
  ## Population parameters

  # Inital age class distribution
  SAD[1:Amax] ~ ddirch(alpha = SAD.alpha[1:Amax])
  
  ## Vital rate averages
  
  # Age-specific maturation rates
  Mu.pMat[1:(min.matA-1)] <- 0
  
  for(a in min.matA:max.matA){
    Mu.pMat[a] ~ dunif(0, 1)
  }

  Mu.pMat[max.matA+1] <- 1
  
  
  # Ovulation rate
  pOvl ~ dunif(0, 1)
  
  
  # Pregnancy rate
  pPrg ~ dunif(0, 1)
  
  
  # Pup survival
  Mu.S_pup ~ dunif(S_pup.lLimit, S_pup.uLimit)
  
  
  # First-year survival & natural mortality
  S_YOY <- S_YOY.fix
  #S_YOY ~ dwhatever()
  
  mN_YOY ~ T(dlnorm(meanlog = mN.logmean, sdlog = mN.logsd), 0, -log(S_YOY))
  
  # Subadult survival & natural mortality
  for(a in 1:max.matA){
    S_SA[a] <- S_SA.fix[a]
    #S_SA[a] ~ dwhatever()
    
    mN_SA[a] ~ T(dlnorm(meanlog = mN.logmean, sdlog = mN.logsd), 0, -log(S_SA[a]))
  }
  
  
  # Adult survival & natural mortality
  S_MA <- S_MA.fix
  #S_MA ~ dwhatever()
  
  mN_MA ~ T(dlnorm(meanlog = mN.logmean, sdlog = mN.logsd), 0, -log(S_MA))
  
  
  ## Fixed effects
  
  ## Random year variation
  
  # Annual variation in maturation rate 
  sigmaY.pMat ~ dunif(0, 5)
  
  # Annual variation in pup mortality
  sigmaY.m_pup ~ dunif(0, 10)
  
})


###########################
# INITAL VALUE SIMULATION #
###########################

## Function for simulating initial values
source(paste0(CodePath, 'Seal_IPM_InitValSim_ext2.R'))

## Function for sampling and checking consistency of initial values
initValSet <- function(data, constants){
  
  Tmin <- constants$sim_Tmin
  Tmax <- constants$sim_Tmax-1
  
  good.inits <- FALSE
  iter <- 0
  
  while(!good.inits){
    
    # Increase counter
    iter <- iter + 1
    
    # Sample initial values
    Inits <- initValSim(data = data, constants = constants)
  
    # Check harvest data consistency
    # class7_check --> not necessary because all unknown class individuals are assigned to class 8 in simulation
    class8_check <- all(Inits$H[8,Tmin:Tmax]*(data$no.ACaH[Tmin:Tmax]/colSums(Inits$H[,Tmin:Tmax])) >= colSums(data$uMat[,Tmin:Tmax]))
    class8_check
    
    good.inits <- class8_check  
  }
  
  message(paste0('Iterations to get good inits: ', iter))
  return(Inits)
}


## Sample initial values
#Inits <- list(initValSim(data = seal.data, constants = seal.constants))
Inits <- list(initValSet(data = seal.data, constants = seal.constants),
              initValSet(data = seal.data, constants = seal.constants),
              initValSet(data = seal.data, constants = seal.constants))

## Change pup survival parameter bounds in constants
seal.constants$S_pup.lLimit <- 0
seal.constants$S_pup.uLimit <- 1

############
# TEST RUN #
############

## Set parameters to monitor
params <- c('SAD', 
            'Ntot', 'lambda_real', 
            'Mu.pMat', 'sigmaY.pMat',
            'pMat', 
            'pOvl', 'pPrg',
            'Mu.S_pup', 'S_pup', 'sigmaY.m_pup',
            'estN.2002', 
            'YOY', 'SubA', 'nMatA', 'MatA',
            'mN_YOY', 'mN_SA', 'mN_MA',
            'mH_YOY', 'mH_SA', 'mH_MA',
            'S_YOY', 'S_SA', 'S_MA',
            'alpha', 'D', 'H', 'extra_ACaH', 'u_nMat_H' 
            )

## MCMC specs
niter <- 10
nburnin <- 0
nthin <- 1
nchains <- 3
#nchains <- 1

#niter <- 100000
#nburnin <- 30000
#nthin <- 10
#nchains <- 3

## Testrun
testRun <- nimbleMCMC(code = seal.IPM, 
                      constants = seal.constants, 
                      data = seal.data, 
                      inits = Inits, 
                      monitors = params, 
                      niter = niter, nburnin = nburnin, nchains = nchains, thin = nthin,
                      samplesAsCodaMCMC = TRUE, setSeed = mySeed)


setwd('/data/P-Prosjekter/41201625_sustainable_harvesting_of_seals_in_svalbard/SealIPM')
save(testRun, file = '220406_Seal_IPM_noPeriodEff_ext2_noSADconstraint_pupTvar.RData')

pdf('220406_IPMtest_noPeriodEff_ext2_noSADconstraint_pupTvar_Traces.pdf', height = 8, width = 11)
plot(testRun)
dev.off()
