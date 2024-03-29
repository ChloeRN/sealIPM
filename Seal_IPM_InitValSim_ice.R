#-----------------------------------------------------#
# FUNCTION FOR SIMULATING INITIAL VALUES FOR SEAL IPM #
#-----------------------------------------------------#

initValSim <- function(data, constants){
  
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
  
  ## Age-specific maturation rates
  Mu.pMat <- c(0, 0, runif(3, 0.3, 0.7), 1)
  
  ## Ovulation rate
  pOvl_exp <- sum(data$noOvl)/sum(data$noMat)
  pOvl <- runif(1, pOvl_exp*0.75, ifelse(pOvl_exp*1.25 < 1, pOvl_exp*1.25, 1))
  
  ## Pregnancy rate
  pPrg_exp <- data$noPrg.P3/data$noMat.P3
  pPrg <- runif(1, pPrg_exp*0.75, ifelse(pPrg_exp*1.25 < 1, pPrg_exp*1.25, 1))
  
  ## Pup survival
  S_pup.ideal <- constants$Mu.S_pup.ideal
  
  ## First-year survival and natural & harvest mortality 
  # Estimation period
  S_YOY <- c(data$S_YOY.fix, NA)
  mN_YOY <- runif(1, exp(data$mN.logmean-2*data$mN.logsd), -log(S_YOY))
  mH_YOY <- -log(S_YOY) - mN_YOY
  
  # Simulation period
  mH_YOY[2] <- mH_YOY[1]*data$mH.fac
  S_YOY[2] <- exp(-(mH_YOY[2] + mN_YOY))
  
  ## Subadult survival & natural mortality
  # Estimation period
  S_SA <- matrix(c(data$S_SA.fix, rep(NA, constants$SA_Amax)), ncol = 2)
  mN_SA <- rep(NA, nrow(S_SA))
  for(a in 1:nrow(S_SA)){
    mN_SA[a] <- runif(1, exp(data$mN.logmean-2*data$mN.logsd), -log(S_SA[a,1]))
  }
  mH_SA <- -log(S_SA) - mN_SA
  
  # Simulation period
  mH_SA[,2] <- mH_SA[,1]*data$mH.fac
  S_SA[,2] <- exp(-(mH_SA[,2] + mN_SA))
  
  
  ## Adult survival & natural mortality
  # Estimation period
  S_MA <- c(data$S_MA.fix, NA)
  mN_MA <- runif(1, exp(data$mN.logmean-2*data$mN.logsd), -log(S_MA))
  mH_MA <- -log(S_MA) - mN_MA
  
  # Simulation period
  mH_MA[2] <- mH_MA[1]*data$mH.fac
  S_MA[2] <- exp(-(mH_MA[2] + mN_MA))
  
  ## Proportions dying due to harvest
  alpha <- matrix(NA, nrow = constants$Amax-1, ncol = 2)
  for(i in 1:2){
    alpha[1,i] <- mH_YOY[i]/(mH_YOY[i]+mN_YOY) 
    alpha[2:(constants$max.matA+1),i] <- mH_SA[1:constants$max.matA,i]/(mH_SA[1:constants$max.matA,i]+mN_SA[1:constants$max.matA]) 
    alpha[constants$Amax-1,i] <- mH_MA[i]/(mH_MA[i]+mN_MA)
  }
  
  # 1.2) Vital rate standard deviations and random effects
  
  ## Maturation rates
  sigmaY.pMat <- runif(1, 0.01, 0.1)
  epsilonY.pMat <- rep(0, constants$sim_Tmax)
  
  
  # 1.3) Time-dependent vital rates
  
  ## Maturation rates 
  pMat <- matrix(0, nrow = constants$max.matA+1, ncol = constants$sim_Tmax)
  
  for(t in 1:constants$sim_Tmax){
    for(a in constants$min.matA:constants$max.matA){
      pMat[a, t] <- icloglog(cloglog(Mu.pMat[a]) + epsilonY.pMat[t])
    }
  }
  
  pMat[constants$max.matA+1, 1:constants$sim_Tmax] <- 1
  
  ## Pup survival
  S_pup <- rep(NA, constants$sim_Tmax)
  env.ideal <- constants$Mu.env.ideal
  
  for(t in 1:constants$sim_Tmax){
    S_pup[t] <- calc.S_pup(S_pup.ideal = S_pup.ideal, env.ideal = env.ideal, env = env.cov[t])
  }
  
  #-----------------------------------------------------------------------
  
  # 2) Initialize population model #
  #--------------------------------#
  
  # 2.1) Draw population size from aerial survey estimation
  estN.2002 <- rnorm(1, mean = data$estN.mean, sd = data$estN.sd)
  
  # 2.2) Extract stable ageclass distribution from projection matrix
  projMat <- make.sealPM(S_YOY = S_YOY[1], S_SA = S_SA[1:constants$SA_Amax,1], S_MA = S_MA[1], 
                         pMat = Mu.pMat[1:constants$SA_Amax], 
                         pOvl = pOvl, pPrg = pPrg, S_pup = S_pup[constants$sim_Tmin+1])
  SAD <- as.numeric(eigen(projMat)$vectors[,1]/sum(eigen(projMat)$vectors[,1]))
  
  # 2.3) Set inital numbers of individuals in each class
  YOY <- rep(NA, constants$sim_Tmax)
  YOY[constants$sim_Tmin] <- round(estN.2002*SAD[1]*0.5)
  
  SubA <- matrix(NA, nrow = constants$SA_Amax, ncol = constants$sim_Tmax)
  SubA[,constants$sim_Tmin] <- round(estN.2002*SAD[2:(constants$SA_Amax+1)]*0.5)
  
  nMatA <- rep(NA, constants$sim_Tmax)
  nMatA[constants$sim_Tmin] <- round(estN.2002*SAD[constants$Amax-1]*0.5)
  
  MatA <- rep(NA, constants$sim_Tmax)
  MatA[constants$sim_Tmin] <- round(estN.2002*SAD[constants$Amax]*0.5)
  
  # 2.4) Prepare additional vectors/matrixes for storing results
  Surv_SubA <- Mat_SubA <- matrix(NA, nrow = constants$SA_Amax+1, ncol = constants$sim_Tmax)
  Surv_SubA[1,] <- Mat_SubA[1,] <- 0
  Surv_nMatA <- Surv_MatA <- R_Mat <- R_nMat <- Ntot <- rep(NA, constants$sim_Tmax)
  lambda_real <- rep(NA, constants$sim_Tmax-1)
  D <- H <- matrix(NA, nrow = constants$Amax-1, ncol = constants$sim_Tmax-1)
  
  # 2.5) Set values prior to simulation start to 0
  YOY[1:(constants$sim_Tmin-1)] <- 0
  SubA[1:constants$SA_Amax, 1:(constants$sim_Tmin-1)] <- 0
  nMatA[1:(constants$sim_Tmin-1)] <- 0
  MatA[1:(constants$sim_Tmin-1)] <- 0
  
  Surv_SubA[, 1:constants$sim_Tmin] <- 0
  Mat_SubA[, 1:constants$sim_Tmin] <- 0
  Surv_nMatA[1:constants$sim_Tmin] <- 0
  Surv_MatA[1:constants$sim_Tmin] <- 0
  R_Mat[1:constants$sim_Tmin] <- 0
  R_nMat[1:constants$sim_Tmin] <- 0
  Ntot[1:(constants$sim_Tmin-1)] <- 0
  lambda_real[1:(constants$sim_Tmin-1)] <- 0
  
  D[, 1:(constants$sim_Tmin-1)] <- 0
  H[, 1:(constants$sim_Tmin-1)] <- 0
  
  #-----------------------------------------------------------------------
  
  # 3) Simulate population dynamics over time #
  #-------------------------------------------#
  
  for(t in constants$sim_Tmin:(constants$sim_Tmax-1)){
    
    # 3.1) Survival to the next year 
    
    ## Young-of-the-year --> age 1 subadults
    SubA[1, t+1] <- rbinom(1, YOY[t], S_YOY)
    
    ## Age 1-5 subadults
    for(a in 1:constants$SA_Amax){
      Surv_SubA[a+1, t+1] <- rbinom(1, SubA[a, t], S_SA[a, constants$p.idx[t]])
    }
    
    ## (Newly) mature adults --> mature adults
    Surv_nMatA[t+1] <- rbinom(1, nMatA[t], S_MA[constants$p.idx[t]])
    Surv_MatA[t+1] <- rbinom(1, MatA[t], S_MA[constants$p.idx[t]])
    
    MatA[t+1] <- Surv_nMatA[t+1] + Surv_MatA[t+1]
    
    
    # 3.2) Deaths due to harvest #
    
    ## Number that died in each age class
    #  Age classes re-defined to match AaH data: 
    #  1 = YOY, 2-6 = Sub[1]-Sub[5], 7 = nMatA+MatA
    D[1,t] <- YOY[t] - SubA[1, t+1]
    
    for(a in 1:constants$SA_Amax){
      D[1+a,t] <- SubA[a, t] - Surv_SubA[a+1, t+1]
    }
    
    D[constants$Amax-1,t] <- nMatA[t] + MatA[t] -  MatA[t+1]
    
    ## Number that died due to harvest
    for(a in 1:(constants$Amax-1)){
      
      if(t > constants$Tmax & HarvestScen == 'None'){
        H[a,t] <- 0
      }else{
        H[a,t] <- extraDistr::rtbinom(1, D[a,t], alpha[a, constants$p.idx[t]], a = 0, b = Inf)
        # NOTE: Using a truncated binomial here to avoid sampling H = 0 when D > 0
      }
    }
    
    
    # 3.3) Maturation to the next year #
    
    ## Surviving age 1-5 (2-6) subadults maturing
    Mat_SubA[1, t+1] <- 0
    
    for(a in 1:constants$SA_Amax){
      Mat_SubA[a+1, t+1] <- rbinom(1, Surv_SubA[a+1, t+1], pMat[a+1, t+1])
    }
    
    nMatA[t+1] <- sum(Mat_SubA[1:(constants$SA_Amax+1), t+1])

    ## Surviving age 1-5 (2-6) subadults remaining immature
    SubA[2:constants$SA_Amax, t+1] <- Surv_SubA[2:constants$SA_Amax, t+1] - Mat_SubA[2:constants$SA_Amax, t+1]
    
    
    # 3.4) Reproduction prior to next year's census #
    
    ## Expected reproduction from surviving newly mature adults
    #R_nMat[t+1] <- Surv_nMatA[t+1] * pPrg * 0.5 * S_pup[t+1]
    R_nMat[t+1] <- Surv_nMatA[t+1] * pPrg * 0.5 * ifelse(S_pup[t+1] > 0.5, S_pup[t+1], 0.6)
    
    ## Expected reproduction from surviving mature adults
    #R_Mat[t+1] <- Surv_MatA[t+1] * pOvl * pPrg * 0.5 * S_pup[t+1]
    R_Mat[t+1] <- Surv_MatA[t+1] * pOvl * pPrg * 0.5 * ifelse(S_pup[t+1] > 0.5, S_pup[t+1], 0.6)
    
    # NOTE: I'm constraining S_pup used for simulations to be larger than 0.5 
    #       because inconsistencies between data and simulation are much more 
    #       likely to happen when simulating a declining population. 
    
    ## Realized reproduction
    YOY[t+1] <- rpois(1, R_nMat[t+1] + R_Mat[t+1])
    
  }
  
  
  # 3.5) Calculation of population size and growth #
  
  for(t in constants$sim_Tmin:constants$sim_Tmax){
    
    ## Total population size at census
    Ntot[t] <- YOY[t] + sum(SubA[1:constants$SA_Amax, t]) + nMatA[t] + MatA[t]
    
    ## Realized population growth rate
    if(t > 1){
      lambda_real[t-1] <- Ntot[t] / Ntot[t-1]
    }
  }
  
  #-----------------------------------------------------------------------
  
  # 4) Calculate additional data sampling parameters #
  #--------------------------------------------------#
  
  # 4.1) Scaling factor for age-class at harvest data
  pACaH <- rep(0, constants$sim_Tmin-1)
  
  for(t in constants$sim_Tmin:(constants$sim_Tmax-1)){
    pACaH[t] <- data$no.ACaH[t]/sum(H[,t])
  }
  
  #-----------------------------------------------------------------------
  
  # 7) Assemble and return initial values #
  #---------------------------------------#
  
  InitVals <- list(
    Mu.pMat = Mu.pMat, pOvl = pOvl, pPrg = pPrg, 
    S_pup.ideal = S_pup.ideal, m_pup.ideal = -log(S_pup.ideal),
    S_YOY = S_YOY, mN_YOY = mN_YOY, mH_YOY = mH_YOY, m_YOY = -log(S_YOY[1]),
    S_SA = S_SA, mN_SA = mN_SA, mH_SA = mH_SA, m_SA = -log(S_SA[,1]),
    S_MA = S_MA, mN_MA = mN_MA, mH_MA = mH_MA, m_MA = -log(S_MA[1]),
    alpha = alpha,
    
    sigmaY.pMat = sigmaY.pMat, epsilonY.pMat = epsilonY.pMat,
    
    pMat = pMat, S_pup = S_pup, env.ideal = unname(env.ideal),
    
    estN.2002 = estN.2002, SAD = SAD, lambda_asym = as.numeric(eigen(projMat)$values[1]),
    
    YOY = YOY, SubA = SubA, nMatA = nMatA, MatA = MatA,
    Surv_SubA = Surv_SubA, Mat_SubA = Mat_SubA, Surv_nMatA = Surv_nMatA, Surv_MatA = Surv_MatA,
    R_nMat = R_nMat, R_Mat = R_Mat, 
    
    Ntot = Ntot, lambda_real = lambda_real,
    
    D = D, H = H, 
    
    pACaH = pACaH,
    
    env = env,
    Mu.env = Mu.env, beta.env = beta.env, sigmaY.env = sigmaY.env
  )
  
  return(InitVals)
}
