library(coda)
library(ggplot2)
library(viridis)
library(gridExtra)
library(magrittr)
library(nimble)

####################################
# SETTING UP POSTERIOR PREDICTIONS #
####################################

## Load posterior samples from model run
testRun <- readRDS('IPMtest_fSAD&eHAD_iceSim.rds')
out.mat <- as.matrix(testRun)

## Set time-span to consider
Tmin <- 22
Tmax <- 40

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

## Assemble posteriors for Nmin and various Rmax
Ntot_min <- rep(NA, dim(out.mat)[1])
rlambda_max <- rlambda_mean <- Ntot_min
alambda_max <- alambda_mean <- Ntot_min

for(i in 1:length(Ntot_min)){
  
  # Minimum population size (multiplied by 2 to include both sexes)
  Ntot_min[i] <- min(out.mat[i, paste0('Ntot[', Tmin:Tmax, ']')])*2
  
  # Maximum and mean realized population growth rates
  rlambda_max[i]  <- max(out.mat[i, paste0('lambda_real[', Tmin:Tmax, ']')])
  rlambda_mean[i]  <- mean(out.mat[i, paste0('lambda_real[', Tmin:Tmax, ']')])

  # Maximum and mean asymptotic growth rate assuming no harvest
  alambda_t <- rep(NA, length(Tmin:Tmax))
  for(t in Tmin:Tmax){
    A <- make.sealPM(S_YOY = exp(-out.mat[i, 'mN_YOY']),
                     S_SA = exp(-out.mat[i, paste0('mN_SA[', 1:5, ']')]),
                     S_MA = exp(-out.mat[i, 'mN_MA']),
                     pMat = out.mat[i, paste0('pMat[', 1:6, ', ', t, ']')],
                     pOvl = out.mat[i, 'pOvl'],
                     pPrg = out.mat[i, 'pPrg'],
                     S_pup = out.mat[i, paste0('S_pup[', t+1, ']')])
    alambda_t[t-Tmin+1] <- as.numeric(eigen(A)$values[1])
  }
  alambda_max[i] <- max(alambda_t)
  alambda_mean[i] <- mean(alambda_t)
}


# Calculate minimum population size for 2002 as used in the rudimentary PBR
# NOTE: The "ideal" way to do this (according to refs in Nelson et al. 2019)
# is to calculate it as the "20th percentile of a log-normal distribution 
# based on an estimate of the number of animals in a stock
N.est <- data.frame(
  FjordID = c(5:11),
  EstMean = c(1142, 577, 0, 268, 1082, 137, 36),
  EstSD = c(187, 99, 0, 52, 149, 29, 7)
)

Isfj.Nsim <- function(n, N.est){
  
  # Vector to store samples
  simN <- rep(NA, n)
  
  # Sample from distributions for each fjord and sum
  for(i in 1:n){
    simN[i] <- sum(rnorm(7, mean = N.est$EstMean, sd = N.est$EstSD))
  }
  
  # Return results
  return(simN)
}

simN <- Isfj.Nsim(n = 1000000, N.est = N.est)
Nmin.2002 <- unname(exp(quantile(log(simN), probs = 0.2)))


## Make vector of different values of Fr
Fr <- seq(0.1, 1.0, by = 0.1)

## Make plotting directory if it does not exist
if(!file.exists("Plots")){
  dir.create("Plots")
}

#######################################
# PBR CALCULATION - Nmin from samples #
#######################################

## General function for PBR calculation
calc.PBR <- function(Nmin, lambda, Fr){
  
  # Make an empty data frame
  PBR.data <- data.frame()
  
  # Calculate Rmax from supplied lambda
  Rmax <- lambda - 1
  
  # Calculate PBR for each Fr and merge with remaining data
  for(i in 1:length(Fr)){
    PBR.data.new <- data.frame(
      Nmin = Nmin,
      Rmax = Rmax, 
      Fr = Fr[i],
      PBR = Nmin*0.5*Rmax*Fr[i]
    )
    PBR.data <- rbind(PBR.data, PBR.data.new)
  }
  return(PBR.data)
}

## Using Rmax based on max(r_lambda)
PBR1 <- calc.PBR(Nmin = Ntot_min, lambda = rlambda_max, Fr = Fr)
PBR1$Rmax_Source <- 'max(r_lambda)'

## Using Rmax based on mean(r_lambda)
PBR2 <- calc.PBR(Nmin = Ntot_min, lambda = rlambda_mean, Fr = Fr)
PBR2$Rmax_Source <- 'mean(r_lambda)'

## Using Rmax based on max(a_lambda)
PBR3 <- calc.PBR(Nmin = Ntot_min, lambda = alambda_max, Fr = Fr)
PBR3$Rmax_Source <- 'max(a_lambda)'

## Using Rmax based on max(a_lambda)
PBR4 <- calc.PBR(Nmin = Ntot_min, lambda = alambda_mean, Fr = Fr)
PBR4$Rmax_Source <- 'mean(a_lambda)'


## Combine data
PBR.dataA <- rbind(PBR1, PBR2, PBR3, PBR4)
PBR.dataA$Rmax_Label <- dplyr::case_when(PBR.dataA$Rmax_Source == 'max(r_lambda)' ~ 'Realized growth rate (max)',
                                         PBR.dataA$Rmax_Source == 'mean(r_lambda)' ~ 'Realized growth rate (mean)',
                                         PBR.dataA$Rmax_Source == 'max(a_lambda)' ~ 'Asymptotic growth rate (max)',
                                         PBR.dataA$Rmax_Source == 'mean(a_lambda)' ~ 'Asymptotic growth rate (mean)',)
PBR.dataA$Rmax_Label <- factor(PBR.dataA$Rmax_Label, levels = c('Realized growth rate (max)', 'Realized growth rate (mean)', 'Asymptotic growth rate (max)', 'Asymptotic growth rate (mean)'))

## Posterior summaries
PBR.summaries <- PBR.dataA %>%
  dplyr::group_by(Rmax_Label, Fr) %>%
  dplyr::summarise(median = median(PBR),
                   lCI = quantile(PBR, probs = 0.025),
                   uCI = quantile(PBR, probs = 0.975), .groups = "keep")
print(PBR.summaries, n = 50, digits = 2)

## Plot
pdf('Plots/IPMPBR_sampleNmin.pdf', width = 6.5, height = 5.5)
ggplot(PBR.dataA, aes(x = PBR)) + 
  geom_density(aes(color = as.factor(Fr), fill = as.factor(Fr)), alpha = 0.15) + 
  scale_fill_viridis(discrete = T, option = 'E', name = 'Fr') + 
  scale_color_viridis(discrete = T, option = 'E', name = 'Fr') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~Rmax_Label, ncol = 2, scale = 'free') + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

####################################################
# PBR CALCULATION - Nmin from a-priori calculation #
####################################################

## Using Rmax based on max(r_lambda)
PBR5 <- calc.PBR(Nmin = Nmin.2002, lambda = rlambda_max, Fr = Fr)
PBR5$Rmax_Source <- 'max(r_lambda)'

## Using Rmax based on mean(r_lambda)
PBR6 <- calc.PBR(Nmin = Nmin.2002, lambda = rlambda_mean, Fr = Fr)
PBR6$Rmax_Source <- 'mean(r_lambda)'

## Using Rmax based on max(a_lambda)
PBR7 <- calc.PBR(Nmin = Nmin.2002, lambda = alambda_max, Fr = Fr)
PBR7$Rmax_Source <- 'max(a_lambda)'

## Using Rmax based on mean(a_lambda)
PBR8 <- calc.PBR(Nmin = Nmin.2002, lambda = alambda_mean, Fr = Fr)
PBR8$Rmax_Source <- 'mean(a_lambda)'


## Combine data
PBR.dataB <- rbind(PBR5, PBR6, PBR7, PBR8)
PBR.dataB$Rmax_Label <- dplyr::case_when(PBR.dataB$Rmax_Source == 'max(r_lambda)' ~ 'Realized growth rate (max)',
                                         PBR.dataB$Rmax_Source == 'mean(r_lambda)' ~ 'Realized growth rate (mean)',
                                         PBR.dataB$Rmax_Source == 'max(a_lambda)' ~ 'Asymptotic growth rate (max)',
                                         PBR.dataB$Rmax_Source == 'mean(a_lambda)' ~ 'Asymptotic growth rate (mean)',)
PBR.dataB$Rmax_Label <- factor(PBR.dataA$Rmax_Label, levels = c('Realized growth rate (max)', 'Realized growth rate (mean)', 'Asymptotic growth rate (max)', 'Asymptotic growth rate (mean)'))

## Plot
pdf('Plots/IPMPBR_aprioriNmin.pdf', width = 6.5, height = 5.5)
ggplot(PBR.dataB, aes(x = PBR)) + 
  geom_density(aes(color = as.factor(Fr), fill = as.factor(Fr)), alpha = 0.15) + 
  scale_fill_viridis(discrete = T, option = 'E', name = 'Fr') + 
  scale_color_viridis(discrete = T, option = 'E', name = 'Fr') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~Rmax_Label, ncol = 2, scale = 'free') + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

#############################################
# PLOTTING POPULATION GROWTH RATES AND NMIN #
#############################################

## Assemble data
Nmin.data <- data.frame(
  Value = Ntot_min
)

lambda.data <- data.frame(
  Parameter = rep(c('Realized growth rate (max)', 'Realized growth rate (mean)', 'Asymptotic growth rate (max)', 'Asymptotic growth rate (mean)'), each = length(rlambda_max)),
  Value = c(rlambda_max, rlambda_mean, alambda_max, alambda_mean)
)

## Plot Nmin
p1 <- ggplot(Nmin.data, aes(x = Value)) + 
  geom_density(color = 'grey40', fill = 'grey40', alpha = 0.5) + 
  geom_vline(aes(xintercept = Nmin.2002), linetype = 'dashed', color = 'hotpink') + 
  ggtitle('A) Minimum population size') + 
  theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

## Plot
p2 <- ggplot(lambda.data, aes(x = Value)) + 
  geom_density(aes(color = Parameter, fill = Parameter), alpha = 0.5) + 
  ggtitle('B) Population growth rates') + 
  geom_vline(xintercept = 1, color = 'grey40', linetype = 'dotted') +
  scale_color_manual(values = c('turquoise4', 'turquoise3', 'slateblue4', 'slateblue2')) +
  scale_fill_manual(values = c('turquoise4', 'turquoise3', 'slateblue4', 'slateblue2')) +
  theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

pdf('Lambda_Nmin_Posteriors.pdf', width = 9, height = 3)
grid.arrange(p1, p2, nrow = 1, widths = c(0.5, 1))
dev.off()