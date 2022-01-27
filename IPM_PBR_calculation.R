library(coda)
library(ggplot2)
library(viridis)
library(gridExtra)
library(magrittr)

## Load posterior samples from model run
load('220125_Seal_IPM_noPeriodEff_Test.RData')
#load('220125_Seal_IPM_noPeriodEff_unknownSpup_Test.RData')

out.mat <- as.matrix(testRun)

## Calculation of (posterior distributions for) Nmin and Rmax

# Set variable names to consider
Tmin <- 22
Tmax <- 40
Ntot.names <- paste0('Ntot[', Tmin:Tmax, ']')
lambda.names <- paste0('lambda_real[', Tmin:(Tmax-1), ']')

# Extract maximum and mean realized lambda for each posterior sample 
rlambda_max <- apply(out.mat[, lambda.names], 1, max)
rlambda_mean <- apply(out.mat[, lambda.names], 1, mean)

hist(rlambda_max)
hist(rlambda_mean)

# Extract asymptotic lambda for each posterior sample
alambda <- out.mat[,'lambda_asym']
hist(alambda)


# Extract the minimum population size across study period for each sample
Nmin <- apply(out.mat[, Ntot.names]*2, 1, min) # Multiply by 2 to get up to full popualtion size

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


## Calculation of overall PBR for different values of Fr
Fr <- seq(0.1, 1.0, by = 0.1)


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
PBR1 <- calc.PBR(Nmin = Nmin, lambda = rlambda_max, Fr = Fr)
PBR1$Rmax_Source <- 'max(r_lambda)'

## Using Rmax based on mean(r_lambda)
PBR2 <- calc.PBR(Nmin = Nmin, lambda = rlambda_mean, Fr = Fr)
PBR2$Rmax_Source <- 'mean(r_lambda)'

## Using Rmax based on a_lambda
PBR3 <- calc.PBR(Nmin = Nmin, lambda = alambda, Fr = Fr)
PBR3$Rmax_Source <- 'a_lambda'

## Combine data
PBR.dataA <- rbind(PBR1, PBR2, PBR3)
PBR.dataA$Rmax_Label <- dplyr::case_when(PBR.dataA$Rmax_Source == 'max(r_lambda)' ~ 'Realized growth rate (maximum)',
                                         PBR.dataA$Rmax_Source == 'mean(r_lambda)' ~ 'Realized growth rate (mean)',
                                         PBR.dataA$Rmax_Source == 'a_lambda' ~ 'Asymptotic growth rate',)
PBR.dataA$Rmax_Label <- factor(PBR.dataA$Rmax_Label, levels = c('Realized growth rate (maximum)', 'Realized growth rate (mean)', 'Asymptotic growth rate'))

## Plot
pdf('220125_IPMPBR_sampleNmin.pdf', width = 5, height = 8)
#pdf('220125_IPMPBR_unknownSpup_sampleNmin.pdf', width = 5, height = 8)
ggplot(PBR.dataA, aes(x = PBR)) + 
  geom_density(aes(color = as.factor(Fr), fill = as.factor(Fr)), alpha = 0.15) + 
  scale_fill_viridis(discrete = T, option = 'E', name = 'Fr') + 
  scale_color_viridis(discrete = T, option = 'E', name = 'Fr') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~Rmax_Label, ncol = 1, scale = 'free') + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

####################################################
# PBR CALCULATION - Nmin from a-priori calculation #
####################################################

## Using Rmax based on max(r_lambda)
PBR4 <- calc.PBR(Nmin = Nmin.2002, lambda = rlambda_max, Fr = Fr)
PBR4$Rmax_Source <- 'max(r_lambda)'

## Using Rmax based on mean(r_lambda)
PBR5 <- calc.PBR(Nmin = Nmin.2002, lambda = rlambda_mean, Fr = Fr)
PBR5$Rmax_Source <- 'mean(r_lambda)'

## Using Rmax based on a_lambda
PBR6 <- calc.PBR(Nmin = Nmin.2002, lambda = alambda, Fr = Fr)
PBR6$Rmax_Source <- 'a_lambda)'

## Combine data
PBR.dataB <- rbind(PBR4, PBR5, PBR6)
PBR.dataB$Rmax_Label <- dplyr::case_when(PBR.dataB$Rmax_Source == 'max(r_lambda)' ~ 'Realized growth rate (maximum)',
                                         PBR.dataB$Rmax_Source == 'mean(r_lambda)' ~ 'Realized growth rate (mean)',
                                         PBR.dataB$Rmax_Source == 'a_lambda' ~ 'Asymptotic growth rate',)
PBR.dataB$Rmax_Label <- factor(PBR.dataB$Rmax_Label, levels = c('Realized growth rate (maximum)', 'Realized growth rate (mean)', 'Asymptotic growth rate'))

## Plot
pdf('220125_IPMPBR_aprioriNmin.pdf', width = 5, height = 8)
#pdf('220125_IPMPBR_unknownSpup_aprioriNmin.pdf', width = 5, height = 8)
ggplot(PBR.dataB, aes(x = PBR)) + 
  geom_density(aes(color = as.factor(Fr), fill = as.factor(Fr)), alpha = 0.15) + 
  scale_fill_viridis(discrete = T, option = 'E', name = 'Fr') + 
  scale_color_viridis(discrete = T, option = 'E', name = 'Fr') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~Rmax_Label, ncol = 1, scale = 'free') + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

#############################################
# PLOTTING POPULATION GROWTH RATES AND NMIN #
#############################################

## Assemble data
Nmin.data <- data.frame(
  Value = Nmin
)

lambda.data <- data.frame(
  Parameter = rep(c('Realized growth rate (maximum)', 'Realized growth rate (mean)', 'Asymptotic growth rate'), each = length(rlambda_max)),
  Value = c(rlambda_max, rlambda_mean, alambda)
)

## Plot Nmin
p1 <- ggplot(Nmin.data, aes(x = Value)) + 
  geom_density(color = 'grey40', fill = 'grey40', alpha = 0.5) + 
  ggtitle('A) Minimum population size') + 
  theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

## Plot
p2 <- ggplot(lambda.data, aes(x = Value)) + 
  geom_density(aes(color = Parameter, fill = Parameter), alpha = 0.5) + 
  ggtitle('B) Population growth rates') + 
  geom_vline(xintercept = 1, color = 'grey40', linetype = 'dashed') +
  scale_color_manual(values = c('#8C085E', '#00A69D', 'orange')) +
  scale_fill_manual(values = c('#8C085E', '#00A69D', 'orange')) +
  theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

pdf('220127_Lambda_Nmin_Posteriors.pdf', width = 9, height = 3)
grid.arrange(p1, p2, nrow = 1, widths = c(0.5, 1))
dev.off()