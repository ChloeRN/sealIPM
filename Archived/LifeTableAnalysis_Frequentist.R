library(magrittr)
library(lubridate)
library(dplyr)
library(ggplot2)

####################
# DATA PREPARATION #
####################

## Set data path
DataPath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/Data/'

## Load collated seal data
load(paste0(DataPath, '210906_DemoData_Combined.RData'))

## Reformat data
SealData <- SealData %>%
  
  # Add additional temporal variables
  dplyr::mutate(SamplingPeriod = case_when(year(Date) %in% c(1981, 1982) ~ '1981-82',
                                           year(Date) %in% c(2002:2004) ~ '2002-04',
                                           TRUE ~ '2012-20'),
                Year = year(Date),
                JulianDay = yday(Date)) %>%
  
  # Discard entries without sex/age/maturation status
  dplyr::filter(!is.na(Sex) & !is.na(Age) & !is.na(Maturity)) %>%

  # Add alternative age which lumps individuals older than 5 years 
  dplyr::mutate(Age2 = ifelse(Age > 5, 6, Age)) %>%
  
  # Add age class as used in the IPM
  dplyr::mutate(AgeClass = ifelse(Maturity == 1, 6, Age2))


###########################
# LIFE TABLE CALCULATIONS #
###########################

## Function to assemble life tables
make.LifeTable <- function(data){
  
  # Make empty life table
  LT <- matrix(NA, nrow = max(data$Age)+2, ncol = 6, 
               dimnames = list(NULL, c('Age', 'nx', 'lx', 'dx', 'qx', 'px')))
  
  # Add ages
  LT[,'Age'] <- 0:(max(data$Age)+1)
  
  # Extract numbers harvested per age from data and enter as dx
  LT[,'dx'] <- c(table(factor(data$Age, levels = 0:max(data$Age))),0)
  
  # Calculate numbers alive at start of interval (nx)
  LT[1,'nx'] <- sum(LT[,'dx'])
  for(x in 1:(max(data$Age))){
    LT[x+1,'nx'] <- LT[x,'nx'] - LT[x,'dx']
  }
  
  # Calculate proportion surviving (lx)
  LT[,'lx'] <- LT[,'nx']/LT[1,'nx']
  
  # Calculate finite rate of mortality (qx)
  LT[,'qx'] <- LT[,'dx']/LT[,'nx']
  
  # Calculate finite rate of survival (px)
  LT[,'px'] <- 1 - LT[,'qx']
  
  # Return life table
  return(LT)
}

## Function for extracting/plotting survival estimates
LT.surv <- function(LT){
  
  # Set up data frame
  surv.data <- data.frame(AgeClass = c(1:5, rep('6+',5)),
                          mean = NA, sd = NA,
                          shape1 = NA, shape2 = NA,
                          dist = c(rep(NA, 5), 'Normal', 'Hazard-LogNormal', 'Beta-mle', 'Beta-mme', 'Beta-mmue'))
  
  # Add point estimates for ages 1-5
  surv.data$mean[1:5] <- LT[which(LT[,'Age']%in%c(1:5)), 'px']
  
  # Add mean and sd for age class 6+
  px6 <- LT[which(LT[,'Age']%in%(6:(max(LT[,'Age'])-2))), 'px']
  surv.data$mean[6] <- mean(px6, na.rm = T)
  surv.data$sd[6] <- sd(px6, na.rm = T)
  
  # Add logmean and logsd for age class 6+ (hazard scale)
  px6.red <- ifelse(px6 < 1, px6, 0.99)
  surv.data$mean[7] <- mean(log(-log(px6.red)), na.rm = T)
  surv.data$sd[7] <- sd(log(-log(px6.red)), na.rm = T)
  
  # Fit beta distributions to the numbers
  surv.data[8,c('shape1', 'shape2')] <- EnvStats::ebeta(px6)$parameter
  surv.data[9,c('shape1', 'shape2')] <- EnvStats::ebeta(px6, method = 'mme')$parameter
  surv.data[10,c('shape1', 'shape2')] <- EnvStats::ebeta(px6, method = 'mmue')$parameter
  
  # Print data
  print(surv.data)
  
  # Simulate values from hazard log-normal distribution
  sim.HlogN <- data.frame(sim = exp(-rlnorm(10000, meanlog = surv.data$mean[7], sdlog = surv.data$sd[7])))
    
  # Simulate values from beta distributions
  sim.beta1 <- data.frame(sim = rbeta(10000, shape1 = surv.data$shape1[8], shape2 = surv.data$shape2[8]))
  sim.beta2 <- data.frame(sim = rbeta(10000, shape1 = surv.data$shape1[9], shape2 = surv.data$shape2[9]))
  sim.beta3 <- data.frame(sim = rbeta(10000, shape1 = surv.data$shape1[10], shape2 = surv.data$shape2[10]))
  
  
  # Plot histogram for distribution of adult survival rates including beta distribution overlay
  ggplot() +
    geom_histogram(data = data.frame(px6), aes(x = px6, y = stat(ncount)), binwidth = 0.025) +
    geom_density(data = data.frame(px6), aes(x = px6, y = stat(density/max(density))), col = 'grey70', linetype = 'dotted') +
    
    geom_density(data = sim.HlogN, aes(x = sim, y = stat(density/max(density))), color = 'orange') +
    
    geom_density(data = sim.beta1, aes(x = sim, y = stat(density/max(density))), color = 'hotpink') +
    geom_density(data = sim.beta2, aes(x = sim, y = stat(density/max(density))), color = 'hotpink', linetype = 'dashed') +
    geom_density(data = sim.beta3, aes(x = sim, y = stat(density/max(density))), color = 'hotpink', linetype = 'dotted') +
     
    ylab('Relative likelihood') + xlab('Adult survival probability') + 
    theme_bw() + theme(panel.grid = element_blank())
} 


## Life table analyses - both sexes

# All sampling periods pooled
LT <- make.LifeTable(data = SealData)
LT
LT.surv(LT)

# Sampling period 81-82
LT <- make.LifeTable(data = subset(SealData, SamplingPeriod == '1981-82'))
LT
LT.surv(LT)

# Sampling period 02-04
LT <- make.LifeTable(data = subset(SealData, SamplingPeriod == '2002-04'))
LT
LT.surv(LT)

# Sampling period 12-20
LT <- make.LifeTable(data = subset(SealData, SamplingPeriod == '2012-20'))
LT
LT.surv(LT)


## Life table analyses - females only

# All sampling periods pooled
LT <- make.LifeTable(data = subset(SealData, Sex == 'F'))
LT
LT.surv(LT)

# Sampling period 81-82
LT <- make.LifeTable(data = subset(SealData, SamplingPeriod == '1981-82' & Sex == 'F'))
LT
LT.surv(LT)

# Sampling period 02-04
LT <- make.LifeTable(data = subset(SealData, SamplingPeriod == '2002-04' & Sex == 'F'))
LT
LT.surv(LT)

# Sampling period 12-20
LT <- make.LifeTable(data = subset(SealData, SamplingPeriod == '2012-20' & Sex == 'F'))
LT
LT.surv(LT)

