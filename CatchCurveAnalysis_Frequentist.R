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


############################
# CATCH CURVE CALCULATIONS #
############################

## Function to assemble age-at-harvest data
make.AaH <- function(data, trail0, plot){
  
  # Make empty data frame
  AaH.data <- data.frame(
    matrix(NA, nrow = max(data$Age)+1+trail0, ncol = 2, 
           dimnames = list(NULL, c('Age', 'Number'))))
  
  # Add ages
  AaH.data$Age <- 0:(max(data$Age)+trail0)
  
  # Add numbers harvested per age from data 
  AaH.data$Number <- c(table(factor(data$Age, levels = 0:max(data$Age))), rep(0, trail0))
  
  # Convert numbers to relative frequencies
  AaH.data$Frequency <- AaH.data$Number / sum(AaH.data$Number)
  
  if(plot){
    # Plot frequencies on natural and log scales
    plot.nat <- ggplot(AaH.data, aes(x = Age, y = Frequency)) + 
      geom_point() + theme_bw()
    
    plot.log <- ggplot(AaH.data, aes(x = Age, y = log(Frequency))) + 
      geom_point() + theme_bw()
    
    gridExtra::grid.arrange(plot.nat, plot.log, ncol = 1) 
  }
  
  # Return data
  return(AaH.data)
}

## Function for extracting/plotting survival estimates
analyse.CC <- function(AaH.data, startAge){
  
  ## Take a subset depending on threshold
  data <- subset(AaH.data, Age >= startAge)
  
  ## Run analysis using log-linear model on frequencies

  # Add a small number to 0 values
  data$Frequency[which(data$Frequency==0)] <- 0.001
    
  # Fit model and print summary
  mod1 <- lm(log(data$Frequency) ~ data$Age)
  
  message('Linear model on log(Frequency):')
  print(summary(mod1)$coefficients)
  message(paste0('R-squared: ', round(summary(mod1)$r.squared, digits = 4)))
  message('')
  
  # Plot data and fit
  plot.fit.lm <- ggplot(data, aes(x = Age, y = log(Frequency))) + 
                 geom_point() + 
                 geom_smooth(method = 'lm', formula = 'y ~ x', color = 'orange', fill = 'orange', alpha = 0.5) + 
                 ggtitle('Linear model on log(Frequency)') + theme_bw()
    
  # Simulate distribution of survival parameter
  h.sim <- rnorm(10000, mean = mod1$coefficients[2], summary(mod1)$coefficients[2,'Std. Error'])
  S.sim <- data.frame(Estimate = exp(h.sim), Model = 'LM')
  
  message('Survival:')
  message(paste0('Mean = ', round(mean(S.sim$Estimate), digits = 4)))
  message(paste0('SD = ', round(sd(S.sim$Estimate), digits = 4)))  
  message('')
  
  ## Run analysis using Poisson GLM on numbers
    
  # Fit model and print summary
  mod2 <- glm(data$Number ~ data$Age, family = poisson(link = 'log'))
  
  message('Poisson GLM on Number:')
  print(summary(mod2)$coefficients)
  message('')
  
  # Plot data and fit
  plot.fit.glm <- ggplot(data, aes(x = Age, y = Number)) + 
    geom_point() + 
    geom_smooth(method = 'glm', method.args = list(family = poisson(link = 'log')), formula = 'y ~ x', color = 'hotpink', fill = 'hotpink', alpha = 0.5) +
    ggtitle('Poisson GLM on Number') + theme_bw()
    
  # Simulate and plot distribution of survival parameter
  h.sim <- rnorm(10000, mean = mod2$coefficients[2], summary(mod2)$coefficients[2,'Std. Error'])
  S.sim <- rbind(S.sim, data.frame(Estimate = exp(h.sim), Model = 'GLM'))
  
  message('Survival:')
  message(paste0('Mean = ', round(mean(S.sim$Estimate), digits = 4)))
  message(paste0('SD = ', round(sd(S.sim$Estimate), digits = 4)))  
  
  ## Visualize fit and predicted survival distribution
  plot.Sdist <- ggplot(S.sim, aes(x = Estimate)) + 
    geom_density(aes(color = Model, fill = Model), alpha = 0.5) + 
    scale_color_manual(values = c('hotpink', 'orange')) + 
    scale_fill_manual(values = c('hotpink', 'orange')) + 
    geom_vline(xintercept = exp(mean(mod1$coefficients[2])), color = 'orange', size = 1) + 
    geom_vline(xintercept = exp(mean(mod2$coefficients[2])), color = 'hotpink', size = 1) + 
    xlab('Adult survival probability') + theme_bw()
  
  gridExtra::grid.arrange(plot.fit.lm, plot.fit.glm, plot.Sdist, ncol = 1)
} 


## Catch curve analyses - both sexes

# All sampling periods pooled
AaH.data <- make.AaH(SealData, trail0 = 0, plot = F)
analyse.CC(AaH.data, startAge = 1)
analyse.CC(AaH.data, startAge = 10)
  
# Sampling period 81-82
AaH.data <- make.AaH(subset(SealData, SamplingPeriod == '1981-82'), trail0 = 0, plot = F)
analyse.CC(AaH.data, startAge = 1)
analyse.CC(AaH.data, startAge = 10)

# Sampling period 02-04
AaH.data <- make.AaH(subset(SealData, SamplingPeriod == '2002-04'), trail0 = 0, plot = F)
analyse.CC(AaH.data, startAge = 1)
analyse.CC(AaH.data, startAge = 10)

# Sampling period 12-20
AaH.data <- make.AaH(subset(SealData, SamplingPeriod == '2012-20'), trail0 = 0, plot = F)
analyse.CC(AaH.data, startAge = 1)
analyse.CC(AaH.data, startAge = 10)


## Life table analyses - females only

# All sampling periods pooled
AaH.data <- make.AaH(subset(SealData, Sex == 'F'), trail0 = 0, plot = F)
analyse.CC(AaH.data, startAge = 1)
analyse.CC(AaH.data, startAge = 10)

# Sampling period 81-82
AaH.data <- make.AaH(subset(SealData, SamplingPeriod == '1981-82' & Sex == 'F'), trail0 = 0, plot = F)
analyse.CC(AaH.data, startAge = 1)
analyse.CC(AaH.data, startAge = 10)

# Sampling period 02-04
AaH.data <- make.AaH(subset(SealData, SamplingPeriod == '2002-04' & Sex == 'F'), trail0 = 0, plot = F)
analyse.CC(AaH.data, startAge = 1)
analyse.CC(AaH.data, startAge = 10)

# Sampling period 12-20
AaH.data <- make.AaH(subset(SealData, SamplingPeriod == '2012-20'  & Sex == 'F'), trail0 = 0, plot = F)
analyse.CC(AaH.data, startAge = 1)
analyse.CC(AaH.data, startAge = 10)

