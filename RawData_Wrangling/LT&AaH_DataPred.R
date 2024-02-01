library(magrittr)
library(lubridate)
library(dplyr)

#####################
# GENERAL DATA PREP #
#####################

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
  
  # Discard entries without sex/age
  dplyr::filter(!is.na(Sex) & !is.na(Age))


#######################
# LIFE TABLE ASSEMBLY #
#######################

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

## Store life tables for different sex-period combinations in a list
LT.data <- list(
  FM.All = make.LifeTable(data = SealData),
  FM.P1 = make.LifeTable(data = subset(SealData, SamplingPeriod == '1981-82')),
  FM.P2 = make.LifeTable(data = subset(SealData, SamplingPeriod == '2002-04')),
  FM.P3 = make.LifeTable(data = subset(SealData, SamplingPeriod == '2012-20')),
  F.All = make.LifeTable(data = subset(SealData, Sex == 'F')),
  F.P1 = make.LifeTable(data = subset(SealData, SamplingPeriod == '1981-82' & Sex == 'F')),
  F.P2 = make.LifeTable(data = subset(SealData, SamplingPeriod == '2002-04' & Sex == 'F')),
  F.P3 = make.LifeTable(data = subset(SealData, SamplingPeriod == '2012-20' & Sex == 'F'))
)

## Save life table data
saveRDS(LT.data, file = 'LifeTable_Analyses/LifeTables.rds')

#####################
# AaH DATA ASSEMBLY #
#####################

## Function to assemble age-at-harvest data
make.AaH <- function(data){
  
  # Make empty data frame
  AaH.data <- data.frame(
    matrix(NA, nrow = max(data$Age)+1, ncol = 2, 
           dimnames = list(NULL, c('Age', 'Number'))))
  
  # Add ages
  AaH.data$Age <- 0:max(data$Age)
  
  # Add numbers harvested per age from data 
  AaH.data$Number <- c(table(factor(data$Age, levels = 0:max(data$Age))))
  
  # Convert numbers to relative frequencies
  AaH.data$Frequency <- AaH.data$Number / sum(AaH.data$Number)
  
  # Return data
  return(AaH.data)
}

## Store Aah data for different sex-period combinations in a list
AaH.data <- list(
  FM.All = make.AaH(data = SealData),
  FM.P1 = make.AaH(data = subset(SealData, SamplingPeriod == '1981-82')),
  FM.P2 = make.AaH(data = subset(SealData, SamplingPeriod == '2002-04')),
  FM.P3 = make.AaH(data = subset(SealData, SamplingPeriod == '2012-20')),
  F.All = make.AaH(data = subset(SealData, Sex == 'F')),
  F.P1 = make.AaH(data = subset(SealData, SamplingPeriod == '1981-82' & Sex == 'F')),
  F.P2 = make.AaH(data = subset(SealData, SamplingPeriod == '2002-04' & Sex == 'F')),
  F.P3 = make.AaH(data = subset(SealData, SamplingPeriod == '2012-20' & Sex == 'F'))
)

## Save age-at-harvest data
saveRDS(AaH.data, file = 'CatchCurve_Analyses/AaHData.rds')
