library(coda)
library(ggplot2)
library(viridis)
library(gridExtra)
library(tidyverse)
library(lubridate)

## Load posterior samples from model run
load('212015_Seal_IPM_Test_noPeriodEff.RData')
out.mat <- as.matrix(testRun)

## Collate SAD posteriors in a data frame
SAD.data <- data.frame(
  Proportion = c(out.mat[,'SAD[1]'], out.mat[,'SAD[2]'], out.mat[,'SAD[3]'], 
                 out.mat[,'SAD[4]'], out.mat[,'SAD[5]'], out.mat[,'SAD[6]'], 
                 out.mat[,'SAD[7]'], out.mat[,'SAD[8]'], out.mat[,'SAD[9]']),
  AgeClass = rep(c('YOY', paste0('SubA[', 1:5, ']'), 'nMatA', 'MatA'), each = dim(out.mat)[1])
)

SAD.data$AgeClass <- factor(SAD.data$AgeClass, levels = c('YOY', paste0('SubA[', 1:5, ']'), 'nMatA', 'MatA'))


## Plot posterior distributions for SAD
ggplot(SAD.data, aes(x = Proportion)) + 
        geom_density(aes(y = ..scaled.., color = AgeClass, fill = AgeClass), alpha = 0.25) + 
        scale_fill_viridis(discrete = T) + 
        scale_color_viridis(discrete = T) + 
        theme_bw() + theme(panel.grid = element_blank())
# --> This is pretty useless...

## Calculate quantiles
round(quantile(out.mat[,'SAD[1]'], probs = c(0.025, 0.5, 0.975)), digits = 3)
round(quantile(out.mat[,'SAD[2]'], probs = c(0.025, 0.5, 0.975)), digits = 3)
round(quantile(out.mat[,'SAD[3]'], probs = c(0.025, 0.5, 0.975)), digits = 3)
round(quantile(out.mat[,'SAD[4]'], probs = c(0.025, 0.5, 0.975)), digits = 3)
round(quantile(out.mat[,'SAD[5]'], probs = c(0.025, 0.5, 0.975)), digits = 3)
round(quantile(out.mat[,'SAD[6]'], probs = c(0.025, 0.5, 0.975)), digits = 3)
round(quantile(out.mat[,'SAD[7]'], probs = c(0.025, 0.5, 0.975)), digits = 3)
round(quantile(out.mat[,'SAD[8]'], probs = c(0.025, 0.5, 0.975)), digits = 3)

## Calculate observed proportions from raw data
DataPath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/Data/'
load(paste0(DataPath, '210906_DemoData_Combined.RData'))

SealData <- SealData %>%
  dplyr::mutate(SamplingPeriod = case_when(year(Date) %in% c(1981, 1982) ~ '1981-82',
                                           year(Date) %in% c(2002:2004) ~ '2002-04',
                                           TRUE ~ '2012-20'),
                Year = year(Date),
                JulianDay = yday(Date))

SealDataF <- subset(SealData, Sex == 'F' & !is.na(Age) & !is.na(Maturity))

SealDataF <- SealDataF %>%
  dplyr::mutate(AgeClass = case_when(Age == 0 ~ 'YOY',
                                     Age %in% c(1:8) & Maturity == 0 ~ paste0('SubA[', Age, ']'),
                                     Age > 2 ~ 'MatA'))


## Add collapsed age class variable
SealDataF$AgeClassColl <- ifelse(SealDataF$AgeClass %in% c(6:8), 6, SealDataF$AgeClass)

## Subset to keep only age classes relevant for fitting of maturation model
SealDataF_sub <- subset(SealDataF, !is.na(AgeClass))

# The remainder of the code seems to be missing here...