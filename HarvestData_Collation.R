library(magrittr)
library(lubridate)
library(dplyr)

options(tibble.width = Inf)

## Set data path
DataPath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/Data/'

##############################
# AGECLASS-AT-HARVEST MATRIX #
##############################

## Load collated seal data
load(paste0(DataPath, '210906_DemoData_Combined.RData'))

## Reformat data
SealData <- SealData %>%
  
  # Add additional temporal variables
  dplyr::mutate(SamplingPeriod = case_when(year(Date) %in% c(1981, 1982) ~ '1981-82',
                                           year(Date) %in% c(2002:2004) ~ '2002-04',
                                           TRUE ~ '2012-20'),
                Year = year(Date),
                SealYear = ifelse(month(Date) < 6, year(Date)-1, year(Date)),
                JulianDay = yday(Date)) %>%
  
  # Add (seal) year index
  dplyr::mutate(YearIdx = Year - 1981 + 1,
                SealYearIdx = SealYear - 1981 + 1) %>%
  
  
  # Keep only entries of females with complete age / reproductive state information
  dplyr::filter(Sex == 'F' & !is.na(Age) & !is.na(Maturity)) %>%

  # Assign age class (index)
  dplyr::mutate(AgeClass = case_when(Age == 0 ~ 'YOY',
                                     Age %in% c(1:5) & Maturity == 0 ~ paste0('SubA[', Age, ']'),
                                     Age > 2 ~ 'MatA'),
                AgeClassIdx = case_when(Age == 0 ~ 1,
                                        Age %in% c(1:5) & Maturity == 0 ~ Age+1,
                                        Age > 2 ~ 7))

## Extract preliminary ageclass-at-harvest matrix
prel.ACaH <- table(SealData$AgeClassIdx, SealData$YearIdx)

## Write complete ageclass-at-harvest matrix (incl. years without sampling)
ACaH <- matrix(NA, nrow = max(SealData$AgeClassIdx), ncol = max(SealData$YearIdx),
               dimnames = list(c('YOY', paste0('SubA[', 1:5, ']'), 'MatA'), min(SealData$Year):max(SealData$Year)))

for(t in 1:dim(ACaH)[2]){
  if(t %in% dimnames(prel.ACaH)[[2]]){
    ACaH[,t] <- prel.ACaH[,which(dimnames(prel.ACaH)[[2]]==t)]
  }  
}

## Count numbers harvested by year and seal year
No.ACaH.yr <- colSums(ACaH)

No.ACaH.syr <- as.vector(table(SealData$SealYear)[paste0(1980:2020)])
names(No.ACaH.syr) <- 1980:2020

################################
# HARVEST COUNT DATA - OVERALL #
################################

## Read collated count data
Hcounts <- read.csv(paste0(DataPath, '220120_NoSealsHarvested_DataSourceComp.csv'), fileEncoding="UTF-8-BOM")

## Add a new column containing the maximum number reported across data sources
Hcounts$No_max <- NA
for(i in 1:nrow(Hcounts)){
  Hcounts$No_max[i] <- max(Hcounts[i,2:4], na.rm = T)
}

## Write vector of harvest numbers
no.H <- c(NA, Hcounts$No_max)


####################################
# HARVEST COUNT DATA - BY LOCATION #
####################################

## Read iNatur data (2019-2021)
iNatur <- read.csv(paste0(DataPath, '200120_iNatur_SealHarvest_2019-21.csv'), fileEncoding="UTF-8-BOM")

## Read location key
Hloc <- read.csv(paste0(DataPath, '200120_SealHarvest_LocationKey.csv'), fileEncoding="UTF-8-BOM")

## Merge location key into data
iNatur <- merge(iNatur, Hloc, by = 'Location')

## Calculate the numbers of ringed seals harvested per year overall and in Isfjorden area only
Years <- c(2019:2021)
count.all <- count.Isfj <- rep(NA, length(Years))

for(i in 1:length(Years)){
  count.all[i] <- sum(subset(iNatur, Year == Years[i])$Ringsel, na.rm = T)
  count.Isfj[i] <- sum(subset(iNatur, Year == Years[i] & IsfjordenArea == 1)$Ringsel, na.rm = T)
}

count.Isfj/count.all

#################
# DATA ASSEMBLY #
#################

## Collect data in a list
HarvestData <- list(
  ACaH = ACaH,
  No.ACaH.yr = No.ACaH.yr,
  No.ACaH.syr = No.ACaH.syr,
  no.H = no.H,
  count.all = count.all,
  count.Isfj = count.Isfj
)

## Save data
saveRDS(HarvestData, file = '220124_HarvestCountData.rds')
