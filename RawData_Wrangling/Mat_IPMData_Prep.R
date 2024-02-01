#-----------------#
# MATURATION DATA #
#-----------------#

## Set data path
DataPath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/Data/'

## Load collated seal data
load(paste0(DataPath, '210906_DemoData_Combined.RData'))

## Add additional temporal variables
SealData <- SealData %>%
  dplyr::mutate(SamplingPeriod = case_when(year(Date) %in% c(1981, 1982) ~ '1981-82',
                                           year(Date) %in% c(2002:2004) ~ '2002-04',
                                           TRUE ~ '2012-20'),
                Year = year(Date),
                JulianDay = yday(Date))

## Take a subset of females only
SealDataF <- subset(SealData, Sex == 'F' & !is.na(Age) & !is.na(Maturity))

## Add age class variables relevant for maturation
SealDataF <- SealDataF %>%
  dplyr::mutate(AgeClass = ifelse(Age %in% c(3:5), Age, NA))
# This assumes age classes 0, 1, 2, 3, 4, 5, 6+,
# with maturation rates of 0 for age classes 0-2, and 1 for age class 6+

## Subset to keep only age classes relevant for fitting of maturation model
SealDataF_mat <- subset(SealDataF, !is.na(AgeClass))


#-----------#
# COLLATION #
#-----------#

save(SealDataF_mat, file = '220114_SealIPM_Data.RData')
