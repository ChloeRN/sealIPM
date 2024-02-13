library(dplyr)
library(tidyverse)

## Set path to raw data
DataPath <- "RawData/"


################################################################################
# SEAL MATURATION DATA #
########################

## Read in collated seal data
SealData <- readRDS(paste0(DataPath, "DemoData_Combined.rds"))

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

## Drop unneccessary columns
SealDataF_mat <- SealDataF_mat %>%
  dplyr::select(IndvID, Year, AgeClass, Maturity)

## Save in input data folder
saveRDS(SealDataF_mat, file = "InputData/SealIPM_DemoData.rds")



################################################################################
# SEAL HARVEST DATA #
#####################

## Read in collated seal data
SealData <- readRDS(paste0(DataPath, "DemoData_Combined.rds"))

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


#-------------------------------------#
# Complete Ageclass-at-harvest matrix #
#-------------------------------------#

## Extract preliminary ageclass-at-harvest matrix
prel.ACaH <- table(SealData$AgeClassIdx, SealData$YearIdx)

## Write complete ageclass-at-harvest matrix (incl. years without sampling)
ACaH <- matrix(0, nrow = max(SealData$AgeClassIdx), ncol = max(SealData$YearIdx),
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
No.ACaH.syr[which(is.na(No.ACaH.syr))] <- 0


#---------------------------------------------#
# Ageclass-at-harvest matrix - Isfjorden Area #
#---------------------------------------------#

## Assign dummy coordinates to location names (https://placenames.npolar.no/)
LocationNames <- sort(unique(SealData$Location[which(!is.na(SealData$Location))]))
DummyCoord <- matrix(rep(NA, length(LocationNames)*2), ncol = 2, dimnames =  list(LocationNames, c('DummyLongitude', 'DummyLatitude')))

DummyCoord[1,] <- c(15.553468, 78.258484) # Adventfjorden
DummyCoord[2:3,] <- rep(c(16.306963, 78.53041), each = 2) # Billefjorden (x2)
DummyCoord[4,] <- c(11.80264, 78.836876) # Engelsbukta
DummyCoord[5,] <- c(11.6226015, 78.60231) # Forlandssundet
DummyCoord[6,] <- c(11.835841, 78.69465) # Hornbækbukta
DummyCoord[7,] <- c(14.78042, 78.249435) # Isfjorden
DummyCoord[8,] <- c(12.33183, 78.87152) # Kongsfjorden 
DummyCoord[9,] <- c(14.751617, 78.504364) # Nordfjorden
DummyCoord[10,] <- c(14.658182, 77.51814) # Recherchefjorden
DummyCoord[11,] <- c(11.432147, 78.70017) # Revflaket
DummyCoord[12,] <- c(16.422525, 78.376656) # Sassenfjorden
DummyCoord[13,] <- c(12.7511425, 78.515015) # St. Johnsfjorden
DummyCoord[14,] <- c(12.061151, 79.20379) # Tinayrebukta
DummyCoord[15,] <- c(15.42829, 77.55777) # Van Keulenfjorden
DummyCoord[16:18,] <- rep(c(15.450018, 79.62894), each = 3) # Wijdefjorden (x3)

DummyCoord <- data.frame(cbind(LocationNames, DummyCoord), row.names = NULL) %>%
  dplyr::rename(Location = 'LocationNames') %>%
  dplyr::mutate(DummyLongitude = as.numeric(DummyLongitude),
                DummyLatitude = as.numeric(DummyLatitude))

## Add dummy coordinates for locations missing GPS information & drop entries outside Isfjorden area
SealDataIsfj <- merge(SealData, DummyCoord, by = 'Location', all.x = T) %>%
  dplyr::mutate(Longitude2 = ifelse(is.na(Longitude), DummyLongitude, Longitude),
                Latitude2 = ifelse(is.na(Latitude), DummyLatitude, Latitude)) %>%
  dplyr::filter(Longitude2 > 13.5 & between(Latitude2, 78.2 ,78.885))

## Extract preliminary ageclass-at-harvest matrix and insert missing age classes 5 and 6 (age 4 and 5 subadults)
prel.ACaH.Isfj <- table(SealDataIsfj$AgeClassIdx, SealDataIsfj$YearIdx)
prel.ACaH.Isfj <- rbind(prel.ACaH.Isfj[1:4,], rep(0, dim(prel.ACaH.Isfj)[2]), rep(0, dim(prel.ACaH.Isfj)[2]), prel.ACaH.Isfj[5,])

## Write complete ageclass-at-harvest matrix (incl. years without sampling)
ACaH.Isfj <- matrix(0, nrow = max(SealDataIsfj$AgeClassIdx), ncol = max(SealDataIsfj$YearIdx),
                    dimnames = list(c('YOY', paste0('SubA[', 1:5, ']'), 'MatA'), min(SealDataIsfj$Year):max(SealDataIsfj$Year)))

for(t in 1:dim(ACaH.Isfj)[2]){
  if(t %in% dimnames(prel.ACaH.Isfj)[[2]]){
    ACaH.Isfj[,t] <- prel.ACaH.Isfj[,which(dimnames(prel.ACaH.Isfj)[[2]]==t)]
  }  
}

## Count numbers harvested by year and seal year
No.ACaH.Isfj.yr <- colSums(ACaH.Isfj)

No.ACaH.Isfj.syr <- as.vector(table(SealDataIsfj$SealYear)[paste0(1980:2020)])
names(No.ACaH.Isfj.syr) <- 1980:2020
No.ACaH.Isfj.syr[which(is.na(No.ACaH.Isfj.syr))] <- 0


#------------------------------#
# HARVEST COUNT DATA - OVERALL #
#------------------------------#

## Read collated count data
Hcounts <- read.csv(paste0(DataPath, "NoSealsHarvested_DataSourceComp.csv"), fileEncoding="UTF-8-BOM")

## Add a new column containing the maximum number reported across data sources
Hcounts$No_max <- NA
for(i in 1:nrow(Hcounts)){
  Hcounts$No_max[i] <- max(Hcounts[i,2:4], na.rm = T)
}

## Write vector of harvest numbers
no.H <- c(rep(NA, 2002-1981+1), Hcounts$No_max)
names(no.H) <- 1981:2021


#----------------------------------#
# HARVEST COUNT DATA - BY LOCATION #
#----------------------------------#

## Read iNatur data (2019-2021)
iNatur <- read.csv(paste0(DataPath, 'iNatur_SealHarvest_2019-21.csv'), fileEncoding = "UTF-8-BOM")

## Read location key
Hloc <- read.csv(paste0(DataPath, 'SealHarvest_LocationKey.csv'), fileEncoding = "UTF-8-BOM")

## Merge location key into data
iNatur <- merge(iNatur, Hloc, by = 'Location')

## Calculate the numbers of ringed seals harvested per year overall and in Isfjorden area only
Years <- c(2019:2021)
count.all <- count.Isfj <- rep(NA, length(Years))
names(count.all) <- names(count.Isfj) <- Years

for(i in 1:length(Years)){
  count.all[i] <- sum(subset(iNatur, Year == Years[i])$Ringsel, na.rm = T)
  count.Isfj[i] <- sum(subset(iNatur, Year == Years[i] & IsfjordenArea == 1)$Ringsel, na.rm = T)
}


#----------------#
# REPORTING RATE #
#----------------#

# NOTE: These numbers are based on "guesstimates" from Kit Kovacs and Christian Lydersen

r <- c(rep(0.75, 2009-1981+1), rep(0.95, 2021-2010+1))
names(r) <- c(1981:2021)


#---------------#
# DATA ASSEMBLY #
#---------------#

## Collect data in a list
HarvestData <- list(
  ACaH = ACaH,
  No.ACaH.yr = No.ACaH.yr,
  No.ACaH.syr = No.ACaH.syr,
  ACaH.Isfj = ACaH.Isfj,
  No.ACaH.Isfj.yr = No.ACaH.Isfj.yr,
  No.ACaH.Isfj.syr = No.ACaH.Isfj.syr,
  no.H = no.H,
  count.all = count.all,
  count.Isfj = count.Isfj,
  r = r
)

## Save data
saveRDS(HarvestData, file = "InputData/HarvestCountData.rds")


################################################################################
# SEA ICE DATA #
################

## Load ice data from Glen Liston
ice.data.GL <- read.delim(paste0(DataPath, "isfjorden_yearly_averages.dat"), header = FALSE, sep =  "")[,2:4]

colnames(ice.data.GL) <- c('year', 'IcePerc', 'LairHabitat')

## Reformat data 
ice.GL.sum <- reshape2::melt(ice.data.GL, id.vars = 'year')[,c(1,3,2)]
colnames(ice.GL.sum) <- c('year', 'Value', 'Variable')

## Add missing early years to data
ice.data.missing <- data.frame(year = 1981:(min(ice.data.GL$year)-1),
                               IcePerc = NA, LairHabitat = NA)
ice.data.GL <- rbind(ice.data.missing, ice.data.GL)

## Arrange data into vectors
ice.GL.mean <- subset(ice.data.GL, year %in% c(1981:2021))$IcePerc
lairH.GL.mean <- subset(ice.data.GL, year %in% c(1981:2021))$LairHabitat

ice.years <- 1981:2021

ice.period <- c(rep(1, 10), rep(2, 15), rep(3, 16))


## Collect vectors into list
ice.IPMdata <- list(
  ice.GL.mean = ice.GL.mean,
  lairH.GL.mean = lairH.GL.mean,
  ice.years = ice.years, ice.period = ice.period
)

## Save data list
saveRDS(ice.IPMdata, file = "InputData/SealIPM_IceData.rds")


