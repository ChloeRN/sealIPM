library(magrittr)
library(lubridate)
library(dplyr)

options(tibble.width = Inf)

## Set data path
DataPath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/Data/'

###############################################
# AGECLASS-AT-HARVEST MATRIX - WHOLE SVALBARD #
###############################################

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

  # Assign conservative age class (index)
  dplyr::mutate(AgeClass = case_when(Age == 0 ~ 'YOY',
                                     Age %in% c(1:5) & Maturity == 0 ~ paste0('SubA[', Age, ']'),
                                     Age > 2 ~ 'MatA'),
                AgeClassIdx = case_when(Age == 0 ~ 1,
                                        Age %in% c(1:5) & Maturity == 0 ~ Age+1,
                                        Age > 2 ~ 7)) %>%
  
  # Assign alternative age class (index) - making max. use of available data
  dplyr::mutate(AgeClass2 = case_when(Age == 0 ~ 'YOY',
                                      Age %in% c(1:5) & Maturity == 0 ~ paste0('SubA[', Age, ']'),
                                      Age == 3 & Maturity == 1 ~ 'nMatA',
                                      Age > 6 ~ 'MatA',
                                      TRUE ~ paste0('uMatA[', Age, ']')),
                AgeClass2Idx = case_when(Age == 0 ~ 1,
                                         Age %in% c(1:5) & Maturity == 0 ~ Age+1,
                                         Age == 3 & Maturity == 1 ~ 7,
                                         Age > 6 ~ 8,
                                         TRUE ~ Age+1))
                
## Make data subset containing only individuals with certain age class
SealData_k <- subset(SealData, !(AgeClass2 %in% paste0('uMatA[', 4:6, ']')))

## Make data subset containing only individuals with uncertain age class
SealData_u <- subset(SealData, AgeClass2 %in% paste0('uMatA[', 4:6, ']'))

## Extract preliminary ageclass-at-harvest matrices
prel.ACaH <- table(SealData$AgeClassIdx, SealData$YearIdx)
prel.ACaH_k <- table(SealData_k$AgeClass2Idx, SealData_k$YearIdx)
prel.uMatA <- table(SealData_u$AgeClass2Idx, SealData_u$YearIdx)

## Write complete ageclass-at-harvest matrices (incl. years without sampling)
ACaH <- matrix(0, nrow = max(SealData$AgeClassIdx), ncol = max(SealData$YearIdx),
               dimnames = list(c('YOY', paste0('SubA[', 1:5, ']'), 'MatA'), min(SealData$Year):max(SealData$Year)))
ACaH_k <- matrix(0, nrow = max(SealData_k$AgeClass2Idx, na.rm = T), ncol = max(SealData_k$YearIdx),
               dimnames = list(c('YOY', paste0('SubA[', 1:5, ']'),'nMatA', 'MatA'), min(SealData_k$Year):max(SealData_k$Year)))
uMatA <- matrix(0, nrow = 3, ncol = max(SealData$YearIdx),
               dimnames = list(paste0('uMatA[', 4:6, ']'), min(SealData$Year):max(SealData$Year)))

for(t in 1:dim(ACaH)[2]){
  if(t %in% dimnames(prel.ACaH)[[2]]){
    ACaH[,t] <- prel.ACaH[,which(dimnames(prel.ACaH)[[2]]==t)]
  }
  if(t %in% dimnames(prel.ACaH2)[[2]]){
    ACaH_k[,t] <- prel.ACaH_k[,which(dimnames(prel.ACaH2)[[2]]==t)]
  }
  if(t %in% dimnames(prel.uMatA)[[2]]){
    uMatA[,t] <- prel.uMatA[,which(dimnames(prel.uMatA)[[2]]==t)]
  }
}

## Count numbers harvested by year and seal year
No.ACaH.yr <- colSums(ACaH)

No.ACaH.syr <- as.vector(table(SealData$SealYear)[paste0(1980:2020)])
names(No.ACaH.syr) <- 1980:2020
No.ACaH.syr[which(is.na(No.ACaH.syr))] <- 0


###############################################
# AGECLASS-AT-HARVEST MATRIX - ISFJORDEN AREA #
###############################################

## Assign dummy coordinates to location names (https://placenames.npolar.no/)
LocationNames <- sort(unique(SealData$Location[which(!is.na(SealData$Location))]))
DummyCoord <- matrix(rep(NA, length(LocationNames)*2), ncol = 2, dimnames =  list(LocationNames, c('DummyLongitude', 'DummyLatitude')))

DummyCoord[1,] <- c(15.553468, 78.258484) # Adventfjorden
DummyCoord[2:3,] <- rep(c(16.306963, 78.53041), each = 2) # Billefjorden (x2)
DummyCoord[4,] <- c(11.80264, 78.836876) # Engelsbukta
DummyCoord[5,] <- c(11.6226015, 78.60231) # Forlandssundet
DummyCoord[6,] <- c(11.835841, 78.69465) # Hornb?kbukta
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
no.H <- c(rep(NA, 2002-1981+1), Hcounts$No_max)
names(no.H) <- 1981:2021


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
names(count.all) <- names(count.Isfj) <- Years

for(i in 1:length(Years)){
  count.all[i] <- sum(subset(iNatur, Year == Years[i])$Ringsel, na.rm = T)
  count.Isfj[i] <- sum(subset(iNatur, Year == Years[i] & IsfjordenArea == 1)$Ringsel, na.rm = T)
}

count.Isfj/count.all


##################
# REPORTING RATE #
##################

# NOTE: These numbers are based on "guesstimates" from Kit Kovacs and Christian Lydersen

r <- c(rep(0.75, 2009-1981+1), rep(0.95, 2021-2010+1))
names(r) <- c(1981:2021)


#################
# DATA ASSEMBLY #
#################

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
  r = r,
  
  ACaH_k = ACaH_k,
  uMatA = uMatA
)

## Save data
saveRDS(HarvestData, file = 'HarvestCountData.rds')
