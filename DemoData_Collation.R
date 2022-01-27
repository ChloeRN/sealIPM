library(tidyverse)
library(ggplot2)
library(readxl)

options(tibble.width = Inf)

## Set data path
DataPath <- 'C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/SealHarvest/Data/'

## List strings to identify as NA
NAstrings <- c('', '*', '?', 'na')

## Load and reformat raw data (1981 - 1982)
raw1 <- read_xlsx(paste0(DataPath, 'Svalbard 1981-82.xlsx'), 
                  sheet = 1, range = 'B1:J285', na = NAstrings) %>%
  dplyr::rename('IndvID' = 'Animal ID',
                'Date' = 'Date (DMY)',
                'Length' = 'Standard L',
                'Girth' = 'Ax.Girth',
                'Mass' = 'Body mass') %>%
  dplyr::mutate(IndvID = paste0(year(Date), '-', IndvID),
                Maturity = case_when(Maturity == 'I' ~ 0,
                                     Maturity == 'M' ~ 1))

## Load and reformat raw data (2002)
raw2a <- read_xls(paste0(DataPath, 'Svalbard 2002-04.xls'), 
                  sheet = 1, range = 'A1:K73', na = NAstrings) %>%
  dplyr::rename('IndvID' = 'Anim. No',
                'Location' = 'Area',
                'Latitude' = 'GPS positon N',
                'Longitude' = 'GPS positon E',
                'Length' = 'St. length',
                'Mass' = 'Body mass') %>%
  dplyr::mutate(IndvID = paste0(year(Date), '-', IndvID),
                Maturity = case_when(Maturity == 'I' ~ 0,
                                     Maturity == 'M' ~ 1),
                Latitude = ifelse(as.numeric(Latitude) > 100, as.numeric(Latitude)/10e4, as.numeric(Latitude)),
                Longitude = as.numeric(Longitude))

## Load and reformat raw data (2003)
raw2b <- read_xls(paste0(DataPath, 'Svalbard 2002-04.xls'), 
                 sheet = 2, range = 'B1:L94', na = NAstrings)[-1,] %>%
  dplyr::rename('IndvID' = 'Anim. No',
                'Date' = 'Dato',
                'Age' = 'age',
                'Location' = 'Area',
                'Latitude' = 'GPS positon N',
                'Longitude' = 'GPS positon E',
                'Length' = 'St. Length',
                'Mass' = 'Body mass') %>%
  dplyr::mutate(IndvID = paste0(year(Date), '-', IndvID),
                Maturity = case_when(Maturity == 'I' ~ 0,
                                     Maturity == 'M' ~ 1),
                Sex = case_when(Sex == 'm' ~ 'M',
                                Sex == 'f' ~ 'F'),
                Latitude = Latitude/10e4,
                Longitude = Longitude/10e4)

## Load and reformat raw data (2004)
raw2c <- read_xls(paste0(DataPath, 'Svalbard 2002-04.xls'), 
                  sheet = 3, range = 'A1:K115', na = NAstrings)[-1,] %>%
  dplyr::rename('IndvID' = 'Anim. No',
                'Location' = 'area',
                'Latitude' = 'GPS positon N',
                'Longitude' = 'GPS positon E',
                'Length' = 'St. Length',
                'Mass' = 'Body mass') %>%
  dplyr::mutate(IndvID = paste0(year(Date), '-', IndvID),
                Maturity = case_when(Maturity == 'I' ~ 0,
                                     Maturity == 'M' ~ 1),
                Latitude = Latitude/10e4,
                Longitude = Longitude/10e4)

## Load and reformat raw data (2012 - 2020)
raw3 <- read_xlsx(paste0(DataPath, 'Svalbard 2012-20.xlsx'), 
                  sheet = 1, range = 'A1:N266', na = NAstrings) %>%
  dplyr::rename('IndvID' = 'Seal ID',
                'Latitude' = 'Lat',
                'Longitude' = 'Long',
                'Age' = 'Age (CL)',
                'Length' = 'Standard length (cm)',
                'Girth' = 'Girth (cm)',
                'Mass' = 'Body mass',
                'Blubber_back' = 'Blubber thickness (on back)',
                'Blubber_sternum' = 'Blubber thickness (over sternum)',
                'MaturityM' = 'Male sex mature',
                'MaturityF' = 'Female sex mature') %>%
  dplyr::mutate(IndvID = ifelse(year(Date) > 2016, paste0(year(Date), '-', IndvID), IndvID),
                Maturity = case_when(MaturityM == 0 | MaturityF == 0 ~ 0,
                                     MaturityM == 1 | MaturityF == 1 ~ 1),
                Age = ifelse(Age == '1+', NA_integer_, as.integer(Age))) %>%
  dplyr::select(-c(MaturityM, MaturityF))

## Bind together dataframes
SealData <- dplyr::bind_rows(raw1, raw2a, raw2b, raw2c, raw3)

## Save as .RData and as .csv
save(SealData, file = '210906_DemoData_Combined.RData')
write.csv(SealData, '210906_DemoData_Combined.csv', row.names = FALSE)
