library(ggplot2)
library(viridis)
library(fields)

## Write data for Fjord-specific population size estimates (from Krafft et al. 2006)
N.est <- data.frame(
  FjordID = c(5:11),
  EstMean = c(1142, 577, 0, 268, 1082, 137, 36),
  EstSD = c(187, 99, 0, 52, 149, 29, 7)
)

## Function for simulating Isfjorden population size
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

## Simulate Isfjorden population size
simN <-Isfj.Nsim(n = 1000000, N.est = N.est)
hist(simN, breaks = 100)

## Extract Nmin for PBR calculations
#  The recommendation (e.g. in Nelson et al. 2019) is to use the 20th percentile
#  of a log-normal distribution of estimates of population size
Nmin <- exp(quantile(log(simN), probs = 0.2))

## Set realistic range for Rmax
Rmax_default <- 0.12 # This is the suggested default for pinnipeds
Rmax_llim <- 0
Rmax_ulim <- 2*Rmax_default 

## Set the range for Fr
Fr_llim <- 0
Fr_ulim <- 1

## Function for calculating PBR across ranges of Rmax and Fr
calc.PBR <- function(Nmin, Rmax, Fr){
  
  # Prepare matrix for storing results
  PBR <- matrix(NA, nrow = length(Fr), ncol = length(Rmax),
                dimnames = list(round(Fr, digits = 2), round(Rmax, digits = 2)))
  
  # Calculate PBR for each Rmax-Fr combination
  for(i in 1:length(Fr)){
    for(j in 1:length(Rmax)){
      PBR[i,j] <- Nmin*0.5*Rmax[j]*Fr[i]
    }
  }
  
  return(PBR)
}

## Determine ranges for Rmax and Fr
Rmax_n <- 101
Fr_n <- 101
Rmax <- seq(Rmax_llim, Rmax_ulim, length.out = Rmax_n)
Fr <- seq(Fr_llim, Fr_ulim, length.out = Fr_n)

## Calculate PBR value matrix
PBR.mat <- calc.PBR(Nmin, Rmax, Fr)

## Make plotting directory if it does not exist
if(!file.exists("Plots")){
  dir.create("Plots")
}

## Plot PBR heatmap
pdf('Plots/RudimentaryPBR.pdf', width = 7, height = 6)
#par(mar=c(5.1,5.1,4.1,2.1))
image.plot(PBR.mat, col = cividis(100), main = 'PBR', xlab = 'Fr', ylab = "", xaxt= "n", yaxt= "n", cex.lab = 1.3, cex.main = 1.5, axis.args = list(cex.axis = 1.2))
axis(1, at = seq(0, 1, length.out = 6), labels = seq(0, 1, length.out = 6), cex.axis = 1.2)
axis(2, at = seq(0, 1, length.out = 6), labels = round(seq(0, Rmax_ulim, length.out = 6), digits = 2), cex.axis = 1.2)
contour(PBR.mat, add = TRUE, levels = c(10, 20, 30, 50, 100, 200, 300), col = 'white', labcex = 1, method = 'edge', lty = 1, crt = 90, drawlabels = F)
abline(h = 0.5, col = 'skyblue1', lty = 3, lwd = 1.5)
title(ylab= expression('R'[max]), line = 2.5, cex.lab = 1.3)
dev.off()