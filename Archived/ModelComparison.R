library(coda)
library(reshape2)
library(ggplot2)
library(viridis)

#################################
# COLLATING POSTERIOR ESTIMATES #
#################################

## Load posterior samples from models of interest
#sam.modA <- readRDS('220429_IPMtest_fSAD_ice.rds')
#sam.modB <- readRDS('220429_IPMtest_eHAD_ice.rds')
#sam.modC <- readRDS('220429_IPMtest_eHAD_extraHdata_ice.rds')
sam.modA <- readRDS('220429_IPMtest_fSAD_ice.rds')
sam.modB <- readRDS('220502_IPMtest_fSAD_ice.rds')
sam.modC <- readRDS('220502_IPMtest_fSAD_ice_2.rds')
sam.modD <- readRDS('220503_IPMtest_fSAD_ice_3.rds')

## Re-format posterior samples
data.modA <- melt(as.matrix(sam.modA))
data.modB <- melt(as.matrix(sam.modB))
data.modC <- melt(as.matrix(sam.modC))
data.modD <- melt(as.matrix(sam.modD))

## Add model identifier
#data.modA$Model <- 'fSAD_ice'
#data.modB$Model <- 'eHAD_ice'
#data.modC$Model <- 'eHAD_extraHdata_ice'
data.modA$Model <- 'base'
data.modB$Model <- '+ ice var.'
data.modC$Model <- '+ surv. var.'
data.modD$Model <- '+ ice & surv. var.'

## Combine dataframes
#data.comb <- rbind(data.modA, data.modB, data.modC)
data.comb <- rbind(data.modA, data.modB, data.modC, data.modD)
colnames(data.comb) <- c('Sample', 'Parameter', 'Estimate', 'Model')

##################################################
# PLOTTING POSTERIOR DISTRIBUTIONS ACROSS MODELS #
##################################################

## Defining sets of parameters to plot
Surv.params <- c('S_YOY', 'mN_YOY', 'mH_YOY',
                 paste0('S_SA[', 1:5, ']'), paste0('mN_SA[', 1:5, ']'), paste0('mH_SA[', 1:5, ']'),
                 'S_MA', 'mN_MA', 'mH_MA', paste0('alpha[', 1:5, ']'))
Rep.params <- c(paste0('Mu.pMat[', 3:5, ']'), 'sigmaY.pMat', 'pOvl', 'pPrg', 'S_pup.ideal', 'ice.ideal')
Spup.params <- c(paste0('S_pup[', 1:40, ']'))
Ntot.params <- c('estN.2002', paste0('Ntot[', 22:40, ']'))
SAD.params <- paste0('SAD[', 1:8, ']')
lam.params <- c('lambda_asym', paste0('lambda_real[', 22:39, ']'))

## Plot comparisons to pdf
pdf('220503_ModelComparison_Uncertainty_eHAD_extraHdata.pdf', width = 11, height = 8)

# Survival parameters
ggplot(subset(data.comb, Parameter%in%Surv.params), aes(x = Estimate)) + 
  geom_density(aes(color = Model, fill = Model), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Survival parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Reproduction
ggplot(subset(data.comb, Parameter%in%Rep.params), aes(x = Estimate)) + 
  geom_density(aes(color = Model, fill = Model), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Reproduction parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Pup survival
ggplot(subset(data.comb, Parameter%in%Spup.params), aes(x = Estimate)) + 
  geom_density(aes(color = Model, fill = Model), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Pup survival parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Population size parameters
ggplot(subset(data.comb, Parameter%in%Ntot.params), aes(x = Estimate)) + 
  geom_density(aes(color = Model, fill = Model), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Population size parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Population structure parameters
ggplot(subset(data.comb, Parameter%in%SAD.params), aes(x = Estimate)) + 
  geom_density(aes(color = Model, fill = Model), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Population structure parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Population growth parameters
ggplot(subset(data.comb, Parameter%in%lam.params), aes(x = Estimate)) + 
  geom_density(aes(color = Model, fill = Model), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Population growth parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

dev.off()


##################################
# PLOTTING ESTIMATES ACROSS TIME #
##################################

## Collect posterior samples in a matrix
out.mat <- list(
  out.matA = as.matrix(sam.modA),
  out.matB = as.matrix(sam.modB),
  out.matC = as.matrix(sam.modC),
  out.matD = as.matrix(sam.modD)
)

## Set up
sim.Tmin <- 22
sim.Tmax <- 40
ModelType <- c('base', '+ ice var.', '+ surv. var', '+ice & surv. var.')
ann.est <- data.frame()

## Collate annual posterior summaries for relevant parameters
for(m in 1:length(out.mat)){
  
  post.mat <- out.mat[[m]]
  Amax <- 7
  
  for(t in 1:sim.Tmax){
    
    Spup.sum <- quantile(post.mat[,paste0('S_pup[',t,']')], probs = c(0.025, 0.5, 0.975))
    
    if(t >= sim.Tmin){
      Ntot.sum <- quantile(post.mat[,paste0('Ntot[',t,']')], probs = c(0.025, 0.5, 0.975))
    }else{
      Ntot.sum <- rep(NA, 3)
    }
    
    if(t >= sim.Tmin & t < sim.Tmax){
      H.sum <- quantile(
        rowSums(post.mat[,paste0('H[',1:Amax, ', ',t,']')]), 
        probs = c(0.025, 0.5, 0.975))
    }else{
      H.sum <- rep(NA, 3)
    }
    
    post.data <- data.frame(rbind(Ntot.sum, H.sum, Spup.sum))
    colnames(post.data) <- c('lCI', 'Median', 'uCI')
    rownames(post.data) <- NULL
    post.data$Year = t+1980
    post.data$Parameter = c('Ntot', 'H', 'S_pup')
    post.data$Model <- ModelType[m]
    
    ann.est <- rbind(ann.est, post.data)
  }
}

## Plot comparisons to pdf
pdf('220503_ModelComparisonTime_Uncertainty_eHAD_extraHdata.pdf', width = 11, height = 6)
ggplot(ann.est) + 
  geom_line(aes(x = Year, y = Median, color = Model)) + 
  geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 0.1) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  facet_wrap(~Parameter, scales = 'free_y', ncol = 1) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'top')
dev.off()