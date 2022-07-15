library(coda)
library(reshape2)
library(ggplot2)
library(viridis)

#################################
# COLLATING POSTERIOR ESTIMATES #
#################################

## Load posterior samples from models of interest
sam.modA <- as.matrix(readRDS('220520_IPMtest_eHAD_iceSim.rds'))
sam.modB <- as.matrix(readRDS('220524_IPMtest_eHAD_iceSim_halfH.rds'))
sam.modC <- as.matrix(readRDS('220524_IPMtest_eHAD_iceSim_noH.rds'))

## Rename parameters where necessary
param.old <- c('mH_YOY',  paste0('mH_SA[', 1:5, ']'), 'mH_MA',
               'S_YOY',  paste0('S_SA[', 1:5, ']'), 'S_MA',
               paste0('alpha[', 1:7, ']'))
colnames(sam.modA)[which(colnames(sam.modA) %in% param.old)]

param.new <- c('S_MA[1]', paste0('S_SA[', 1:5, ', 1]'), 'S_YOY[1]', 
               paste0('alpha[', 1:7, ', 1]'),
               'mH_MA[1]',paste0('mH_SA[', 1:5, ', 1]'),   'mH_YOY[1]')

colnames(sam.modA)[which(colnames(sam.modA) %in% param.old)] <- param.new


## Retain only parameters required for plotting
params1 <- c(
  paste0('SAD[', 1:8, ']'),
  paste0('Ntot[', 22:70, ']'), paste0('lambda_real[', 1:69,']'),
  paste0('Mu.pMat[', 3:5,']'), 'sigmaY.pMat',
  'pOvl', 'pPrg',
  'mN_YOY', paste0('mN_SA[', 1:5, ']'), 'mN_MA',
  'mH_YOY[1]', paste0('mH_SA[', 1:5, ', 1]'), 'mH_MA[1]',
  'S_YOY[1]', paste0('S_SA[', 1:5, ', 1]'), 'S_MA[1]',
  paste0('alpha[', 1:7, ', 1]'),
  'S_pup.ideal', paste0('S_pup[', 1:70,']'), 'env.ideal',
  paste0('env[', 1:70, ']'), 'Mu.env', 'beta.env', 'sigmaY.env'
)

params2 <- c(
  paste0('SAD[', 1:8, ']'),
  paste0('Ntot[', 22:70, ']'), paste0('lambda_real[', 1:69,']'),
  paste0('Mu.pMat[', 3:5,']'), 'sigmaY.pMat',
  'pOvl', 'pPrg',
  'mN_YOY', paste0('mN_SA[', 1:5, ']'), 'mN_MA',
  'mH_YOY[1]', paste0('mH_SA[', 1:5, ', 1]'), 'mH_MA[1]',
  'mH_YOY[2]', paste0('mH_SA[', 1:5, ', 2]'), 'mH_MA[2]',
  'S_YOY[1]', paste0('S_SA[', 1:5, ', 1]'), 'S_MA[1]',
  'S_YOY[2]', paste0('S_SA[', 1:5, ', 2]'), 'S_MA[2]',
  paste0('alpha[', 1:7, ', 1]'), paste0('alpha[', 1:7, ', 2]'),
  'S_pup.ideal', paste0('S_pup[', 1:70,']'), 'env.ideal',
  paste0('env[', 1:70, ']'), 'Mu.env', 'beta.env', 'sigmaY.env'
)

## Re-format posterior samples
data.modA <- melt(sam.modA[, params1])
data.modB <- melt(sam.modB[, params2])
data.modC <- melt(sam.modC[, params2])

## Add information on harvest scenario
data.modA$Harvest <- 'Unchanged'
data.modB$Harvest <- 'Half'
data.modC$Harvest <- 'None'

## Add information on environment scenario
data.modA$Environment <- 'Stable'
data.modB$Environment <- 'Stable'
data.modC$Environment <- 'Stable'

## Combine dataframes
#data.comb <- rbind(data.modA, data.modB, data.modC)
data.comb <- rbind(data.modA, data.modB, data.modC)
colnames(data.comb) <- c('Sample', 'Parameter', 'Estimate', 'Harvest', 'Environment')

#####################################################
# PLOTTING POSTERIOR DISTRIBUTIONS ACROSS SCENARIOS #
#####################################################

## Defining sets of parameters to plot
Surv.params <- c('mN_YOY', paste0('mN_SA[', 1:5, ']'), 'mN_MA',
                 'mH_YOY[1]', paste0('mH_SA[', 1:5, ', 1]'), 'mH_MA[1]',
                 'mH_YOY[2]', paste0('mH_SA[', 1:5, ', 2]'), 'mH_MA[2]',
                 'S_YOY[1]', paste0('S_SA[', 1:5, ', 1]'), 'S_MA[1]',
                 'S_YOY[2]', paste0('S_SA[', 1:5, ', 2]'), 'S_MA[2]',
                 paste0('alpha[', 1:7, ', 1]'), paste0('alpha[', 1:7, ', 2]'))
Rep.params <- c(paste0('Mu.pMat[', 3:5, ']'), 'sigmaY.pMat', 'pOvl', 'pPrg', 'S_pup.ideal', 'ice.ideal')
SpupEnv.params <- c('S_pup.ideal', 'env.ideal', 'Mu.env', 'beta.env', 'sigmaY.env')
SAD.params <- paste0('SAD[', 1:8, ']')

## Plot comparisons to pdf
pdf('220525_ModelComparison_Scenarios.pdf', width = 11, height = 8)

# Survival parameters
ggplot(subset(data.comb, Parameter%in%Surv.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Survival parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Reproduction
ggplot(subset(data.comb, Parameter%in%Rep.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Reproduction parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Pup survival 6 environment
ggplot(subset(data.comb, Parameter%in%SpupEnv.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Pup survival & environment parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Population structure parameters
ggplot(subset(data.comb, Parameter%in%SAD.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  ggtitle('Population structure parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

dev.off()


##########################################################
# PLOTTING ESTIMATES ACROSS TIME FOR DIFFERENT SCENARIOS #
##########################################################

## Collect posterior samples in a matrix
out.mat <- list(
  out.matA = as.matrix(sam.modA),
  out.matB = as.matrix(sam.modB),
  out.matC = as.matrix(sam.modC)
)

## Set up
sim.Tmin <- 22
sim.Tmax <- 70
ModelType <- c('Unchanged', 'Halved', 'None')
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

## Re-arrange factor levels
ann.est$Model <- factor(ann.est$Model, levels = c('Unchanged', 'Halved', 'None'))

## Plot comparisons to pdf
pdf('220715_ModelComparisonTime_Scenarios.pdf', width = 11, height = 6)
ggplot(ann.est) + 
  geom_line(aes(x = Year, y = Median, color = Model)) + 
  geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 0.1) + 
  scale_color_viridis(discrete = T, name = 'Harvest') + 
  scale_fill_viridis(discrete = T, name = 'Harvest') + 
  facet_wrap(~Parameter, scales = 'free_y', ncol = 1) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'top')
dev.off()
