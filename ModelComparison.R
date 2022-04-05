library(coda)
library(reshape2)
library(ggplot2)
library(viridis)

#################################
# COLLATING POSTERIOR ESTIMATES #
#################################

## Load posterior samples from models of interest
#load('220221_Seal_IPM_noPeriodEff_ext_Test.RData')
#load('220222_Seal_IPM_noPeriodEff_ext_FullRangePupS.RData')
load('220223_Seal_IPM_noPeriodEff_ext_noSADconstraint.RData')
sam.modA <- testRun

#load('220222_Seal_IPM_noPeriodEff_ext_FullRangePupS.RData')
#load('220222_Seal_IPM_noPeriodEff_ext_mHest.RData')
#load('220222_Seal_IPM_noPeriodEff_ext_noSADconstraint.RData')
load('220223_Seal_IPM_noPeriodEff_ext_mHest.RData')
sam.modB <- testRun

#load('220223_Seal_IPM_noPeriodEff_ext_noSADconstraint.RData')
#load('220223_Seal_IPM_noPeriodEff_ext_mHest.RData')
load('220223_Seal_IPM_noPeriodEff_ext_mHest&noSADconstraint.RData')
sam.modC <- testRun

## Re-format posterior samples
data.modA <- melt(as.matrix(sam.modA))
#data.modA$Model <- 'Ext-base'
#data.modA$Model <- 'Ext-FullRangePupS'
data.modA$Model <- 'Ext-FullRangePupS-noSADconstraint'

data.modB <- melt(as.matrix(sam.modB))
#data.modB$Model <- 'Ext-FullRangePupS'
#data.modB$Model <- 'Ext-mHest'
#data.modB$Model <- 'Ext-noSADconstraint'
data.modB$Model <- 'Ext-FullRangePupS-mHest'

data.modC <- melt(as.matrix(sam.modC))
#data.modC$Model <- 'Ext-FullRangePupS-noSADconstraint'
#data.modC$Model <- 'Ext-FullRangePupS-mHest'
data.modC$Model <- 'Ext-FullRangePupS-mHest-noSADconstraint'

## Combine dataframes
data.comb <- rbind(data.modA, data.modB, data.modC)
colnames(data.comb) <- c('Sample', 'Parameter', 'Estimate', 'Model')

#################################################
# PLOTTING POSERIOR DISTRIBUTIONS ACROSS MODELS #
#################################################

## Defining sets of parameters to plot
Surv.params <- c('S_YOY', 'mN_YOY', 'mH_YOY',
                 paste0('S_SA[', 1:5, ']'), paste0('mN_SA[', 1:5, ']'), paste0('mH_SA[', 1:5, ']'),
                 'S_MA', 'mN_MA', 'mH_MA')
Rep.params <- c(paste0('Mu.pMat[', 3:5, ']'), 'sigmaY.pMat', 'pOvl', 'pPrg', 'S_pup[1]')
Ntot.params <- c('estN.2002', paste0('Ntot[', 22:40, ']'))
SAD.params <- paste0('SAD[', 1:8, ']')
lam.params <- c('lambda_asym', paste0('lambda_real[', 22:39, ']'))

## Plot comparisons to pdf
pdf('220223_ModelComparison.pdf', width = 11, height = 8)

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

