library(coda)
library(ggplot2)
library(reshape2)
library(viridis)

## Load MCMC samples from separately run survival models
sam.LTA <- readRDS('LifeTable_Analyses/LifeTableAnalysis_MCMC.rds')
sam.CCA <- readRDS('CatchCurve_Analyses/CatchCurveAnalysis_Truncated_MCMC.rds')
sam.HOE <- readRDS('HoeningModel_Sim/mO_HoeningSim.rds')

## Load MCMC samples from integrated survival model
sam.ISM <- readRDS('220119_ISM_Exp2.rds')


#######################################
# COMPARISON OF INDEPENDENT ESTIMATES #
#######################################

## Reformat data from Life Table and Catch Curve Analyses
data.LTA <- data.CCA <- data.frame()

for(i in 1:8){
  
  sub.LTA <- melt(as.matrix(sam.LTA[[i]]))
  colnames(sub.LTA) <- c('Sample', 'Parameter', 'Value')
  sub.LTA$Data <- names(sam.LTA)[[i]]
  
  sub.CCA <- melt(as.matrix(sam.CCA[[i]]))
  colnames(sub.CCA) <- c('Sample', 'Parameter', 'Value')
  sub.CCA$Data <- names(sam.CCA)[[i]]
  
  data.LTA <- rbind(data.LTA, sub.LTA)
  data.CCA <- rbind(data.CCA, sub.CCA)
}

data.LTA$Method <- 'Life table'
data.CCA$Method <- 'Catch curve'

## Reformat data from Hoening model simulation
mat.HOE <- exp(-as.matrix(sam.HOE))
data.HOE <- melt(mat.HOE)
colnames(data.HOE) <- c('Sample', 'Parameter', 'Value')
data.HOE$Data <- dplyr::case_when(data.HOE$Parameter == 'mO[1]' ~ 'F.P1',
                                  data.HOE$Parameter == 'mO[2]' ~ 'F.P2',
                                  data.HOE$Parameter == 'mO[3]' ~ 'F.P3')
data.HOE$Parameter <- 'S_MA'
data.HOE$Method <- 'Hoening model'

## Assemble data
SurvPriors1 <- rbind(data.LTA, data.CCA, data.HOE)
SurvPriors1 <- subset(SurvPriors1, !(Parameter == 'Mu'))

## Plot posterior distributions for all survival rates
pdf('Plots/SurvPriors_All_MethodComp.pdf', width = 10, height = 7.5)
ggplot(SurvPriors1, aes(x = Value)) + 
  geom_density(aes(color = Method, fill = Method, linetype = Data), alpha = 0.1) + 
  facet_wrap(~ Parameter, scales = 'free') + 
  scale_color_manual(values = c('#8C085E', '#00A69D', 'orange')) + 
  scale_fill_manual(values = c('#8C085E', '#00A69D', 'orange')) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

## "Zoom in" on posterior distributions for adults survival rate
pdf('Plots/SurvPriors_MA_MethodComp.pdf', width = 8, height = 5)
ggplot(subset(SurvPriors1, Parameter == 'S_MA'), aes(x = Value)) + 
  geom_density(aes(color = Method, fill = Method, linetype = Data), alpha = 0.1) + 
  xlim(0.775, 1) + 
  scale_color_manual(values = c('#8C085E', '#00A69D', 'orange')) + 
  scale_fill_manual(values = c('#8C085E', '#00A69D', 'orange')) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

## "Zoom in" on posterior distributions for adults survival rate - Female only
pdf('Plots/SurvPriors_F_MA_MethodComp.pdf', width = 8, height = 5)
ggplot(subset(SurvPriors1, Parameter == 'S_MA' & Data%in%c('F.All', 'F.P1', 'F.P2', 'F.P3')), aes(x = Value)) + 
  geom_density(aes(color = Method, fill = Method, linetype = Data), alpha = 0.1) + 
  xlim(0.775, 1) + 
  scale_color_manual(values = c('#8C085E', '#00A69D', 'orange')) + 
  scale_fill_manual(values = c('#8C085E', '#00A69D', 'orange')) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()


######################################################
# COMPARISON OF INDEPENDENT AND INTEGRATED ESTIMATES #
######################################################

## Reformat data from Life Table and Catch Curve Analyses
data.LTA <- data.CCA <- data.frame()

i <- 5
sub.LTA <- melt(as.matrix(sam.LTA[[i]]))
colnames(sub.LTA) <- c('Sample', 'Parameter', 'Value')
  
sub.CCA <- melt(as.matrix(sam.CCA[[i]]))
colnames(sub.CCA) <- c('Sample', 'Parameter', 'Value')
  
data.LTA <- rbind(data.LTA, sub.LTA)
data.CCA <- rbind(data.CCA, sub.CCA)

data.LTA$Method <- 'Life table'
data.CCA$Method <- 'Catch curve'

## Reformat data from Hoening model simulation
mat.HOE <- as.matrix(sam.HOE)
data.HOE <- melt(mat.HOE)
colnames(data.HOE) <- c('Sample', 'Parameter', 'Value')
data.HOE <- subset(data.HOE, Parameter == 'mO[1]')
data.HOE$Parameter <- 'mO_MA'
data.HOE$Method <- 'Hoening model'

## Reformat data from integrated model
data.ISM <- melt(as.matrix(sam.ISM))
colnames(data.ISM) <- c('Sample', 'Parameter', 'Value')
data.ISM$Method <- 'Integrated'

## Assemble data
SurvPriors2 <- rbind(data.LTA, data.CCA, data.HOE, data.ISM)
SurvPriors2 <- subset(SurvPriors2, !(Parameter == 'Mu'))

## Plot posterior distributions for all survival rates
pdf('Plots/SurvPriors_All_IndepVSInt2.pdf', width = 10, height = 10)
ggplot(SurvPriors2, aes(x = Value)) + 
  geom_density(aes(color = Method, fill = Method), alpha = 0.1) + 
  facet_wrap(~ Parameter, scales = 'free') + 
  scale_color_manual(values = c('#8C085E', '#00A69D', 'grey40', 'orange')) + 
  scale_fill_manual(values = c('#8C085E', '#00A69D', 'grey40', 'orange')) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

## "Zoom in" on posterior distributions for adults survival rate
pdf('Plots/SurvPriors_MA_IndepVSInt2.pdf', width = 8, height = 5)
ggplot(subset(SurvPriors2, Parameter == 'S_MA'), aes(x = Value)) + 
  geom_density(aes(color = Method, fill = Method), alpha = 0.1) + 
  scale_color_manual(values = c('#8C085E', 'grey40', 'orange')) + 
  scale_fill_manual(values = c('#8C085E', 'grey40', 'orange')) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()
