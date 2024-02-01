library(coda)
library(reshape2)
library(ggplot2)
library(viridis)
library(ggstance)

#################################
# COLLATING POSTERIOR ESTIMATES #
#################################

## Load posterior samples from models of interest
sam.mod_uH <- as.matrix(readRDS('IPMtest_fSAD&eHAD_iceSim.rds'))
sam.mod_hH <- as.matrix(readRDS('IPMtest_fSAD&eHAD_halfH_iceSim.rds'))
sam.mod_nH <- as.matrix(readRDS('IPMtest_fSAD&eHAD_noH_iceSim.rds'))

sam.mod_uH_T <- as.matrix(readRDS('IPMtest_fSAD&eHAD_iceSimTrend.rds'))
sam.mod_hH_T <- as.matrix(readRDS('IPMtest_fSAD&eHAD_halfH_iceSimTrend.rds'))
sam.mod_nH_T <- as.matrix(readRDS('IPMtest_fSAD&eHAD_noH_iceSimTrend.rds'))


## Retain only parameters required for plotting
params <- c(
  paste0('SAD[', 1:8, ']'),
  paste0('Ntot[', 22:70, ']'), paste0('lambda_real[', 1:69,']'),
  paste0('Mu.pMat[', 3:5,']'), 'sigmaY.pMat',
  'pOvl', 'pPrg',
  'mN_YOY', paste0('mN_SA[', 1:5, ']'), 'mN_MA',
  'mH_YOY[1]', paste0('mH_SA[', 1:5, ', 1]'), 'mH_MA[1]',
  'mH_YOY[2]', paste0('mH_SA[', 1:5, ', 2]'), 'mH_MA[2]',
  'S_YOY[1]', paste0('S_SA[', 1:5, ', 1]'), 'S_MA[1]',
  'S_YOY[2]', paste0('S_SA[', 1:5, ', 2]'), 'S_MA[2]',
  #paste0('alpha[', 1:7, ', 1]'), paste0('alpha[', 1:7, ', 2]'),
  'S_pup.ideal', paste0('S_pup[', 1:70,']'), 'env.ideal',
  paste0('env[', 1:70, ']'), 'Mu.env', 'beta.env', 'sigmaY.env'
)

## Re-format posterior samples
data.mod_uH <- melt(sam.mod_uH[, params])
data.mod_hH <- melt(sam.mod_hH[, params])
data.mod_nH <- melt(sam.mod_nH[, params])

data.mod_uH_T <- melt(sam.mod_uH_T[, params])
data.mod_hH_T <- melt(sam.mod_hH_T[, params])
data.mod_nH_T <- melt(sam.mod_nH_T[, params])

## Add information on harvest scenario
data.mod_uH$Harvest <- 'Unchanged'
data.mod_hH$Harvest <- 'Half'
data.mod_nH$Harvest <- 'None'

data.mod_uH_T$Harvest <- 'Unchanged'
data.mod_hH_T$Harvest <- 'Half'
data.mod_nH_T$Harvest <- 'None'

## Add information on environment scenario
data.mod_uH$Environment <- 'Stable'
data.mod_hH$Environment <- 'Stable'
data.mod_nH$Environment <- 'Stable'

data.mod_uH_T$Environment <- 'Trend'
data.mod_hH_T$Environment <- 'Trend'
data.mod_nH_T$Environment <- 'Trend'

## Combine dataframes
#data.comb <- rbind(data.modA, data.modB, data.modC)
data.comb <- rbind(data.mod_uH, data.mod_hH, data.mod_nH, data.mod_uH_T, data.mod_hH_T, data.mod_nH_T)
colnames(data.comb) <- c('Sample', 'Parameter', 'Estimate', 'Harvest', 'Environment')

## Sort factor levels
data.comb$Harvest <- factor(data.comb$Harvest, levels = c('Unchanged', 'Half', 'None'))

#####################################################
# PLOTTING POSTERIOR DISTRIBUTIONS ACROSS SCENARIOS #
#####################################################

## Defining sets of parameters to plot
Surv.params <- c('mN_YOY', paste0('mN_SA[', 1:5, ']'), 'mN_MA',
                 'mH_YOY[1]', paste0('mH_SA[', 1:5, ', 1]'), 'mH_MA[1]',
                 'mH_YOY[2]', paste0('mH_SA[', 1:5, ', 2]'), 'mH_MA[2]',
                 'S_YOY[1]', paste0('S_SA[', 1:5, ', 1]'), 'S_MA[1]',
                 'S_YOY[2]', paste0('S_SA[', 1:5, ', 2]'), 'S_MA[2]'#,
                 #paste0('alpha[', 1:7, ', 1]'), paste0('alpha[', 1:7, ', 2]')
                 )
Rep.params <- c(paste0('Mu.pMat[', 3:5, ']'), 'sigmaY.pMat', 'pOvl', 'pPrg', 'S_pup.ideal', 'ice.ideal')
SpupEnv.params <- c('S_pup.ideal', 'env.ideal', 'Mu.env', 'beta.env', 'sigmaY.env')
SAD.params <- paste0('SAD[', 1:8, ']')

## Determine colors for plotting
#plotColors <- rev(c('#7CCBA2', '#FCDE9C', '#F0746E'))
plotColors <- inferno(10)[c(2,5,8)]

## Make plotting directory if it does not exist
if(!file.exists("Plots")){
  dir.create("Plots")
}

## Plot comparisons to pdf
pdf('Plots/ModelComparison_Scenarios_CC0.pdf', width = 11, height = 8)

# Survival parameters
ggplot(subset(data.comb, Environment == 'Stable' & Parameter%in%Surv.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_manual(values = plotColors) +
  scale_fill_manual(values = plotColors) +
  ggtitle('Survival parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Reproduction
ggplot(subset(data.comb, Environment == 'Stable' & Parameter%in%Rep.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_manual(values = plotColors) +
  scale_fill_manual(values = plotColors) +
  ggtitle('Reproduction parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Pup survival 6 environment
ggplot(subset(data.comb, Environment == 'Stable' & Parameter%in%SpupEnv.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_manual(values = plotColors) +
  scale_fill_manual(values = plotColors) +
  ggtitle('Pup survival & environment parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Population structure parameters
ggplot(subset(data.comb, Environment == 'Stable' & Parameter%in%SAD.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_manual(values = plotColors) +
  scale_fill_manual(values = plotColors) +
  ggtitle('Population structure parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

dev.off()


pdf('Plots/ModelComparison_Scenarios_CCT.pdf', width = 11, height = 8)

# Survival parameters
ggplot(subset(data.comb, Environment != 'Stable' & Parameter%in%Surv.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_manual(values = plotColors) +
  scale_fill_manual(values = plotColors) +
  ggtitle('Survival parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Reproduction
ggplot(subset(data.comb, Environment != 'Stable' & Parameter%in%Rep.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_manual(values = plotColors) +
  scale_fill_manual(values = plotColors) +
  ggtitle('Reproduction parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Pup survival 6 environment
ggplot(subset(data.comb, Environment != 'Stable' & Parameter%in%SpupEnv.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_manual(values = plotColors) +
  scale_fill_manual(values = plotColors) +
  ggtitle('Pup survival & environment parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

# Population structure parameters
ggplot(subset(data.comb, Environment != 'Stable' & Parameter%in%SAD.params), aes(x = Estimate)) + 
  geom_density(aes(color = Harvest, fill = Harvest), alpha = 0.5) + 
  scale_color_manual(values = plotColors) +
  scale_fill_manual(values = plotColors) +
  ggtitle('Population structure parameters') + 
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid = element_blank())

dev.off()

##########################################################
# PLOTTING ESTIMATES ACROSS TIME FOR DIFFERENT SCENARIOS #
##########################################################

## Collect posterior samples in a matrix
out.mat <- list(
  out.mat_uH = as.matrix(sam.mod_uH),
  out.mat_hH = as.matrix(sam.mod_hH),
  out.mat_nH = as.matrix(sam.mod_nH),
  out.mat_uH_T = as.matrix(sam.mod_uH_T),
  out.mat_hH_T = as.matrix(sam.mod_hH_T),
  out.mat_nH_T = as.matrix(sam.mod_nH_T)
)

## Set up
sim.Tmin <- 22
sim.Tmax <- 70
ModelType <- rep(c('Unchanged', 'Half', 'None'), 2)
Environment <- rep(c('Stable', 'Trend'), each = 3)
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
      lambda.sum <- quantile(post.mat[,paste0('lambda_real[',t,']')], probs = c(0.025, 0.5, 0.975))
    }else{
      H.sum <- rep(NA, 3)
      lambda.sum <- rep(NA, 3)
    }
    
    post.data <- data.frame(rbind(Ntot.sum, lambda.sum, H.sum, Spup.sum))
    colnames(post.data) <- c('lCI', 'Median', 'uCI')
    rownames(post.data) <- NULL
    post.data$Year = t+1980
    post.data$Parameter = c('Ntot', 'lambda', 'H', 'S_pup')
    post.data$Model <- ModelType[m]
    post.data$Environment <- Environment[m]
    
    ann.est <- rbind(ann.est, post.data)
  }
}

## Re-arrange factor levels
ann.est$Model <- factor(ann.est$Model, levels = c('Unchanged', 'Half', 'None'))

## Add full-text Parameter label
ann.est$Label <- dplyr::case_when(ann.est$Parameter == "Ntot" ~ "Female population size",
                                  ann.est$Parameter == "lambda" ~ "Realized population growth rate",
                                  ann.est$Parameter == "H" ~ "Number of females harvested",
                                  ann.est$Parameter == "S_pup" ~ "Pup survival")

## Plot lambda comparisons to pdf
pdf('Plots/ModelComparisonTime_Scenarios_lambda.pdf', width = 9, height = 3)
ggplot(subset(ann.est, Year >= 2002 & Parameter == "lambda")) + 
  geom_line(aes(x = Year, y = Median, color = Model, linetype = Environment)) + 
  geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model, linetype = Environment), alpha = 0.05) + 
  geom_vline(aes(xintercept = 2020), color = 'grey60', linetype = 'dotted') + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_fill_manual(values = plotColors, name = 'Harvest') + 
  ylab('Estimate') + 
  theme_bw() + theme(panel.grid = element_blank(),
                     strip.text = element_text(face = "bold", size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 11),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 11))
dev.off()


## Plot comparisons (Ntot, H, S_pup) to pdf
pdf('Plots/ModelComparisonTime_Scenarios_CC0.pdf', width = 9, height = 6)
ggplot(subset(ann.est, Environment == 'Stable')) + 
  geom_line(aes(x = Year, y = Median, color = Model)) + 
  geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 0.1) + 
  geom_vline(aes(xintercept = 2020), color = 'grey60', linetype = 'dotted') + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_fill_manual(values = plotColors, name = 'Harvest') + 
  facet_wrap(~Label, scales = 'free_y', ncol = 1) +
  ylab('Estimate') + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = 'top',
                     strip.text = element_text(face = "bold", size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 11),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 11))

ggplot(subset(ann.est, Environment == 'Stable' & Year >= 2002)) + 
  geom_line(aes(x = Year, y = Median, color = Model)) + 
  geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 0.1) + 
  geom_vline(aes(xintercept = 2020), color = 'grey60', linetype = 'dotted') + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_fill_manual(values = plotColors, name = 'Harvest') + 
  facet_wrap(~Label, scales = 'free_y', ncol = 1) +
  ylab('Estimate') + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = 'top',
                     strip.text = element_text(face = "bold", size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 11),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 11))
dev.off()

pdf('Plots/ModelComparisonTime_Scenarios_CCT.pdf', width = 9, height = 6)
ggplot(subset(ann.est, Environment != 'Stable')) + 
  geom_line(aes(x = Year, y = Median, color = Model)) + 
  geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 0.1) + 
  geom_vline(aes(xintercept = 2020), color = 'grey60', linetype = 'dotted') + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_fill_manual(values = plotColors, name = 'Harvest') + 
  facet_wrap(~Label, scales = 'free_y', ncol = 1) +
  ylab('Estimate') + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = 'top',
                     strip.text = element_text(face = "bold", size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 11),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 11))

ggplot(subset(ann.est, Environment != 'Stable' & Year >= 2002)) + 
  geom_line(aes(x = Year, y = Median, color = Model)) + 
  geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 0.1) + 
  geom_vline(aes(xintercept = 2020), color = 'grey60', linetype = 'dotted') + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_fill_manual(values = plotColors, name = 'Harvest') + 
  facet_wrap(~Label, scales = 'free_y', ncol = 1) +
  ylab('Estimate') + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = 'top',
                     strip.text = element_text(face = "bold", size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 11),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 11))

dev.off()


pdf('Plots/ModelComparisonTime_Scenarios_all.pdf', width = 9, height = 6)
ggplot(subset(ann.est, Year >= 2002)) + 
  geom_line(aes(x = Year, y = Median, color = Model, linetype = Environment)) + 
  geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model, linetype = Environment), alpha = 0.05) + 
  geom_vline(aes(xintercept = 2020), color = 'grey60', linetype = 'dotted') + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_fill_manual(values = plotColors, name = 'Harvest') + 
  facet_wrap(~Label, scales = 'free_y', ncol = 1) +
  ylab('Estimate') + 
  theme_bw() + theme(panel.grid = element_blank(),
                     strip.text = element_text(face = "bold", size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 11),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 11))
dev.off()

###########################################################
# PLOTTING ESTIMATES ACROSS SCENARIOS FOR DIFFERENT YEARS #
###########################################################

## Subset to contain only years of choice
ann.sub <- subset(ann.est, Year %in% c(2020, 2030, 2040, 2050))

## Forestplot for population size
pdf('Plots/CombScenarios_N.pdf', width = 6, height = 3)
ggplot(subset(ann.sub, Parameter == 'Ntot')) + 
  geom_pointrange(aes(x = Model, y = Median, ymin = lCI, ymax = uCI, colour = Model, shape = Environment), size = 0.3, fatten = 4, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_shape_manual(values = c(19,1)) + 
  scale_y_continuous(breaks = seq(100, 1400, by = 200)) + 
  facet_wrap(~Year, nrow = 1) + 
  ylab('N') + 
  ggtitle('Predicted population size (females)') + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'right', panel.spacing.y = unit(1, "lines"))
dev.off()

## Forestplot for pup survival
pdf('Plots/CombScenarios_S_pup.pdf', width = 6, height = 3)
ggplot(subset(ann.sub, Parameter == 'S_pup')) + 
  geom_pointrange(aes(x = Model, y = Median, ymin = lCI, ymax = uCI, colour = Model, shape = Environment), size = 0.3, fatten = 4, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_shape_manual(values = c(19,1)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
  facet_wrap(~Year, nrow = 1) + 
  ylab('S_pup') + 
  ggtitle('Predicted pup survival') + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'right', panel.spacing.y = unit(1, "lines"))
dev.off()

## Forestplot for harvest
pdf('Plots/CombScenarios_H.pdf', width = 6, height = 3)
ggplot(subset(ann.sub, Parameter == 'H')) + 
  geom_pointrange(aes(x = Model, y = Median, ymin = lCI, ymax = uCI, colour = Model, shape = Environment), size = 0.3, fatten = 4, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = plotColors, name = 'Harvest') + 
  scale_shape_manual(values = c(19,1)) + 
  scale_y_continuous(breaks = seq(0, 30, by = 5)) + 
  facet_wrap(~Year, nrow = 1) + 
  ylab('H') + 
  ggtitle('Predicted harvest (females)') + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'right', panel.spacing.y = unit(1, "lines"))
dev.off()