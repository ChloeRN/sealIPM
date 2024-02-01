# sealIPM

## About this repository

This repository contains code for running an integrated population model (IPM) & harvest sustainability analysis for ringed seals on Svalbard. The project directory contains primary code for a simple PBR (= Potential Biological Removal) analysis, a more advanced PBR that uses outputs from the IPM, the IPM itself including PVA ( = population viability analysis) extension, as well as for visualizing different model aspects and results. Additionally, folders within the project directory contain further code for models that were part of assessing suitability of informative survival priors.

The files constituting the different parts of the analyses are described below.

## File contents

The following table gives an overview over each script in the repository, including a brief description of its contents/function.

| Script                                                      | Description                                                                                                                                                                                                                                              |
|-------------------------|-----------------------------------------------|
| RudimentaryPBR_Exploration.R                                | A simple PBR analysis employing traditionally used, fixed values.                                                                                                                                                                                        |
| Pup_S&SeaIce_Plots.R                                        | Visualization of the relationship between relative sea ice availability and pup survival as implemented in the IPM.                                                                                                                                      |
| Seal_IPM_fSAD&eHAD_ice.R                                    | Workflow for formatting data for and implementing and running the IPM (including PVA scenarios) using NIMBLE. A range of "toggles" at the start of the script allow to control what kind of scenario (duration, harvest change, sea ice change) is run.  |
| Seal_IPM_InitValSim_ice.R                                   | Function for simulating inital value for the IPM. Dependency of Seal_IPM_fSAD&eHAD_ice.R.                                                                                                                                                                |
| IPM_PBR_calculation.R                                       | A PBR analysis using estimates from the IPM as input.                                                                                                                                                                                                    |
| ModelComparison_Scenarios.R                                 | Graphical comparison of IPMs run under different PVA scenarios.                                                                                                                                                                                          |
| SurvivalPriors_Comparison.R                                 | Comparison of informative survival priors used in the model vs. equivalent estimates obtained from traditional life table and catch curve analyses, as well as the Hoening Model.                                                                        |
| LifeTable_Analyses/LifeTableAnalysis_Bayesian.R             | Bayesian implementation of a classical life table analysis for the ringed seal data.                                                                                                                                                                     |
| CatchCurve_Analyses/CatchCurveAnalysis_Bayesian.R           | Bayesian implementation of a classical catch curve analysis for the ringed seal data (no age truncation).                                                                                                                                                |
| CatchCurve_Analyses/CatchCurveAnalysis_Truncated_Bayesian.R | Bayesian implementation of a classical catch curve analysis for the ringed seal data (including age truncation).                                                                                                                                         |
| HoeningModel_Sim/HoeningMod_PriorPred.R                     | Simulation of prior distributions for natural morality using the cross-species Hoening Model developed by Tom Porteus.                                                                                                                                   |

## Underlying data

The analyses require a set of input data (RDS files) for running. These will be made publicly available through OSF once the manuscript associated with this work has been accepted for publication. In the meantime, the input data can be requested directly by emailing me: [chloe.nater\@nina.no](mailto:chloe.nater@nina.no){.email}

The data files underlying the main analyses are:

-   *SealIPM_DemoData.rds:* a tibble containing demographic data (age class, maturation status) for female seals harvested in Svalbard.

-   *mO_HoeningParameters_Seal.rds:* a list containing log mean and log sd's for priors for natural mortality hazard rates calculated using the Hoening model by Tom Porteus.

-   *HarvestCountData.rds:* a list containing age-at-harvest and harvest count data for female ringed seals from Svalbard.

-   *SealIPM_IceData.rds:* a list containing covariate sea ice data that can be used as a predictor of (relative) seal pup mortality in the IPM.

The auxiliary analyses (life table analysis, catch curve analysis, and Hoening model simulation) also require on input file each (*LifeTables.rds*, *AaHData.rds*, and *Hoenig_Posteriors_NoTitle.txt*, respectively).

## Contact

For questions, inquiries, etc. contact: [chloe.nater\@nina.no](mailto:chloe.nater@nina.no){.email}
