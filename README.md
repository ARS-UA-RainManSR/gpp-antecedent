# Antecedent modeling of tent GPP at RainMan

This repo contains scripts to run the Stochastic Antecedent Model ([Ogle et al. 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12399)) on tent gas exchange data from 2021-2022 at the Rainfall Manipulation site on the Santa Rita Experimental Range, AZ. 

## Installation instructions
This repo uses `renv` to reproduce the R environment. Type `renv::restore()` into the console to install the correct versions of each package in your session. 

The flux and environmental data are stored in a `data/` folder in the root directory. Please contact [Fangyue Zhang](mailto:fangyuezhang@arizona.edu) for the raw data files. 

In addition to R and RStudio, this repo requires the installation of [JAGS 4.3.1](https://sourceforge.net/projects/mcmc-jags/files/) on your computer. 

## Model formulations
Models are saved in the `scripts/` directory. 
 - `model-v1` contains a simple version with only concurrent PAR and antecedent daytime mean VPD for 10 days at 1 day intervals. No VWC or random effects are included yet. 
 - `model-v2` explores the addition of VWC at two depths, which are indexed by treatment. Shallow VWC seems to require a longer antecedent period. Which two-way interactions are the most interesting? Random effects are not included yet. 
 - `model-v3` separate out VWC antecedent effects by treatment and slims down the number of interaction terms. Random effects are not included yet. 
 - `model-v4` adds random effects with post sweeping to improve identifiability with the intercept terms. Includes code for replicated data and model evaluation