# Prescribed fires alter the spatial pattern of wildland fuels

## Abstract
We apply Bayesian hierarchical spatial models with Gaussian process spatial 
random effects to describe the spatial pattern of litter and fine woody debris 
both before and prescribed fires in a Sierra Nevada mixed conifer forest. The 
analysis reveals that prescribed fire alters not only mean fuel loads, but also 
the fine-scale spatial pattern of biomass. Prescribed fire increased the relative 
strength of the fine-scale spatial pattern for litter, 1-hour fine woody debris, 
and 10-hour fine woody debris. The burns also increased the length scale of the 
spatial pattern (the distance over which autocorrelation occurs) for litter and 
1-hour fine woody debris. Finally, the Gaussian process noise parameter describing 
very fine-scale autocorrelation increased for 1-hour fuels after the prescribed 
burns. Changes to the fine-scale spatial pattern of litter and fine woody debris 
are likely to impact the behavior of future fires and the ecological function of 
these forests. 


## Short summary for non-experts
This work introduces a novel application of spatial statistics to quantify the 
fine scale (sub-meter to tens of meters) spatial pattern of wildland fuels (the 
fallen litter, sticks, and duff that provide the fuel for wildland fires). The 
method is applied to quantify the impacts of prescribed fire on the spatial 
pattern of these surface fuels. We find that prescribed fire significantly 
alters the spatial pattern of surface fuels, with implications for the fire 
ecology of forests managed with prescribed fire.

## This repo

- **00-R/** contains generic functions for using Brown's transects data to estimate 
fuel loads
- **01-context/** contains some background and planning documents
- **02-data/** contains the raw data and intermediate data products fed to the 
stan models. The actual model results are deposited in 02-data/03-results, but 
are too large to commit to git. See notes below to reproduce them.
- **03-analysis/** contains the core analysis scripts for ingesting the data, 
defining the models, estimating the parameters, and analyzing the results.
  - R scripts are numbered according to their order in the data pipeline
  - xx- and zz- files are work branches that did not wind up in the final paper
  - Stan models (both stan code and compiled executables) are labeled accordingly. 
  See `03-analysis/05-estimate_parameters.R` for code pointing to the specific 
  stan code used for the paper.
- **04-communication/** contains figures, tables, and the manuscript.

## Reproducing fitted model results

The stan model objects and draws .csvs are produced by `03-analysis/05-estimate_parameters.R`, 
which loads the data objects and runs the stan models to estimate parameters. 
These stan model objects and the draws .csvs are too large to commit to github.
Scripts downstream of `03-analysis/05-estimate_parameters.R` require the stan 
model objects. To reproduce the figures and tables included in the manuscript, 
you will need to re-run `03-analysis/05-estimate_parameters.R` to produce those 
files locally. 
