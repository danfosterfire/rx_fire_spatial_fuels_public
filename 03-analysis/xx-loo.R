library(here)
library(cmdstanr)
library(loo)
library(bayesplot)


litterduff_fit = 
  readRDS(here::here('02-data', '03-results', 'litterduff_fit.rds'))

loo_output = 
  litterduff_fit$loo()

loo_output

plot(loo_output)
