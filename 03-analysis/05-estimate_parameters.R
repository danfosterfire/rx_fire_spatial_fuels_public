
#### setup #####################################################################

library(here)
library(cmdstanr)
library(tidyverse)
library(posterior)
library(bayesplot)


litter_data = readRDS(here::here('02-data',
                                     '02-for_analysis',
                                     'litter_training_data.rds'))

duff_data = readRDS(here::here('02-data',
                               '02-for_analysis',
                               'duff_training_data.rds'))

fwd1h_data = readRDS(here::here('02-data',
                                     '02-for_analysis',
                                     'fwd1h_training_data.rds'))

fwd10h_data = readRDS(here::here('02-data',
                                     '02-for_analysis',
                                     'fwd10h_training_data.rds'))

fwd100h_data = readRDS(here::here('02-data',
                                  '02-for_analysis',
                                  'fwd100h_training_data.rds'))

stan_model.litterduff = cmdstan_model(here::here('03-analysis',
                                      'litterduff_depths_2.stan'))

stan_model.fwd = cmdstan_model(here::here('03-analysis',
                                          'fwd_tallies.stan'))



#### litter ####################################################################

litter_fit = stan_model.litterduff$sample(data = litter_data,
                                   parallel_chains = 4,
                                   iter_warmup = 1000,
                                   iter_sampling = 1000,
                                   output_dir = 
                                     here::here('02-data', '03-results'),
                                   output_basename = 'litter_testing',
                                   seed = 110819,
                                   adapt_delta = 0.80)

litter_fit$cmdstan_diagnose()


litter_fit$save_object(here::here('02-data', '06-results', 'real_fits', 
                                  'litter_fit.rds'))

litter_fit = readRDS(here::here('02-data', '06-results', 'real_fit', 'litter_fit.rds'))


litter_fit$summary(c('beta', 'beta_alpha', 'beta_rho',
                     'sigmaPlot', 'beta_kappa')) %>%
  print(n = Inf)

litter_pairs_1 = 
  mcmc_pairs(litter_fit$draws(),
             pars = c('beta[1]', 'alpha[1]', 'rho[1]', 
                      'sigmaPlot[1]', 'kappa[1]'),
             grid_args = list(top = 'Litter'))

litter_pairs_1

ggsave(litter_pairs_1,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'litter_pairs_1.png'),
       height = 6.5, width = 6.5, units = 'in')

litter_pairs_2 = 
  mcmc_pairs(litter_fit$draws(),
             pars = c('beta[2]', 'alpha[2]', 'rho[2]', 
                      'sigmaPlot[2]', 'kappa[2]'),
             grid_args = list(top = 'Litter'))

litter_pairs_2

ggsave(litter_pairs_2,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'litter_pairs_2.png'),
       height = 6.5, width = 6.5, units = 'in')


litter_pairs_3 = 
  mcmc_pairs(litter_fit$draws(),
             pars = c('beta[3]', 'alpha[3]', 'rho[3]', 
                      'sigmaPlot[3]', 'kappa[3]'),
             grid_args = list(top = 'Litter'))

litter_pairs_3

ggsave(litter_pairs_3,
       filename = 
         here::here('04-communication', 
                    'figures',
                    'manuscript',
                    'litter_pairs_3.png'),
       height = 6.5, width = 6.5, units = 'in')

mcmc_pairs(litter_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]'))

mcmc_dens_overlay(litter_fit$draws(variables = c('alpha[1]', 'rho[1]','beta[1]',
                                                 'alpha[2]', 'rho[2]','beta[2]',
                                                 'alpha[3]', 'rho[3]','beta[3]',
                                                 'sigmaPlot[1]',
                                                 'sigmaPlot[2]', 
                                                 'sigmaPlot[3]', 'kappa[1]',
                                                 'kappa[2]', 'kappa[3]')))





#### duff ######################################################################

duff_fit = stan_model.litterduff$sample(data = duff_data,
                                   parallel_chains = 4,
                                   output_dir = 
                                     here::here('02-data', '03-results'),
                                   output_basename = 'duff',
                                   seed = 110819,
                                   adapt_delta = 0.9)

duff_fit$cmdstan_diagnose()

duff_fit$summary(c('beta', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post',
                         'sigmaPlot', 'kappa')) %>%
  print(n = Inf)

mcmc_pairs(duff_fit$draws(),
           pars = c('beta[1]', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post', 
                    'sigmaPlot', 'kappa'))

mcmc_pairs(duff_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]', 'beta[5]', 'beta[6]'))

mcmc_dens_overlay(duff_fit$draws(variables = c('alpha_pre', 'alpha_post',
                                                     'rho_pre', 'rho_post',
                                                     'sigmaPlot', 'kappa')))

duff_fit$save_object(here::here('02-data', '03-results', 
                                      'duff_fit.rds'))

#### litterduff ################################################################

# NOTE: this section was used for my dissertation, but for the version submitted 
# to IJWF I separated litter and duff into separate models

litterduff_fit = stan_model.litterduff$sample(data = litterduff_data,
                                   parallel_chains = 4,
                                   output_dir = 
                                     here::here('02-data', '03-results'),
                                   output_basename = 'litterduff',
                                   seed = 110819,
                                   adapt_delta = 0.8,
                                   iter_sampling = 10000)

litterduff_fit$save_object(here::here('02-data', '03-results',
                                      'litterduff_fit.rds'))

litterduff_fit = readRDS(here::here('02-data', '03-results', 'litterduff_fit.rds'))

litterduff_fit$cmdstan_diagnose()

litterduff_fit$summary(c('beta', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post',
                         'sigmaPlot', 'kappa')) %>%
  print(n = Inf)

mcmc_pairs(litterduff_fit$draws(),
           pars = c('beta[1]', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post', 
                    'sigmaPlot', 'kappa'))

mcmc_pairs(litterduff_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]', 'beta[5]', 'beta[6]'))

mcmc_dens_overlay(litterduff_fit$draws(variables = c('alpha_pre', 'alpha_post',
                                                     'rho_pre', 'rho_post',
                                                     'sigmaPlot', 'kappa')))




litterduff_df = 
  litterduff_data$X %>%
  as_tibble() %>%
  mutate(Y = litterduff_data$Y,
         plot_id = litterduff_data$plot_id,
         timestep = 
           ifelse(timePost==0,
                  'pre',
                  'post'),
         comp = 
           ifelse(compB==1,
                  'B',
                  ifelse(compC==1,
                         'C',
                         'A'))) %>%
  mutate(timestep = factor(timestep, levels = c('pre', 'post')),
         comp = factor(comp, levels = c('A', 'B', 'C')))


library(brms)

litterduff_fit_weibull = 
  brm(data = litterduff_df,
      Y ~ comp * timestep + (1|plot_id),
      family = hurdle_lognormal(),
      backend = 'cmdstanr',
      cores = 4)

stancode(litterduff_fit_weibull)

weibull_model = cmdstan_model(here::here('03-analysis', 'weibull.stan'))

weibull_fit = 
  weibull_model$sample(data = 
                         list(N = litterduff_data$N,
                              y = litterduff_data$Y+0.01),
                       seed = 110819,
                       parallel_chains = 4,
                       output_dir = here::here('02-data', '03-results'),
                       output_basename = 'weibull',
                       init = list(list('alpha' = 1,
                                   'sigma' = 1),
                                   list('alpha' = 1,
                                   'sigma' = 1),
                                   list('alpha' = 1,
                                   'sigma' = 1),
                                   list('alpha' = 1,
                                   'sigma' = 1)))

weibull_fit$summary()  


#### fwd 1h ####################################################################
fwd1h_fit = stan_model.fwd$sample(data = fwd1h_data,
                                   parallel_chains = 4,
                                  iter_warmup = 1000,
                                  iter_sampling = 9000,
                                   output_dir = 
                                     here::here('02-data', '03-results'),
                                   output_basename = 'fwd1h',
                                   seed = 110819)

fwd1h_fit$save_object(here::here('02-data', '03-results', 
                                      'fwd1h_fit.rds'))

fwd1h_fit = readRDS(here::here('02-data', '03-results', 'fwd1h_fit.rds'))

fwd1h_fit$cmdstan_diagnose()

fwd1h_fit$summary(c('beta', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post',
                    'tau_pre', 'tau_post',
                         'sigmaPlot', 'kappa')) %>%
  print(n = Inf)

mcmc_pairs(fwd1h_fit$draws(),
           pars = c('beta[1]', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post', 
                    'sigmaPlot', 'kappa'))

mcmc_pairs(fwd1h_fit$draws(),
           pars = c('kappa', 'tau_pre', 'tau_post', 'alpha_pre', 'alpha_post'))

mcmc_pairs(fwd1h_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]', 'beta[5]', 'beta[6]'))

mcmc_dens_overlay(fwd1h_fit$draws(variables = c('alpha_pre', 'alpha_post',
                                                     'rho_pre', 'rho_post',
                                                'tau_pre', 'tau_post',
                                                     'sigmaPlot', 'kappa')))

fwd1h_fit$save_object(here::here('02-data', '03-results', 
                                      'fwd1h_fit.rds'))


#### fwd 10h ###################################################################

fwd10h_fit = stan_model.fwd$sample(data = fwd10h_data,
                                   parallel_chains = 4,
                                   iter_warmup = 1000,
                                   iter_sampling = 9000,
                                   output_dir = 
                                     here::here('02-data', '03-results'),
                                   output_basename = 'fwd10h',
                                   seed = 112188,
                                   adapt_delta = 0.9)


fwd10h_fit$save_object(here::here('02-data', '03-results',
                                      'fwd10h_fit.rds'))

fwd10h_fit = readRDS(here::here('02-data', '03-results', 'fwd10h_fit.rds'))

fwd10h_fit$cmdstan_diagnose()


fwd10h_fit$summary(c('beta', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post',
                    'tau_pre', 'tau_post',
                         'sigmaPlot', 'kappa')) %>%
  print(n = Inf)

mcmc_pairs(fwd10h_fit$draws(),
           pars = c('beta[1]', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post', 
                    'sigmaPlot', 'kappa'))

mcmc_pairs(fwd10h_fit$draws(),
           pars = c('kappa', 'tau_pre', 'tau_post', 'alpha_pre', 'alpha_post'))

mcmc_pairs(fwd10h_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]', 'beta[5]', 'beta[6]'))

mcmc_dens_overlay(fwd10h_fit$draws(variables = c('alpha_pre', 'alpha_post',
                                                     'rho_pre', 'rho_post',
                                                'tau_pre', 'tau_post',
                                                     'sigmaPlot', 'kappa')))



#### fwd 100h ###################################################################

fwd100h_fit = stan_model.fwd$sample(data = fwd100h_data,
                                   parallel_chains = 4,
                                   iter_warmup = 1000,
                                   iter_sampling = 9000,
                                   output_dir = 
                                     here::here('02-data', '03-results'),
                                   output_basename = 'fwd100h',
                                   seed = 020190)


fwd100h_fit$save_object(here::here('02-data', '03-results',
                                      'fwd100h_fit.rds'))

fwd100h_fit = readRDS(here::here('02-data', '03-results', 'fwd100h_fit.rds'))


fwd100h_fit$cmdstan_diagnose()


fwd100h_fit$summary(c('beta', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post',
                    'tau_pre', 'tau_post',
                         'sigmaPlot', 'kappa')) %>%
  print(n = Inf)

mcmc_pairs(fwd100h_fit$draws(),
           pars = c('beta[1]', 'alpha_pre', 'alpha_post', 'rho_pre', 'rho_post', 
                    'sigmaPlot', 'kappa'))

mcmc_pairs(fwd100h_fit$draws(),
           pars = c('kappa', 'tau_pre', 'tau_post', 'alpha_pre', 'alpha_post'))

mcmc_pairs(fwd100h_fit$draws(),
           pars = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]', 'beta[5]', 'beta[6]'))

mcmc_dens_overlay(fwd100h_fit$draws(variables = c('alpha_pre', 'alpha_post',
                                                     'rho_pre', 'rho_post',
                                                'tau_pre', 'tau_post',
                                                     'sigmaPlot', 'kappa')))


