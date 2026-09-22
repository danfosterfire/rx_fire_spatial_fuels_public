
#### setup #####################################################################

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)


#### function defs #############################################################
rank_parameters = 
  function(id){
    
    simulated_data = readRDS(here::here('02-data',
                                        '03-results',
                                        'sim_fits',
                                        paste0('simdata_',id,'.rds')))
    
    sim_fit = readRDS(here::here('02-data',
                                 '03-results',
                                 'sim_fits',
                                 paste0('simfit_',id,'.rds')))
    # compare posterior samples against true parameter values
    parameter_ranks = 
      
      # start with the (thinned) posterior samples
      as_draws_df(sim_fit$draws(variables = 
                                  c('beta', 'rho_pre', 'rho_post', 
                                    'alpha_pre', 'alpha_post', 'sigmaPlot', 
                                    'kappa'))) %>%
      
      # create indicator columns comparing posterior samples to true values
      mutate(
        Rbeta01 = as.integer(`beta[1]` < simulated_data$true_params$beta[1]),
        Rbeta02 = as.integer(`beta[2]` < simulated_data$true_params$beta[2]),
        Rbeta03 = as.integer(`beta[3]` < simulated_data$true_params$beta[3]),
        Rbeta04 = as.integer(`beta[4]` < simulated_data$true_params$beta[4]),
        Rbeta05 = as.integer(`beta[5]` < simulated_data$true_params$beta[5]),
        Rbeta06 = as.integer(`beta[6]` < simulated_data$true_params$beta[6]),
        
        Ralphapre = as.integer(alpha_pre < simulated_data$true_params$alpha_pre),
        Ralphapost = as.integer(alpha_post < simulated_data$true_params$alpha_post),
        Rrhopre = as.integer(rho_pre < simulated_data$true_params$rho_pre),
        Rrhopost = as.integer(rho_post < simulated_data$true_params$rho_post),
        
        Rsigmaplot = as.integer(sigmaPlot < simulated_data$true_params$sigmaPlot),
        
        Rkappa = as.integer(kappa < simulated_data$true_params$kappa)) %>%
      
      select(contains(c('Rbeta', 'Ralphapre', 'Ralphapost', 'Rrhopre', 
                        'Rrhopost', 'Rsigmaplot', 'Rkappa'))) %>%
      
      # sum the indicators to get the rank for each variable
      summarise_all(sum)
    
    return(parameter_ranks)
    
  }


#### check model diagnostics ###################################################

model_diagnostics = 
  readRDS(here::here('02-data', '03-results', 'sim_fits', 'model_diagnostics.rds'))

# check these to see which sets of parameters were causing the trouble
bad_batches = 
  c(1:100)[!str_detect(model_diagnostics, pattern = 'Processing complete, no problems detected.')]

# great
bad_batches

#### extract parameter ranks and check uniformity ##############################

batch_param_ranks = 
  do.call('bind_rows',
          lapply(X = c(1:100),
                 FUN = function(i){
                   tryCatch({
                     rank_parameters(i)
                   },
                   error = function(cond){
                     return(NULL)
                   }
                   )
                 }))


batch_param_ranks %>%
  rowid_to_column('run') %>%
  pivot_longer(cols = c(-run), names_to = 'param', values_to = 'rank') %>%
  mutate(rank_sim = sample(1:1000, size = nrow(.), replace = TRUE)) %>%
  mutate(rank_bin = cut(rank, breaks = seq(from = 0, to = 1001, by = 200), include.lowest = TRUE)) %>%
  ggplot(aes(x = rank_bin))+
  geom_bar()+
  geom_hline(yintercept = qbinom(0.005, size = 100, prob = 2/10))+
  geom_hline(yintercept = qbinom(0.995, size = 100, prob = 2/10))+
  facet_wrap(~param, scales = 'free')+
  theme_minimal()

# alpha_pre looks a little questionable but probalby OK, given that I only ran 100
# simulations. beta03 looks iffy as well, but not really a parameter of interest 
# here in any case.
