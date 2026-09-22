
#### setup #####################################################################

library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(foreach)
library(doParallel)

#### function definitions ######################################################
simulate_data = 
  function(id){
    
    # load in the observed X data
    obs_data = readRDS(here::here('02-data',
                                  '02-for_analysis',
                                  'litterduff_data.rds'))
    
    # pull parameter values from the prior distribution
    # NOTE: because we're interested in the GP parameters for this study,
    # not as concerned with doing SBC for the fixed effects and plot SD; 
    # just setting these rather than drawing them from the prior facilitates 
    # creating realistic datasets to do SBC with
    #beta = rnorm(n = 6, mean = 0, sd = 2)
    #sigmaPlot = truncnorm::rtruncnorm(n = 1, a = 0, mean = 0, sd =0.5)
    beta = c(2.5, 0.5, -0.5, -1, 0.25, -0.25)
    sigmaPlot = 0.25
    kappa = rcauchy(n = 1, location = 0, scale = 5)
    while (kappa <= 0){
      kappa = rcauchy(n = 1, location = 0, scale = 5)
    }
    alpha_pre = truncnorm::rtruncnorm(n = 1, a = 0, mean = 0, sd = 0.5)
    alpha_post = truncnorm::rtruncnorm(n = 1, a = 0, mean = 0, sd = 0.5)
    rho_pre = invgamma::rinvgamma(n = 1, shape = 5, rate = 40)
    rho_post = invgamma::rinvgamma(n = 1, shape = 5, rate = 40)
    
    # fixed effects
    XB = obs_data$X %*% beta
  
    # "pre" GP realization
    D = dist(obs_data$coords) %>% as.matrix()
    Sigma.pre = 
      alpha_pre**2 * exp(-(1/(2*(rho_pre**2)))*(D**2))+
      diag(nrow = nrow(D), ncol = ncol(D), x = 1e-9)
    KL.pre = t(chol(Sigma.pre))
    zGP.pre = 
      matrix(nrow = obs_data$L, ncol = obs_data$P, byrow = FALSE,
             data = sapply(X = 1:obs_data$P,
                           FUN = function(plot){
                             rnorm(n = obs_data$L, mean = 0, sd = 1)
                           }))
    GP.pre = matrix(nrow = obs_data$L, ncol = obs_data$P, byrow = FALSE,
                    data = sapply(X = 1:obs_data$P,
                                  FUN = function(plot){
                                    KL.pre %*% zGP.pre[,plot]
                                  }))
    
    # "post" GP realization
    Sigma.post = 
      alpha_post**2 * exp(-(1/(2*(rho_post**2)))*(D**2))+
      diag(nrow = nrow(D), ncol = ncol(D), x = 1e-9)
    KL.post = t(chol(Sigma.post))
    zGP.post = 
      matrix(nrow = obs_data$L, ncol = obs_data$P, byrow = FALSE,
             data = sapply(X = 1:obs_data$P,
                           FUN = function(plot){
                             rnorm(n = obs_data$L, mean = 0, sd = 1)
                           }))
    GP.post = matrix(nrow = obs_data$L, ncol = obs_data$P, byrow = FALSE,
                    data = sapply(X = 1:obs_data$P,
                                  FUN = function(plot){
                                    KL.post %*% zGP.post[,plot]
                                  }))
    
    # plot random effect realization
    plotEffects = rnorm(n = obs_data$P, mean = 0, sd = sigmaPlot)
    
    # linear predictor for each observation
    logmu = 
      sapply(X = 1:obs_data$N,
             FUN = function(i){
               XB[i] +
                 GP.pre[obs_data$location_id[i],obs_data$plot_id[i]] * obs_data$pre[i] +
                 GP.post[obs_data$location_id[i], obs_data$plot_id[i]] * obs_data$post[i] +
                 plotEffects[obs_data$plot_id[i]]
             })
    
    # draw responses
    Y = rnbinom(n = obs_data$N, mu = exp(logmu), size = kappa)
    
    sim_data = obs_data
    sim_data$Y = Y
    
    # save resulting data
    results = 
      list('id' = id,
           'true_params' = list('beta' = beta,
                                'alpha_pre' = alpha_pre,
                                'alpha_post' = alpha_post,
                                'rho_pre' = rho_pre,
                                'rho_post' = rho_post,
                                'sigmaPlot' = sigmaPlot,
                                'kappa' = kappa),
           'sim_data' = sim_data)
    
    saveRDS(results,
            here::here('02-data',
                       '03-results',
                       'sim_fits',
                       paste0('simdata_', id, '.rds')))
  }


estimate_parameters = 
  function(id){
    library(here)
    library(cmdstanr)
    
    # load the data and the stan model
    stan_model = 
      cmdstanr::cmdstan_model(here::here('03-analysis', 'fuel_loading.stan'))
    
    # fit the model
    sim_fit = 
      stan_model$sample(
        data = readRDS(here::here('02-data', '03-results', 'sim_fits', 
                         paste0('simdata_',id,'.rds')))$sim_data,
        parallel_chains = 1,
        output_dir = here::here('02-data', '03-results', 'sim_fits'),
        output_basename = paste0('simfit_', id),
        thin = 4,
        iter_warmup = 1000,
        iter_sampling = 1000)
    
    sim_fit$save_object(here::here('02-data', '03-results', 'sim_fits',
                                   paste0('simfit_',id,'.rds')))
    
    print(paste0('successful run ', id))
    return(paste0('Successful run ', id))
  }

#### simulate data and prior predictive checks #################################

set.seed(110819)
simulated_data = 
  lapply(X = 1:100,
         FUN = function(i){
           simulate_data(i)
         })


sim_data = 
  do.call('bind_rows',
          lapply(X = 1:100,
                 FUN = function(id){
                   
                   sim_data = readRDS(here::here('02-data',
                                                 '03-results',
                                                 'sim_fits',
                                                 paste0('simdata_',id,'.rds')))
                   
                   results = 
                     sim_data$sim_data$X %>%as_tibble()
                   results$Y = sim_data$sim_data$Y
                   results$beta1 = sim_data$true_params$beta[1]
                   results$beta2 = sim_data$true_params$beta[2]
                   results$beta3 = sim_data$true_params$beta[3]
                   results$beta4 = sim_data$true_params$beta[4]
                   results$beta5 = sim_data$true_params$beta[5]
                   results$beta6 = sim_data$true_params$beta[6]
                   results$alpha_pre = sim_data$true_params$alpha_pre
                   results$alpha_post = sim_data$true_params$alpha_post
                   results$rho_pre = sim_data$true_params$rho_pre
                   results$rho_post = sim_data$true_params$rho_post
                   results$sigmaPlot = sim_data$true_params$sigmaPlot
                   results$kappa = sim_data$true_params$kappa
                   results$id = id
                   results$comp = 
                     ifelse(results$compB, 'B',
                            ifelse(results$compC, 'C', 'A'))
                   results$timestep = ifelse(results$timePost, 'post', 'pre')
                   return(results)
                 })
          )


ggplot(data = sim_data,
       aes(x = Y))+
  geom_histogram()+
  facet_grid(timestep~comp)+
  scale_x_continuous(limits = c(-2, 45))

ggplot(data = sim_data %>%
         mutate(kappa_bin = cut(kappa, 
                                breaks = quantile(kappa, c(0, 0.25, 0.5, 0.75, 1)),
                                include.lowest = TRUE)),
       aes(x = Y))+
  geom_histogram()+
  facet_wrap(~kappa_bin, scales = 'free')+
  scale_x_continuous(limits = c(-2, 25))

ggplot(data = sim_data %>%
         mutate(beta1_bin = cut(beta1, 
                                breaks = quantile(beta1, c(0, 0.25, 0.5, 0.75, 1)),
                                include.lowest = TRUE),
                kappa_bin = cut(kappa, 
                                breaks = quantile(kappa, c(0, 0.25, 0.5, 0.75, 1)),
                                include.lowest = TRUE)) %>%
         filter(comp=='A'&timestep=='pre'),
       aes(x = Y))+
  geom_histogram()+
  facet_grid(kappa_bin~beta1_bin, scales = 'free')+
  scale_x_continuous(limits = c(-2, 40))




#### estimate parameters from simulated data ##################################

registerDoParallel()


foreach(i = 1:100) %dopar% {
  estimate_parameters(i)
}

stopImplicitCluster()

#### extract model diagnostics #################################################

# only works properly on savio cluster, some filepaths in the fitted model 
# object are full and point to locations on the cluster; could maybe reach 
# into the objects and modify the filepaths but simpler to just run this 
# part on the cluster and save the results
model_diagnostics = 
  sapply(X = c(1:100),
         FUN = function(i){
           output = 
             tryCatch(
               {
                 sim_fit = readRDS(here::here('02-data', '03-results', 'sim_fits', 
                                              paste0('simfit_', i, '.rds')))
                 
                 diagnostics = sim_fit$cmdstan_diagnose()$stdout
                 
                 return(diagnostics)
               },
               error = function(cond){
                 message(paste0('Sim ', i,' was not completed'))
                 return(paste0('sim ',i, ' not completed'))
               }
             )
         })

saveRDS(model_diagnostics,
        here::here('02-data', '03-results', 'sim_fits', 'model_diagnostics.rds'))



