
library(here)
library(tidyverse)
library(posterior)
library(cmdstanr)

set.seed(110819)

#### litter ####################################################################

litter_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'litter_validation_data.rds'))

litter_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'litter_fit.rds'))

litter_posterior = 
  as_draws_df(litter_fit$draws())

litter_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw ', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     litter_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha_pre = litter_posterior$alpha_pre[draw]
                   alpha_post = litter_posterior$alpha_post[draw]
                   
                   rho_pre = litter_posterior$rho_pre[draw]
                   rho_post = litter_posterior$rho_post[draw]
                   
                   plotEffects = 
                     litter_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP_pre = 
                     matrix(nrow = litter_data$L, ncol = litter_data$P,
                            byrow = FALSE, 
                            data = 
                              litter_posterior %>%
                              select(contains('zGP_pre')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   zGP_post = 
                     matrix(nrow = litter_data$L, ncol = litter_data$P,
                            byrow = FALSE, 
                            data = 
                              litter_posterior %>%
                              select(contains('zGP_post')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   kappa = litter_posterior$kappa[draw]
                   
                   dists = 
                     litter_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   K_pre = alpha_pre * exp(-(dists**2)/(2*(rho_pre**2)))
                   K_pre = K_pre + diag(nrow = nrow(K_pre),
                                        ncol = ncol(K_pre),
                                        x = 1e-9)
                   KL_pre = t(chol(K_pre))
                   
                   GP_pre = 
                     matrix(nrow = litter_data$L,
                            ncol = litter_data$P,
                            data = 
                              sapply(X = 1:litter_data$P,
                                     FUN = function(p){
                                       KL_pre %*% zGP_pre[,p]
                                     }))
                   
                   K_post = alpha_post * exp(-(dists**2)/(2*(rho_post**2)))
                   K_post = K_post + diag(nrow = nrow(K_post),
                                          ncol = ncol(K_post),
                                          x = 1e-9)
                   KL_post = t(chol(K_post))
                   
                   GP_post = 
                     matrix(nrow = litter_data$L,
                            ncol = litter_data$P,
                            data = 
                              sapply(X = 1:litter_data$P,
                                     FUN = function(p){
                                       KL_post %*% zGP_post[,p]
                                     }))
                   
                   logmu = 
                     sapply(X = 1:litter_data$N,
                            FUN = function(i){
                              as.numeric(litter_data$X[i,] %*% beta)+
                                plotEffects[litter_data$plot_id[i]]+
                                (GP_pre[litter_data$location_id[i],litter_data$plot_id[i]] *
                                   litter_data$pre[i]) +
                                (GP_post[litter_data$location_id[i],
                                         litter_data$plot_id[i]] *
                                   litter_data$post[i])
                            })
                   
                   Y = rnbinom(n = litter_data$N,
                               size = kappa,
                               mu = exp(logmu))
                   
                   results = 
                     litter_data$X %>% as_tibble()
                   results$draw = draw
                   results$Y_sim = Y
                   results$logmu = logmu
                   results$Y_obs = litter_data$Y
                   results$sample = 1:litter_data$N
                   
                   return(results)
                 }))



litter_mae.validation = 
  litter_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(mean_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-mean_sim)) %>%
  pull(ae) %>%
  mean()

litter_mean.validation = 
  mean(litter_data$Y)

litter_wMAPE.validation = litter_mae.validation/litter_mean.validation


litter_prediction1 = 
  litter_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(y.50 = median(Y_sim),
            y.mean = mean(Y_sim),
            y.05 = quantile(Y_sim, probs = 0.05),
            y.95 = quantile(Y_sim, probs = 0.95),
            logmu = mean(logmu)) %>%
  ungroup() %>% 
  ggplot()+
  geom_ribbon(aes(x = logmu, ymin = y.05, ymax = y.95), alpha = 0.5, color = 'grey')+
  geom_point(aes(x =logmu, y = Y_obs), size = 1)+
  theme_minimal()+
  labs(y = 'Observed depth (cm)', x = 'Mean log predicted depth (log cm)')
litter_prediction1

litter_prediction3 = 
  litter_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(y.50 = median(Y_sim),
            y.mean = mean(Y_sim),
            y.05 = quantile(Y_sim, probs = 0.05),
            y.95 = quantile(Y_sim, probs = 0.95),
            logmu = mean(logmu)) %>%
  ungroup() %>%
  ggplot(aes(x = y.mean, y = Y_obs))+
  geom_point(size = 0)+
  geom_abline(intercept= 0, slope = 1, color = 'red')+
  geom_smooth(method = 'lm')+
  coord_fixed()+
  theme_minimal()+
  labs(x = 'Mean simulated depth (cm)', y = 'Observed depth (cm)')

litter_prediction3

litter_prediction2 = 
  litter_sim %>%
  pivot_longer(cols = c('Y_sim', 'Y_obs'),
               names_to = 'source', values_to = 'value') %>%
  mutate(timestep = ifelse(timePost==1,'post','pre'),
         comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
  ggplot(aes(x = value, fill = source))+
  geom_bar(position = position_dodge())+
  scale_x_continuous(limits = c(-2, 40))+
  facet_grid(timestep~comp, scale = 'free_x')+
  theme_minimal()

litter_prediction2

litter_prediction4 = 
  litter_sim %>%
  pivot_longer(cols = c('Y_sim', 'Y_obs'),
               names_to = 'source', values_to = 'value') %>%
  mutate(timestep = ifelse(timePost==1,'post','pre'),
         comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
  group_by(timestep, comp, value, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    litter_sim %>%
      pivot_longer(cols = c('Y_sim', 'Y_obs'),
                   names_to = 'source', values_to = 'value') %>%
      
      mutate(timestep = ifelse(timePost==1,'post','pre'),
             comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
      group_by(timestep, comp, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  ggplot(aes(x = value, y = proportion, fill = source))+
  geom_col(position = position_dodge())+
  scale_x_continuous(limits = c(-2, 40))+
  facet_grid(timestep~comp, scale = 'free_x')+
  theme_minimal()

litter_prediction4  


litter_prediction5 = 
  litter_sim %>%
  filter(draw == 1) %>%
  group_by(draw, Y_obs) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  ggplot()+
  geom_col(aes(x = Y_obs, y = count), 
           fill = 'red', col = NA, alpha = 0.5)+
  theme_minimal()+
  geom_errorbar(
    data = 
      litter_sim %>%
      group_by(draw, Y_sim) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      complete(draw,
               Y_sim = 1:max(.$Y_sim)) %>%
      mutate(count = ifelse(is.na(count), 0, count)) %>%
      group_by(Y_sim) %>%
      summarise(count.025 = quantile(count, 0.025),
                count.975 = quantile(count, 0.975),
                count.min = min(count),
                count.max = max(count)) %>%
      ungroup() %>%
      mutate(source = 'sim'),
    aes(x = Y_sim, ymin = count.025, ymax = count.975),
    width = 1
  )+
  labs(x = 'Litter depth (cm)', y = 'N observations')

litter_prediction5

ggsave(litter_prediction1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction1.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(litter_prediction2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction2.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(litter_prediction3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction3.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(litter_prediction4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction4.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(litter_prediction5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_prediction5.png'),
       height = 4, width = 6.5, units = 'in')


rm(litter_sim)



#### fwd 1h ####################################################################
fwd1h_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'fwd1h_validation_data.rds'))

fwd1h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'fwd1h_fit.rds'))

fwd1h_posterior = 
  as_draws_df(fwd1h_fit$draws())

fwd1h_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     fwd1h_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha_pre = fwd1h_posterior$alpha_pre[draw]
                   alpha_post = fwd1h_posterior$alpha_post[draw]
                   
                   rho_pre = fwd1h_posterior$rho_pre[draw]
                   rho_post = fwd1h_posterior$rho_post[draw]
                   
                   tau_pre = fwd1h_posterior$tau_pre[draw]
                   tau_post = fwd1h_posterior$tau_post[draw]
                   
                   plotEffects = 
                     fwd1h_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP_pre = 
                     matrix(nrow = fwd1h_data$L, ncol = fwd1h_data$P,
                            byrow = FALSE, 
                            data = 
                              fwd1h_posterior %>%
                              select(contains('zGP_pre')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   zGP_post = 
                     matrix(nrow = fwd1h_data$L, ncol = fwd1h_data$P,
                            byrow = FALSE, 
                            data = 
                              fwd1h_posterior %>%
                              select(contains('zGP_post')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   kappa = fwd1h_posterior$kappa[draw]
                   
                   dists = 
                     fwd1h_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   K_pre = alpha_pre * exp(-(dists**2)/(2*(rho_pre**2)))
                   K_pre = K_pre + diag(nrow = nrow(K_pre),
                                        ncol = ncol(K_pre),
                                        x = tau_pre)
                   KL_pre = t(chol(K_pre))
                   
                   GP_pre = 
                     matrix(nrow = fwd1h_data$L,
                            ncol = fwd1h_data$P,
                            data = 
                              sapply(X = 1:fwd1h_data$P,
                                     FUN = function(p){
                                       KL_pre %*% zGP_pre[,p]
                                     }))
                   
                   K_post = alpha_post * exp(-(dists**2)/(2*(rho_post**2)))
                   K_post = K_post + diag(nrow = nrow(K_post),
                                          ncol = ncol(K_post),
                                          x = tau_post)
                   KL_post = t(chol(K_post))
                   
                   GP_post = 
                     matrix(nrow = fwd1h_data$L,
                            ncol = fwd1h_data$P,
                            data = 
                              sapply(X = 1:fwd1h_data$P,
                                     FUN = function(p){
                                       KL_post %*% zGP_post[,p]
                                     }))
                   
                   logmu = 
                     sapply(X = 1:fwd1h_data$N,
                            FUN = function(i){
                              as.numeric(fwd1h_data$X[i,] %*% beta)+
                                plotEffects[fwd1h_data$plot_id[i]]+
                                (GP_pre[fwd1h_data$location_id[i],fwd1h_data$plot_id[i]] *
                                   fwd1h_data$pre[i]) +
                                (GP_post[fwd1h_data$location_id[i],
                                         fwd1h_data$plot_id[i]] *
                                   fwd1h_data$post[i])
                            })
                   
                   Y = rnbinom(n = fwd1h_data$N,
                               size = kappa,
                               mu = exp(logmu))
                   
                   results = 
                     fwd1h_data$X %>% as_tibble()
                   results$draw = draw
                   results$Y_sim = Y
                   results$logmu = logmu
                   results$Y_obs = fwd1h_data$Y
                   results$sample = 1:fwd1h_data$N
                   
                   return(results)
                 }))


fwd1h_mae.validation = 
  fwd1h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(mean_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-mean_sim)) %>%
  pull(ae) %>%
  mean()

fwd1h_mean.validation = 
  mean(fwd1h_data$Y)

fwd1h_wMAPE.validation = fwd1h_mae.validation/fwd1h_mean.validation


fwd1h_prediction1 = 
  fwd1h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(y.50 = median(Y_sim),
            y.mean = mean(Y_sim),
            y.05 = quantile(Y_sim, probs = 0.05),
            y.95 = quantile(Y_sim, probs = 0.95),
            logmu = mean(logmu)) %>%
  ungroup() %>% 
  ggplot()+
  geom_ribbon(aes(x = logmu, ymin = y.05, ymax = y.95), alpha = 0.5, color = 'grey')+
  geom_point(aes(x =logmu, y = Y_obs), size = 1)+
  theme_minimal()
fwd1h_prediction1

fwd1h_prediction3 = 
  fwd1h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(y.50 = median(Y_sim),
            y.mean = mean(Y_sim),
            y.05 = quantile(Y_sim, probs = 0.05),
            y.95 = quantile(Y_sim, probs = 0.95),
            logmu = mean(logmu)) %>%
  ungroup() %>%
  ggplot(aes(x = y.mean, y = Y_obs))+
  geom_point(size = 0)+
  geom_abline(intercept= 0, slope = 1, color = 'red')+
  geom_smooth(method = 'lm')+
  coord_fixed()+
  theme_minimal()+
  labs(x = 'Mean simulated tally', y = 'Observed tally')

fwd1h_prediction3

fwd1h_prediction2 = 
  fwd1h_sim %>%
  pivot_longer(cols = c('Y_sim', 'Y_obs'),
               names_to = 'source', values_to = 'value') %>%
  mutate(timestep = ifelse(timePost==1,'post','pre'),
         comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
  ggplot(aes(x = value, fill = source))+
  geom_bar(position = position_dodge())+
  scale_x_continuous(limits = c(-2, 50))+
  facet_grid(timestep~comp)+
  theme_minimal()

fwd1h_prediction2

fwd1h_prediction4 = 
  fwd1h_sim %>%
  pivot_longer(cols = c('Y_sim', 'Y_obs'),
               names_to = 'source', values_to = 'value') %>%
  mutate(timestep = ifelse(timePost==1,'post','pre'),
         comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
  group_by(timestep, comp, value, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    fwd1h_sim %>%
      pivot_longer(cols = c('Y_sim', 'Y_obs'),
                   names_to = 'source', values_to = 'value') %>%
      
      mutate(timestep = ifelse(timePost==1,'post','pre'),
             comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
      group_by(timestep, comp, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  ggplot(aes(x = value, y = proportion, fill = source))+
  geom_col(position = position_dodge())+
  scale_x_continuous(limits = c(-2, 50))+
  facet_grid(timestep~comp, scale = 'free_x')+
  theme_minimal()

fwd1h_prediction4

ggsave(fwd1h_prediction1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd1h_prediction1.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd1h_prediction2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd1h_prediction2.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd1h_prediction3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd1h_prediction3.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd1h_prediction4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd1h_predictions4.png'))

fwd1h_prediction5 = 
  fwd1h_sim %>%
  filter(draw == 1) %>%
  group_by(draw, Y_obs) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  ggplot()+
  geom_col(aes(x = Y_obs, y = count), 
           fill = 'red', col = NA, alpha = 0.5)+
  theme_minimal()+
  geom_errorbar(
    data = 
      fwd1h_sim %>%
      group_by(draw, Y_sim) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      complete(draw,
               Y_sim = 1:max(.$Y_sim)) %>%
      mutate(count = ifelse(is.na(count), 0, count)) %>%
      group_by(Y_sim) %>%
      summarise(count.025 = quantile(count, 0.025),
                count.975 = quantile(count, 0.975),
                count.min = min(count),
                count.max = max(count)) %>%
      ungroup() %>%
      mutate(source = 'sim'),
    aes(x = Y_sim, ymin = count.025, ymax = count.975),
    width = 1
  )+
  labs(x = '1-hour FWD tally', y = 'N observations')

ggsave(fwd1h_prediction5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd1h_prediction5.png'),
       height = 4, width = 6.5, units = 'in')


rm(fwd1h_sim)

#### fwd 10h ###################################################################

fwd10h_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'fwd10h_validation_data.rds'))

fwd10h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'fwd10h_fit.rds'))

fwd10h_posterior = 
  as_draws_df(fwd10h_fit$draws())

fwd10h_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     fwd10h_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha_pre = fwd10h_posterior$alpha_pre[draw]
                   alpha_post = fwd10h_posterior$alpha_post[draw]
                   
                   rho_pre = fwd10h_posterior$rho_pre[draw]
                   rho_post = fwd10h_posterior$rho_post[draw]
                   
                   tau_pre = fwd10h_posterior$tau_pre[draw]
                   tau_post = fwd10h_posterior$tau_post[draw]
                   
                   plotEffects = 
                     fwd10h_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP_pre = 
                     matrix(nrow = fwd10h_data$L, ncol = fwd10h_data$P,
                            byrow = FALSE, 
                            data = 
                              fwd10h_posterior %>%
                              select(contains('zGP_pre')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   zGP_post = 
                     matrix(nrow = fwd10h_data$L, ncol = fwd10h_data$P,
                            byrow = FALSE, 
                            data = 
                              fwd10h_posterior %>%
                              select(contains('zGP_post')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   kappa = fwd10h_posterior$kappa[draw]
                   
                   dists = 
                     fwd10h_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   K_pre = alpha_pre * exp(-(dists**2)/(2*(rho_pre**2)))
                   K_pre = K_pre + diag(nrow = nrow(K_pre),
                                        ncol = ncol(K_pre),
                                        x = tau_pre)
                   KL_pre = t(chol(K_pre))
                   
                   GP_pre = 
                     matrix(nrow = fwd10h_data$L,
                            ncol = fwd10h_data$P,
                            data = 
                              sapply(X = 1:fwd10h_data$P,
                                     FUN = function(p){
                                       KL_pre %*% zGP_pre[,p]
                                     }))
                   
                   K_post = alpha_post * exp(-(dists**2)/(2*(rho_post**2)))
                   K_post = K_post + diag(nrow = nrow(K_post),
                                          ncol = ncol(K_post),
                                          x = tau_post)
                   KL_post = t(chol(K_post))
                   
                   GP_post = 
                     matrix(nrow = fwd10h_data$L,
                            ncol = fwd10h_data$P,
                            data = 
                              sapply(X = 1:fwd10h_data$P,
                                     FUN = function(p){
                                       KL_post %*% zGP_post[,p]
                                     }))
                   
                   logmu = 
                     sapply(X = 1:fwd10h_data$N,
                            FUN = function(i){
                              as.numeric(fwd10h_data$X[i,] %*% beta)+
                                plotEffects[fwd10h_data$plot_id[i]]+
                                (GP_pre[fwd10h_data$location_id[i],fwd10h_data$plot_id[i]] *
                                   fwd10h_data$pre[i]) +
                                (GP_post[fwd10h_data$location_id[i],
                                         fwd10h_data$plot_id[i]] *
                                   fwd10h_data$post[i])
                            })
                   
                   Y = rnbinom(n = fwd10h_data$N,
                               size = kappa,
                               mu = exp(logmu))
                   
                   results = 
                     fwd10h_data$X %>% as_tibble()
                   results$draw = draw
                   results$Y_sim = Y
                   results$logmu = logmu
                   results$Y_obs = fwd10h_data$Y
                   results$sample = 1:fwd10h_data$N
                   
                   return(results)
                 }))


fwd10h_mae.validation = 
  fwd10h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(mean_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-mean_sim)) %>%
  pull(ae) %>%
  mean()

fwd10h_mean.validation = 
  mean(fwd10h_data$Y)

fwd10h_wMAPE.validation = fwd10h_mae.validation/fwd10h_mean.validation


fwd10h_prediction1 = 
  fwd10h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(y.50 = median(Y_sim),
            y.mean = mean(Y_sim),
            y.05 = quantile(Y_sim, probs = 0.05),
            y.95 = quantile(Y_sim, probs = 0.95),
            logmu = mean(logmu)) %>%
  ungroup() %>% 
  ggplot()+
  geom_ribbon(aes(x = logmu, ymin = y.05, ymax = y.95), alpha = 0.5, color = 'grey')+
  geom_point(aes(x =logmu, y = Y_obs), size = 1)+
  theme_minimal()
fwd10h_prediction1

fwd10h_prediction3 = 
  fwd10h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(y.50 = median(Y_sim),
            y.mean = mean(Y_sim),
            y.05 = quantile(Y_sim, probs = 0.05),
            y.95 = quantile(Y_sim, probs = 0.95),
            logmu = mean(logmu)) %>%
  ungroup() %>%
  ggplot(aes(x = y.mean, y = Y_obs))+
  geom_point(size = 0)+
  geom_abline(intercept= 0, slope = 1, color = 'red')+
  geom_smooth(method = 'lm')+
  coord_fixed()+
  theme_minimal()+
  labs(x = 'Mean simulated tally', y = 'Observed tally')

fwd10h_prediction3

fwd10h_prediction2 = 
  fwd10h_sim %>%
  pivot_longer(cols = c('Y_sim', 'Y_obs'),
               names_to = 'source', values_to = 'value') %>%
  mutate(timestep = ifelse(timePost==1,'post','pre'),
         comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
  ggplot(aes(x = value, fill = source))+
  geom_bar(position = position_dodge())+
  scale_x_continuous(limits = c(-2, 50))+
  facet_grid(timestep~comp)+
  theme_minimal()

fwd10h_prediction2

fwd10h_prediction4 = 
  fwd10h_sim %>%
  pivot_longer(cols = c('Y_sim', 'Y_obs'),
               names_to = 'source', values_to = 'value') %>%
  mutate(timestep = ifelse(timePost==1,'post','pre'),
         comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
  group_by(timestep, comp, value, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    fwd10h_sim %>%
      pivot_longer(cols = c('Y_sim', 'Y_obs'),
                   names_to = 'source', values_to = 'value') %>%
      
      mutate(timestep = ifelse(timePost==1,'post','pre'),
             comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
      group_by(timestep, comp, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  ggplot(aes(x = value, y = proportion, fill = source))+
  geom_col(position = position_dodge())+
  scale_x_continuous(limits = c(-2, 50))+
  facet_grid(timestep~comp, scale = 'free_x')+
  theme_minimal()

fwd10h_prediction4

ggsave(fwd10h_prediction1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd10h_prediction1.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd10h_prediction2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd10h_prediction2.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd10h_prediction3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd10h_prediction3.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd10h_prediction4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd10h_predictions4.png'))

fwd10h_prediction5 = 
  fwd10h_sim %>%
  filter(draw == 1) %>%
  group_by(draw, Y_obs) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  ggplot()+
  geom_col(aes(x = Y_obs, y = count), 
           fill = 'red', col = NA, alpha = 0.5)+
  theme_minimal()+
  geom_errorbar(
    data = 
      fwd10h_sim %>%
      group_by(draw, Y_sim) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      complete(draw,
               Y_sim = 1:max(.$Y_sim)) %>%
      mutate(count = ifelse(is.na(count), 0, count)) %>%
      group_by(Y_sim) %>%
      summarise(count.025 = quantile(count, 0.025),
                count.975 = quantile(count, 0.975),
                count.min = min(count),
                count.max = max(count)) %>%
      ungroup() %>%
      mutate(source = 'sim'),
    aes(x = Y_sim, ymin = count.025, ymax = count.975),
    width = 1
  )+
  labs(x = '10-hour FWD tally', y = 'N observations')
ggsave(fwd10h_prediction5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd10h_prediction5.png'),
       height = 4, width = 6.5, units = 'in')


rm(fwd10h_sim)

#### fwd 100h ##################################################################

fwd100h_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'fwd100h_validation_data.rds'))

fwd100h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'fwd100h_fit.rds'))

fwd100h_posterior = 
  as_draws_df(fwd100h_fit$draws())

fwd100h_sim = 
  do.call('bind_rows',
          lapply(X = 1:4000,
                 FUN = function(draw){
                   
                   print(paste0('Working on draw', draw))
                   
                   # extract parameters for this draw
                   beta = 
                     fwd100h_posterior %>%
                     select(contains('beta')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   alpha_pre = fwd100h_posterior$alpha_pre[draw]
                   alpha_post = fwd100h_posterior$alpha_post[draw]
                   
                   rho_pre = fwd100h_posterior$rho_pre[draw]
                   rho_post = fwd100h_posterior$rho_post[draw]
                   
                   tau_pre = fwd100h_posterior$tau_pre[draw]
                   tau_post = fwd100h_posterior$tau_post[draw]
                   
                   plotEffects = 
                     fwd100h_posterior %>%
                     select(contains('plot_effect')) %>%
                     slice(draw) %>%
                     as.data.frame() %>%
                     as.numeric()
                   
                   zGP_pre = 
                     matrix(nrow = fwd100h_data$L, ncol = fwd100h_data$P,
                            byrow = FALSE, 
                            data = 
                              fwd100h_posterior %>%
                              select(contains('zGP_pre')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   zGP_post = 
                     matrix(nrow = fwd100h_data$L, ncol = fwd100h_data$P,
                            byrow = FALSE, 
                            data = 
                              fwd100h_posterior %>%
                              select(contains('zGP_post')) %>%
                              slice(draw) %>%
                              as.data.frame() %>%
                              as.numeric())
                   
                   kappa = fwd100h_posterior$kappa[draw]
                   
                   dists = 
                     fwd100h_data$coords %>%
                     dist() %>%
                     as.matrix()
                   
                   K_pre = alpha_pre * exp(-(dists**2)/(2*(rho_pre**2)))
                   K_pre = K_pre + diag(nrow = nrow(K_pre),
                                        ncol = ncol(K_pre),
                                        x = tau_pre)
                   KL_pre = t(chol(K_pre))
                   
                   GP_pre = 
                     matrix(nrow = fwd100h_data$L,
                            ncol = fwd100h_data$P,
                            data = 
                              sapply(X = 1:fwd100h_data$P,
                                     FUN = function(p){
                                       KL_pre %*% zGP_pre[,p]
                                     }))
                   
                   K_post = alpha_post * exp(-(dists**2)/(2*(rho_post**2)))
                   K_post = K_post + diag(nrow = nrow(K_post),
                                          ncol = ncol(K_post),
                                          x = tau_post)
                   KL_post = t(chol(K_post))
                   
                   GP_post = 
                     matrix(nrow = fwd100h_data$L,
                            ncol = fwd100h_data$P,
                            data = 
                              sapply(X = 1:fwd100h_data$P,
                                     FUN = function(p){
                                       KL_post %*% zGP_post[,p]
                                     }))
                   
                   logmu = 
                     sapply(X = 1:fwd100h_data$N,
                            FUN = function(i){
                              as.numeric(fwd100h_data$X[i,] %*% beta)+
                                plotEffects[fwd100h_data$plot_id[i]]+
                                (GP_pre[fwd100h_data$location_id[i],fwd100h_data$plot_id[i]] *
                                   fwd100h_data$pre[i]) +
                                (GP_post[fwd100h_data$location_id[i],
                                         fwd100h_data$plot_id[i]] *
                                   fwd100h_data$post[i])
                            })
                   
                   Y = rnbinom(n = fwd100h_data$N,
                               size = kappa,
                               mu = exp(logmu))
                   
                   results = 
                     fwd100h_data$X %>% as_tibble()
                   results$draw = draw
                   results$Y_sim = Y
                   results$logmu = logmu
                   results$Y_obs = fwd100h_data$Y
                   results$sample = 1:fwd100h_data$N
                   
                   return(results)
                 }))


fwd100h_mae.validation = 
  fwd100h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(mean_sim = mean(Y_sim)) %>%
  ungroup() %>%
  mutate(ae = abs(Y_obs-mean_sim)) %>%
  pull(ae) %>%
  mean()

fwd100h_mean.validation = 
  mean(fwd100h_data$Y)

fwd100h_wMAPE.validation = fwd100h_mae.validation/fwd100h_mean.validation


fwd100h_prediction1 = 
  fwd100h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(y.50 = median(Y_sim),
            y.mean = mean(Y_sim),
            y.05 = quantile(Y_sim, probs = 0.05),
            y.95 = quantile(Y_sim, probs = 0.95),
            logmu = mean(logmu)) %>%
  ungroup() %>% 
  ggplot()+
  geom_ribbon(aes(x = logmu, ymin = y.05, ymax = y.95), alpha = 0.5, color = 'grey')+
  geom_point(aes(x =logmu, y = Y_obs), size = 1)+
  theme_minimal()
fwd100h_prediction1

fwd100h_prediction3 = 
  fwd100h_sim %>%
  group_by(sample, Y_obs) %>%
  summarise(y.50 = median(Y_sim),
            y.mean = mean(Y_sim),
            y.05 = quantile(Y_sim, probs = 0.05),
            y.95 = quantile(Y_sim, probs = 0.95),
            logmu = mean(logmu)) %>%
  ungroup() %>%
  ggplot(aes(x = y.mean, y = Y_obs))+
  geom_point(size = 0)+
  geom_abline(intercept= 0, slope = 1, color = 'red')+
  geom_smooth(method = 'lm')+
  coord_fixed()+
  theme_minimal()+
  labs(x = 'Mean simulated tally', y = 'Observed tally')

fwd100h_prediction3

fwd100h_prediction2 = 
  fwd100h_sim %>%
  pivot_longer(cols = c('Y_sim', 'Y_obs'),
               names_to = 'source', values_to = 'value') %>%
  mutate(timestep = ifelse(timePost==1,'post','pre'),
         comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
  ggplot(aes(x = value, fill = source))+
  geom_bar(position = position_dodge())+
  scale_x_continuous(limits = c(-2, 50))+
  facet_grid(timestep~comp)+
  theme_minimal()

fwd100h_prediction2

fwd100h_prediction4 = 
  fwd100h_sim %>%
  pivot_longer(cols = c('Y_sim', 'Y_obs'),
               names_to = 'source', values_to = 'value') %>%
  mutate(timestep = ifelse(timePost==1,'post','pre'),
         comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
  group_by(timestep, comp, value, source) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(
    fwd100h_sim %>%
      pivot_longer(cols = c('Y_sim', 'Y_obs'),
                   names_to = 'source', values_to = 'value') %>%
      
      mutate(timestep = ifelse(timePost==1,'post','pre'),
             comp = ifelse(compB==1,'B',ifelse(compC==1,'C','A'))) %>%
      group_by(timestep, comp, source) %>%
      summarise(total = n()) %>%
      ungroup() 
  ) %>%
  mutate(proportion = count / total) %>%
  ggplot(aes(x = value, y = proportion, fill = source))+
  geom_col(position = position_dodge())+
  scale_x_continuous(limits = c(-2, 50))+
  facet_grid(timestep~comp, scale = 'free_x')+
  theme_minimal()

fwd100h_prediction4

ggsave(fwd100h_prediction1,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd100h_prediction1.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd100h_prediction2,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd100h_prediction2.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd100h_prediction3,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd100h_prediction3.png'),
       height = 4, width = 6.5, units = 'in')

ggsave(fwd100h_prediction4,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd100h_predictions4.png'))

fwd100h_prediction5 = 
  fwd100h_sim %>%
  filter(draw == 1) %>%
  group_by(draw, Y_obs) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(source = 'obs') %>%
  ggplot()+
  geom_col(aes(x = Y_obs, y = count), 
           fill = 'red', col = NA, alpha = 0.5)+
  theme_minimal()+
  geom_errorbar(
    data = 
      fwd100h_sim %>%
      group_by(draw, Y_sim) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      complete(draw,
               Y_sim = 1:max(.$Y_sim)) %>%
      mutate(count = ifelse(is.na(count), 0, count)) %>%
      group_by(Y_sim) %>%
      summarise(count.025 = quantile(count, 0.025),
                count.975 = quantile(count, 0.975),
                count.min = min(count),
                count.max = max(count)) %>%
      ungroup() %>%
      mutate(source = 'sim'),
    aes(x = Y_sim, ymin = count.025, ymax = count.975),
    width = 1
  )+
  labs(x = '100-hour FWD tally', y = 'N observations')

ggsave(fwd100h_prediction5,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd100h_prediction5.png'),
       height = 4, width = 6.5, units = 'in')


rm(fwd100h_sim)

#### MAEs and wMAPEs ###########################################################

prediction_summary_stats = 
  c(litter_mae.validation,
    litter_mean.validation,
    litter_wMAPE.validation,
    fwd1h_mae.validation,
    fwd1h_mean.validation,
    fwd1h_wMAPE.validation,
    fwd10h_mae.validation,
    fwd10h_mean.validation,
    fwd10h_wMAPE.validation,
    fwd100h_mae.validation,
    fwd100h_mean.validation,
    fwd100h_wMAPE.validation)

saveRDS(prediction_summary_stats,
        here::here('02-data', '03-results', 'prediction_summary_stats.rds'))
