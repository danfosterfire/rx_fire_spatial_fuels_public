# now the fun part


library(here)
library(tidyverse)
library(posterior)
library(bayesplot)
library(cowplot)
library(bayestestR)



#### litter prior v post ###################################################


litter_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'litter_fit.rds'))

litter_posterior = as_draws_df(litter_fit$draws(variables = c('beta','sigmaPlot','kappa',
                                               'alpha_pre', 'alpha_post',
                                               'rho_pre', 'rho_post')))

litter_summary = 
  litter_fit$summary(variables = c('beta', 'sigmaPlot', 
                                   'alpha_pre', 'alpha_post',
                                   'rho_pre', 'rho_post', 'kappa'))

litter_summary

write.csv(litter_summary,
          here::here('04-communication',
                     'tables',
                     'litter_summary.csv'))


litter_posterior_long = 
  litter_posterior %>%
  pivot_longer(cols = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]',
                        'beta[5]', 'beta[6]', 'alpha_pre', 'alpha_post',
                        'sigmaPlot', 'rho_pre', 'rho_post', 'kappa'),
               names_to = 'parameter',
               values_to = 'value') %>%
  left_join(data.frame(parameter = c('beta[1]', 'beta[2]', 'beta[3]', 
                                     'beta[4]', 'beta[5]', 'beta[6]',
                                     'alpha_pre', 'alpha_post', 'sigmaPlot',
                                     'rho_pre', 'rho_post', 'kappa'),
                       name = c('Intercept', 'compB', 'compC', 'timePost',
                                'compB:timePost', 'compC:timePost',
                                'GP magnitude (pre)', 'GP magnitude (post)', 'SD Plot',
                                'GP length scale (pre)', 'GP length scale (post)', 'NB Dispersion')) %>%
              mutate(name = factor(name, 
                                   levels = c('Intercept', 'compB', 'compC', 'timePost',
                                              'compB:timePost', 'compC:timePost',
                                              'GP magnitude (pre)', 'GP magnitude (post)', 'SD Plot',
                                              'GP length scale (pre)', 'GP length scale (post)', 
                                              'NB Dispersion'))))


litter_alpha_priorvpost = 
  litter_posterior %>%
  select(alpha_pre, alpha_post) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c('alpha_pre', 'alpha_post', 'prior'),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

litter_alpha_priorvpost

litter_beta_priorvpost = 
  litter_posterior %>%
  select(contains('beta')) %>%
  mutate(prior = rnorm(n = 36000, mean = 0, sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

litter_beta_priorvpost

litter_rho_priorvpost = 
  litter_posterior %>%
  select(rho_pre, rho_post) %>%
  mutate(prior = invgamma::rinvgamma(n = 36000,
                                     shape = 5,
                                     rate = 40)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

litter_rho_priorvpost

litter_sigmaplot_priorvpost = 
  litter_posterior %>%
  select(sigmaPlot) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

litter_sigmaplot_priorvpost

half_cauchy_samples = 
  sapply(X = 1:36000,
         FUN = function(d){
           value = rcauchy(n = 1, location = 0, scale = 5)
           while(value<=0){
             value = rcauchy(n = 1, location = 0, scale = 5)
           }
           return(value)
         })

litter_kappa_priorvpost = 
  litter_posterior %>%
  select(kappa) %>%
  mutate(prior = half_cauchy_samples) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'kappa',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = kappa))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 500))


litter_kappa_priorvpost

litter_priorvpost = 
  cowplot::plot_grid(litter_alpha_priorvpost,
                     litter_beta_priorvpost,
                     litter_rho_priorvpost,
                     litter_sigmaplot_priorvpost,
                     litter_kappa_priorvpost,
                     nrow = 5)

litter_priorvpost

ggsave(litter_priorvpost,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'litter_priorvpost.png'),
       height = 9, width = 6.5, units = 'in')

#### fwd1h prior v post ########################################################

fwd1h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'fwd1h_fit.rds'))

fwd1h_posterior = as_draws_df(fwd1h_fit$draws(variables = c('beta','sigmaPlot','kappa',
                                               'alpha_pre', 'alpha_post',
                                               'rho_pre', 'rho_post',
                                               'tau_pre', 'tau_post')))

fwd1h_summary = 
  fwd1h_fit$summary(variables = c('beta', 'sigmaPlot', 
                                   'alpha_pre', 'alpha_post',
                                   'rho_pre', 'rho_post', 
                                  'tau_pre', 'tau_post', 'kappa'))

fwd1h_summary

write.csv(fwd1h_summary,
          here::here('04-communication',
                     'tables',
                     'fwd1h_summary.csv'))


fwd1h_alpha_priorvpost = 
  fwd1h_posterior %>%
  select(alpha_pre, alpha_post) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c('alpha_pre', 'alpha_post', 'prior'),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd1h_alpha_priorvpost

fwd1h_beta_priorvpost = 
  fwd1h_posterior %>%
  select(contains('beta')) %>%
  mutate(prior = rnorm(n = 36000, mean = 0, sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd1h_beta_priorvpost

fwd1h_rho_priorvpost = 
  fwd1h_posterior %>%
  select(rho_pre, rho_post) %>%
  mutate(prior = invgamma::rinvgamma(n = 36000,
                                     shape = 5,
                                     rate = 40)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd1h_rho_priorvpost

fwd1h_tau_priorvpost = 
  fwd1h_posterior %>%
  select(tau_pre, tau_post) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c('tau_pre', 'tau_post', 'prior'),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd1h_tau_priorvpost



fwd1h_sigmaplot_priorvpost = 
  fwd1h_posterior %>%
  select(sigmaPlot) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()




fwd1h_sigmaplot_priorvpost

half_cauchy_samples = 
  sapply(X = 1:36000,
         FUN = function(d){
           value = rcauchy(n = 1, location = 0, scale = 5)
           while(value<=0){
             value = rcauchy(n = 1, location = 0, scale = 5)
           }
           return(value)
         })

fwd1h_kappa_priorvpost = 
  fwd1h_posterior %>%
  select(kappa) %>%
  mutate(prior = half_cauchy_samples) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 2000))

fwd1h_kappa_priorvpost

fwd1h_priorvpost = 
  cowplot::plot_grid(fwd1h_alpha_priorvpost,
                     fwd1h_beta_priorvpost,
                     fwd1h_rho_priorvpost,
                     fwd1h_tau_priorvpost,
                     fwd1h_sigmaplot_priorvpost,
                     fwd1h_kappa_priorvpost,
                     nrow = 6)

fwd1h_priorvpost

ggsave(fwd1h_priorvpost,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'fwd1h_priorvpost.png'),
       height = 9, width = 6.5, units = 'in')

#### fwd 10h prior v post ######################################################

fwd10h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'fwd10h_fit.rds'))

fwd10h_posterior = as_draws_df(fwd10h_fit$draws(variables = c('beta','sigmaPlot','kappa',
                                               'alpha_pre', 'alpha_post',
                                               'rho_pre', 'rho_post',
                                               'tau_pre', 'tau_post')))

fwd10h_summary = 
  fwd10h_fit$summary(variables = c('beta', 'sigmaPlot', 
                                   'alpha_pre', 'alpha_post',
                                   'rho_pre', 'rho_post', 
                                  'tau_pre', 'tau_post', 'kappa'))

fwd10h_summary

write.csv(fwd10h_summary,
          here::here('04-communication',
                     'tables',
                     'fwd10h_summary.csv'))


fwd10h_alpha_priorvpost = 
  fwd10h_posterior %>%
  select(alpha_pre, alpha_post) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c('alpha_pre', 'alpha_post', 'prior'),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd10h_alpha_priorvpost

fwd10h_beta_priorvpost = 
  fwd10h_posterior %>%
  select(contains('beta')) %>%
  mutate(prior = rnorm(n = 36000, mean = 0, sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd10h_beta_priorvpost

fwd10h_rho_priorvpost = 
  fwd10h_posterior %>%
  select(rho_pre, rho_post) %>%
  mutate(prior = invgamma::rinvgamma(n = 36000,
                                     shape = 5,
                                     rate = 40)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd10h_rho_priorvpost

fwd10h_tau_priorvpost = 
  fwd10h_posterior %>%
  select(tau_pre, tau_post) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c('tau_pre', 'tau_post', 'prior'),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd10h_tau_priorvpost



fwd10h_sigmaplot_priorvpost = 
  fwd10h_posterior %>%
  select(sigmaPlot) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()




fwd10h_sigmaplot_priorvpost

half_cauchy_samples = 
  sapply(X = 1:36000,
         FUN = function(d){
           value = rcauchy(n = 1, location = 0, scale = 5)
           while(value<=0){
             value = rcauchy(n = 1, location = 0, scale = 5)
           }
           return(value)
         })

fwd10h_kappa_priorvpost = 
  fwd10h_posterior %>%
  select(kappa) %>%
  mutate(prior = half_cauchy_samples) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 2000))

fwd10h_kappa_priorvpost

fwd10h_priorvpost = 
  cowplot::plot_grid(fwd10h_alpha_priorvpost,
                     fwd10h_beta_priorvpost,
                     fwd10h_rho_priorvpost,
                     fwd10h_tau_priorvpost,
                     fwd10h_sigmaplot_priorvpost,
                     fwd10h_kappa_priorvpost,
                     nrow = 6)

fwd10h_priorvpost

ggsave(fwd10h_priorvpost,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'fwd10h_priorvpost.png'),
       height = 9, width = 6.5, units = 'in')


#### fwd100h prior v post ######################################################

fwd100h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'fwd100h_fit.rds'))

fwd100h_posterior = as_draws_df(fwd100h_fit$draws(variables = c('beta','sigmaPlot','kappa',
                                               'alpha_pre', 'alpha_post',
                                               'rho_pre', 'rho_post',
                                               'tau_pre', 'tau_post')))

fwd100h_summary = 
  fwd100h_fit$summary(variables = c('beta', 'sigmaPlot', 
                                   'alpha_pre', 'alpha_post',
                                   'rho_pre', 'rho_post', 
                                  'tau_pre', 'tau_post', 'kappa'))

fwd100h_summary

write.csv(fwd100h_summary,
          here::here('04-communication',
                     'tables',
                     'fwd100h_summary.csv'))


fwd100h_alpha_priorvpost = 
  fwd100h_posterior %>%
  select(alpha_pre, alpha_post) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c('alpha_pre', 'alpha_post', 'prior'),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd100h_alpha_priorvpost

fwd100h_beta_priorvpost = 
  fwd100h_posterior %>%
  select(contains('beta')) %>%
  mutate(prior = rnorm(n = 36000, mean = 0, sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd100h_beta_priorvpost

fwd100h_rho_priorvpost = 
  fwd100h_posterior %>%
  select(rho_pre, rho_post) %>%
  mutate(prior = invgamma::rinvgamma(n = 36000,
                                     shape = 5,
                                     rate = 40)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd100h_rho_priorvpost

fwd100h_tau_priorvpost = 
  fwd100h_posterior %>%
  select(tau_pre, tau_post) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c('tau_pre', 'tau_post', 'prior'),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()

fwd100h_tau_priorvpost



fwd100h_sigmaplot_priorvpost = 
  fwd100h_posterior %>%
  select(sigmaPlot) %>%
  mutate(prior = truncnorm::rtruncnorm(n = 36000,
                                       a = 0,
                                       mean = 0,
                                       sd = 5)) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()




fwd100h_sigmaplot_priorvpost

half_cauchy_samples = 
  sapply(X = 1:36000,
         FUN = function(d){
           value = rcauchy(n = 1, location = 0, scale = 5)
           while(value<=0){
             value = rcauchy(n = 1, location = 0, scale = 5)
           }
           return(value)
         })

fwd100h_kappa_priorvpost = 
  fwd100h_posterior %>%
  select(kappa) %>%
  mutate(prior = half_cauchy_samples) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot(aes(x = value, fill = parameter))+
  geom_density(alpha = 0.6, color = NA)+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 2000))

fwd100h_kappa_priorvpost

fwd100h_priorvpost = 
  cowplot::plot_grid(fwd100h_alpha_priorvpost,
                     fwd100h_beta_priorvpost,
                     fwd100h_rho_priorvpost,
                     fwd100h_tau_priorvpost,
                     fwd100h_sigmaplot_priorvpost,
                     fwd100h_kappa_priorvpost,
                     nrow = 6)

fwd100h_priorvpost

ggsave(fwd100h_priorvpost,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'fwd100h_priorvpost.png'),
       height = 9, width = 6.5, units = 'in')



#### alphas ####################################################################

build_density_slice_df = 
  function(post, param, comp, time){
    result = 
      data.frame(
        x = 
          post %>% 
          filter(timestep==param&component==comp) %>%
          pull(value) %>%
          density() %>%
          .$x,
        y =
          post %>%
          filter(timestep==param & component==comp) %>%
          pull(value) %>%
          density() %>%
          .$y,
        xmin = 
          post %>%
          filter(timestep == param & component == comp) %>%
          pull(value) %>%
          #quantile(., 0.025),
          bayestestR::ci(x = ., ci = 0.95, method = 'ETI') %>%
          .$CI_low,
        xmax = 
          post %>%
          filter(timestep == param & component == comp) %>%
          pull(value) %>%
          #quantile(., 0.975),
          bayestestR::ci(x = ., ci = 0.95, method = 'ETI') %>%
          .$CI_high,
        Time = as.character(time),
        component = factor(comp, levels = c('Litter', '1-hour', '10-hour'))) %>%
      filter(x >= xmin & x <= xmax)
    return(result)
  }

long_alphas_posteriors = 
  litter_posterior %>%
  select(alpha_pre, alpha_post) %>%
  rowid_to_column('draw') %>%
  mutate(component = 'Litter') %>%
  bind_rows(
    fwd1h_posterior %>%
      select(alpha_pre, alpha_post) %>%
      rowid_to_column('draw') %>%
      mutate(component = '1-hour')
  ) %>%
  bind_rows(
    fwd10h_posterior %>%
      select(alpha_pre, alpha_post) %>%
      rowid_to_column('draw') %>%
      mutate(component = '10-hour')
  ) %>%
  pivot_longer(cols = c('alpha_pre', 'alpha_post'),
               names_to = 'timestep',
               values_to = 'value') %>%
  mutate(Time = 
           ifelse(timestep=='alpha_pre',
                  'Pre-burn',
                  'Post-burn')) %>%
  mutate(#Time = factor(Time, levels = c('Pre-burn', 'Post-burn')),
         component = factor(component, 
                            levels = c('Litter', '1-hour', '10-hour')))


alphas_posteriors_plot = 
  ggplot()+
  geom_line(data = long_alphas_posteriors, 
            aes(x = value, color = Time), stat = 'density', alpha = 0.75, lwd = 1)+
  geom_ribbon(
    data = build_density_slice_df(post = long_alphas_posteriors,
                                  param = 'alpha_pre',
                                  comp = 'Litter',
                                  time = 'Pre-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_alphas_posteriors,
                                  param = 'alpha_post',
                                  comp = 'Litter',
                                  time = 'Post-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_alphas_posteriors,
                                  param = 'alpha_pre',
                                  comp = '1-hour',
                                  time = 'Pre-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_alphas_posteriors,
                                  param = 'alpha_post',
                                  comp = '1-hour',
                                  time = 'Post-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_alphas_posteriors,
                                  param = 'alpha_pre',
                                  comp = '10-hour',
                                  time = 'Pre-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_alphas_posteriors,
                                  param = 'alpha_post',
                                  comp = '10-hour',
                                  time = 'Post-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_line(
    data = 
      data.frame(x = seq(long_alphas_posteriors %>% pull(value) %>% min(),
                         long_alphas_posteriors %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 5, a = 0, b= Inf)),
    aes(x = x, y = y, lty = 'Prior'),
    color = 'red'
  )+
  facet_grid(component~., scales = 'free_y')+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'D',
                        breaks = c('Pre-burn', 'Post-burn'))+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D',
                       breaks = c('Pre-burn', 'Post-burn'))+
  scale_linetype_manual(values = c('Prior' = 2),
                        name = NULL)+
  labs(y = 'Posterior density', x = 'GP magnitude')+
  theme(axis.text.y = element_blank())

alphas_posteriors_plot

ggsave(
  alphas_posteriors_plot+
  theme(text = element_text(size = 16)),
  filename = 
    here::here('04-communication',
               'figures',
               'powerpoint',
               'alphas_posteriors.png'),
  height = 6.5, width = 6.5, units = 'in')

ggsave(alphas_posteriors_plot,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'alphas_posteriors.png'),
       height = 6.5, width = 6.5, units = 'in')


#### rhos ######################################################################

long_rhos_posteriors = 
  litter_posterior %>%
  select(rho_pre, rho_post) %>%
  rowid_to_column('draw') %>%
  mutate(component = 'Litter') %>%
  bind_rows(
    fwd1h_posterior %>%
      select(rho_pre, rho_post) %>%
      rowid_to_column('draw') %>%
      mutate(component = '1-hour')
  ) %>%
  bind_rows(
    fwd10h_posterior %>%
      select(rho_pre, rho_post) %>%
      rowid_to_column('draw') %>%
      mutate(component = '10-hour')
  ) %>%
  pivot_longer(cols = c('rho_pre', 'rho_post'),
               names_to = 'timestep',
               values_to = 'value') %>%
  mutate(Time = 
           ifelse(timestep=='rho_pre',
                  'Pre-burn',
                  'Post-burn')) %>%
  mutate(#Time = factor(Time, levels = c('Pre-burn', 'Post-burn')),
    component = factor(component, 
                       levels = c('Litter', '1-hour', '10-hour')))


rhos_posteriors_plot = 
  ggplot()+
  geom_line(data = long_rhos_posteriors, 
            aes(x = value, color = Time), stat = 'density', alpha = 0.75, lwd = 1)+
  geom_ribbon(
    data = build_density_slice_df(post = long_rhos_posteriors,
                                  param = 'rho_pre',
                                  comp = 'Litter',
                                  time = 'Pre-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_rhos_posteriors,
                                  param = 'rho_post',
                                  comp = 'Litter',
                                  time = 'Post-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_rhos_posteriors,
                                  param = 'rho_pre',
                                  comp = '1-hour',
                                  time = 'Pre-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_rhos_posteriors,
                                  param = 'rho_post',
                                  comp = '1-hour',
                                  time = 'Post-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_rhos_posteriors,
                                  param = 'rho_pre',
                                  comp = '10-hour',
                                  time = 'Pre-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_rhos_posteriors,
                                  param = 'rho_post',
                                  comp = '10-hour',
                                  time = 'Post-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_line(
    data = 
      data.frame(x = seq(0,
                         long_rhos_posteriors %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = invgamma::dinvgamma(x = x, shape = 5, rate = 40)),
    aes(x = x, y = y, lty = 'Prior'),
    color = 'red'
  )+
  facet_grid(component~., scales = 'free_y')+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'D',
                        breaks = c('Pre-burn', 'Post-burn'))+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D',
                       breaks = c('Pre-burn', 'Post-burn'))+
  scale_linetype_manual(values = c('Prior' = 2),
                        name = NULL)+
  labs(y = 'Posterior density', x = 'GP length scale')+
  theme(axis.text.y = element_blank())+
  scale_x_continuous(limits = c(0, 12))

rhos_posteriors_plot

ggsave(
  rhos_posteriors_plot+
    theme(text = element_text(size = 16)),
  filename = 
    here::here('04-communication',
               'figures',
               'powerpoint',
               'rhos_posteriors.png'),
  height = 6.5, width = 6.5, units = 'in')

ggsave(rhos_posteriors_plot,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'rhos_posteriors.png'),
       height = 6.5, width = 6.5, units = 'in')



#### taus ######################################################################
long_taus_posteriors = 
  bind_rows(
    fwd1h_posterior %>%
      select(tau_pre, tau_post) %>%
      rowid_to_column('draw') %>%
      mutate(component = '1-hour')
  ) %>%
  bind_rows(
    fwd10h_posterior %>%
      select(tau_pre, tau_post) %>%
      rowid_to_column('draw') %>%
      mutate(component = '10-hour')
  ) %>%
  pivot_longer(cols = c('tau_pre', 'tau_post'),
               names_to = 'timestep',
               values_to = 'value') %>%
  mutate(Time = 
           ifelse(timestep=='tau_pre',
                  'Pre-burn',
                  'Post-burn')) %>%
  mutate(#Time = factor(Time, levels = c('Pre-burn', 'Post-burn')),
    component = factor(component, 
                       levels = c('1-hour', '10-hour')))


taus_posteriors_plot = 
  ggplot()+
  geom_line(data = long_taus_posteriors, 
            aes(x = value, color = Time), stat = 'density', alpha = 0.75, lwd = 1)+
  geom_ribbon(
    data = build_density_slice_df(post = long_taus_posteriors,
                                  param = 'tau_pre',
                                  comp = '1-hour',
                                  time = 'Pre-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_taus_posteriors,
                                  param = 'tau_post',
                                  comp = '1-hour',
                                  time = 'Post-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_taus_posteriors,
                                  param = 'tau_pre',
                                  comp = '10-hour',
                                  time = 'Pre-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_ribbon(
    data = build_density_slice_df(post = long_taus_posteriors,
                                  param = 'tau_post',
                                  comp = '10-hour',
                                  time = 'Post-burn'),
    aes(ymax = y, ymin = 0, x = x, fill = Time),
    alpha = 0.75
  )+
  geom_line(
    data = 
      data.frame(x = seq(0,
                         long_taus_posteriors %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 5, a = 0, b = Inf)),
    aes(x = x, y = y, lty = 'Prior'),
    color = 'red'
  )+
  facet_grid(component~., scales = 'free_y')+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'D',
                        breaks = c('Pre-burn', 'Post-burn'))+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D',
                       breaks = c('Pre-burn', 'Post-burn'))+
  scale_linetype_manual(values = c('Prior' = 2),
                        name = NULL)+
  labs(y = 'Posterior density', x = 'GP nugget')+
  theme(axis.text.y = element_blank())

taus_posteriors_plot

ggsave(
  taus_posteriors_plot+
    theme(text = element_text(size = 16)),
  filename = 
    here::here('04-communication',
               'figures',
               'powerpoint',
               'taus_posteriors.png'),
  height = 6.5, width = 6.5, units = 'in')

ggsave(taus_posteriors_plot,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'taus_posteriors.png'),
       height = 6.5, width = 6.5, units = 'in')


#### summary table #############################################################

summary_table = 
  litter_posterior %>%
  mutate(component = 'Litter') %>%
  bind_rows(fwd1h_posterior %>%
              mutate(component = '1-hour FWD')) %>%
  bind_rows(fwd10h_posterior %>%
              mutate(component = '10-hour FWD')) %>%
  as_tibble() %>%
  select(component,
         intercept = `beta[1]`,
         compB = `beta[2]`,
         compC = `beta[3]`,
         timePost = `beta[4]`,
         compBtimePost = `beta[5]`,
         compCtimePost = `beta[6]`,
         sigma_plot = sigmaPlot,
         alpha_pre,
         alpha_post,
         rho_pre,
         rho_post,
         tau_pre,
         tau_post,
         kappa,
         draw = .draw) %>%
  pivot_longer(cols = c(intercept, compB, compC, timePost,
                        compBtimePost, compCtimePost, sigma_plot,
                        alpha_pre, alpha_post, rho_pre, rho_post, tau_pre, tau_post,
                        kappa),
               names_to = 'parameter',
               values_to = 'value') %>%
  filter(!is.na(value)) %>%
  group_by(component, parameter) %>%
  summarise(val_mean = mean(value),
            val_median = median(value),
            q025 = quantile(value, p = c(0.025), na.rm = TRUE),
            q975 = quantile(value, p = c(0.975), na.rm = TRUE),
            cilow = bayestestR::ci(x = value, ci = 0.95, method = 'ETI')$CI_low,
            cihigh = bayestestR::ci(x = value, ci = 0.95, method = 'ETI')$CI_high) %>%
  ungroup() %>%
  mutate(Estimate = paste0(round(val_median, 2),
                           ' (',
                           round(cilow, 2),
                           ', ',
                           round(cihigh, 2),
                           ')')) %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'component',
              values_from = 'Estimate') %>%
  mutate(parameter = 
           factor(parameter, levels = c('intercept', 'compB', 'compC',
                                        'timePost', 'compBtimePost', 
                                        'compCtimePost', 'sigma_plot',
                                        'alpha_pre', 'alpha_post',
                                        'rho_pre', 'rho_post',
                                        'tau_pre', 'tau_post',
                                        'kappa'))) %>%
  arrange(parameter) %>%
  select(parameter, `Litter`, `1-hour FWD`, `10-hour FWD`)

summary_table

write.csv(summary_table,
          here::here('04-communication',
                     'tables',
                     'summary_table.csv'),
          row.names = FALSE)

#### scratch ###################################################################

litter_posterior %>%
  select(alpha_pre, alpha_post, rho_pre, rho_post) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = c('parameter', 'timestep'),
               names_sep = '_',
               values_to = 'value') %>%
  left_join(data.frame(parameter = c('alpha', 'rho'),
                       param_name = c('GP std. dev.', 'GP range'))) %>%
  left_join(data.frame(timestep = c('pre', 'post'),
                       Time = c('Preburn', 'Postburn')) %>%
              mutate(Time = factor(Time,
                                            levels = c('Preburn', 'Postburn')))) %>%
  ggplot(aes(x = value, fill = Time))+
  geom_density(alpha = 0.6)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  facet_wrap(~param_name, scales = 'free')+
  theme_minimal()+
  labs(y = 'Posterior density', x = 'Parameter value',
       title = 'Litter')+
  theme(plot.title = element_text(size = 16))














#### scratch ###################################################################


litter_posteriors_alpha = 
  litter_posterior_long %>%
  filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
  ggplot(aes(fill = name))+
  geom_density(aes(x = value), alpha = 0.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  theme_minimal()+
  geom_errorbarh(
    data = 
      litter_posterior_long %>%
      filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
      group_by(name, parameter) %>%
      summarise(value.med = median(value),
                value.25 = quantile(value, probs = 0.25),
                value.75 = quantile(value, probs = 0.75),
                value.025 = quantile(value, probs = 0.025),
                value.975 = quantile(value, probs = 0.975)) %>%
      ungroup(),
    aes(y = 11, xmin = value.025, xmax = value.975))+
  geom_point(
    data = 
      litter_posterior_long %>%
        filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
        group_by(name, parameter) %>%
        summarise(value.med = median(value),
                  value.25 = quantile(value, probs = 0.25),
                  value.75 = quantile(value, probs = 0.75),
                  value.025 = quantile(value, probs = 0.025),
                  value.975 = quantile(value, probs = 0.975)) %>%
        ungroup(),
    aes(y = 11, x = value.med), size = 5, pch = 21)+
  labs(y = 'Posterior Probability Density', x = 'Value', fill = 'Parameter')

litter_posteriors_alpha

ggsave(litter_posteriors_alpha,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_posteriors_alpha.png'),
       width = 6.5, height = 4, units = 'in')

litter_posteriors_rho = 
  litter_posterior_long %>%
  filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
  ggplot(aes(fill = name))+
  geom_density(aes(x = value), alpha = 0.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  theme_minimal()+
  geom_errorbarh(
    data = 
      litter_posterior_long %>%
      filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
      group_by(name, parameter) %>%
      summarise(value.med = median(value),
                value.25 = quantile(value, probs = 0.25),
                value.75 = quantile(value, probs = 0.75),
                value.025 = quantile(value, probs = 0.025),
                value.975 = quantile(value, probs = 0.975)) %>%
      ungroup(),
    aes(y = 7, xmin = value.025, xmax = value.975))+
  geom_point(
    data = 
      litter_posterior_long %>%
        filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
        group_by(name, parameter) %>%
        summarise(value.med = median(value),
                  value.25 = quantile(value, probs = 0.25),
                  value.75 = quantile(value, probs = 0.75),
                  value.025 = quantile(value, probs = 0.025),
                  value.975 = quantile(value, probs = 0.975)) %>%
        ungroup(),
    aes(y = 7, x = value.med), size = 5, pch = 21)+
  labs(y = 'Posterior Probability Density', x = 'Value', fill = 'Parameter')

litter_posteriors_rho

ggsave(litter_posteriors_rho,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_posteriors_rho.png'),
       height = 4, width = 6.5, units = 'in')



litter_posteriors_corrfn = 
  litter_posterior_long %>%
  filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
  group_by(name, parameter) %>%
  summarise(value.025 = quantile(value, probs = 0.025),
            value.med = median(value),
            value.975 = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  pivot_longer(cols = c('value.025', 'value.med', 'value.975'),
               names_to = 'quantile',
               values_to = 'rho') %>%
  expand(nesting(name, parameter, quantile, rho),
         dist_m = seq(from = 0, to = 8, by = 0.1)) %>%
  mutate(timestep = gsub(x = parameter,
                         pattern = '^.*_',
                         replacement = ''),
         quantile = gsub(x = quantile,
                         pattern = 'value',
                         replacement = '')) %>%
  select(-parameter,-name) %>%
  mutate(corr = exp(-(dist_m**2)/(2*rho**2))) %>%
  select(-rho) %>%
  pivot_wider(names_from = c('quantile'),
              values_from = 'corr',
              names_sep = '_') %>%
  ggplot(aes(x = dist_m))+
  geom_line(aes(y = .med, color = timestep), lwd = 1)+
  geom_ribbon(aes(ymin = `.025`, ymax = `.975`, fill = timestep),
              color = NA, alpha = 0.5)+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  labs(x = 'Distance (m)', y = 'Correlation', fill = 'Litter/duff', color = 'Litter/duff')


litter_posteriors_corrfn

ggsave(litter_posteriors_corrfn,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'litter_posteriors_corrfn.png'),
       height = 4, width = 6.5, units = 'in')


litter_posteriors_beta = 
  litter_posterior_long %>%
  filter(is.element(parameter, c('beta[1]', 'beta[2]', 'beta[3]', 
                                 'beta[4]', 'beta[5]', 'beta[6]'))) %>%
  group_by(name, parameter) %>%
  summarise(value.med = median(value),
            value.25 = quantile(value, probs = 0.25),
            value.75 = quantile(value, probs = 0.75),
            value.025 = quantile(value, probs = 0.025),
            value.975 = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  mutate(name = factor(as.character(name),
                       levels = c('compC:timePost','compB:timePost','timePost',
                                  'compC', 'compB', 'Intercept'))) %>%
  ggplot(aes(y = name))+
  geom_errorbarh(aes(xmin = value.025, xmax = value.975), height = 0.3)+
  geom_point(aes(x = value.med), size = 5)+
  labs(y = 'Fixed EFfect Parameter', x = 'Value')+
  theme_minimal()+
  geom_vline(xintercept = 0, lwd = 1, lty = 2)

litter_posteriors_beta

ggsave(litter_posteriors_beta,
       filename = here::here('04-communication',
                  'figures',
                  'manuscript',
                  'litter_posteriors_beta.png'),
       height = 4, width = 6.5, units = 'in')


litter_summary = 
  litter_fit$summary(variables = c('beta', 'alpha_pre', 'alpha_post',
                                     'rho_pre', 'rho_post', 'sigmaPlot', 'kappa')) %>%
  left_join(data.frame(parameter = c('beta[1]', 'beta[2]', 'beta[3]', 
                                     'beta[4]', 'beta[5]', 'beta[6]',
                                     'alpha_pre', 'alpha_post', 'sigmaPlot',
                                     'rho_pre', 'rho_post', 'kappa'),
                       name = c('Intercept', 'compB', 'compC', 'timePost',
                                'compB:timePost', 'compC:timePost',
                                'GP SD (pre)', 'GP SD (post)', 'SD Plot',
                                'GP Range (pre)', 'GP Range (post)', 'NB Dispersion')),
            by = c('variable' = 'parameter')) %>%
  select(name, mean, median, sd, mad, q5, q95, rhat, ess_bulk, ess_tail)

write.csv(litter_summary,
          here::here('04-communication',
                     'tables',
                     'litter_gp_summary.csv'),
          row.names = FALSE)


set.seed(112188)
litter_posterior_med = 
  litter_posterior %>%
  summarise_all(median)


# simulate grids of pre and post
litter_coords = 
  expand.grid(x_coord = seq(from = 0, to = 10, by = 0.15),
              y_coord = seq(from = 0, to = 10, by = 0.15))

litter_D = 
  litter_coords %>%
  dist() %>%
  as.matrix()

litter_Sigma_pre = 
  (litter_posterior_med$alpha_pre**2)*
  exp(-(1/(2*(litter_posterior_med$rho_pre**2)))*(litter_D**2))+
  diag(x = 1e-9, nrow = nrow(litter_D), ncol = ncol(litter_D))

litter_Sigma_post = 
  (litter_posterior_med$alpha_post**2)*
  exp(-(1/(2*(litter_posterior_med$rho_post**2)))*(litter_D**2))+
  diag(x = 1e-9, nrow = nrow(litter_D), ncol = ncol(litter_D))
  
litter_GP_pre = 
  t(chol(litter_Sigma_pre)) %*% 
  rnorm(n = nrow(litter_D), mean = 0, sd = 1)

litter_GP_post = 
  t(chol(litter_Sigma_post)) %*% 
  rnorm(n = nrow(litter_D), mean = 0, sd = 1)


litter_sim_grid = 
  litter_coords %>%
  mutate(timestep = 'pre',
         GP = as.numeric(litter_GP_pre),
         intercept = 1,
         post = 0) %>%
  bind_rows(
    litter_coords %>%
      mutate(timestep='post',
             GP = as.numeric(litter_GP_post),
             intercept = 1,
             post = 1)
  ) %>%
  mutate(XB = 
           as.numeric((intercept*litter_posterior_med$`beta[1]`)+
                        (post*litter_posterior_med$`beta[4]`))) %>%
  mutate(logmu = XB+GP,
         mu = exp(logmu),
         timestep = factor(timestep, levels = c('pre', 'post')))

litter_sim_grid$depth_cm = 
  rnbinom(n = nrow(litter_sim_grid),
          mu = litter_sim_grid$mu,
          size = litter_posterior_med$kappa)


ggplot(data = litter_sim_grid,
       aes(x = x_coord, y = y_coord, fill = depth_cm))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  facet_wrap(~timestep)+
  coord_fixed()


#### fwd 1h ####################################################################

fwd1h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'fwd1h_fit.rds'))

fwd1h_posterior = as_draws_df(fwd1h_fit$draws(variables = c('beta','sigmaPlot','kappa',
                                               'alpha_pre', 'alpha_post',
                                               'rho_pre', 'rho_post', 'tau_pre',
                                               'tau_post')))

fwd1h_posterior_long = 
  fwd1h_posterior %>%
  pivot_longer(cols = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]',
                        'beta[5]', 'beta[6]', 'alpha_pre', 'alpha_post',
                        'sigmaPlot', 'rho_pre', 'rho_post', 'kappa',
                        'tau_pre', 'tau_post'),
               names_to = 'parameter',
               values_to = 'value') %>%
  left_join(data.frame(parameter = c('beta[1]', 'beta[2]', 'beta[3]', 
                                     'beta[4]', 'beta[5]', 'beta[6]',
                                     'alpha_pre', 'alpha_post', 'sigmaPlot',
                                     'rho_pre', 'rho_post', 'tau_pre',
                                     'tau_post', 'kappa'),
                       name = c('Intercept', 'compB', 'compC', 'timePost',
                                'compB:timePost', 'compC:timePost',
                                'GP SD (pre)', 'GP SD (post)', 'SD Plot',
                                'GP Range (pre)', 'GP Range (post)', 
                                'GP Noise (pre)', 'GP Noise (post)', 'NB Dispersion')) %>%
              mutate(name = factor(name, 
                                   levels = c('Intercept', 'compB', 'compC', 'timePost',
                                              'compB:timePost', 'compC:timePost',
                                              'GP SD (pre)', 'GP SD (post)', 'SD Plot',
                                              'GP Range (pre)', 'GP Range (post)', 
                                              'GP Noise (pre)', 'GP Noise (post)',
                                              'NB Dispersion'))))


fwd1h_posteriors_alpha = 
  fwd1h_posterior_long %>%
  filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
  ggplot(aes(fill = name))+
  geom_density(aes(x = value), alpha = 0.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  theme_minimal()+
  geom_errorbarh(
    data = 
      fwd1h_posterior_long %>%
      filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
      group_by(name, parameter) %>%
      summarise(value.med = median(value),
                value.25 = quantile(value, probs = 0.25),
                value.75 = quantile(value, probs = 0.75),
                value.025 = quantile(value, probs = 0.025),
                value.975 = quantile(value, probs = 0.975)) %>%
      ungroup(),
    aes(y = 15, xmin = value.025, xmax = value.975))+
  geom_point(
    data = 
      fwd1h_posterior_long %>%
        filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
        group_by(name, parameter) %>%
        summarise(value.med = median(value),
                  value.25 = quantile(value, probs = 0.25),
                  value.75 = quantile(value, probs = 0.75),
                  value.025 = quantile(value, probs = 0.025),
                  value.975 = quantile(value, probs = 0.975)) %>%
        ungroup(),
    aes(y = 15, x = value.med), size = 5, pch = 21)+
  labs(y = 'Posterior Probability Density', x = 'Value', fill = 'Parameter')

fwd1h_posteriors_alpha

ggsave(fwd1h_posteriors_alpha,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd1h_posteriors_alpha.png'),
       width = 6.5, height = 4, units = 'in')

fwd1h_posteriors_rho = 
  fwd1h_posterior_long %>%
  filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
  ggplot(aes(fill = name))+
  geom_density(aes(x = value), alpha = 0.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  theme_minimal()+
  geom_errorbarh(
    data = 
      fwd1h_posterior_long %>%
      filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
      group_by(name, parameter) %>%
      summarise(value.med = median(value),
                value.25 = quantile(value, probs = 0.25),
                value.75 = quantile(value, probs = 0.75),
                value.025 = quantile(value, probs = 0.025),
                value.975 = quantile(value, probs = 0.975)) %>%
      ungroup(),
    aes(y = 3.5, xmin = value.025, xmax = value.975), height = 0.25)+
  geom_point(
    data = 
      fwd1h_posterior_long %>%
        filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
        group_by(name, parameter) %>%
        summarise(value.med = median(value),
                  value.25 = quantile(value, probs = 0.25),
                  value.75 = quantile(value, probs = 0.75),
                  value.025 = quantile(value, probs = 0.025),
                  value.975 = quantile(value, probs = 0.975)) %>%
        ungroup(),
    aes(y = 3.5, x = value.med), size = 5, pch = 21)+
  labs(y = 'Posterior Probability Density', x = 'Value', fill = 'Parameter')

fwd1h_posteriors_rho

ggsave(fwd1h_posteriors_rho,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd1h_posteriors_rho.png'),
       height = 4, width = 6.5, units = 'in')


fwd1h_posteriors_corrfn = 
  fwd1h_posterior_long %>%
  filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
  group_by(name, parameter) %>%
  summarise(value.025 = quantile(value, probs = 0.025),
            value.med = median(value),
            value.975 = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  pivot_longer(cols = c('value.025', 'value.med', 'value.975'),
               names_to = 'quantile',
               values_to = 'rho') %>%
  expand(nesting(name, parameter, quantile, rho),
         dist_m = seq(from = 0, to = 25, by = 0.1)) %>%
  mutate(timestep = gsub(x = parameter,
                         pattern = '^.*_',
                         replacement = ''),
         quantile = gsub(x = quantile,
                         pattern = 'value',
                         replacement = '')) %>%
  select(-parameter,-name) %>%
  mutate(corr = exp(-(dist_m**2)/(2*rho**2))) %>%
  select(-rho) %>%
  pivot_wider(names_from = c('quantile'),
              values_from = 'corr',
              names_sep = '_') %>%
  ggplot(aes(x = dist_m))+
  geom_line(aes(y = .med, color = timestep), lwd = 1)+
  geom_ribbon(aes(ymin = `.025`, ymax = `.975`, fill = timestep),
              color = NA, alpha = 0.5)+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  labs(x = 'Distance (m)', y = 'Correlation', fill = 'FWD 1h', color = 'FWD 1h')


fwd1h_posteriors_corrfn

ggsave(fwd1h_posteriors_corrfn,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'fwd1h_posteriors_corrfn.png'),
       height = 4, width = 6.5, units = 'in')

fwd1h_posteriors_tau = 
  fwd1h_posterior_long %>%
  filter(is.element(parameter, c('tau_pre', 'tau_post'))) %>%
  ggplot(aes(fill = name))+
  geom_density(aes(x = value), alpha = 0.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  theme_minimal()+
  geom_errorbarh(
    data = 
      fwd1h_posterior_long %>%
      filter(is.element(parameter, c('tau_pre', 'tau_post'))) %>%
      group_by(name, parameter) %>%
      summarise(value.med = median(value),
                value.25 = quantile(value, probs = 0.25),
                value.75 = quantile(value, probs = 0.75),
                value.025 = quantile(value, probs = 0.025),
                value.975 = quantile(value, probs = 0.975)) %>%
      ungroup(),
    aes(y = 30, xmin = value.025, xmax = value.975), height = 0.25)+
  geom_point(
    data = 
      fwd1h_posterior_long %>%
        filter(is.element(parameter, c('tau_pre', 'tau_post'))) %>%
        group_by(name, parameter) %>%
        summarise(value.med = median(value),
                  value.25 = quantile(value, probs = 0.25),
                  value.75 = quantile(value, probs = 0.75),
                  value.025 = quantile(value, probs = 0.025),
                  value.975 = quantile(value, probs = 0.975)) %>%
        ungroup(),
    aes(y = 30, x = value.med), size = 5, pch = 21)+
  labs(y = 'Posterior Probability Density', x = 'Value', fill = 'Parameter')

fwd1h_posteriors_tau

ggsave(fwd1h_posteriors_tau,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd1h_posteriors_tau.png'),
       height = 4, width = 6.5, units = 'in')




fwd1h_posteriors_beta = 
  fwd1h_posterior_long %>%
  filter(is.element(parameter, c('beta[1]', 'beta[2]', 'beta[3]', 
                                 'beta[4]', 'beta[5]', 'beta[6]'))) %>%
  group_by(name, parameter) %>%
  summarise(value.med = median(value),
            value.25 = quantile(value, probs = 0.25),
            value.75 = quantile(value, probs = 0.75),
            value.025 = quantile(value, probs = 0.025),
            value.975 = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  mutate(name = factor(as.character(name),
                       levels = c('compC:timePost','compB:timePost','timePost',
                                  'compC', 'compB', 'Intercept'))) %>%
  ggplot(aes(y = name))+
  geom_errorbarh(aes(xmin = value.025, xmax = value.975), height = 0.3)+
  geom_point(aes(x = value.med), size = 5)+
  labs(y = 'Fixed EFfect Parameter', x = 'Value')+
  theme_minimal()+
  geom_vline(xintercept = 0, lwd = 1, lty = 2)

fwd1h_posteriors_beta

ggsave(fwd1h_posteriors_beta,
       filename = here::here('04-communication',
                  'figures',
                  'manuscript',
                  'fwd1h_posteriors_beta.png'),
       height = 4, width = 6.5, units = 'in')


fwd1h_summary = 
  fwd1h_fit$summary(variables = c('beta', 'alpha_pre', 'alpha_post',
                                     'rho_pre', 'rho_post', 'sigmaPlot', 'kappa')) %>%
  left_join(data.frame(parameter = c('beta[1]', 'beta[2]', 'beta[3]', 
                                     'beta[4]', 'beta[5]', 'beta[6]',
                                     'alpha_pre', 'alpha_post', 'sigmaPlot',
                                     'rho_pre', 'rho_post', 'tau_pre', 'tau_post', 'kappa'),
                       name = c('Intercept', 'compB', 'compC', 'timePost',
                                'compB:timePost', 'compC:timePost',
                                'GP SD (pre)', 'GP SD (post)', 'SD Plot',
                                'GP Range (pre)', 'GP Range (post)', 
                                'GP Noise (pre)', 'GP Noise (post)', 'NB Dispersion')),
            by = c('variable' = 'parameter')) %>%
  select(name, mean, median, sd, mad, q5, q95, rhat, ess_bulk, ess_tail)

write.csv(fwd1h_summary,
          here::here('04-communication',
                     'tables',
                     'fwd1h_gp_summary.csv'),
          row.names = FALSE)

#### fwd 10h ###################################################################

fwd10h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'fwd10h_fit.rds'))

fwd10h_posterior = as_draws_df(fwd10h_fit$draws(variables = c('beta','sigmaPlot','kappa',
                                               'alpha_pre', 'alpha_post',
                                               'rho_pre', 'rho_post', 'tau_pre',
                                               'tau_post')))

fwd10h_posterior_long = 
  fwd10h_posterior %>%
  pivot_longer(cols = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]',
                        'beta[5]', 'beta[6]', 'alpha_pre', 'alpha_post',
                        'sigmaPlot', 'rho_pre', 'rho_post', 'kappa',
                        'tau_pre', 'tau_post'),
               names_to = 'parameter',
               values_to = 'value') %>%
  left_join(data.frame(parameter = c('beta[1]', 'beta[2]', 'beta[3]', 
                                     'beta[4]', 'beta[5]', 'beta[6]',
                                     'alpha_pre', 'alpha_post', 'sigmaPlot',
                                     'rho_pre', 'rho_post', 'tau_pre',
                                     'tau_post', 'kappa'),
                       name = c('Intercept', 'compB', 'compC', 'timePost',
                                'compB:timePost', 'compC:timePost',
                                'GP SD (pre)', 'GP SD (post)', 'SD Plot',
                                'GP Range (pre)', 'GP Range (post)', 
                                'GP Noise (pre)', 'GP Noise (post)', 'NB Dispersion')) %>%
              mutate(name = factor(name, 
                                   levels = c('Intercept', 'compB', 'compC', 'timePost',
                                              'compB:timePost', 'compC:timePost',
                                              'GP SD (pre)', 'GP SD (post)', 'SD Plot',
                                              'GP Range (pre)', 'GP Range (post)', 
                                              'GP Noise (pre)', 'GP Noise (post)',
                                              'NB Dispersion'))))


fwd10h_posteriors_alpha = 
  fwd10h_posterior_long %>%
  filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
  ggplot(aes(fill = name))+
  geom_density(aes(x = value), alpha = 0.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  theme_minimal()+
  geom_errorbarh(
    data = 
      fwd10h_posterior_long %>%
      filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
      group_by(name, parameter) %>%
      summarise(value.med = median(value),
                value.25 = quantile(value, probs = 0.25),
                value.75 = quantile(value, probs = 0.75),
                value.025 = quantile(value, probs = 0.025),
                value.975 = quantile(value, probs = 0.975)) %>%
      ungroup(),
    aes(y = 15, xmin = value.025, xmax = value.975))+
  geom_point(
    data = 
      fwd10h_posterior_long %>%
        filter(is.element(parameter, c('alpha_pre', 'alpha_post'))) %>%
        group_by(name, parameter) %>%
        summarise(value.med = median(value),
                  value.25 = quantile(value, probs = 0.25),
                  value.75 = quantile(value, probs = 0.75),
                  value.025 = quantile(value, probs = 0.025),
                  value.975 = quantile(value, probs = 0.975)) %>%
        ungroup(),
    aes(y = 15, x = value.med), size = 5, pch = 21)+
  labs(y = 'Posterior Probability Density', x = 'Value', fill = 'Parameter')

fwd10h_posteriors_alpha

ggsave(fwd10h_posteriors_alpha,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd10h_posteriors_alpha.png'),
       width = 6.5, height = 4, units = 'in')

fwd10h_posteriors_rho = 
  fwd10h_posterior_long %>%
  filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
  ggplot(aes(fill = name))+
  geom_density(aes(x = value), alpha = 0.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  theme_minimal()+
  geom_errorbarh(
    data = 
      fwd10h_posterior_long %>%
      filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
      group_by(name, parameter) %>%
      summarise(value.med = median(value),
                value.25 = quantile(value, probs = 0.25),
                value.75 = quantile(value, probs = 0.75),
                value.025 = quantile(value, probs = 0.025),
                value.975 = quantile(value, probs = 0.975)) %>%
      ungroup(),
    aes(y = 3.5, xmin = value.025, xmax = value.975), height = 0.25)+
  geom_point(
    data = 
      fwd10h_posterior_long %>%
        filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
        group_by(name, parameter) %>%
        summarise(value.med = median(value),
                  value.25 = quantile(value, probs = 0.25),
                  value.75 = quantile(value, probs = 0.75),
                  value.025 = quantile(value, probs = 0.025),
                  value.975 = quantile(value, probs = 0.975)) %>%
        ungroup(),
    aes(y = 3.5, x = value.med), size = 5, pch = 21)+
  labs(y = 'Posterior Probability Density', x = 'Value', fill = 'Parameter')

fwd10h_posteriors_rho

ggsave(fwd10h_posteriors_rho,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd10h_posteriors_rho.png'),
       height = 4, width = 6.5, units = 'in')


fwd10h_posteriors_corrfn = 
  fwd10h_posterior_long %>%
  filter(is.element(parameter, c('rho_pre', 'rho_post'))) %>%
  group_by(name, parameter) %>%
  summarise(value.025 = quantile(value, probs = 0.025),
            value.med = median(value),
            value.975 = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  pivot_longer(cols = c('value.025', 'value.med', 'value.975'),
               names_to = 'quantile',
               values_to = 'rho') %>%
  expand(nesting(name, parameter, quantile, rho),
         dist_m = seq(from = 0, to = 15, by = 0.1)) %>%
  mutate(timestep = gsub(x = parameter,
                         pattern = '^.*_',
                         replacement = ''),
         quantile = gsub(x = quantile,
                         pattern = 'value',
                         replacement = '')) %>%
  select(-parameter,-name) %>%
  mutate(corr = exp(-(dist_m**2)/(2*rho**2))) %>%
  select(-rho) %>%
  pivot_wider(names_from = c('quantile'),
              values_from = 'corr',
              names_sep = '_') %>%
  ggplot(aes(x = dist_m))+
  geom_line(aes(y = .med, color = timestep), lwd = 1)+
  geom_ribbon(aes(ymin = `.025`, ymax = `.975`, fill = timestep),
              color = NA, alpha = 0.5)+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  labs(x = 'Distance (m)', y = 'Correlation', fill = 'FWD 10h', color = 'FWD 10h')


fwd10h_posteriors_corrfn

ggsave(fwd10h_posteriors_corrfn,
       filename = 
         here::here('04-communication',
                    'figures',
                    'manuscript',
                    'fwd10h_posteriors_corrfn.png'),
       height = 4, width = 6.5, units = 'in')


fwd10h_posteriors_tau = 
  fwd10h_posterior_long %>%
  filter(is.element(parameter, c('tau_pre', 'tau_post'))) %>%
  ggplot(aes(fill = name))+
  geom_density(aes(x = value), alpha = 0.5)+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'D')+
  theme_minimal()+
  geom_errorbarh(
    data = 
      fwd10h_posterior_long %>%
      filter(is.element(parameter, c('tau_pre', 'tau_post'))) %>%
      group_by(name, parameter) %>%
      summarise(value.med = median(value),
                value.25 = quantile(value, probs = 0.25),
                value.75 = quantile(value, probs = 0.75),
                value.025 = quantile(value, probs = 0.025),
                value.975 = quantile(value, probs = 0.975)) %>%
      ungroup(),
    aes(y = 15, xmin = value.025, xmax = value.975), height = 0.25)+
  geom_point(
    data = 
      fwd10h_posterior_long %>%
        filter(is.element(parameter, c('tau_pre', 'tau_post'))) %>%
        group_by(name, parameter) %>%
        summarise(value.med = median(value),
                  value.25 = quantile(value, probs = 0.25),
                  value.75 = quantile(value, probs = 0.75),
                  value.025 = quantile(value, probs = 0.025),
                  value.975 = quantile(value, probs = 0.975)) %>%
        ungroup(),
    aes(y = 15, x = value.med), size = 5, pch = 21)+
  labs(y = 'Posterior Probability Density', x = 'Value', fill = 'Parameter')

fwd10h_posteriors_tau

ggsave(fwd10h_posteriors_tau,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwd10h_posteriors_tau.png'),
       height = 4, width = 6.5, units = 'in')




fwd10h_posteriors_beta = 
  fwd10h_posterior_long %>%
  filter(is.element(parameter, c('beta[1]', 'beta[2]', 'beta[3]', 
                                 'beta[4]', 'beta[5]', 'beta[6]'))) %>%
  group_by(name, parameter) %>%
  summarise(value.med = median(value),
            value.25 = quantile(value, probs = 0.25),
            value.75 = quantile(value, probs = 0.75),
            value.025 = quantile(value, probs = 0.025),
            value.975 = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  mutate(name = factor(as.character(name),
                       levels = c('compC:timePost','compB:timePost','timePost',
                                  'compC', 'compB', 'Intercept'))) %>%
  ggplot(aes(y = name))+
  geom_errorbarh(aes(xmin = value.025, xmax = value.975), height = 0.3)+
  geom_point(aes(x = value.med), size = 5)+
  labs(y = 'Fixed EFfect Parameter', x = 'Value')+
  theme_minimal()+
  geom_vline(xintercept = 0, lwd = 1, lty = 2)

fwd10h_posteriors_beta

ggsave(fwd10h_posteriors_beta,
       filename = here::here('04-communication',
                  'figures',
                  'manuscript',
                  'fwd10h_posteriors_beta.png'),
       height = 4, width = 6.5, units = 'in')


fwd10h_summary = 
  fwd10h_fit$summary(variables = c('beta', 'alpha_pre', 'alpha_post',
                                     'rho_pre', 'rho_post', 'sigmaPlot', 'kappa')) %>%
  left_join(data.frame(parameter = c('beta[1]', 'beta[2]', 'beta[3]', 
                                     'beta[4]', 'beta[5]', 'beta[6]',
                                     'alpha_pre', 'alpha_post', 'sigmaPlot',
                                     'rho_pre', 'rho_post', 'tau_pre', 'tau_post', 'kappa'),
                       name = c('Intercept', 'compB', 'compC', 'timePost',
                                'compB:timePost', 'compC:timePost',
                                'GP SD (pre)', 'GP SD (post)', 'SD Plot',
                                'GP Range (pre)', 'GP Range (post)', 
                                'GP Noise (pre)', 'GP Noise (post)', 'NB Dispersion')),
            by = c('variable' = 'parameter')) %>%
  select(name, mean, median, sd, mad, q5, q95, rhat, ess_bulk, ess_tail)

write.csv(fwd10h_summary,
          here::here('04-communication',
                     'tables',
                     'fwd10h_gp_summary.csv'),
          row.names = FALSE)

