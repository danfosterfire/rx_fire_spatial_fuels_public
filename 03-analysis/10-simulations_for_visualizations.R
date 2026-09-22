
library(here)
library(tidyverse)
library(posterior)
library(cmdstanr)

set.seed(110819)

#### setup #####################################################################

# load constants tables
species_codes = readRDS(here::here('02-data',
                                   '02-intermediate',
                                   'van_wagtendonk',
                                   'species_codes.rds'))

kvals = readRDS(here::here('02-data',
                           '02-intermediate',
                           'van_wagtendonk',
                           'kvals.rds'))


litter_coeffs = readRDS(here::here('02-data',
                                       '02-intermediate',
                                       'van_wagtendonk',
                                       'litter_coeffs.rds'))


QMDcm = readRDS(here::here('02-data',
                           '02-intermediate',
                           'van_wagtendonk',
                           'QMDcm.rds'))

head(QMDcm)

SEC = readRDS(here::here('02-data',
                         '02-intermediate',
                         'van_wagtendonk',
                         'SEC.rds')) 

SG = readRDS(here::here('02-data',
                        '02-intermediate',
                        'van_wagtendonk',
                        'SG.rds'))

sg_1000r = readRDS(here::here('02-data',
                              '02-intermediate',
                              'van_wagtendonk',
                              'vw96_sg_1000r.rds'))

overstory_list = 
  read_csv(here::here('02-data',
                      '01-cleaned',
                      'trees.csv'))

head(overstory_list)

overstory_composition = 
  overstory_list %>%
  filter(COMP == 'A') %>%
  
  # map species wihout coefficients from van wagtendonk's work to 'OTHER' 
  # category, will get the 'All species' coefficient
  mutate(SPECIES = ifelse(!is.element(SPECIES, as.character(species_codes$species_code)),
                      'OTHER',
                      SPECIES)) %>%
  
  mutate(ba_m2 = pi*((dbh_cm/100)/2), # basal area of each tree in m2
         # scaled by 0.05 ha plot size (10x30m + 10x10m + 10x10m)
         ba_m2ha = ba_m2*(1/0.05)) %>%
  
  group_by(COMP, PLOT, SPECIES) %>%
  summarise(ba_m2ha = sum(ba_m2ha)) %>%
  ungroup()  %>%
  
  # fill in 0s for species which werent present on a plot
  complete(nesting(COMP, PLOT), SPECIES) %>%
  mutate(ba_m2ha = ifelse(is.na(ba_m2ha), 0, ba_m2ha)) %>%
  
  # get the total BA/ha on each plot
  left_join(x = .,
            y = 
              group_by(., COMP, PLOT) %>% 
              summarise(total_ba_m2ha = sum(ba_m2ha)) %>%
              ungroup()) %>%
  
  # get the proportion of plot total basal area occoupied by each species
  mutate(pBA = ba_m2ha / total_ba_m2ha) %>%
  
  # aggregate to mortality-class average composition
  group_by(COMP, SPECIES) %>%
  summarise(pBA = mean(pBA)) %>%
  ungroup()


head(overstory_composition)

litter_coeff = 
  overstory_composition %>%
  left_join(litter_coeffs,
            by = c('SPECIES' = 'spp')) %>%
  mutate(weighted = pBA * litter_coeff) %>%
  group_by(COMP) %>%
  summarise(weighted_coeff = sum(weighted)) %>%
  ungroup()

# get BA-weighted SEC (secant of acute angle) for each timelag class
SEC = 
  overstory_composition %>%
  left_join(SEC,
            by = c('SPECIES' = 'spp')) %>%
  select(COMP, SPECIES, pBA, x1h, x10h, x100h, x1000s = x1000h, x1000r = x1000h) %>%
  pivot_longer(c(x1h, x10h, x100h, x1000s, x1000r),
               names_to = 'timelag_class',
               values_to = 'sec',
               names_prefix = 'x') %>%
  mutate(weighted = sec*pBA) %>%
  group_by(COMP, timelag_class) %>%
  summarise(weighted_sec = sum(weighted)) %>%
  ungroup()


# get the BA-weighted average specific gravity for each timelag class
SG = 
  overstory_composition %>%
  left_join(SG,
            by = c('SPECIES' = 'spp')) %>%
  select(COMP, SPECIES, pBA, x1h, x10h, x100h, x1000s) %>%
  pivot_longer(c(x1h, x10h, x100h, x1000s),
               names_to = 'timelag_class',
               values_to = 'sg',
               names_prefix = 'x') %>%
  mutate(weighted = sg*pBA) %>%
  group_by(COMP, timelag_class) %>%
  summarise(weighted_sg = sum(weighted)) %>%
  ungroup() %>%
  
  # and add an entry for 1000h rotten, which doesn't vary by species
  rbind(.,
        data.frame(COMP = c('A'),
                   timelag_class = c('1000r'),
                   weighted_sg = rep(sg_1000r, times = 1)) %>%
          mutate(COMP = as.character(COMP),
                 timelag_class = as.character(timelag_class)))

# get BA-weighted average QMD (cm2) for each FWD timelag class
QMDcm = 
  overstory_composition %>%
  left_join(QMDcm,
            by = c('SPECIES' = 'spp')) %>%
  select(COMP, SPECIES, pBA, x1h, x10h, x100h) %>%
  pivot_longer(c(x1h, x10h, x100h),
               names_to = 'timelag_class',
               values_to = 'qmd_cm',
               names_prefix = 'x') %>%
  mutate(weighted = qmd_cm*pBA) %>%
  group_by(COMP, timelag_class) %>%
  summarise(weighted_qmd = sum(weighted)) %>%
  ungroup()



#### litter ####################################################################

litter_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'litter_training_data.rds'))

litter_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'litter_fit.rds'))

litter_posterior = 
  as_draws_df(litter_fit$draws())

litter_medians = 
  litter_posterior %>%
  select(-lp__, -contains('zGP_'),-contains('z_plot'), -contains('plot_effect'), -`.iteration`, -`.chain`, -`.draw`) %>%
  as_tibble() %>%
  summarise(across(everything(), ~median(.x)))

litter_medians

sim_points = 
  expand_grid(x_coord = seq(0, 60, 1),
              y_coord = seq(0, 60, 1)) %>%
  rowid_to_column('location_id') %>%
  crossing(expand_grid(stand_id = c('a', 'b', 'c'),
                  pre_post = c('pre', 'post'))) %>%
  mutate(pre_post = factor(pre_post, levels = c('pre', 'post')),
         stand_id = factor(stand_id, levels = c('a', 'b', 'c'))) %>%
  mutate(intercept = 1,
         stand_b = as.numeric(stand_id == 'b'),
         stand_c = as.numeric(stand_id == 'c'),
         pre = as.numeric(pre_post == 'pre'),
         post = as.numeric(pre_post == 'post'),
         stand_b_post = as.numeric(stand_id == 'b' & pre_post == 'post'),
         stand_c_post = as.numeric(stand_id == 'c' & pre_post == 'post'))

litter_sim_results = 
  do.call('bind_rows',
          lapply(X = c('a', 'b', 'c'),
                 FUN = function(stand_id){
                   
                   
                      distance_matrix = 
                        sim_points %>%
                        filter(stand_id == !!stand_id & pre_post == 'pre') %>%
                        arrange(location_id) %>%
                        select(x_coord, y_coord) %>%
                        as.matrix() %>%
                        dist() %>%
                        as.matrix()
                      
                      K_pre = 
                       as.numeric(litter_medians$alpha_pre) * 
                       exp(-(distance_matrix**2)/(2*(as.numeric(litter_medians$rho_pre)**2)))
                      
                      
                      K_pre = K_pre + diag(nrow = nrow(K_pre),
                                          ncol = ncol(K_pre),
                                          x = 1e-9)
                      KL_pre = t(chol(K_pre))
                      
                      GP_pre = KL_pre %*% rnorm(mean = 0, sd = 1, n = nrow(distance_matrix)) %>%
                        as.numeric()
                      
                      
                      
                      K_post = 
                       as.numeric(litter_medians$alpha_post) * 
                       exp(-(distance_matrix**2)/(2*(as.numeric(litter_medians$rho_post)**2)))
                      
                      
                      K_post = K_post + diag(nrow = nrow(K_post),
                                          ncol = ncol(K_post),
                                          x = 1e-9)
                      KL_post = t(chol(K_post))
                      
                      GP_post = KL_post %*% rnorm(mean = 0, sd = 1, n = nrow(distance_matrix)) %>%
                        as.numeric()
                      
                      GP_combined = c(GP_pre, GP_post)
                      
                      results = 
                        sim_points %>%
                        filter(stand_id == !!stand_id) %>%
                        mutate(
                         XB = litter_medians$`beta[1]`*intercept + 
                           litter_medians$`beta[2]`*stand_b +
                           litter_medians$`beta[3]`*stand_c +
                           litter_medians$`beta[4]`*post +
                           litter_medians$`beta[5]`*stand_b_post +
                           litter_medians$`beta[6]`*stand_c_post) %>%
                        arrange(pre_post, location_id) %>%
                        bind_cols(
                          tibble(GP_combined)
                        ) %>%
                        mutate(logmu = XB+GP_combined)
                      
                      Y = rnbinom(n = nrow(results),
                                 size = as.numeric(litter_medians$kappa),
                                 mu = exp(results$logmu))
                      
                      results = results %>% bind_cols(tibble(Y))
                      
                      return(results)
                   
                   
                 }))

#### fwd 1h ####################################################################


fwd1h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'fwd1h_fit.rds'))

fwd1h_posterior = 
  as_draws_df(fwd1h_fit$draws())

fwd1h_medians = 
  fwd1h_posterior %>%
  select(-lp__, -contains('zGP_'),-contains('z_plot'), -contains('plot_effect'), -`.iteration`, -`.chain`, -`.draw`) %>%
  as_tibble() %>%
  summarise(across(everything(), ~median(.x)))

fwd1h_medians


fwd1h_sim_results = 
  do.call('bind_rows',
          lapply(X = c('a', 'b', 'c'),
                 FUN = function(stand_id){
                   
                   
                      distance_matrix = 
                        sim_points %>%
                        filter(stand_id == !!stand_id & pre_post == 'pre') %>%
                        arrange(location_id) %>%
                        select(x_coord, y_coord) %>%
                        as.matrix() %>%
                        dist() %>%
                        as.matrix()
                      
                      K_pre = 
                       as.numeric(fwd1h_medians$alpha_pre) * 
                       exp(-(distance_matrix**2)/(2*(as.numeric(fwd1h_medians$rho_pre)**2)))
                      
                      
                      K_pre = K_pre + diag(nrow = nrow(K_pre),
                                          ncol = ncol(K_pre),
                                          x = fwd1h_medians$tau_pre)
                      KL_pre = t(chol(K_pre))
                      
                      GP_pre = KL_pre %*% rnorm(mean = 0, sd = 1, n = nrow(distance_matrix)) %>%
                        as.numeric()
                      
                      
                      
                      K_post = 
                       as.numeric(fwd1h_medians$alpha_post) * 
                       exp(-(distance_matrix**2)/(2*(as.numeric(fwd1h_medians$rho_post)**2)))
                      
                      
                      K_post = K_post + diag(nrow = nrow(K_post),
                                          ncol = ncol(K_post),
                                          x = fwd1h_medians$tau_post)
                      KL_post = t(chol(K_post))
                      
                      GP_post = KL_post %*% rnorm(mean = 0, sd = 1, n = nrow(distance_matrix)) %>%
                        as.numeric()
                      
                      GP_combined = c(GP_pre, GP_post)
                      
                      results = 
                        sim_points %>%
                        filter(stand_id == !!stand_id) %>%
                        mutate(
                         XB = fwd1h_medians$`beta[1]`*intercept + 
                           fwd1h_medians$`beta[2]`*stand_b +
                           fwd1h_medians$`beta[3]`*stand_c +
                           fwd1h_medians$`beta[4]`*post +
                           fwd1h_medians$`beta[5]`*stand_b_post +
                           fwd1h_medians$`beta[6]`*stand_c_post) %>%
                        arrange(pre_post, location_id) %>%
                        bind_cols(
                          tibble(GP_combined)
                        ) %>%
                        mutate(logmu = XB+GP_combined)
                      
                      Y = rnbinom(n = nrow(results),
                                 size = as.numeric(fwd1h_medians$kappa),
                                 mu = exp(results$logmu))
                      
                      results = results %>% bind_cols(tibble(Y))
                      
                      return(results)
                   
                   
                 }))

#### fwd 10h ###################################################################



fwd10h_fit = 
  readRDS(here::here('02-data',
                     '03-results',
                     'real_fits',
                     'fwd10h_fit.rds'))

fwd10h_posterior = 
  as_draws_df(fwd10h_fit$draws())

fwd10h_medians = 
  fwd10h_posterior %>%
  select(-lp__, -contains('zGP_'),-contains('z_plot'), -contains('plot_effect'), -`.iteration`, -`.chain`, -`.draw`) %>%
  as_tibble() %>%
  summarise(across(everything(), ~median(.x)))

fwd10h_medians


fwd10h_sim_results = 
  do.call('bind_rows',
          lapply(X = c('a', 'b', 'c'),
                 FUN = function(stand_id){
                   
                   
                      distance_matrix = 
                        sim_points %>%
                        filter(stand_id == !!stand_id & pre_post == 'pre') %>%
                        arrange(location_id) %>%
                        select(x_coord, y_coord) %>%
                        as.matrix() %>%
                        dist() %>%
                        as.matrix()
                      
                      K_pre = 
                       as.numeric(fwd10h_medians$alpha_pre) * 
                       exp(-(distance_matrix**2)/(2*(as.numeric(fwd10h_medians$rho_pre)**2)))
                      
                      
                      K_pre = K_pre + diag(nrow = nrow(K_pre),
                                          ncol = ncol(K_pre),
                                          x = fwd10h_medians$tau_pre)
                      KL_pre = t(chol(K_pre))
                      
                      GP_pre = KL_pre %*% rnorm(mean = 0, sd = 1, n = nrow(distance_matrix)) %>%
                        as.numeric()
                      
                      
                      
                      K_post = 
                       as.numeric(fwd10h_medians$alpha_post) * 
                       exp(-(distance_matrix**2)/(2*(as.numeric(fwd10h_medians$rho_post)**2)))
                      
                      
                      K_post = K_post + diag(nrow = nrow(K_post),
                                          ncol = ncol(K_post),
                                          x = fwd10h_medians$tau_post)
                      KL_post = t(chol(K_post))
                      
                      GP_post = KL_post %*% rnorm(mean = 0, sd = 1, n = nrow(distance_matrix)) %>%
                        as.numeric()
                      
                      GP_combined = c(GP_pre, GP_post)
                      
                      results = 
                        sim_points %>%
                        filter(stand_id == !!stand_id) %>%
                        mutate(
                         XB = fwd10h_medians$`beta[1]`*intercept + 
                           fwd10h_medians$`beta[2]`*stand_b +
                           fwd10h_medians$`beta[3]`*stand_c +
                           fwd10h_medians$`beta[4]`*post +
                           fwd10h_medians$`beta[5]`*stand_b_post +
                           fwd10h_medians$`beta[6]`*stand_c_post) %>%
                        arrange(pre_post, location_id) %>%
                        bind_cols(
                          tibble(GP_combined)
                        ) %>%
                        mutate(logmu = XB+GP_combined)
                      
                      Y = rnbinom(n = nrow(results),
                                 size = as.numeric(fwd10h_medians$kappa),
                                 mu = exp(results$logmu))
                      
                      results = results %>% bind_cols(tibble(Y))
                      
                      return(results)
                   
                   
                 }))





#### plotting sim stands #######################################################

head(litter_sim_results)

head(litter_coeff)

litter_sim_results_loading = 
  litter_sim_results %>%
  mutate(COMP = 'A') %>%
  left_join(litter_coeff,
            by = c('COMP' = 'COMP')) %>%
  # coefficient goes from depth in cm to biomass in kg/m2, so multiply by 10 
  # to get biomass in Mg/ha
  mutate(litter_kgm2 = Y * weighted_coeff)


litter_sim_plot = 
  ggplot(
    data = litter_sim_results_loading %>%
      filter(stand_id == 'a') %>%
      mutate(kgm2_binned = 
               cut(litter_kgm2, 
                   c(0,0.5, 1, 2, 4, 8, 16, Inf), 
                   include.lowest = TRUE,
                   labels = c('0-0.5', '0.5-1', '1-2', '2-4', '4-8', '8-16', '>16'))) %>%
      mutate(
        stand_id = 
          case_when(stand_id == 'a' ~ 'A',
                    stand_id == 'b' ~ 'B',
                    stand_id == 'c' ~ 'C'),
        pre_post = 
          case_when(pre_post == 'pre' ~ 'Preburn',
                    pre_post == 'post' ~ 'Postburn')
      ) %>%
      mutate(pre_post = factor(pre_post, levels = c('Preburn', 'Postburn'))),
    aes(x = x_coord, y = y_coord, fill = kgm2_binned)
  )+
  geom_tile()+
  scale_fill_viridis_d()+
  facet_grid(.~pre_post)+
  coord_fixed()+
  theme_bw()+
  labs(x = 'Easting (m)',
       y = 'Northing (m)',
       fill = expression(Litter~fuel~load~(kg/m^2)))


litter_sim_plot

head(QMDcm)
head(SG)
head(SEC)
fwd1h_sim_results_loading = 
  fwd1h_sim_results %>%
  mutate(COMP = 'A',
         timelag_class = '1h') %>%
  # add columns for slope, QMD, SEC, SG, and transect length, and the conversion 
  # constant k
  left_join(., 
            y = QMDcm) %>%
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>%
  mutate(slp_c = 1) %>%
  mutate(transect_length_m = 1,
         k = 1.234) %>%
  
  # use brown's equations to estimate fuel load
  mutate(mgha = 
           (k*weighted_qmd*weighted_sec*slp_c*weighted_sg*Y)/
           transect_length_m,
         kgm2 = mgha/10) 

hist(fwd1h_sim_results_loading$kgm2)

fwd1h_sim_plot = 
  ggplot(
    data = fwd1h_sim_results_loading %>%
      filter(stand_id == 'a') %>%
      mutate(kgm2_binned = 
               cut(kgm2, 
                   c(0, 0.25, 0.5, 1, 2, 4, Inf), 
                   include.lowest = TRUE,
                   labels = c('0-0.25', '0.25-0.5', '0.5-1', '1-2', '2-4',  '>4'))) %>%
      mutate(
        stand_id = 
          case_when(stand_id == 'a' ~ 'A',
                    stand_id == 'b' ~ 'B',
                    stand_id == 'c' ~ 'C'),
        pre_post = 
          case_when(pre_post == 'pre' ~ 'Preburn',
                    pre_post == 'post' ~ 'Postburn')
      ) %>%
      mutate(pre_post = factor(pre_post, levels = c('Preburn', 'Postburn'))),
    aes(x = x_coord, y = y_coord, fill = kgm2_binned)
  )+
  geom_tile()+
  scale_fill_viridis_d(option = 'B', drop = FALSE)+
  facet_grid(.~pre_post)+
  coord_fixed()+
  theme_bw()+
  labs(x = 'Easting (m)',
       y = 'Northing (m)',
       fill = expression(FWD1h~fuel~load~(kg/m^2)))

fwd1h_sim_plot

fwd10h_sim_results_loading = 
  fwd10h_sim_results %>%
  mutate(COMP = 'A',
         timelag_class = '10h') %>%
  # add columns for slope, QMD, SEC, SG, and transect length, and the conversion 
  # constant k
  left_join(., 
            y = QMDcm) %>%
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>%
  mutate(slp_c = 1) %>%
  mutate(transect_length_m = 1,
         k = 1.234) %>%
  
  # use brown's equations to estimate fuel load
  mutate(mgha = 
           (k*weighted_qmd*weighted_sec*slp_c*weighted_sg*Y)/
           transect_length_m,
         kgm2 = mgha/10) 


fwd10h_sim_plot = 
  ggplot(
    data = fwd10h_sim_results_loading %>%
      filter(stand_id == 'a') %>%
      mutate(kgm2_binned = 
               cut(kgm2, 
                   c(0, 0.25, 0.5, 1, 2, 4, Inf), 
                   include.lowest = TRUE,
                   labels = c('0-0.25','0.25-0.5', '0.5-1', '1-2', '2-4', '>4'))) %>%
      mutate(
        stand_id = 
          case_when(stand_id == 'a' ~ 'A',
                    stand_id == 'b' ~ 'B',
                    stand_id == 'c' ~ 'C'),
        pre_post = 
          case_when(pre_post == 'pre' ~ 'Preburn',
                    pre_post == 'post' ~ 'Postburn')
      ) %>%
      mutate(pre_post = factor(pre_post, levels = c('Preburn', 'Postburn'))),
    aes(x = x_coord, y = y_coord, fill = kgm2_binned)
  )+
  geom_tile()+
  scale_fill_viridis_d(option = 'B', drop = FALSE)+
  facet_grid(.~pre_post)+
  coord_fixed()+
  theme_bw()+
  labs(x = 'Easting (m)',
       y = 'Northing (m)',
       fill = expression(FWD10h~fuel~load~(kg/m^2)))

fwd10h_sim_plot

library(cowplot)

all_sims_plot = 
  cowplot::plot_grid(
    litter_sim_plot, fwd1h_sim_plot, fwd10h_sim_plot,
    nrow = 3,
    align = 'v'
  )


all_sims_plot


ggsave(
  all_sims_plot,
  filename = here::here('04-communication',
             'figures',
             'simulations_plot.png'),
  height = 10, width = 7, units = 'in'
)


#### error statistics ##########################################################

prediction_summary_stats = 
  readRDS(
        here::here('02-data', '03-results', 'prediction_summary_stats.rds'))


retrodiction_summary_stats = 
  readRDS(
        here::here('02-data', '03-results', 'retrodiction_summary_stats.rds'))

head(prediction_summary_stats)
head(retrodiction_summary_stats)

error_stats = 
  tibble(
    statistic = 
      c('Litter_MAE_Validation',
        'Litter_Mean_Validation',
        'Litter_wMAPE_Validation',
        'FWD1h_MAE_Validation',
        'FWD1h_Mean_Validation',
        'FWD1h_wMAPE_Validation',
        'FWD10h_MAE_Validation',
        'FWD10h_Mean_Validation',
        'FWD10h_wMAPE_Validation',
        'FWD100h_MAE_Validation',
        'FWD100h_Mean_Validation',
        'FWD100h_wMAPE_Validation',
        'Litter_MAE_Training',
        'Litter_Mean_Training',
        'Litter_wMAPE_Training',
        'FWD1h_MAE_Training',
        'FWD1h_Mean_Training',
        'FWD1h_wMAPE_Training',
        'FWD10h_MAE_Training',
        'FWD10h_Mean_Training',
        'FWD10h_wMAPE_Training',
        'FWD100h_MAE_Training',
        'FWD100h_Mean_Training',
        'FWD100h_wMAPE_Training'
        ),
    value = 
      c(prediction_summary_stats, retrodiction_summary_stats)
  ) %>%
  
  separate_wider_delim(statistic, 
                       names = c('Fuel Component', 'Statistic', 'Dataset'),
                       delim = '_') %>%
  select(Dataset, `Fuel Component`, Statistic, value)


error_stats

error_stats %>%
  write_csv(here::here('04-communication','tables', 'error_stats.csv'))


