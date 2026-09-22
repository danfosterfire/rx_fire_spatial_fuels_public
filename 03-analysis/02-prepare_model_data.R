

#### setup #####################################################################


library(here)
library(tidyverse)

set.seed(110819)

unique_plots = 
  read_csv(here::here('02-data', '01-cleaned', 'litter_duff.csv')) %>%
  select(comp = COMP, plot_id = PLOT) %>%
  group_by(comp, plot_id) %>%
  summarise() %>%
  ungroup() %>%
  crossing(tibble(timestep = c('pre', 'post'))) %>%
  rowid_to_column('plot_id.i') %>%
  mutate(comp_id.i = as.integer(as.factor(comp)))

unique_plots$comp_id.i


#### litter and duff ###########################################################

# load  the data, do some basic munging to get pre and post in one table
litterduff = 
  read_csv(here::here('02-data', '01-cleaned', 'litter_duff.csv')) %>%
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT,
         az = AZ, location_m = Displacement_m, litter_cm = LITTER_cm, duff_cm = DUFF_cm) %>% 
  mutate(timestep = 'pre') %>%
  bind_rows(read_csv(here::here('02-data', '01-cleaned', 'litter_duff_post.csv')) %>%
              mutate(timestep = 'post')) %>%
  mutate(litterduff_cm = litter_cm+duff_cm) %>%
  
  # make a version clamped to the nearest cm
  mutate(litterduff_cm0 = round(litterduff_cm, 0),
         litter_cm0 = round(litter_cm, 0),
         duff_cm0 = round(duff_cm, 0)) %>%
  # get relative spatial coordinates
  mutate(x_rel = round(location_m*sin(az*pi/180),2),
         y_rel = round(location_m*cos(az*pi/180),2)) %>%
  
  # drop NA observations
  filter(!is.na(litterduff_cm0) & !is.na(litter_cm0) & !is.na(duff_cm0)) %>%
  
  # compXtime = group
  mutate(group_id = factor(paste0(comp,'_',timestep),
                              levels = c('A_pre','B_pre', 'C_pre', 'A_post',
                                         'B_post', 'C_post')),
         group_id.i = as.integer(group_id),
         timestep_pre = as.integer(timestep == 'pre'),
         timestep_post = as.integer(timestep == 'post'),
         comp_A = as.integer(comp == 'A'),
         comp_B = as.integer(comp == 'B'),
         comp_C = as.integer(comp == 'C'),
         pre_A = timestep_pre * comp_A,
         pre_B = timestep_pre * comp_B,
         pre_C = timestep_pre * comp_C,
         post_A = timestep_post*comp_A,
         post_B = timestep_post*comp_B,
         post_C = timestep_post*comp_C,
         intercept = 1)
  
# unique relative coordinates within a plot  
coords_litterduff = 
  litterduff %>%
  select(x_rel,y_rel) %>%
  group_by(x_rel, y_rel) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('location_id')

# add in the location IDs, timestep binary flags, and fixed effect covariates
litterduff = 
  litterduff %>%
  left_join(coords_litterduff) %>%
  left_join(unique_plots) %>%
  rowid_to_column('obs_id')

# split into training and validation
set.seed(110819)
litterduff_training = 
  litterduff %>%
  sample_frac(replace = FALSE, size = 0.9)

litterduff_validation = 
  litterduff %>%
  filter(!is.element(obs_id,
                     litterduff_training$obs_id))


fixeff_params = c('intercept', 'timestep_pre', 'timestep_post',
                  'comp_A', 'comp_B', 'comp_C', 'pre_A', 'pre_B', 'pre_C',
                  'post_A', 'post_B', 'post_C')

duff_training_data = 
  list(N = nrow(litterduff_training),
       P = nrow(unique_plots),
       L = nrow(coords_litterduff),
       Y = round(litterduff_training$duff_cm, 0),
       Q = length(fixeff_params),
       plot_id = litterduff_training$plot_id.i,
       X = litterduff_training %>%
         select(all_of(c('plot_id.i', fixeff_params))) %>%
         distinct() %>%
         arrange(plot_id.i) %>%
         select(-plot_id.i) %>%
         as.matrix(),
       location_id = litterduff_training$location_id,
       coords = coords_litterduff[,c('x_rel', 'y_rel')] %>% as.matrix())

litter_training_data = 
  list(N = nrow(litterduff_training),
       P = nrow(unique_plots),
       L = nrow(coords_litterduff),
       Y = round(litterduff_training$litter_cm, 0),
       Q = length(fixeff_params),
       plot_id = litterduff_training$plot_id.i,
       X = litterduff_training %>%
         select(all_of(c('plot_id.i', fixeff_params))) %>%
         distinct() %>%
         arrange(plot_id.i) %>%
         select(-plot_id.i) %>%
         as.matrix(),
       location_id = litterduff_training$location_id,
       coords = coords_litterduff[,c('x_rel', 'y_rel')] %>% as.matrix())


duff_validation_data = 
  list(N = nrow(litterduff_validation),
       P = nrow(unique_plots),
       L = nrow(coords_litterduff),
       Y = round(litterduff_validation$duff_cm, 0),
       Q = length(fixeff_params),
       plot_id = litterduff_validation$plot_id.i,
       X = litterduff_validation %>%
         select(all_of(c('plot_id.i', fixeff_params))) %>%
         distinct() %>%
         arrange(plot_id.i) %>%
         select(-plot_id.i) %>%
         as.matrix(),
       location_id = litterduff_validation$location_id,
       coords = coords_litterduff[,c('x_rel', 'y_rel')] %>% as.matrix())

litter_validation_data = 
  list(N = nrow(litterduff_validation),
       P = nrow(unique_plots),
       L = nrow(coords_litterduff),
       Y = round(litterduff_validation$litter_cm, 0),
       Q = length(fixeff_params),
       plot_id = litterduff_validation$plot_id.i,
       X = litterduff_validation %>%
         select(all_of(c('plot_id.i', fixeff_params))) %>%
         distinct() %>%
         arrange(plot_id.i) %>%
         select(-plot_id.i) %>%
         as.matrix(),
       location_id = litterduff_validation$location_id,
       coords = coords_litterduff[,c('x_rel', 'y_rel')] %>% as.matrix())




#### fine woody debris #########################################################

fwd = 
  read_csv(here::here('02-data', '01-cleaned', 'fwd_tallies.csv')) %>%
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT,
         az = AZ, location_m = Displacement_m, a1h, a10h, a100h, b1h, b10h, b100h) %>%
  mutate(timestep = 'pre') %>%
  bind_rows(read_csv(here::here('02-data', '01-cleaned', 'fwd_tallies_post.csv')) %>%
              mutate(timestep = 'post')) %>%
  
  # get relative spatial coordinates
  mutate(x_rel = round(location_m*sin(az*pi/180),2),
         y_rel = round(location_m*cos(az*pi/180),2)) %>%
  
  # reformat tallies to longwise
  pivot_longer(cols = c(a1h, a10h, a100h, b1h, b10h, b100h),
               names_to = 'subsample_id',
               values_to = 'count') %>%
  mutate(timelag_class = gsub(x = subsample_id,
                              pattern = 'a|b',
                              replacement = ''),
         subsample = gsub(x = subsample_id,
                          pattern = '1h|10h|100h',
                          replacement = ''))

coords_fwd = 
  fwd %>%
  group_by(x_rel, y_rel) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('location_id')

fwd = 
  fwd %>%
  left_join(coords_fwd) %>%
  left_join(unique_plots) %>%
  mutate(pre = as.integer(timestep=='pre'),
         post = as.integer(timestep=='post'),
         intercept = 1,
         compB = as.integer(comp=='B'),
         compC = as.integer(comp=='C'),
         timePost = as.integer(timestep=='post'), 
         compBtimePost = as.integer(timestep=='post'&comp=='B'),
         compCtimePost = as.integer(timestep=='post'&comp=='C'))

head(fwd)

fwd_1h = 
  fwd %>%
  filter(timelag_class == '1h' & !is.na(count)) %>%
  rowid_to_column('obs_id')

fwd_10h = 
  fwd %>%
  filter(timelag_class == '10h' & !is.na(count)) %>%
  rowid_to_column('obs_id')

fwd_100h = 
  fwd %>%
  filter(timelag_class == '100h' & !is.na(count)) %>%
  rowid_to_column('obs_id')


fwd_1h_training = 
  fwd_1h %>%
  sample_frac(replace = FALSE,
              size = 0.9)

fwd_1h_validation = 
  fwd_1h %>%
  filter(!is.element(obs_id,
                     fwd_1h_training$obs_id))


fwd_10h_training = 
  fwd_10h %>%
  sample_frac(replace = FALSE,
              size = 0.9)

fwd_10h_validation = 
  fwd_10h %>%
  filter(!is.element(obs_id,
                     fwd_10h_training$obs_id))


fwd_100h_training = 
  fwd_100h %>%
  sample_frac(replace = FALSE,
              size = 0.9)

fwd_100h_validation = 
  fwd_100h %>%
  filter(!is.element(obs_id,
                     fwd_100h_training$obs_id))


ggplot(data = fwd_1h_training,
       aes(x = count))+
  geom_bar()

ggplot(data = fwd_1h_training,
       aes(x = count))+
  geom_bar()+
  facet_grid(timestep~comp)

ggplot(data = fwd_10h_training,
       aes(x = count))+
  geom_bar()

ggplot(data = fwd_10h_training,
       aes(x = count))+
  geom_bar()+
  facet_grid(timestep~comp)

ggplot(data = fwd_100h_training,
       aes(x = count))+
  geom_bar()

ggplot(data = fwd_100h_training,
       aes(x = count))+
  geom_bar()+
  facet_grid(timestep~comp)+
  scale_x_continuous(limits = c(0,20))


fwd1h_training_data = 
  list(N = nrow(fwd_1h_training),
       P = nrow(unique_plots),
       L = nrow(coords_fwd),
       Y = fwd_1h_training$count,
       plot_id = fwd_1h_training$plot_id.i,
       location_id = fwd_1h_training$location_id,
       X = fwd_1h_training[,c('intercept', 'compB', 'compC', 'timePost',
                         'compBtimePost', 'compCtimePost')] %>%
         as.matrix(),
       pre = fwd_1h_training$pre,
       post = fwd_1h_training$post,
       coords = coords_fwd[,c('x_rel', 'y_rel')] %>% as.matrix())


fwd1h_validation_data = 
  list(N = nrow(fwd_1h_validation),
       P = nrow(unique_plots),
       L = nrow(coords_fwd),
       Y = fwd_1h_validation$count,
       plot_id = fwd_1h_validation$plot_id.i,
       location_id = fwd_1h_validation$location_id,
       X = fwd_1h_validation[,c('intercept', 'compB', 'compC', 'timePost',
                         'compBtimePost', 'compCtimePost')] %>%
         as.matrix(),
       pre = fwd_1h_validation$pre,
       post = fwd_1h_validation$post,
       coords = coords_fwd[,c('x_rel', 'y_rel')] %>% as.matrix())

fwd10h_training_data = 
  list(N = nrow(fwd_10h_training),
       P = nrow(unique_plots),
       L = nrow(coords_fwd),
       Y = fwd_10h_training$count,
       plot_id = fwd_10h_training$plot_id.i,
       location_id = fwd_10h_training$location_id,
       X = fwd_10h_training[,c('intercept', 'compB', 'compC', 'timePost',
                         'compBtimePost', 'compCtimePost')] %>%
         as.matrix(),
       pre = fwd_10h_training$pre,
       post = fwd_10h_training$post,
       coords = coords_fwd[,c('x_rel', 'y_rel')] %>% as.matrix())


fwd10h_validation_data = 
  list(N = nrow(fwd_10h_validation),
       P = nrow(unique_plots),
       L = nrow(coords_fwd),
       Y = fwd_10h_validation$count,
       plot_id = fwd_10h_validation$plot_id.i,
       location_id = fwd_10h_validation$location_id,
       X = fwd_10h_validation[,c('intercept', 'compB', 'compC', 'timePost',
                         'compBtimePost', 'compCtimePost')] %>%
         as.matrix(),
       pre = fwd_10h_validation$pre,
       post = fwd_10h_validation$post,
       coords = coords_fwd[,c('x_rel', 'y_rel')] %>% as.matrix())

fwd100h_training_data = 
  list(N = nrow(fwd_100h_training),
       P = nrow(unique_plots),
       L = nrow(coords_fwd),
       Y = fwd_100h_training$count,
       plot_id = fwd_100h_training$plot_id.i,
       location_id = fwd_100h_training$location_id,
       X = fwd_100h_training[,c('intercept', 'compB', 'compC', 'timePost',
                         'compBtimePost', 'compCtimePost')] %>%
         as.matrix(),
       pre = fwd_100h_training$pre,
       post = fwd_100h_training$post,
       coords = coords_fwd[,c('x_rel', 'y_rel')] %>% as.matrix())


fwd100h_validation_data = 
  list(N = nrow(fwd_100h_validation),
       P = nrow(unique_plots),
       L = nrow(coords_fwd),
       Y = fwd_100h_validation$count,
       plot_id = fwd_100h_validation$plot_id.i,
       location_id = fwd_100h_validation$location_id,
       X = fwd_100h_validation[,c('intercept', 'compB', 'compC', 'timePost',
                         'compBtimePost', 'compCtimePost')] %>%
         as.matrix(),
       pre = fwd_100h_validation$pre,
       post = fwd_100h_validation$post,
       coords = coords_fwd[,c('x_rel', 'y_rel')] %>% as.matrix())



#### particle diameters ########################################################

particle_diams = 
  read.csv(here::here('02-data',
                      '01-cleaned',
                      'fwd_diams.csv')) %>%
  select(comp = COMP,
         plot_id = PLOT,
         transect_id = TRANSECT,
         az = AZ,
         diam_mm = DIAM_mm) %>%
  mutate(timestep = 'pre') %>%
  bind_rows(
    read.csv(here::here('02-data',
                        '01-cleaned',
                        'fwd_diams_post.csv')) %>%
      mutate(timestep = 'post')
  ) %>%
  mutate(component = 
           ifelse(diam_mm < 6.4,
                  '1-hour',
                  ifelse(diam_mm >= 6.4 & diam_mm < 25.4,
                         '10-hour',
                         ifelse(diam_mm >= 25.4 & diam_mm < 76.2,
                                '100-hour',
                                NA))))

particle_diams %>%
  ggplot(aes(x = diam_mm, fill = timestep))+
  geom_density(color = NA, alpha = 0.6)+
  facet_wrap(~component, scales = 'free')+
  theme_minimal()

particle_diams %>%
  ggplot(aes(x = diam_mm, y = timestep))+
  geom_boxplot()+
  facet_grid(.~component, scales = 'free')+
  theme_minimal()

particle_diams %>%
  ggplot(aes(y = -diam_mm, x = timestep))+
  geom_violin()+
  geom_jitter(width = 0.1, size = 0)+
  facet_grid(component ~., scales = 'free')

particle_diams %>%
  ggplot(aes(x = diam_mm, fill = timestep))+
  geom_density(color = NA, alpha = 0.6)+
  geom_vline(lty = 2, xintercept = 6.4)+
  geom_vline(lty = 2, xintercept = 25.4)+
  geom_vline(lty = 2, xintercept = 76.2)+
  theme_minimal()


head(particle_diams)

head(read.csv(here::here('02-data',
                         '01-cleaned',
                         'fwd_diams_post.csv')))



#### write results #############################################################

saveRDS(litter_training_data,
        here::here('02-data', '02-for_analysis', 'litter_training_data.rds'))
saveRDS(duff_training_data,
        here::here('02-data', '02-for_analysis', 'duff_training_data.rds'))
saveRDS(fwd1h_training_data,
        here::here('02-data', '02-for_analysis', 'fwd1h_training_data.rds'))
saveRDS(fwd10h_training_data,
        here::here('02-data', '02-for_analysis', 'fwd10h_training_data.rds'))
saveRDS(fwd100h_training_data,
        here::here('02-data', '02-for_analysis', 'fwd100h_training_data.rds'))



saveRDS(litter_validation_data,
        here::here('02-data', '02-for_analysis', 'litter_validation_data.rds'))
saveRDS(duff_validation_data,
        here::here('02-data', '02-for_analysis', 'duff_validation_data.rds'))
saveRDS(fwd1h_validation_data,
        here::here('02-data', '02-for_analysis', 'fwd1h_validation_data.rds'))
saveRDS(fwd10h_validation_data,
        here::here('02-data', '02-for_analysis', 'fwd10h_validation_data.rds'))
saveRDS(fwd100h_validation_data,
        here::here('02-data', '02-for_analysis', 'fwd100h_validation_data.rds'))


