
#### setup #####################################################################

library(here)
library(tidyverse)

litter_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'litter_training_data.rds'))

duff_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'duff_training_data.rds'))

litterduff_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'litterduff_training_data.rds'))

fwd1h_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'fwd1h_training_data.rds'))

fwd10h_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'fwd10h_training_data.rds'))

fwd100h_data = 
  readRDS(here::here('02-data',
                     '02-for_analysis',
                     'fwd100h_training_data.rds'))


#### litter ####################################################################

names(litter_data)

litter_data$N

litter_data$P

litter_data$L

litter_df = 
  litter_data$X %>%
  as.data.frame() %>%
  mutate(y = litter_data$Y,
         plot_id = litter_data$plot_id,
         location_id = litter_data$location_id,
         pre = litter_data$pre,
         post = litter_data$post,
         comp = 
           ifelse(compB==1,
                  'B',
                  ifelse(compC==1,
                         'C',
                         'A')),
         timestep = 
           ifelse(timePost=='1',
                  'post',
                  'pre')) %>%
  left_join(
    litter_data$coords %>%
      as.data.frame() %>%
      rowid_to_column('location_id')
  )

# X distributions
head(litter_df)

ggplot(litter_df,
       aes(x = comp))+
  geom_bar()

ggplot(litter_df,
       aes(x = timestep))+
  geom_bar()

ggplot(litter_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(litter_df,
       aes(x = location_id))+
  geom_bar()

ggplot(litter_df,
       aes(x = x_rel))+
  geom_histogram()

ggplot(litter_df,
       aes(x = y_rel))+
  geom_histogram()

# Y distribution
ggplot(litter_df,
       aes(x = y))+
  geom_bar()

# XX distribution
litter_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = comp_timestep))+
  geom_bar()

litter_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = location_id))+
  geom_bar()+
  facet_wrap(~comp_timestep)

# XY distribution
litter_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(comp~.)

litter_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~.)

litter_df %>%
  ggplot(aes(y = y, x = as.factor(location_id)))+
  geom_boxplot()

litter_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~comp)

litter_df %>%
  mutate(logy = log(y+min(y[y>0]))) %>%
  ggplot(aes(x = logy))+
  geom_bar()+
  facet_grid(timestep~comp)


#### duff ######################################################################


names(duff_data)

duff_data$N

duff_data$P

duff_data$L

duff_df = 
  duff_data$X %>%
  as.data.frame() %>%
  mutate(y = duff_data$Y,
         plot_id = duff_data$plot_id,
         location_id = duff_data$location_id,
         pre = duff_data$pre,
         post = duff_data$post,
         comp = 
           ifelse(compB==1,
                  'B',
                  ifelse(compC==1,
                         'C',
                         'A')),
         timestep = 
           ifelse(timePost=='1',
                  'post',
                  'pre')) %>%
  left_join(
    duff_data$coords %>%
      as.data.frame() %>%
      rowid_to_column('location_id')
  )

# X distributions
head(duff_df)

ggplot(duff_df,
       aes(x = comp))+
  geom_bar()

ggplot(duff_df,
       aes(x = timestep))+
  geom_bar()

ggplot(duff_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(duff_df,
       aes(x = location_id))+
  geom_bar()

ggplot(duff_df,
       aes(x = x_rel))+
  geom_histogram()

ggplot(duff_df,
       aes(x = y_rel))+
  geom_histogram()

# Y distribution
ggplot(duff_df,
       aes(x = y))+
  geom_bar()

# XX distribution
duff_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = comp_timestep))+
  geom_bar()

duff_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = location_id))+
  geom_bar()+
  facet_wrap(~comp_timestep)

# XY distribution
duff_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(comp~.)

duff_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~.)

duff_df %>%
  ggplot(aes(y = y, x = as.factor(location_id)))+
  geom_boxplot()

duff_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~comp)


#### litterduff ################################################################

names(litterduff_data)

litterduff_data$N

litterduff_data$P

litterduff_data$L

litterduff_df = 
  litterduff_data$X %>%
  as.data.frame() %>%
  mutate(y = litterduff_data$Y,
         plot_id = litterduff_data$plot_id,
         location_id = litterduff_data$location_id,
         pre = litterduff_data$pre,
         post = litterduff_data$post,
         comp = 
           ifelse(compB==1,
                  'B',
                  ifelse(compC==1,
                         'C',
                         'A')),
         timestep = 
           ifelse(timePost=='1',
                  'post',
                  'pre')) %>%
  left_join(
    litterduff_data$coords %>%
      as.data.frame() %>%
      rowid_to_column('location_id')
  )

# X distributions
head(litterduff_df)

ggplot(litterduff_df,
       aes(x = comp))+
  geom_bar()

ggplot(litterduff_df,
       aes(x = timestep))+
  geom_bar()

ggplot(litterduff_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(litterduff_df,
       aes(x = location_id))+
  geom_bar()

ggplot(litterduff_df,
       aes(x = x_rel))+
  geom_histogram()

ggplot(litterduff_df,
       aes(x = y_rel))+
  geom_histogram()

# Y distribution
ggplot(litterduff_df,
       aes(x = y))+
  geom_bar()

# XX distribution
litterduff_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = comp_timestep))+
  geom_bar()

litterduff_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = location_id))+
  geom_bar()+
  facet_wrap(~comp_timestep)

# XY distribution
litterduff_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(comp~.)

litterduff_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~.)

litterduff_df %>%
  ggplot(aes(y = y, x = as.factor(location_id)))+
  geom_boxplot()

litterduff_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~comp)

litterduff_df %>%
  mutate(logy = log(y+min(y[y>0]))) %>%
  ggplot(aes(x = logy))+
  geom_bar()+
  facet_grid(timestep~comp)



#### fwd 1h ####################################################################


names(fwd1h_data)

fwd1h_data$N

fwd1h_data$P

fwd1h_data$L

fwd1h_df = 
  fwd1h_data$X %>%
  as.data.frame() %>%
  mutate(y = fwd1h_data$Y,
         plot_id = fwd1h_data$plot_id,
         location_id = fwd1h_data$location_id,
         pre = fwd1h_data$pre,
         post = fwd1h_data$post,
         comp = 
           ifelse(compB==1,
                  'B',
                  ifelse(compC==1,
                         'C',
                         'A')),
         timestep = 
           ifelse(timePost=='1',
                  'post',
                  'pre')) %>%
  left_join(
    fwd1h_data$coords %>%
      as.data.frame() %>%
      rowid_to_column('location_id')
  )

# X distributions
head(fwd1h_df)

ggplot(fwd1h_df,
       aes(x = comp))+
  geom_bar()

ggplot(fwd1h_df,
       aes(x = timestep))+
  geom_bar()

ggplot(fwd1h_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(fwd1h_df,
       aes(x = location_id))+
  geom_bar()

ggplot(fwd1h_df,
       aes(x = x_rel))+
  geom_histogram()

ggplot(fwd1h_df,
       aes(x = y_rel))+
  geom_histogram()

# Y distribution
ggplot(fwd1h_df,
       aes(x = y))+
  geom_bar()

# XX distribution
fwd1h_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = comp_timestep))+
  geom_bar()

fwd1h_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = location_id))+
  geom_bar()+
  facet_wrap(~comp_timestep)

# XY distribution
fwd1h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(comp~.)

fwd1h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~.)

fwd1h_df %>%
  ggplot(aes(y = y, x = as.factor(location_id)))+
  geom_boxplot()

fwd1h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~comp)

#### fwd 10h ###################################################################

names(fwd10h_data)

fwd10h_data$N

fwd10h_data$P

fwd10h_data$L

fwd10h_df = 
  fwd10h_data$X %>%
  as.data.frame() %>%
  mutate(y = fwd10h_data$Y,
         plot_id = fwd10h_data$plot_id,
         location_id = fwd10h_data$location_id,
         pre = fwd10h_data$pre,
         post = fwd10h_data$post,
         comp = 
           ifelse(compB==1,
                  'B',
                  ifelse(compC==1,
                         'C',
                         'A')),
         timestep = 
           ifelse(timePost=='1',
                  'post',
                  'pre')) %>%
  left_join(
    fwd10h_data$coords %>%
      as.data.frame() %>%
      rowid_to_column('location_id')
  )

# X distributions
head(fwd10h_df)

ggplot(fwd10h_df,
       aes(x = comp))+
  geom_bar()

ggplot(fwd10h_df,
       aes(x = timestep))+
  geom_bar()

ggplot(fwd10h_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(fwd10h_df,
       aes(x = location_id))+
  geom_bar()

ggplot(fwd10h_df,
       aes(x = x_rel))+
  geom_histogram()

ggplot(fwd10h_df,
       aes(x = y_rel))+
  geom_histogram()

# Y distribution
ggplot(fwd10h_df,
       aes(x = y))+
  geom_bar()

# XX distribution
fwd10h_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = comp_timestep))+
  geom_bar()

fwd10h_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = location_id))+
  geom_bar()+
  facet_wrap(~comp_timestep)

# XY distribution
fwd10h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(comp~.)

fwd10h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~.)

fwd10h_df %>%
  ggplot(aes(y = y, x = as.factor(location_id)))+
  geom_boxplot()

fwd10h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~comp)

#### fwd 100h ##################################################################

names(fwd100h_data)

fwd100h_data$N

fwd100h_data$P

fwd100h_data$L

fwd100h_df = 
  fwd100h_data$X %>%
  as.data.frame() %>%
  mutate(y = fwd100h_data$Y,
         plot_id = fwd100h_data$plot_id,
         location_id = fwd100h_data$location_id,
         pre = fwd100h_data$pre,
         post = fwd100h_data$post,
         comp = 
           ifelse(compB==1,
                  'B',
                  ifelse(compC==1,
                         'C',
                         'A')),
         timestep = 
           ifelse(timePost=='1',
                  'post',
                  'pre')) %>%
  left_join(
    fwd100h_data$coords %>%
      as.data.frame() %>%
      rowid_to_column('location_id')
  )

# X distributions
head(fwd100h_df)

ggplot(fwd100h_df,
       aes(x = comp))+
  geom_bar()

ggplot(fwd100h_df,
       aes(x = timestep))+
  geom_bar()

ggplot(fwd100h_df,
       aes(x = plot_id))+
  geom_bar()

ggplot(fwd100h_df,
       aes(x = location_id))+
  geom_bar()

ggplot(fwd100h_df,
       aes(x = x_rel))+
  geom_histogram()

ggplot(fwd100h_df,
       aes(x = y_rel))+
  geom_histogram()

# Y distribution
ggplot(fwd100h_df,
       aes(x = y))+
  geom_bar()

# XX distribution
fwd100h_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = comp_timestep))+
  geom_bar()

fwd100h_df %>%
  mutate(comp_timestep = 
           paste0(comp, '-', timestep)) %>%
  ggplot(aes(x = location_id))+
  geom_bar()+
  facet_wrap(~comp_timestep)

# XY distribution
fwd100h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(comp~.)

fwd100h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~.)

fwd100h_df %>%
  ggplot(aes(y = y, x = as.factor(location_id)))+
  geom_boxplot()

fwd100h_df %>%
  ggplot(aes(x = y))+
  geom_bar()+
  facet_grid(timestep~comp)

