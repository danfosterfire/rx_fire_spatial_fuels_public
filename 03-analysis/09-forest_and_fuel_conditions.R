
library(here)
library(tidyverse)


#### setup #####################################################################

set.seed(110819)

unique_plots = 
  read_csv(here::here('02-data', '01-cleaned', 'litter_duff.csv')) %>%
  select(comp = COMP, plot_id = PLOT) %>%
  group_by(comp, plot_id) %>%
  summarise() %>%
  ungroup() %>%
  rowid_to_column('plot_id.i')


# get biomass estimates for each category in each timestep on each plot, 
# calculate the total deltas, convert into a % of total consumption for each category
# on each plot

# store the constants tables reported from the van wagtendonk papers and 
# dump them into some RDS files; this is a modified version of the RFuels 
# script (which dumps data into RDS files instead of constructing it 
# implicitly as part of the R package)
source(here::here('00-R', 'build_constants_tables.R'))
rm(list = ls())

# load data

# contains both the pre- and post-burn trees data
trees = read_csv(here::here('02-data',
                            '01-cleaned',
                            'trees_post.csv'))

litter_duff = 
  read_csv(here::here('02-data',
                      '01-cleaned',
                      'litter_duff.csv')) %>%
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT, 
         az = AZ, location_m = Displacement_m,
         litter_cm = LITTER_cm, duff_cm = DUFF_cm) %>%
  mutate(timestep = 'pre') %>%
  bind_rows(read_csv(here::here('02-data',
                                '01-cleaned',
                                'litter_duff_post.csv')) %>%
              select(-comment) %>%
              mutate(timestep = 'post'))

head(litter_duff)

fwd_tallies = 
  read_csv(here::here('02-data','01-cleaned',
                      'fwd_tallies.csv')) %>%
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT,
         az = AZ, location_m = Displacement_m, 
         a1h, a10h, a100h, b1h, b10h, b100h) %>%
  mutate(timestep = 'pre') %>%
  bind_rows(
    read_csv(here::here('02-data', '01-cleaned', 
                        'fwd_tallies_post.csv')) %>%
      select(-comment) %>%
      mutate(timestep = 'post')
  )

head(fwd_tallies)

cwd = 
  read_csv(here::here('02-data','01-cleaned','cwd.csv')) %>%
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT, az = AZ,
         spp = SPECIES, decay = DECAY, location_m = LOCATION_m, 
         diam_i_cm = DIAM_I_cm, diam_s_cm = DIAM_S_cm, diam_l_cm = DIAM_L_cm,
         length_m = LENGTH_m) %>%
  mutate(timestep = 'pre') %>%
  bind_rows(
    read_csv(here::here('02-data', '01-cleaned', 
                        'cwd_post.csv')) %>%
      select(-comment) %>%
      mutate(timestep = 'post')
  ) %>%
  mutate(decay = as.character(decay))


plot_metadata = 
  read_csv(here::here('02-data','01-cleaned',
                      'plot_metadata_post.csv')) %>%
  left_join(
    read_csv(here::here('02-data', '01-cleaned', 
                        'plot_metadata.csv')) %>%
      select(transect_id = TRANSECT,
             inv_date = Date,
             herb_pcent = HerbCover),
    by = c('transect_id' = 'transect_id'),
    suffix = c('_post', '_pre')
  )


shrubs = 
  read_csv(here::here('02-data', '01-cleaned', 
                      'shrubs.csv')) %>%
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT,
         az = AZ, spp = SPECIES, status = STATUS, location_m = LOCATION_m,
         distance_m = DISTANCE_m, maj_diam_m = MAJ_DIAM_m, 
         min_diam_m = MINOR_DIAM_m,
         height_m = HEIGHT_m) %>%
  mutate(timestep = 'pre') %>%
  bind_rows(
    read_csv(here::here('02-data', '01-cleaned', 
                        'shrubs_post.csv')) %>%
      mutate(timestep = 'post') %>%
      select(-comment)
  )

head(shrubs)

# load constants tables
species_codes = readRDS(here::here('02-data',
                                   '00-source',
                                   'van_wagtendonk',
                                   'species_codes.rds'))

kvals = readRDS(here::here('02-data',
                           '00-source',
                           'van_wagtendonk',
                           'kvals.rds'))


litterduff_coeffs = readRDS(here::here('02-data',
                                       '00-source',
                                       'van_wagtendonk',
                                       'litterduff_coeffs.rds'))


QMDcm = readRDS(here::here('02-data',
                           '00-source',
                           'van_wagtendonk',
                           'QMDcm.rds'))

SEC = readRDS(here::here('02-data',
                         '00-source',
                         'van_wagtendonk',
                         'SEC.rds'))

SG = readRDS(here::here('02-data',
                        '00-source',
                        'van_wagtendonk',
                        'SG.rds'))

sg_1000r = readRDS(here::here('02-data',
                              '00-source',
                              'van_wagtendonk',
                              'vw96_sg_1000r.rds'))

#### make the mcginnis table ###################################################

# full version of the table, because i'm sick of customizing it
# McGinnis, T. W., Shook, C. D., & Keeley, J. E. (2010). 
# Estimating aboveground biomass for broadleaf woody plants and young 
# conifers in Sierra Nevada, California, forests. Western Journal of Applied 
# Forestry, 25(mm), 203–209. Table 1: "Regression equatinos for aboveground 
# biomass of broadleaf shrubs and trees in northern, central, and southern 
# sierra nevada western slope collection sites". According to the text,
# "The following simple linear regression formula was used to 
# predict biomass (Y) in grams, for individual species and multiple-species 
# groups: Y = exp(a+b*ln(x)) where x is the independent variables. The multiple
# linear regression formula (where diameter and height were used to predict 
# biomass) was in the following form: Y = exp(a+b_1*ln(crown_diameter)+b_2*ln(heiht))
# Predicted biomass should be multiplied by the correction factor (CF) to 
# correct for log bias. The initial generic set of equations includes 
# all the species listed in the first column, plus: Ceanothus cuneatus, Chamaebatia 
# fotiolosa, Eriodictyon californicum, Heteromeles arbutifolia, Nama lobbii, 
# Purshia tridentata, and Toxicodendron diversilobum
mcginnis_t1 = 
  data.frame(
    
    spp = 
      rep(c('generic', 'ARPA', 'PARVI', 'CECO', 'CEA-', 
            'CEIN', 'CHSE', 'CHR-', 'PREM', 'QUE-', 'QUKE',
            'RIB-', 'SYM-'),
          each = 3),
    component = 
      rep(c('<0.64cm Live', 'Live foliage', 'Total biomass'), times = 13),
    
    crown_diam_cm.a = 
      c(-4.254, -4.168, -4.305, 
        -3.844, -6.008, -4.218, 
        -4.250, -3.080, -4.315, 
        -4.912, -7.693, -6.025, 
        -3.406, -4.539, -3.326, 
        -2.194, -4.541, -3.000,
        -3.657, -2.918, -4.965, 
        -3.502, -3.069, -3.871, 
        -3.086, -10.083, -3.370, 
        -3.896, -4.035, -3.669, 
        -3.537, -5.602, -3.393, 
        -4.318, -6.576, -4.859, 
        -6.141, -8.378, -6.840),
    crown_diam_cm.b = 
      c(2.121, 1.889, 2.312, 
        2.042, 2.611, 2.383, 
        2.125, 1.914, 2.484, 
        2.295, 2.579, 2.653, 
        2.016, 2.115, 2.120, 
        1.720, 1.849, 2.051, 
        1.963, 1.811, 2.533, 
        1.954, 1.637, 2.263, 
        1.780, 3.006, 1.970, 
        1.995, 1.960, 2.216, 
        1.829, 2.184, 2.109, 
        2.094, 2.293, 2.344, 
        2.548, 2.648, 2.784), 
    crown_diam_cm.cf = 
      c(1.25, 1.69, 1.30, 
        1.05, 1.07, 1.05, 
        1.26, 1.44, 1.24, 
        1.17, 1.25, 1.19, 
        1.15, 1.22, 1.13, 
        1.31, 1.50, 1.27, 
        1.03, 1.11, 1.05, 
        1.22, 1.32, 1.63, 
        1.07, 1.99, 1.10, 
        1.40, 1.33, 1.14, 
        1.04, 1.02, 1.02, 
        1.33, 1.40, 1.44, 
        1.49, 1.91, 1.57), 
    
    plant_height_cm.a = 
      c(-1.597, -2.076, -2.062, 
        -0.966, -2.208, -1.002, 
        -5.706, -4.093, -5.584, 
        -4.506, -6.401, -5.557, 
        0.780, -1.039, 0.803, 
        -0.840, -4.185, -2.509, 
        -3.940, 0.047, -6.019, 
        -7.069, -6.404, -8.690, 
        -5.538, -20.045, -6.683, 
        0.421, -0.317, -0.238, 
        -1.608, -3.662, -1.936, 
        -7.644, -9.646, -8.709, 
        -0.482, -1.586, -0.442), 
    plant_height_cm.b = 
      c(1.557, 1.461, 1.863, 
        1.368, 1.714, 1.638, 
        2.374, 2.069, 2.697, 
        2.611, 2.728, 3.019, 
        1.212, 1.540, 1.356, 
        1.471, 1.812, 1.989, 
        2.146, 1.147, 2.960, 
        2.681, 2.327, 3.265, 
        2.308, 5.162, 2.691, 
        1.059, 1.154, 1.472, 
        1.342, 1.677, 1.705, 
        3.015, 3.187, 3.409, 
        1.445, 1.235, 1.515), 
    plant_height_cm.cf = 
      c(1.37, 1.42, 1.37, 
        0.87, 1.12, 0.97, 
        0.89, 1.12, 1.12, 
        1.00, 1.29, 1.12, 
        0.85, 0.69, 0.79, 
        1.03, 1.03, 0.91, 
        0.79, 0.71, 0.87, 
        1.12, 1.01, 1.31, 
        0.92, 1.57, 0.98, 
        1.81, 1.69, 1.65, 
        0.79, 0.78, 0.58, 
        1.10, 1.08, 1.21, 
        0.99, 1.30, 1.08), 
    
    crown_vol_m3.a = 
      c(6.362, 5.297, 7.278, 
        6.147, 6.753, 7.439, 
        6.188, 6.331, 7.868, 
        6.644, 5.292, 7.334, 
        6.591, 6.051, 7.233, 
        6.365, 4.661, 7.200, 
        6.470, 6.394, 8.099, 
        6.350, 5.207, 7.563, 
        5.693, 4.810, 6.375, 
        6.008, 5.703, 7.340, 
        5.508, 5.209, 7.059, 
        6.291, 5.013, 7.023, 
        7.132, 4.952, 7.550), 
    crown_vol_m3.b = 
      c(0.739, 0.662, 0.812, 
        0.644, 0.819, 0.750, 
        0.741, 0.671, 0.858, 
        0.812, 0.919, 0.939, 
        0.606, 0.671, 0.654, 
        0.618, 0.670, 0.725, 
        0.746, 0.673, 0.961, 
        0.721, 0.612, 0.844, 
        0.594, 1.082, 0.672, 
        0.688, 0.684, 0.774, 
        0.575, 0.697, 0.685, 
        0.766, 0.822, 0.585, 
        0.947, 0.820, 0.995), 
    crown_vol_m3.cf = 
      c(1.25, 1.65, 1.27, 
        1.07, 1.12, 1.08, 
        1.16, 1.32, 1.20, 
        1.16, 1.23, 1.18, 
        1.21, 1.20, 1.16, 
        1.22, 1.37, 1.20, 
        1.04, 1.08, 1.09, 
        1.24, 1.29, 1.57, 
        1.11, 1.52, 1.11, 
        1.46, 1.33, 1.14, 
        1.12, 1.08, 1.03, 
        1.35, 1.45, 1.46, 
        1.28, 1.92, 1.36), 
    
    diam_height_cm.a = 
      c(-4.250, -4.363, -4.658, 
        -3.949, -6.095, -4.396, 
        -4.073, -1.832, -2.892, 
        -5.909, -8.228, -7.179, 
        -3.416, -3.717, -3.083, 
        -2.282, -5.323, -3.807, 
        -3.747, -5.774, -5.807, 
        -4.497, -4.459, -6.128, 
        -1.018, -11.353, -2.176, 
        -2.928, -3.360, -3.313, 
        -3.593, -5.623, -3.362, 
        -6.079, -8.476, -6.998, 
        -6.593, -8.677, -7.294), 
    diam_height_cm.b1 = 
      c(2.124, 1.754, 2.079, 
        1.937, 2.524, 2.204, 
        2.305, 3.192, 3.941, 
        1.903, 2.370, 2.199, 
        2.026, 1.293, 1.877, 
        1.659, 1.309, 1.494, 
        1.918, 1.713, 2.112, 
        1.705, 1.289, 1.698, 
        2.190, 2.760, 2.215, 
        2.693, 2.447, 2.473, 
        2.454, 2.425, 1.763, 
        1.640, 1.672, 1.794, 
        1.931, 2.241, 2.165), 
    diam_height_cm.b2 = 
      c(-0.004, 0.193, 0.336, 
        0.143, 0.118, 0.224, 
        -0.214, -1.514, -1.726, 
        0.715, 0.382, 0.828, 
        -0.009, 0.761, 0.225, 
        0.081, 0.715, 0.738, 
        0.072, 0.860, 0.677, 
        0.467, 0.653, 1.060, 
        -0.860, 0.525, -0.512, 
        -0.908, -0.633, -0.334, 
        -0.583, -0.225, 0.322, 
        0.923, 1.140, 1.121, 
        0.883, 0.583, 0.885), 
    diam_height_cm.cf = 
      c(1.25, 1.69, 1.28, 
        1.05, 1.07, 1.05, 
        1.33, 1.48, 1.20, 
        1.15, 1.26, 1.15, 
        1.18, 1.19, 1.16, 
        1.35, 1.51, 1.25, 
        1.03, 1.12, 1.04, 
        1.25, 1.35, 1.65, 
        1.06, 2.78, 1.12, 
        1.19, 1.25, 1.13, 
        1.02, 1.02, 1.02, 
        1.29, 1.34, 1.38, 
        1.43, 1.99, 1.50)
    
  )


#### get overstory composition #################################################

# need plot-level species composition expressed as propotion of total BA
head(trees)

# first aggregate to plot level by summing individual trees
# note that we lump live and dead trees together for the purpose of estimating 
# fuels, because both contribute (or have contributed) to the existing fuel bed
overstory_composition = 
  trees %>%
  
  # keep only big trees (>11.4cm DBH)
  filter(dbh_cm >= 11.4) %>%
  
  # map species wihout coefficients from van wagtendonk's work to 'OTHER' 
  # category, will get the 'All species' coefficient
  mutate(spp = ifelse(!is.element(spp, as.character(species_codes$species_code)),
                      'OTHER',
                      spp)) %>%
  
  mutate(ba_m2 = pi*((dbh_cm/100)/2), # basal area of each tree in m2
         # scaled by 0.05ha plot size
         ba_m2ha = ba_m2*(1/0.05)) %>%
  
  group_by(comp, plot_id, spp) %>%
  summarise(ba_m2ha = sum(ba_m2ha)) %>%
  ungroup() %>%
  
  # fill in 0s for species which werent present on a plot
  complete(nesting(comp, plot_id), spp) %>%
  mutate(ba_m2ha = ifelse(is.na(ba_m2ha), 0, ba_m2ha)) %>%
  
  # get the total BA/ha on each plot
  left_join(x = .,
            y = 
              group_by(., comp, plot_id) %>% 
              summarise(total_ba_m2ha = sum(ba_m2ha)) %>%
              ungroup()) %>%
  
  # get the proportion of plot total basal area occoupied by each species
  mutate(pba = ba_m2ha / total_ba_m2ha) %>%
  
  select(comp, plot_id, spp, pba)


overstory_composition

#### calculate composition-weighted coefficients ###############################

# get BA-weighted average coefficient for fuel load (kg/m2) as a function of 
# litterduff_depth (cm)
litterduff_coeffs = 
  overstory_composition %>%
  left_join(litterduff_coeffs,
            by = c('spp' = 'spp')) %>%
  mutate(weighted_litter = pba * litter_coeff,
         weighted_duff = pba * duff_coeff) %>%
  group_by(comp, plot_id) %>%
  summarise(litter_coeff = sum(weighted_litter),
            duff_coeff = sum(weighted_duff)) %>%
  ungroup()


# get BA-weighted SEC (secant of acute angle) for each timelag class
SEC = 
  overstory_composition %>%
  left_join(SEC,
            by = c('spp' = 'spp')) %>%
  select(comp, plot_id, spp, pba, x1h, x10h, x100h, 
         x1000s = x1000h, x1000r = x1000h) %>%
  pivot_longer(c(x1h, x10h, x100h, x1000s, x1000r),
               names_to = 'timelag_class',
               values_to = 'sec',
               names_prefix = 'x') %>%
  mutate(weighted = sec*pba) %>%
  group_by(comp, plot_id, timelag_class) %>%
  summarise(weighted_sec = sum(weighted)) %>%
  ungroup()

SEC

# get the BA-weighted average specific gravity for each timelag class
SG = 
  overstory_composition %>%
  left_join(SG,
            by = c('spp' = 'spp')) %>%
  select(comp, plot_id, spp, pba, x1h, x10h, x100h, x1000s) %>%
  pivot_longer(c(x1h, x10h, x100h, x1000s),
               names_to = 'timelag_class',
               values_to = 'sg',
               names_prefix = 'x') %>%
  mutate(weighted = sg*pba) %>%
  group_by(comp, plot_id, timelag_class) %>%
  summarise(weighted_sg = sum(weighted)) %>%
  ungroup() %>%
  
  # SG of 1000r doesn't vary by species, so make a table with every comp and plot_id 
  # where timelag class is 1000r and the weighted SG is the value from van 
  # wagtendonk's paper
  bind_rows(plot_metadata %>%
              select(comp, plot_id) %>%
              group_by(comp, plot_id) %>%
              summarise() %>%
              ungroup() %>%
              mutate(timelag_class = '1000r',
                     weighted_sg = sg_1000r))


# get BA-weighted average QMD (cm2) for each FWD timelag class
QMDcm = 
  overstory_composition %>%
  left_join(QMDcm,
            by = c('spp' = 'spp')) %>%
  select(comp, plot_id, spp, pba, x1h, x10h, x100h) %>%
  pivot_longer(c(x1h, x10h, x100h),
               names_to = 'timelag_class',
               values_to = 'qmd_cm',
               names_prefix = 'x') %>%
  mutate(weighted = qmd_cm*pba) %>%
  group_by(comp, plot_id, timelag_class) %>%
  summarise(weighted_qmd = sum(weighted)) %>%
  ungroup()


#### estimate litter + duff loads ##############################################

#' Litter and duff are measured as depths at specific points along a sampling
#' transect. Van WAgtendonk et al. (1998) developed regressions for litter,
#' duff, and combined-litter-and-duff loading (kg / m^2) as a function of
#' depth (cm) for 19 different Sierra Nevada conifer species. See vignette
#' for details

head(litter_duff)

# get the estimated litter and duff loads for each observation, keep the duff 
# depth, which FOFEM wants
litter_duff = 
  litter_duff %>%
  left_join(litterduff_coeffs, by = c('comp' = 'comp',
                                      'plot_id' = 'plot_id')) %>%
  mutate(litter_mgha = 10 * litter_cm * litter_coeff,
         duff_mgha = 10*duff_cm*duff_coeff) %>%
  
  # standardize columns
  select(timestep, comp, plot_id, transect_id,  az = az, location_m, litter_cm, duff_cm,
         litter_mgha, duff_mgha)

head(litter_duff)

summary(litter_duff)

# duff depths are lower and litter loads higher than I would have expected, 
# based on the FFS controls (which had ~20-25 mg/ha of litter and 40-60mg/ha of 
# duff, compared to ~34 mg/ha of litter here and ~37 mg/ha of duff here). 
# thinking about the stand structure, I guess I'm not surprised to have so many 
# 0 duff observations. There are many skid trails in these stands with some 
# litter but no duff. I also think the distinction between duff and litter is 
# not obvious, especially in situations where the fuelbed has been disturbed by 
# heavy machinery. 

#### FWD loading ###############################################################

fwd = 
  fwd_tallies %>%
  
  # reformat tallies to longwise
  pivot_longer(cols = c(a1h, a10h, a100h, b1h, b10h, b100h),
               names_to = 'subsample_id',
               values_to = 'count') %>%
  mutate(timelag_class = gsub(x = subsample_id,
                              pattern = 'a|b',
                              replacement = ''),
         subsample = gsub(x = subsample_id,
                          pattern = '1h|10h|100h',
                          replacement = '')) %>%
  
  # add columns for slope, QMD, SEC, SG, and transect length, and the conversion 
  # constant k
  left_join(., 
            y = QMDcm) %>%
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>%
  left_join(.,
            y = 
              plot_metadata %>% 
              select(transect_id, slp_pcent) %>%
              mutate(slp_c = sqrt(1+(slp_pcent/100)^2)) %>%
              select(transect_id, slp_c)) %>%
  mutate(transect_length_m = 1,
         k = 1.234) %>%
  
  # use brown's equations to estimate fuel load
  mutate(mgha = 
           (k*weighted_qmd*weighted_sec*slp_c*weighted_sg*count)/
           transect_length_m)  %>%
  
  # standardize columns
  select(timestep, comp, plot_id, transect_id, 
         az, location_m,
         subsample,
         timelag_class,
         count,
         mgha)


head(fwd)

#### CWD loading ###############################################################

head(cwd)

cwd = 
  
  cwd %>%
  
  # bin into sound and rotten
  mutate(decay_class = ifelse(decay<=2,'S','R'),
         diam_cm = diam_i_cm) %>%
  
  # get the sum of squared diameters for each subtransect and timelag class
  group_by(timestep, comp, plot_id, transect_id, az, decay_class) %>%
  summarize(ssd_cm2 = sum(diam_cm^2)) %>%
  ungroup() %>%
  
  # left join with the plot metadata to get the length of each transect
  right_join(plot_metadata %>%
               select(comp, plot_id, transect_id, az, 
                      transect_length_m = max_dist_m) %>%
               expand_grid(timestep = c('pre', 'post'))) %>%
  
  # left join with the complete combination of transect and timelag class to 
  # get the structural missings, fill them with 0s
  complete(nesting(timestep, comp, plot_id, transect_id, az, transect_length_m), 
           decay_class, 
           fill = list(ssd_cm2 = 0)) %>%
  
  # ditch the NAs we just created (there were a few plots with no CWD, so 
  # NA timelag class after the right join above)
  filter(!is.na(decay_class)) %>%
  
  # add a 'timelag class' (1000s or 1000r) for matching up with the van 
  #wagtendonk tables
  mutate(timelag_class = ifelse(decay_class=='S',
                                '1000s',
                                '1000r')) %>%
  
  # and join in the coefficients 
  left_join(.,
            y = SEC) %>%
  left_join(.,
            y = SG) %>% 
  left_join(.,
            y = 
              plot_metadata %>% 
              select(transect_id, slp_pcent) %>%
              mutate(slp_c = sqrt(1+(slp_pcent/100)^2)) %>%
              select(transect_id, slp_c)) %>%
  mutate(k = 1.234) %>%
  
  # and use brown's equation (with the new k const) to estimate fuel load
  mutate(mgha = 
           (k*ssd_cm2*weighted_sec*slp_c*weighted_sg)/
           transect_length_m) %>%
  
  # drop some intermediate columns
  select(timestep, comp, plot_id, transect_id, az, timelag_class, mgha) %>%
  
  # pivot wider
  mutate(timelag_class = paste0('x', timelag_class, '_mgha')) %>%
  pivot_wider(id_cols = 
                c('timestep', 'comp', 'plot_id', 'transect_id', 'az'),
              names_from = timelag_class, values_from = mgha)




#### trees #####################################################################
head(trees)


trees_summary = 
  trees %>%
  filter(dbh_cm >= 11.3 & status_pre == 'L') %>%
  mutate(tph = 1/0.05,
         ba_m2ha = pi*((dbh_cm/200)**2) * tph) %>%
  group_by(comp, plot_id) %>%
  summarise(ba_m2ha = sum(ba_m2ha),
            tph = sum(tph)) %>%
  ungroup() %>%
  right_join(plot_metadata %>%
               group_by(comp, plot_id) %>%
               summarise()) %>%
  mutate(ba_m2ha = ifelse(is.na(ba_m2ha),0,ba_m2ha),
         tph = ifelse(is.na(tph), 0, tph)) %>%
  summarise(ba_m2ha.mean = mean(ba_m2ha),
            ba_m2ha.sd = sd(ba_m2ha),
            tph.mean = mean(tph),
            tph.sd = sd(tph)) %>%
  mutate(timestep = 'pre') %>%
  bind_rows(
    trees %>%
      filter(dbh_cm >= 11.3 & status_post == 'L') %>%
      mutate(tph = 1/0.05,
             ba_m2ha = pi*((dbh_cm/200)**2) * tph) %>%
      group_by(comp, plot_id) %>%
      summarise(ba_m2ha = sum(ba_m2ha),
                tph = sum(tph)) %>%
      ungroup() %>%
      right_join(plot_metadata %>%
                   group_by(comp, plot_id) %>%
                   summarise()) %>%
      mutate(ba_m2ha = ifelse(is.na(ba_m2ha),0,ba_m2ha),
             tph = ifelse(is.na(tph), 0, tph)) %>%
      summarise(ba_m2ha.mean = mean(ba_m2ha),
                ba_m2ha.sd = sd(ba_m2ha),
                tph.mean = mean(tph),
                tph.sd = sd(tph)) %>%
      mutate(timestep = 'post')
  )

#### shrubs ####################################################################


shrubs = 
  shrubs %>%
  
  # height values were rounded to the nearest 0.25m; this results in some 
  # individuals being assigned height = 0, which causes problems with the 
  # biomass estimations below; make these a small number instead
  mutate(height_m = ifelse(height_m==0, 0.1, height_m)) %>%
  
  # get the transect length for the relevant transect
  left_join(plot_metadata %>%
              select(transect_id, transect_length_m = max_dist_m)) %>%
  
  # estimate the crown perimeter of each shrub; 
  # https://www.mathsisfun.com/geometry/ellipse-perimeter.html
  mutate(crown_perim_m = 
           pi*
           (maj_diam_m+min_diam_m)*
           (1+
              0.25*(((maj_diam_m-min_diam_m)**2)/
                      ((maj_diam_m+min_diam_m)**2))+
              (1/64)*((((maj_diam_m-min_diam_m)**2)/
                         ((maj_diam_m+min_diam_m)**2))**2)+
              (1/256)*((((maj_diam_m-min_diam_m)**2)/
                          ((maj_diam_m+min_diam_m)**2))**3)+
              (25/16384)*((((maj_diam_m-min_diam_m)**2)/
                             ((maj_diam_m+min_diam_m)**2))**4))) %>%
  # using the perimeter, calculate the effective diameter (the average diameter 
  # of the ellipse across all possible intersection angles); also 
  # calculate the mean crown diameter (mean of major and minor), which is 
  # used in the mcginnis biomass calcs
  mutate(eff_diam_m = crown_perim_m / pi,
         mean_diam_m = (maj_diam_m+min_diam_m)/2) %>%
  
  # calculate the quantity of interest, which is the individual's biomass
  ## first, left-join in the relevant coefficients; this triples the number of 
  ## rows because we have 3 biomass categories (total, live foliage, and live fine)
  mutate(spp.match = 
           ifelse(spp=='RIRO',
                  'RIB-',
                  ifelse(spp=='SYMO',
                         'SYM-',
                         ifelse(is.element(spp, mcginnis_t1$spp),
                                spp,
                                'generic')))) %>%
  left_join(mcginnis_t1 %>%
              select(spp, component, 
                     diam_height_cm.a, diam_height_cm.b1, diam_height_cm.b2,
                     diam_height_cm.cf),
            by = c('spp.match' = 'spp')) %>%
  
  ## use the equation given in mcginnis to calculate biomass in grams
  mutate(biomass_g = 
           exp(diam_height_cm.a+
                 (diam_height_cm.b1*log(mean_diam_m*100))+
                 (diam_height_cm.b2*log(height_m*100))) *
           diam_height_cm.cf) %>%
  
  ## where relevant, discount the biomass by the % torch
  mutate(
    biomass_g = 
      ifelse((component=='Live foliage'|component=='<0.64cm Live')&timestep=='post',
             biomass_g*(1-(torch/100)),
             biomass_g)) %>%
  
  # now use equation from battles 1996 to estimate a quantity on a per-area 
  # basis
  ## biomass per meter of intersection (calculated for each individual)
  mutate(biomass_gm = biomass_g / eff_diam_m) %>%
  
  ## aggregate by component:transect:timestep and get the 'mean' transect 
  ## length (theres only 1 value) and the sum of the biomass per meter of intersection
  group_by(timestep, comp, plot_id, transect_id, component) %>%
  summarise(transect_length_m = mean(transect_length_m), 
            biomass_gm = sum(biomass_gm)) %>%
  ungroup() %>%
  mutate(biomass_gm2 = biomass_gm / transect_length_m,
         biomass_mgha = biomass_gm2 * (0.001)*(0.001)*10000,
         biomass_tac = biomass_mgha * 0.44609) %>%
  select(timestep, comp, plot_id, transect_id, component, biomass_mgha) %>%
  pivot_wider(names_from = 'component', values_from = 'biomass_mgha') %>%
  rename(live_fine_biomass_mgha = `<0.64cm Live`,
         live_foliage_bioamss_mgha = `Live foliage`,
         total_biomass_mgha = `Total biomass`) %>%
  
  # for fofem, just want the total
  select(timestep, comp, plot_id, transect_id, total_biomass_mgha)



#### combine estimates #########################################################

head(trees_summary)

head(litter_duff)

summary_stats = 
  trees_summary %>%
  pivot_longer(cols = c('ba_m2ha.mean', 'ba_m2ha.sd', 'tph.mean', 'tph.sd'),
               names_to = c('metric', 'stat'),
               names_sep = '\\.') %>%
  bind_rows(
    litter_duff %>%
      pivot_longer(cols = c('litter_mgha', 'duff_mgha', 'litter_cm' ,'duff_cm'),
                   names_to = 'metric', values_to = 'value') %>%
      group_by(timestep, comp, plot_id, metric) %>%
      summarise(value = mean(value, na.rm = TRUE))%>%
      ungroup() %>%
      group_by(timestep, metric) %>%
      summarise(mean = mean(value),
                sd = sd(value)) %>%
      ungroup() %>%
      pivot_longer(cols = c('mean', 'sd'),
                   names_to = 'stat',
                   values_to = 'value')
  ) %>%
  bind_rows(
    fwd %>%
      pivot_longer(cols = c('count', 'mgha'),
                   names_to = 'metric', values_to = 'value') %>%
      mutate(metric = paste0(timelag_class,'_',metric)) %>%
      group_by(timestep, comp, plot_id, metric) %>%
      summarise(value = mean(value, na.rm = TRUE)) %>%
      ungroup() %>%
      group_by(timestep, metric) %>%
      summarise(mean = mean(value),
                sd = sd(value)) %>%
      ungroup() %>%
      pivot_longer(cols = c('mean', 'sd'),
                   names_to = 'stat',
                   values_to = 'value')
  ) %>%
  bind_rows(
    cwd %>%
      mutate(x1000h_mgha = x1000r_mgha+x1000s_mgha) %>%
      group_by(timestep, comp, plot_id) %>%
      summarise(value = mean(x1000h_mgha, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(metric = '1000h_mgha') %>%
      group_by(timestep, metric)%>%
      summarise(mean = mean(value),
                sd = sd(value)) %>%
      ungroup() %>%
      pivot_longer(cols = c('mean', 'sd'),
                   names_to = 'stat',
                   values_to = 'value')
  ) %>%
  bind_rows(
    shrubs %>%
      mutate(shrub_mgha = total_biomass_mgha) %>%
      group_by(timestep, comp, plot_id)%>%
      summarise(value = mean(shrub_mgha, na.rm = TRUE))%>%
      ungroup() %>%
      mutate(metric = 'shrub_mgha') %>%
      group_by(timestep, metric) %>%
      summarise(mean = mean(value),
                sd = sd(value))%>%
      ungroup() %>%
      pivot_longer(cols = c('mean', 'sd'),
                   names_to = 'stat',
                   values_to = 'value')
  ) %>%
  pivot_wider(id_cols = c('timestep', 'metric'),
              names_from = 'stat',
              values_from = 'value') %>%
  mutate(mean_sd = paste0(round(mean, 1),' (', round(sd, 1),')')) %>%
  select(timestep, metric, mean_sd) %>%
  pivot_wider(id_cols = 'metric',
              names_from = 'timestep',
              values_from = 'mean_sd')%>%
  mutate(metric = factor(metric, levels = c('tph', 'ba_m2ha','duff_cm', 'duff_mgha',
                                            'litter_cm', 'litter_mgha', '1h_count', 
                                            '1h_mgha', '10h_count', '10h_mgha',
                                            '100h_count', '100h_mgha', '1000h_mgha',
                                            'shrub_mgha'))) %>%
  arrange(metric) %>%
  select(metric, pre, post)
  

summary_stats

write.csv(summary_stats,
          here::here('04-communication',
                     'tables',
                     'summary_stats.csv'),
          row.names = FALSE)

#### interplot distances #######################################################
library(sf)
plot_locations = 
  st_read(here::here('02-data', '00-source', 'bfrs', 'CARB_BFRSPlots.shp')) %>%
  mutate(x_coord = sf::st_coordinates(.)[,'X'],
        y_coord = sf::st_coordinates(.)[,'Y']) %>%
  select(COMP, PLOT, x_coord, y_coord) %>%
  as.data.frame() %>%
  select(-geometry) %>%
  as.tibble() 
head(plot_locations)

dist(as.matrix(plot_locations[,c('x_coord', 'y_coord')])) %>% summary()


#### modeling stand structure and fire effects #################################

head(trees)

library(glmmTMB)

hist(trees$torch_ht)
hist(trees$scorch_ht)
hist(trees$char_ht)
hist(trees$char_p)

summary(trees$char_ht)

trees_imputed = 
  trees %>%
  mutate(char_ht_imputed = ifelse(is.na(char_ht), 
                                  min(char_ht[char_ht>0], na.rm = TRUE), 
                                  char_ht+min(char_ht[char_ht>0], na.rm = TRUE)),
         scorch_ht_imputed = ifelse(is.na(scorch_ht),
                                    min(scorch_ht[scorch_ht > 0], na.rm=TRUE),
                                    scorch_ht+min(scorch_ht[scorch_ht>0], na.rm = TRUE)))

summary(trees_imputed$char_ht_imputed)

char_fit = 
  glmmTMB(data = trees_imputed,
          family = gaussian(),
          log(char_ht_imputed) ~ comp + (1|plot_id))

scorch_fit = 
  glmmTMB(data = trees %>%
            filter(!is.na(scorch_ht)) %>%
            filter(dbh_cm >= 11.4),
          family = gaussian(),
          scorch_ht ~ comp + (1|plot_id))

trees %>%
  filter(dbh_cm >= 11.4) %>%
  mutate(scorch_meas = !is.na(scorch_ht),
         char_meas = !is.na(char_ht),
         torch_meas = !is.na(torch_ht)) %>%
  group_by(comp, scorch_meas, char_meas, torch_meas) %>%
  summarise(n = n(),
            med_dbh = median(dbh_cm),
            med_ht = median(height_m))

summary(scorch_fit)
summary(char_fit)

char_fit_2 = 
  glmmTMB(data = trees %>%
            mutate(char_at_all = ifelse(is.na(char_ht)|char_ht == 0, FALSE, TRUE)),
          family = binomial(),
          char_at_all ~ comp + (1|plot_id))

summary(char_fit_2)

summary_stats


ggplot(data = trees_imputed,
       aes(x = char_ht_imputed, color = comp))+
  geom_density()+
  scale_x_log10()+
  theme_minimal()

scorch_fit = 
  glmmTMB(data = trees,
          family = tweedie(),
          scorch_ht ~ comp)

summary(char_fit)

ggplot(data = trees,
       aes(x = scorch_ht, color = comp))+
  geom_density()+
  theme_minimal()


trees %>%
  filter(is.na(char_ht))

trees %>%
  filter(char_ht == 0)

head(litter_duff)

litter_fit = 
  glmmTMB(data = litter_duff %>%
            filter(timestep == 'pre'),
          litter_mgha ~ comp + (1|plot_id),
          family = tweedie())

summary(litter_fit)
