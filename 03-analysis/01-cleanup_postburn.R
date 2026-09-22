
#### setup #####################################################################

library(tidyverse)
library(here)
library(readxl)
library(ggplot2)
library(stringr)

#### load data #################################################################
#need to update filepath when there is new data imported 
filepath = 
  here("02-data", "00-source", "postburn", "bfrs_fuels_090721.xlsx")

plot_metadata <- 
  read_excel(filepath, 
             sheet = "plot_metadata", 
             col_types = c("text", "text", "numeric", "numeric", 
                           "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                           "date", "text", "text")) %>%
  mutate(Date = as.Date(as.character(Date)))


# warnings are about getting NA; fine
cwd <- 
  read_excel(filepath, 
             sheet="cwd_new", 
             col_types = c("text", "text", "text", 
                           "numeric", "numeric", 
                           "numeric", "numeric", "numeric", "numeric", 
                           "text", "text", "text"))


fwd_tallies <- 
  read_excel(filepath, 
             sheet="fwd_tallies", 
             col_types = c("text", "text", "numeric", "numeric", 
                           "numeric", "numeric", "numeric", 
                           "numeric", "numeric", "numeric", 
                           "text")) %>%
  rename(comment = `...11`)
# NAs investigated during cleaning


fwd_diams <- 
  read_excel(filepath, 
             sheet="fwd_diams", 
             col_types = c("text", "text", "numeric", "numeric", "text")) %>%
  rename(comment = `...5`)

litter_duff <- 
  read_excel(filepath, 
             sheet="litterduff", 
             col_types = c("text", "text", "numeric", "numeric", 
                           "numeric", 
                           "numeric", "text")) %>%
  rename(comment = `...7`)

shrubs <- 
  read_excel(filepath, 
             sheet="shrubs_new", 
             col_types = c("text", "text", "numeric", "numeric", "numeric", 
                           "text", "numeric", "numeric", "numeric", "text",
                           "numeric", "numeric", "text")) %>%
  rename(comment = `...13`)

trees <- 
  read_excel(filepath, 
             sheet="trees", 
             col_types = c("text", "text", "text", "numeric", "numeric",
                           'text', 'text', 'text', 'numeric', 'numeric', 'numeric', 
                           'text', 
                           'text', 'text', 'numeric', 'numeric', 'numeric', 
                           'numeric', 'numeric', 'numeric', 'text', 'numeric', 'text'))

consumption = 
  read_excel(filepath,
             sheet = 'consumption',
             col_types = c('text', 'text', 'numeric', 'numeric', 'numeric', 'text')) %>%
  rename(comment = `...6`)

valid_transects <- 
  expand.grid(c("A-001", "A-003", "A-005","A-021", "A-023", "A-025", "A-020",'A-021', "A-024", 
                "A-028", "B-008", "B-013", "B-010", "B-016", "B-026", "C-010", 
                "C-011", "C-025", "C-026", "C-022"), 
              c("000", "090", "180", "270")) %>%
  select(plot = Var1, az = Var2)%>%
  mutate(transect_id = paste0(plot, "-", az)) %>%
  pull(transect_id)

#### cleaning plot metadata ####################################################

head(plot_metadata)

plot_metadata = 
  plot_metadata %>%
  mutate(PLOT = paste0(COMP,'-',str_pad(PLOT, width = 3, side = 'left', pad = 0))) %>%
  mutate(TRANSECT = paste0(PLOT, '-', str_pad(AZ, width = 3, side = 'left', pad = 0))) %>%
  # keep only relevant columns
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT, 
         inv_date = Date, az = AZ, slp_pcent = SLP_pcent, herb_pcent = HerbCover, 
         comment_pre, comment_post) %>%
  
  # add column for out-of-bounds distance (transect extended beyond compartment) 
  mutate(max_dist_m = 30)

# manually update a few out of bounds distances based on reading the datasheet
plot_metadata[plot_metadata$transect_id=='A-023-000','max_dist_m'] = 18.1
plot_metadata[plot_metadata$transect_id=='A-024-090','max_dist_m'] = 17.0
plot_metadata[plot_metadata$transect_id=='A-016-270','max_dist_m'] = 20.6
plot_metadata[plot_metadata$transect_id=='C-022-270','max_dist_m'] = 17.0

#### cleaning litter and duff ##################################################

# criteria:
## depth should always be >= 0
## only valid locations
## NA values ok in depth cols
## no NA in COMP, PLOT, AZ, or Location
## COMP in (a, b, c)
## set number of rows
## specific set of COMP, PLOT combinations 
## 4 az for each plot
## for each comp/plot/az, a specific set of locatinos
## check high outliers for litter and duff depth

head(litter_duff)

summary(litter_duff)

valid_locations = c(2, 2.1, 2.2, 2.3, 2.4, 2.5, 3, 3.5, 4, 4.5, 5.5, 7, 7.5, 9.5, 19.5, 29.5)

litter_duff %>% filter(!is.element(Displacement_m, valid_locations)) # good

litter_duff %>% filter(is.na(litter_cm)|is.na(duff_cm)) %>% print(n = Inf) #fine

litter_duff %>% filter(is.na(COMP)|is.na(PLOT)|is.na(Displacement_m)) # good

litter_duff %>% filter(!is.element(COMP, c('A', 'B', 'C'))) # good

# there are 3 comp with 19 different plots, 4 transects per plot, and 
# 16 locatinos per transect
nrow(litter_duff) == 19*4*16 # good

# create ID versions of comp/plot/transect
litter_duff = 
  litter_duff %>%
  mutate(TRANSECT = paste0(PLOT, '-', str_pad(AZ, width = 3, side = 'left', pad = 0))) %>%
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT, 
         az = AZ, location_m = Displacement_m, litter_cm, duff_cm, comment)

# check for problematic comments
litter_duff %>% filter(!is.na(comment)) %>% pull(comment) # looks fine

#### fwd tallies ###############################################################

head(fwd_tallies)

fwd_tallies = 
  fwd_tallies %>%
  mutate(TRANSECT = paste0(PLOT, '-', str_pad(AZ, width = 3, side = 'left', pad = 0))) %>%
  select(comp = COMP, plot_id = PLOT, transect_id = TRANSECT,
         az = AZ, location_m = Displacement_m, 
         a1h, a10h, a100h, b1h, b10h, b100h, comment)

summary(fwd_tallies)

# tally columns are >= 0 good
# only certain valid locations

valid_locations_fwdtallies = 
  c(2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5, 29.5)

fwd_tallies %>% filter(!is.element(location_m, valid_locations_fwdtallies)) # good

fwd_tallies %>%
  filter(is.na(a1h)|is.na(a10h)|is.na(a100h)|is.na(b1h)|is.na(b10h)|is.na(b100h)) # fine

fwd_tallies %>%
  filter(is.na(comp)|is.na(plot_id)|is.na(location_m)|is.na(az)) # good

fwd_tallies %>% filter(!is.element(comp, c('A', 'B', 'C'))) # good

nrow(fwd_tallies) == 19*4*8 # good

valid_plots = c("A-001", "A-003", "A-005", "A-023", "A-025", "A-020", "A-024", 'A-021',
                "A-028", "B-008", "B-013", "B-010", "B-016", "B-026", "C-010", 
                "C-011", "C-025", "C-026", "C-022")

valid_locations_fwdtallies_full = 
  expand.grid(valid_plots, 
              c("000", "090", "180", "270"), 
              c( 2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 19.5, 29.5)) %>%
  select(PLOT = Var1, TRANSECT = Var2, Displacement_m = Var3) %>%
  mutate(LOCATION = paste0(PLOT, "-", TRANSECT, "-", Displacement_m)) %>%
  pull(LOCATION)

fwd_tallies %>% mutate(location_full = paste0(transect_id,'-', location_m)) %>%
  filter(!is.element(location_full, valid_locations_fwdtallies_full)) # good


#### fwd_diams #################################################################

head(fwd_diams)

fwd_diams = 
  fwd_diams %>%
  mutate(plot_id = paste0(COMP, '-', str_pad(PLOT, width = 3, side = 'left', pad = 0)),
         transect_id = paste0(plot_id, '-', str_pad(AZ, width = 3, side = 'left', pad = 0))) %>%
  select(comp = COMP, plot_id, transect_id, az = AZ, diam_mm = DIAM_mm, comment)

summary(fwd_diams) # all diams are positive, max diam is 75mm, some NAs

fwd_diams %>% filter(is.na(diam_mm)|!is.na(comment)) %>% print(n = Inf) # all cases with 0 fwd;
# drop them from this table
fwd_diams = fwd_diams %>% filter(!is.na(diam_mm))

fwd_diams %>% filter(!is.element(comp, c('A', 'B', 'C')))

#### shrubs ####################################################################

head(shrubs)
shrubs = 
  shrubs %>%
  mutate(plot_id = paste0(COMP, '-', str_pad(PLOT, width = 3, side = 'left', pad = 0)),
         transect_id = paste0(plot_id, '-', str_pad(AZ, width = 3, side = 'left', pad = 0))) %>%
  select(comp = COMP, plot_id, transect_id, az = AZ, location_m = LOCATION_m,
         distance_m = DISTANCE_m, spp = SPECIES, maj_diam_m = MAJ_DIAM_m, 
         min_diam_m = MINOR_DIAM_m, height_m = HEIGHT_m, status = STATUS, 
         scorch = SCORCH, torch = TORCH, comment)

# distance, diameters, should be >= 0; check plausibility of maxima and minima
summary(shrubs) # looks good

# check for high outliers in distance
shrubs %>% filter(distance_m >= 2) # both are big tanoak

# check for NAs in any data column, or non-NA in comment
shrubs %>% filter(is.na(comp)|is.na(plot_id)|is.na(transect_id)|is.na(az)|
                    is.na(location_m)|is.na(distance_m)|is.na(spp)|
                    is.na(maj_diam_m)|is.na(min_diam_m)|is.na(height_m)|
                    is.na(status)|is.na(scorch)|is.na(torch)|!is.na(comment))
# all are transects with no shrubs; drop them
shrubs = 
  shrubs %>%
  filter(is.na(comment))

shrubs %>% filter(!is.element(comp, c('A', 'B', 'C'))) # good

shrubs %>% filter(!is.element(status, c('L', 'D'))) # good

unique(shrubs$spp)[order(unique(shrubs$spp))]

shrubs %>% filter(maj_diam_m < min_diam_m)
# swap these
shrubs = 
  shrubs %>%
  mutate(actual_maj_diam = 
           ifelse(maj_diam_m >= min_diam_m,
                  maj_diam_m,
                  min_diam_m),
         actual_min_diam = 
           ifelse(maj_diam_m >= min_diam_m,
                  min_diam_m,
                  maj_diam_m)) %>%
  select(comp, plot_id, transect_id, az, location_m, distance_m,
         spp, maj_diam_m = actual_maj_diam, min_diam_m = actual_min_diam,
         height_m, status, scorch, torch, comment)

#### cwd #######################################################################

head(cwd)

cwd = 
  cwd %>%
  mutate(plot_id = paste0(COMP, '-', str_pad(PLOT, width = 3, side = 'left', pad = 0)),
         transect_id = paste0(plot_id, '-', str_pad(AZ, width = 3, side = 'left', pad = 0))) %>%
  select(comp = COMP, plot_id, transect_id, az = AZ,
         location_m = LOCATION_m, spp = SPECIES, diam_i_cm = DIAM_I_cm, diam_s_cm = DIAM_S_cm,
         diam_l_cm = DIAM_L_cm, length_m = LENGTH_m, decay = DECAY, comment = NOTES)

summary(cwd) # good range for az, location; some minimums too small for diam _s 
# (min should be 7.6); probably want to drop NAs

cwd %>% filter(is.na(location_m))

cwd = cwd %>% filter(!is.na(location_m))

# only valid combinations of comp:plot:az
cwd %>% filter(!is.element(transect_id, valid_transects))

#check that big piece; i think this is legit i remember a whole down tree
cwd %>% filter(length_m > 30) # yeah its legit
 
# not goint to worry about proper ordering or large/intersection/small diameters;
# according to the crew (pre-treatment) it was pretty common for pieces to be 
# bigger in the middle than the potentially-damaged ends

# also going to not worry about the too-small diam_s_cms; probably crew error 
# and they should have truncated the log earlier according to the protocol, but 
# if they didn't i don't know the true length

# decay = 1:5
unique(cwd$decay)

unique(cwd$spp)[order(unique(cwd$spp))]

cwd %>% filter(spp=='NA')
cwd = 
  cwd %>% 
  mutate(spp = 
           ifelse(spp=='NA',
                  'UNK',
                  ifelse(spp=='UNK HARDWOOD',
                         'HARDWOOD',
                         ifelse(spp=='PINUS',
                                'PINE',
                                ifelse(spp=='UF',
                                       'FIR',
                                       spp)))))

cwd %>% filter(!is.na(comment)) %>% print(n = Inf)
cwd %>% filter(!is.na(comment)) %>% pull(comment)

#### trees #####################################################################

print(trees, width = Inf)

trees = 
  trees %>%
  mutate(tree_id = paste0(PLOT, '-', TS, '-', LOCATION, '-', DISTANCE)) %>%
  select(comp = COMP, plot_id = PLOT, tree_id, ts = TS, location_m = LOCATION,
         distance_m = DISTANCE, tag = TAG, status_pre = STATUS, spp = SPECIES,
         dbh_cm, height_m, htcb_m, decay_pre = decay, 
         status_post, torch_p, scorch_p, char_p, decay_post, height_m_post = height_post,
         ht_meas, torch_ht, scorch_ht, char_ht, comment = CONCERN)

summary(trees)
# these are semi-cleaned pretreatment data? not sure; redo all the cleaning steps
# with some 
# exceptions which should have comments (?)
# bounds look OK for torch_p, scorch_p, char_p;

# only valid combinations of comp:plot
trees %>% filter(!is.element(plot_id, valid_plots))

# valid location for big trees
trees %>%
  filter(
    !(
      (location_m <= 15 & distance_m <=5) | (location_m <= 5 & distance_m <= 15)
    )
  )

# valid location for saplings
trees %>%
  filter(
    !(
      (location_m <= 15 & distance_m <=1) | (location_m <= 1 & distance_m <= 15)
    ) & dbh_cm < 11.4
  )

# valid status
unique(trees$status_pre)
unique(trees$status_post)




# check validity of status_post
unique(trees$status_post)

# if status is live, !is.na(htcb) & (height>htcb) & (htcb >0)
trees %>%
  filter(status_pre=='L' & ((is.na(htcb_m) | !(height_m > htcb_m) | !(htcb_m > 0))))
# if status_post is live, !is.na(torch_p, scorch_p, char_p)
trees %>%
  filter(status_post=='L' & (is.na(torch_p) | is.na(scorch_p) | is.na(char_p)))

# if status is dead, !is.na(decay) and decay in c(1:5)
trees %>%
  filter(status_pre=='D'& 
           (is.na(decay_pre) | !is.element(decay_pre, c('1', '2', '3', '4', '5'))))
# if status_post is dead, !is.na(decay_post) and decay_post in c(1:5)
trees %>% 
  filter(status_post=='D'&
           (is.na(decay_post) | !is.element(decay_post, c('1', '2', '3', '4', '5')))) %>%
  print(n = Inf, width = Inf)

# snags shouldn't have torch_ht or scorch_ht unless ht_meas is true; looks like the 
# crew frequently recorded torch and scorch height for new snags, which means they 
# did not record a decay_post; safe to assume decay_post = 1 for new snags, and that 
# the snag height is the same as the tree height
trees %>% filter(is.na(ht_meas)&status_post=='D'&(!is.na(torch_ht)|!is.na(scorch_ht))) %>% print(n = Inf, width = Inf)
trees %>% filter(status_pre=='L'&status_post=='D'&(decay_post!=1|is.na(decay_post))&is.na(ht_meas)) %>% print(n = Inf, width = Inf)
trees = 
  trees %>%
  mutate(
    decay_post = 
      ifelse(is.na(decay_post)&status_pre=='L'&status_post=='D',
             '1',
             decay_post),
    height_m_post = 
      ifelse(is.na(height_m_post)&status_pre=='L'&status_post=='D',
             height_m,
             height_m_post)
  )

# why do these snags not have a decay_post? looks like the crew just forgot to collect it
# note to self: don't double-purpose columns on paper datasheets
trees %>% filter(is.na(decay_post)&status_pre=='D'&status_post=='D') %>% print(n = Inf, width = Inf)
# just back fill them
trees  = 
  trees %>%
  mutate(decay_post = 
           ifelse(status_pre=='D'&status_post=='D'&is.na(decay_post),
                  decay_pre,
                  decay_post))


# some height_m_post too low; these should have status_post of X because 
# they're no longer in the trees/snags dataset
trees %>% filter(height_m_post < 1.37) %>% print(n = Inf, width = Inf)
trees = 
  trees %>% 
  mutate(status_post = 
           ifelse(height_m_post < 1.37 & !is.na(height_m_post),
                  'X',
                  status_post))

# check valid species
unique(trees$spp)[order(unique(trees$spp))]
valid_tree_spp = c('ABCO', 'CADE', 'QUKE', 'PILA', 'PIPO', 'PSME', 'NODE',
                   'SEGI', 'ARME', 'CHCH', 'CONU', 'UNK')

trees %>% filter(!is.element(spp, valid_tree_spp))
# 'U' is clearly unknown (and all of those are decayed snags). Not sure what 
# 'A' is; probably ABCO but could be ARME? No comment in metadata. Row above 
# it in the source spreadsheet is an ABCO, so either the crew expected it to 
# autofill (and it didnt) or it tried to autofill and they deleted it to 
# code ARME...
# going to assume ABCO
trees[trees$comp=='C' & trees$plot_id=='C-022' & trees$location_m == 9.6 & 
        trees$distance_m==0.2,
      'spp'] = 'ABCO'

trees %>% filter(!is.element(ts, c('NW', 'NE', 'SW', 'SE', 'EN', 'ES', 'WN', 'WS')))

# check for comments
trees %>% filter(!is.na(comment)) %>% print(n = Inf, width = Inf)

# looks like c-011-sw-4.8-5 and c-011-ws-5.2-4.4 are the same tree; remove one
trees = 
  trees %>%
  filter(tree_id != 'C-011-WS-5.2-4.4')

# check validity of decay_post
trees %>% filter(status_post=='D'&is.na(decay_post)) %>% print(n = Inf, width = Inf) # good after backfilling

# check validity of height_m_post; assume no change in height if its missing
trees %>% filter(status_post=='D'&is.na(height_m_post)) %>% print(n = Inf, width = Inf)
trees = 
  trees %>%
  mutate(height_m_post = 
           ifelse(is.na(height_m_post)&status_post=='D',
                  height_m,
                  height_m_post))


# check that torch_p+scorch_p <= 100
trees %>% filter(torch_p + scorch_p > 100) %>% print(n = Inf, width = Inf)
# torch ht is 100% of prefire height; call itt 100
trees[trees$tree_id=='B-026-NE-2-0.7','scorch_p'] = 0

# check NAs for torch_p, scorch_p, char_p, for trees which were alive prefire and 
# not missing postfire
trees %>% 
  filter(status_pre=='L'&!is.element(status_post, c('X', 'F'))&(is.na(scorch_p)|is.na(torch_p)|is.na(char_p))) %>% 
  print(n = Inf, width = Inf)
# A-025-NW-6.9-4.8 has scorch_p 100, torch_p should be 0
trees[trees$tree_id=='A-025-NW-6.9-4.8','torch_p'] = 0
# A-023-NW-10.6-2.3 has a comment indicating it was broken, so NA for scorch/torch makes sense
# A-020-WS-5.1-1.6 not sure what happened; data was entered correctly; crew wrote in '/' 
# in the space for torch/scorch percent, but 0 for char %; marked as live prefire 
# and dead postfire... like they treated it as if it were a snag prefire. 
# leave as NA

# check for usage of torch_ht, scorch_ht, char_ht for trees which werent supposed 
# to get that data (ht_meas = FALSE)
trees %>% 
  filter((!is.na(torch_ht)|!is.na(scorch_ht)|!is.na(char_ht))&is.na(ht_meas)) %>% 
  print(n = Inf, width = Inf) # idk wtf is going on here; don't ever do the double
# purpose columns thing again; most of these are from the first day where we 
# tried recording heights for every tree, and looks like there are a bunch of 
# snags where they recorded char height (without needing to)? I'm going to 
# trust all these

# check that torch_ht <= scorch_ht
trees %>%
  filter(!is.na(torch_ht)&!is.na(scorch_ht)&torch_ht>scorch_ht) %>%
  print(n = Inf, width = Inf) # set the status to X we'll toss it anways

#### consumption ###############################################################

head(consumption)

consumption = 
  consumption %>%
  mutate(plot_id = paste0(COMP, '-', str_pad(PLOT, width = 3, side = 'left', pad = 0)),
         transect_id = paste0(plot_id, '-', str_pad(AZ, width = 3, side = 'left', pad = 0))) %>%
  select(comp = COMP, plot_id, transect_id, az = AZ,
         location_m = LOCATION_m, consumption_cm = CONSUMPTION_cm, comment)

summary(consumption)

consumption %>% filter(is.na(consumption_cm)) %>% print(n = Inf) # these are all cases 
# which got a text code in the original data (either indicating that we couldnt find the pin, 
# or that the pin fell over or was disturbed, etc) and got coerced to NA; if we do want to 
# do something with the 'fell over' values need to reimport as text and then parse it, but 
# for now just leave these all as NA

# also just want to check how many pins we couldn't find...
read_excel(filepath,
           sheet = 'consumption',
           col_types = c('text', 'text', 'numeric', 'numeric', 'text', 'text')) %>%
  filter(CONSUMPTION_cm=='M') %>% 
  print(n = Inf) # 10 pins total, not too shabby; emailed ariel to give her a dheads up and apologize

# check comments
consumption %>% filter(!is.na(comment)) # looks good, NA is appropriate here

#### write results #############################################################

write.csv(consumption, here::here('02-data', '01-cleaned', 'consumption.csv'),
          row.names = FALSE)
write.csv(cwd, here::here('02-data', '01-cleaned', 'cwd_post.csv'),
          row.names = FALSE)
write.csv(fwd_diams, here::here('02-data', '01-cleaned', 'fwd_diams_post.csv'),
          row.names = FALSE)
write.csv(fwd_tallies, here::here('02-data', '01-cleaned', 'fwd_tallies_post.csv'),
          row.names = FALSE)
write.csv(litter_duff, here::here('02-data', '01-cleaned', 'litter_duff_post.csv'),
          row.names = FALSE)
write.csv(plot_metadata, here::here('02-data', '01-cleaned', 'plot_metadata_post.csv'),
          row.names = FALSE)
write.csv(shrubs, here::here('02-data', '01-cleaned', 'shrubs_post.csv'),
          row.names = FALSE)
write.csv(trees, here::here('02-data', '01-cleaned', 'trees_post.csv'),
          row.names = FALSE)
