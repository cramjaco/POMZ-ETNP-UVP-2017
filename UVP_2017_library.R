### Load stuff

library(transformr)
library(tidyverse)
library(mgcv)
library(broom)
library(furrr)
library(lubridate)
pass <- function(x, ...) {x}

ParticleSizeCutoff <- 0.5

### Function to load in initial data and prepare a data frame.

## Define global parameter

C_f = 133 # SmoothsAndFlusRevisited.Rmd, calculated from trap flux
ag = 2.00 # SmoothsAndFlusRevisited
alpha = (ag + 1) / 2 # assuming spherical drag profile
gamma = alpha - 1

C_f_global <- C_f
alpha_global <- alpha
gamma_global <- gamma
ag_global <- ag

# Get the times of the CTD casts
CTD_Unique <- read_csv("data/SKQ201617S_CTD_Profile.csv") %>%
  filter(Station == "sta016") %>%
  mutate(DateJ = mdy(`mon/dd/yyyy`), TimeJ = hms(`hh:mm`))  %>%
  mutate(DateTimeJ = mdy_hms(paste(`mon/dd/yyyy`, `hh:mm`))) %>%
  select(time = DateTimeJ, depth = `Pressure [db]`) %>%
  group_by(time) %>% summarize(MaxDepth = max(depth))

## Times on the UVP are the time it was switched on, not the time the cast started. I will find the next time of the cast
Find_Next_Cast <- function(LTime){
  NextTime <- min(CTD_Unique$time[CTD_Unique$time - LTime > 0]) # Find the next cast
  NextTime
}

Find_Next_Cast_2 <- Vectorize(Find_Next_Cast)

Find_Next_Cast_3 <- function(LTimeVec){
  NextVec <- Find_Next_Cast_2(LTimeVec)
  attributes(NextVec) <- attributes(LTimeVec)
  NextVec
}

bring_in_p2 <- function(){
  
  # bring in metadata that specifies relevant files
  uvpMeta <- read_tsv("data/uvpdata/export_detailed_20190304_23_14_Export_metadata_summary.tsv") %>% rename(time = `yyyy-mm-dd hh:mm`) %>% arrange(time)
  # just take p2 data
  uvpMetaP2 <- uvpMeta %>% filter(Site == "016")
  
  # bring in particle data
  uvp_data_path <- "data/uvpdata"
  particleData <- uvpMetaP2 %>% pull(`Particle filename`) %>%
    map(~read_tsv(file.path(uvp_data_path, .), locale = locale(encoding = "latin1", tz="UTC"))) %>%
    reduce(rbind)
  
  # some initial processing
  particleNumbers <- particleData %>%
    select(profile = Profile, time = `yyyy-mm-dd hh:mm`, depth = `Depth [m]`, vol = `Sampled volume[L]`, `LPM (102-128 µm)[#/L]`:`LPM (>26 mm)[#/L]`) %>%
    gather(key = "sizeclass", value = "nparticles", `LPM (102-128 µm)[#/L]`:`LPM (>26 mm)[#/L]`) %>%
    # convert to central time, originals are in UTC
    #mutate(time = lubridate::with_tz(time, tzone = "Mexico/General"))
    mutate(time = Find_Next_Cast_3(time)) %>%
    pass()
  
  #
  classData <- particleNumbers %>% select(sizeclass) %>% unique() %>%
    mutate(lb0 = as.numeric(str_extract(sizeclass,"(?<=\\().*(?=-)")),
           ub0 = str_extract(sizeclass, "(?<=-).*(?=\\s)") %>% as.numeric(),
           #unit0 = str_extract(sizeclass, "(?<=\\s).(?=m)"),
           ismm = str_detect(sizeclass, "mm"),
           lb1 = if_else(ismm, lb0, lb0 / 1000),
           ub1 = if_else(ismm, ub0, ub0 / 1000),
           lb2 = if_else(str_detect(sizeclass, "\\>26"),26, lb1),
           ub2 = if_else(str_detect(sizeclass, "\\>26"),32, ub1) # arbitrary
    ) %>%
    select(sizeclass, lb = lb2, ub = ub2)
  
  particleNumbers01 <- left_join(particleNumbers, classData, by = "sizeclass") %>%
    mutate(TotalParticles = nparticles * vol) %>%
    mutate(binsize = ub - lb) %>%
    mutate(n_nparticles = nparticles / binsize)
  
  particleNumbers001 <- particleNumbers01 %>% select(profile, time, depth, vol, sizeclass, lb, ub, binsize, TotalParticles, nparticles, n_nparticles)
  
  particleNumbers001
}

### Import P16 station 100

bring_in_p16_s100 <- function(){
  require(tidyverse)
  
  ## initial read in
  s100data00 <- read_tsv("data/p16uvp/export_detailed_20190908_00_44_PAR_p16n_100.tsv", locale = locale(encoding = "latin1"))
  
  ## just the data we need
  s100data01 <- s100data00 %>%
    select(profile = Profile, time = `yyyy-mm-dd hh:mm`, depth = `Depth [m]`, vol = `Sampled volume [L]`, `LPM (102-128 µm) [# l-1]`:`LPM (>26 mm) [# l-1]`) %>%
    gather(key = "sizeclass", value = "nparticles", `LPM (102-128 µm) [# l-1]`:`LPM (>26 mm) [# l-1]`)
  
  ## Convert bin names to values in microns
  renDf00 <- s100data01 %>% select(sizeclass) %>% unique()
  
  renDf01 <- renDf00 %>%
    mutate(
      ##bounds, ignoring microns or mm designation
      lb0 = str_extract(sizeclass, "(?<=\\().*(?=-)") %>% str_extract( ".*(?=-)") %>% as.numeric(),
      ub0 = str_extract(sizeclass, "(?<=-).*(?=\\s)") %>% str_extract(".*(?=\\s)") %>% str_extract(".*(?=\\s)") %>% as.numeric(),
      ## test for mm
      ismm = str_detect(sizeclass, "mm"),
      ## set everthing in microns
      lb1 = if_else(ismm, lb0, lb0 / 1000),
      ub1 = if_else(ismm, ub0, ub0 / 1000),
      ## deal with  biggist size bin
      
      lb2 = if_else(str_detect(sizeclass, "\\>26"),26, lb1),
      ub2 = if_else(str_detect(sizeclass, "\\>26"),32, ub1) # arbitrary
    ) %>%
    select(sizeclass, lb = lb2, ub = ub2)
  
  ## join the names
  s100data02 <- s100data01 %>% left_join(renDf01, by = "sizeclass")
  
  ## merge small bin with no particles wiht the bin one smaller
  s100data03 <- s100data02 %>%
    filter(lb != 0.128) %>% # removing small bin with no particles
    mutate(ub = if_else(ub == 0.128, 0.161, ub))
  
  ## Preliminary math
  s100data04 <- s100data03 %>%
    mutate(
      TotalParticles = nparticles * vol,
      binsize = ub - lb,
      n_nparticles = nparticles/binsize
    )
  
  s100data04
}

combine_projects <- function(proj1, proj2, name1 = "ETNP", name2 = "P16"){
  proj1 <- proj1 %>% mutate(project = name1)
  proj2 <- proj2 %>% mutate(project = name2)
  bothProj = bind_rows(proj1, proj2)
  bothProj
}

## I tend to work with lists of two element. One is the data of each size, and one is the depth summary.
## The following makes the early instance of that list, on whihc I work.

make_twin_df_list <- function(EachSize){
  
  DepthSummary <- EachSize %>% group_by(project, profile, time, depth) %>%
    summarize()
  list(EachSize, DepthSummary)
}


## A function that takes a variable x that can either be a list of two elements, "DepthSummary" and "EachSize" 
## and saves them to the parent environment. Alternatively, if DepthSummary is specified, EachSize is saved to the parent environment, and DepthSummary
## is passed forward
## I stopped using this because the passing data directly to the parent thing, while concise, I think may be confusing to some readers.
## The second approach is more conventional.

parse_jac_input <- function(x, DepthSummary = NULL){
  p <- parent.frame()
  if(!is.null(DepthSummary)){
    p$EachSize = x
    p$DepthSummary = DepthSummary
  }else{
    if(length(x) == 2){
      p$EachSize = x[[1]]
      p$DepthSummary = x[[2]]
    }else{
      stop("Either DepthSummary must be specified, or else the first variable must be a two element list containing EachSize, and DepthSummary elements")}
  }
  
}

# As above, but this just sends out a two element list of EachSize and DepthSummary

parse_jac_input2 <- function(x, DepthSummary){
  if(!is.null(DepthSummary)){
    return(list(ES = x, DS = DepthSummary))
  }else{
    if(length(x) == 2){
      return(x)
    }else{
      stop("Either DepthSummary must be specified, or else the first variable must be a two element list containing EachSize, and DepthSummary elements")}
  }
}
  
## I dont' actually use this for anything,but these twin functions all have common elements, so here's a template
function_template <- function(x, DepthSummary = NULL){
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  # Code does stuff here
  return(list(ES = EachSize, DS = DepthSummary))
}

###  Calculate biovolume, flux, speed, both for size classes and in aggregate
calc_particle_parameters <- function(x, DepthSummary = NULL){
  #Allow passing in either a two elemet list of Eachsize and DepthSummary, or passing in as two variables.

  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  ## Particle fractal dimension estimates
  # alpha = 1.7 # Klips et al. 1994
  # gamma = alpha - 1 # assume spherical drag profile as in Guidi
  # ag_fit = 0.26 # Fitted alpha + gamma parameter; hardcoded but copied from NormalizeUVP_Flux.Rmd
  # C_f_fit = 3.98 # Fitted C_f parameter
  
  # alpha = 0.52 # Alldredge
  # gamma = 0.26 # Alldredge & Gotschalk
  # C_f_fit = 10.51 # Normalize_UVP_Flux.Rmd, nonlinear for now
  # ag_fit = alpha + gamma # Zerod out for now; I'd like to clean this all up soon.
  
  EachSize2 <- EachSize %>% 
    mutate(
      biovolume = nparticles * lb ^ alpha,
      speed = lb ^ gamma,
      flux = biovolume * speed,
      flux_fit = nparticles * C_f_global * lb ^ ag_global
    )
  DepthSummary2 <- EachSize2 %>%
    group_by(project, profile, time, depth) %>%
    summarize(
      vol = first(vol),
      tot_TotParticles = sum(TotalParticles),
      tot_nparticles = sum(nparticles),
      tot_nnparticles = sum(n_nparticles),
      tot_biovolume = sum(biovolume),
      tot_flux = sum(flux),
      tot_flux_fit = sum(flux_fit),
      tot_speed = tot_flux/tot_biovolume
    )
  list(ES = EachSize2, DS = DepthSummary2)
}

## Seperate parameters for small < 530   micron and big > 530 micron particles
calc_small_and_big <- function(x, DepthSummary = NULL){
  #Allow passing in either a two elemet list of Eachsize and DepthSummary, or passing in as two variables.
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  small <- EachSize %>%
    filter(lb < ParticleSizeCutoff)
  
  big <- EachSize %>%
    filter(lb >= ParticleSizeCutoff)
  
  small2 <- small %>%
    group_by(project, profile, time, depth) %>%
    summarise(
      small_TotParticles = sum(TotalParticles),
      small_nparticles = sum(nparticles),
      small_nnparticles = sum(n_nparticles),
      small_biovolume = sum(biovolume),
      small_flux = sum(flux),
      small_flux_fit = sum(flux_fit),
      small_speed = small_flux/small_biovolume
    )
  
  big2 <- big %>%
    group_by(project, profile, time, depth) %>%
    summarise(
      big_TotParticles = sum(TotalParticles),
      big_nparticles = sum(nparticles),
      big_nnparticles = sum(n_nparticles),
      big_biovolume = sum(biovolume),
      big_flux = sum(flux),
      big_flux_fit = sum(flux_fit),
      big_speed = big_flux/big_biovolume
    )
  
  DepthSummary2 <- left_join(DepthSummary, small2, by = c("project","profile", "time", "depth")) %>%
    left_join(big2, by = c("project","profile", "time", "depth"))
  
  list(ES = EachSize, DS = DepthSummary2)
  
}


## Used in a few functions

## Poisson
fit_model = function(df) glm(TotalParticles ~ log(lb), offset = log(binsize * vol), data = df, family = "poisson")
## Neg Bin, doesn't work
fit_nb = function(df) MASS::glm.nb(TotalParticles ~ log(lb) + offset(log(vol * binsize)), data = df)
safe_fit_nb <- safely(fit_nb)

# I switched to gam because it can handle negative binomial. If one does this, one also needs to add the call `parametric = TRUE` to the tidy function
fit_nb_2 <- function(df) gam(TotalParticles ~ log(lb), offset = log(binsize * vol), data = df, family = "nb")

### Calculate particle size distribution, intercept and slope
## Takes EachSize, a data frame of particle size specific stuff, and Depth Summary, which is depth specific stuff
## Returns a list with the above.
calc_psd <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  
  
  psdCalc01 <- EachSize %>% 
    group_by(project, profile, time, depth) %>%
    nest() %>%
    mutate(model = map(data, fit_nb_2)) %>%
    mutate(tidied = map(model, tidy, parametric = TRUE)) %>%
    select(-data, -model) %>%
    unnest(tidied) %>%
    select(project, profile:estimate) %>%
    spread(key = "term", value = "estimate") %>%
    rename(icp = `(Intercept)`, psd = `log(lb)`)
  
  DepthSummary2 <- left_join(DepthSummary, psdCalc01, by = c("project","profile", "time", "depth"))
  
  list(ES = EachSize, DS = DepthSummary2)
  
}

# Get the particle size distribution of small particles
calc_small_psd <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  #fit_model = function(df) glm(TotalParticles ~ log(lb), offset = log(binsize * vol), data = df, family = "poisson")
  
  psdCalc01 <- EachSize %>% 
    filter(lb < ParticleSizeCutoff) %>%
    group_by(project, profile, time, depth) %>%
    nest() %>%
    mutate(model = map(data, fit_nb_2)) %>%
    mutate(tidied = map(model, tidy, parametric = TRUE)) %>%
    select(-data, -model) %>%
    unnest(tidied) %>%
    select(project, profile:estimate) %>%
    spread(key = "term", value = "estimate") %>%
    rename(small_icp = `(Intercept)`, small_psd = `log(lb)`)
  
  DepthSummary2 <- left_join(DepthSummary, psdCalc01, by = c("project","profile", "time", "depth"))
  
  list(ES = EachSize, DS = DepthSummary2)
}


# Get the particle size distribution of large particles
calc_big_psd <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  #fit_model = function(df) glm(TotalParticles ~ log(lb), offset = log(binsize * vol), data = df, family = "poisson")
  
  psdCalc01 <- EachSize %>% 
    filter(lb >= ParticleSizeCutoff) %>%
    group_by(project, profile, time, depth) %>%
    nest() %>%
    mutate(model = map(data, fit_nb_2)) %>%
    mutate(tidied = map(model, tidy, parametric = TRUE)) %>%
    select(-data, -model) %>%
    unnest(tidied) %>%
    select(project, profile:estimate) %>%
    spread(key = "term", value = "estimate") %>%
    rename(big_icp = `(Intercept)`, big_psd = `log(lb)`)
  
  DepthSummary2 <- left_join(DepthSummary, psdCalc01, by = c("project","profile", "time", "depth"))
  
  list(ES = EachSize, DS = DepthSummary2)
}

## Old Way
calc_psd_gam <- function(x, DepthSummary = NULL){
  # this currently slushes everything together. I need to nest and run over everything.
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  
  psd_gam_model <- gam(psd ~ s(depth), data = DepthSummary)
  intercept_gam_model <- gam(icp ~ s(depth), data = DepthSummary)
  
  psd_pred = predict(psd_gam_model, DepthSummary, se.fit = TRUE)
  icp_pred = predict(intercept_gam_model, DepthSummary, se.fit = TRUE)
  
  DepthSummary2 <- bind_cols(DepthSummary,
                                       psd_gam= psd_pred$fit, psd_seg = psd_pred$se.fit,
                                       icp_gam = icp_pred$fit, icp_seg = icp_pred$se.fit) 
  
  list(ES = EachSize, DS = DepthSummary2)
  
}


# Treat each cast seperately when calculating smooded particle size distributions

calc_psd_gam_multiprofile <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  DSN <- DepthSummary %>% group_by(project, profile, time) %>%
    nest()
  
  
  DSN$psd_gam_model <- map(DSN$data, ~ gam(psd ~ s(depth), data = .))
  
  
  DSN$intercept_gam_model <- map(DSN$data, ~ gam(icp ~ s(depth), data = .))
  
  DSN$psd_pred = map2(DSN$psd_gam_model, DSN$data, ~predict(.x, .y, se.fit = TRUE))
  DSN$icp_pred = map2(DSN$intercept_gam_model, DSN$data, ~predict(.x, .y, se.fit = TRUE))
  
  DSN$DepthSummary2 <- pmap(
    .l = list(DSN$data, DSN$psd_pred, DSN$icp_pred),
                            .f = function(DepthSummary, psd_pred, icp_pred){
    bind_cols(DepthSummary,
                                       psd_gam= psd_pred$fit, psd_seg = psd_pred$se.fit,
                                       icp_gam = icp_pred$fit, icp_seg = icp_pred$se.fit) 
  }
  )
  
  DS2 <- DSN %>%
    select(project, profile, time, DepthSummary2) %>% unnest(cols = c(DepthSummary2))
  
  mods = DSN$psd_gam_model
  names(mods) <- DSN$profile
  
  #list(out = list(ES = EachSize, DS = DS2), psd_gam_mod = mods)
  list(ES = EachSize, DS = DS2)
  
  }

# predict total particles for each size class from the gam
pred_tp_gam <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  EachSize2 <- DepthSummary %>% select(project, profile, time, depth, icp_gam, psd_gam) %>%
    right_join(EachSize, by = c("project","profile", "time", "depth")) %>%
    mutate(GamPredictTP = vol * lb * (exp(icp_gam + log(lb) * psd_gam))) %>% 
    select(-icp_gam, psd_gam)
  
  list(ES = EachSize2, DS = DepthSummary)
}
  
  

### Predicted quantiles of particle numbers from gam

## Helper functions
fit_model_secondary = function(df) glm(rp ~ log(lb), offset = log(binsize * vol), data = df, family = "poisson")
pos_fit_model_secondary = quietly(fit_model_secondary)

quantilater <- function(df, niter = 10, q1 = 0.025, q2 = 0.975){
  t2Test00 <- df
  t2Test01 <- t2Test00 %>% mutate(poisDraw = map(GamPredictTP, ~data.frame(iter = 1:niter, rp = rpois(niter, .)))) %>%
    unnest(poisDraw)
  t2Test02 <- t2Test01 %>% group_by(iter) %>% nest()
  t2Test03 <- t2Test02 %>% mutate(modelAndWarnings = map(data, pos_fit_model_secondary)) %>%
    mutate(model = map(modelAndWarnings, ~.[[1]]))
  t2Test04 <- t2Test03 %>% mutate(tidied = map(model, tidy))
  t2Test05 <- t2Test04 %>% select(iter,tidied)%>% unnest(tidied) %>% select(iter, estimate, term) %>% spread(value = "estimate", key = "term")
  t2Test06 <- t2Test05 %>% ungroup() %>% summarise(qpt05 = quantile(`log(lb)`, probs = 0.05), qpt95 = quantile(`log(lb)`, probs = .95))
  t2Test06

  }

## I don't do this analyis in the paper, but it was for deciding whether bins with zero particles are likely from some
## non power law particle size distribution function or not.
tp_quantiles <-  function(x, DepthSummary = NULL,  niter = 10){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  
  t2Calc01 <- EachSize %>% select(project, profile, time, depth, binsize, lb, vol, GamPredictTP) %>%
    mutate(vol2 = vol) %>%
    group_by(project, profile, time, depth) %>%
    nest(data = c(binsize, lb, vol, GamPredictTP))
  
  #t2Calc02 <- t2Calc01 %>% mutate(quantiles = future_map(data, quantilater, niter = niter))
  # the above doesn't work, trying an alternative
  future_map(t2Calc01[["data"]][], quantilater, niter = niter) -> moo
  t2Calc02 <- t2Calc01
  t2Calc02$quantiles = moo
  
  t2Calc03 <- t2Calc02 %>% select(-data) %>% unnest(quantiles)
  
  DepthSummary2 <- DepthSummary %>% left_join(t2Calc03, by = c("project","profile", "time", "depth"))
  
  return(list(ES = EachSize, DS = DepthSummary2))
  
  
}

### Binning
## Daniele's scheme
#unique(c(0:20:100, 100:25:200, 200:50:1000, 1000:100:2000, 2000:200:5000))

BianchiBins <- c(
  seq(from = 0, to = 100, by = 20),
  seq(from = 100, to = 200, by = 25),
  seq(from = 200, to = 1000, by = 50),
  seq(from = 1000, to = 2000, by = 100),
  seq(from = 2000, to = 5600, by = 200)
) %>% unique

## Go from highly resoved, to binned data
bin_depths <- function(x, DepthSummary = NULL, bins = BianchiBins){
  # Input from twin please.
  #Allow passing in either a two elemet list of Eachsize and DepthSummary, or passing in as two variables.

  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  dlb <- bins[1:length(bins)-1]
  dub <- bins[2:length(bins)]
  mids = (dlb + dub)/2
  
  EachSize2 <- EachSize %>%  mutate(DepthBin = cut(depth, bins, labels = mids))
  
  EachSize3 <- EachSize2 %>%
    select(-nparticles, -n_nparticles, -depth) %>%
    group_by(project, profile, time, DepthBin, lb, ub, binsize) %>%
    summarize(vol = sum(vol), TotalParticles = sum(TotalParticles)) %>%
    ungroup()
  
  # Recalculate nparticles and n_nparticles
  EachSize4 <- EachSize3 %>%
    mutate(nparticles = TotalParticles/vol,
           n_nparticles = nparticles/binsize) %>%
    mutate(depth = as.numeric(as.character(DepthBin))) %>%
    select(-DepthBin) %>%
    pass
  
  # Match DepthSummary
  
  DepthSummary2 <- DepthSummary %>%
    mutate(DepthBin = cut(depth, bins, labels = mids)) %>%
    group_by(project, profile, time, DepthBin) %>%
    summarize() %>%
    ungroup() %>%
    mutate(depth = as.numeric(as.character(DepthBin))) %>%
    select(-DepthBin)
    pass
  
  return(list(ES = EachSize4, DS = DepthSummary2))
}

# average profiles by summing TotalParticles and volume and then recalculating
sum_profiles <- function(x, DepthSummary = NULL){
  # Input from twin please.
  #Allow passing in either a two elemet list of Eachsize and DepthSummary, or passing in as two variables.

  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  # dlb <- bins[1:length(bins)-1]
  # dub <- bins[2:length(bins)]
  # mids = (dlb + dub)/2
  
  #EachSize2 <- EachSize %>%  mutate(DepthBin = cut(depth, bins, labels = mids))
  
  # combines the profiiles, loosing npartincles and n_nparticles
  EachSize3 <- EachSize %>%
    select(-nparticles, -n_nparticles) %>%
    group_by(project, depth, lb) %>%
    summarize(vol = sum(vol), TotalParticles = sum(TotalParticles), ub = first(ub), binsize = first(binsize)) %>%
    ungroup()
  
  # Recalculate nparticles and n_nparticles
  EachSize4 <- EachSize3 %>%
    mutate(nparticles = TotalParticles/vol,
           n_nparticles = nparticles/binsize) %>%
    #mutate(depth = as.numeric(as.character(DepthBin))) %>%
    #select(-DepthBin) %>%
    pass
  
  # Match DepthSummary
  
  DepthSummary2 <- DepthSummary %>%
    #mutate(DepthBin = cut(depth, bins, labels = mids)) %>%
    group_by(project, time, depth) %>%
    summarize() %>%
    ungroup() %>%
    mutate(profile = "multiple")
    #mutate(depth = as.numeric(as.character(DepthBin))) %>%
    #select(-DepthBin)
    pass
  
  return(list(ES = EachSize4, DS = DepthSummary))
}



# combine_timesteps <- function(x, DepthSummary = NULL){
#   # Load in data, flexably
#   x2 <- parse_jac_input2(x, DepthSummary)
#   EachSize = x2[[1]]
#   DepthSummary = x2[[2]]
#   
#   # consistency with other expressions
#   EachSize2 <- EachSize
#   
#   EachSize3 <- EachSize2 %>%
#     select(-nparticles, -n_nparticles) %>%
#     group_by(project, profile, time, DepthBin, lb, ub, binsize) %>%
#     summarize(vol = sum(vol), TotalParticles = sum(TotalParticles)) %>%
#     ungroup()
#   
# }

## 24 November 2020

# my_double_gam <- function(df){
#   gam(TotalParticles ~s(log(lb), log(depth)), offset = log(vol * binsize), family = nb(), data = df)
# }
# 
# safe_double_gam <- safely(my_double_gam)

## The following functions are never used as far as I can tell, but if I delete them, they'll probably turn
## out to be required somewhere cryptic and break everything, so here they stay.
## I'm not sure what they do

expand_with_gam <- function(df, mod){
  loc_pred <- predict(mod, type = "link", se.fit = TRUE) %>% as.data.frame %>% mutate(lower = fit - 2 * se.fit, upper = fit + 2 * se.fit) %>% mutate(resp_fit = exp(fit), resp_lower = exp(lower), resp_upper = exp(upper))
  loc_df <- bind_cols(df, loc_pred)
  loc_df
}

gam_size_ggplot <- function(df){
  ggplot(df, aes(x = lb)) + geom_point(aes(y = resp_fit), shape = 1) +
  geom_errorbar(aes(ymin = resp_lower, ymax = resp_upper)) + 
  geom_point(aes(y = n_nparticles)) + scale_x_log10() + scale_y_log10()
}

gam_size_ggplot_2d <- function(df){
  df %>% ggplot(aes(x = resp_fit, y = depth, col = log(lb), group = lb)) +
    scale_y_reverse() + geom_point() + scale_x_log10(limits = c(10^-8, NA)) +
    scale_color_viridis_c() + geom_path() + geom_vline(xintercept = 1) + geom_vline(xintercept = 5)  +
    geom_errorbar(aes(xmin = resp_lower, xmax = resp_upper), width = 10, alpha = 0.5) + theme_bw()
}

nnp_size_ggplot_2d <- function(df){
  df %>% ggplot(aes(x = n_nparticles, y = depth, col = log(lb), group = lb)) +
    scale_y_reverse() + geom_point() + scale_x_log10(limits = c(10^-8, NA)) +
    scale_color_viridis_c() + geom_path() + geom_vline(xintercept = 1) + geom_vline(xintercept = 5)  +
     theme_bw()
}

##Here's a commented out function. I think a newer version is defined below.

# double_gam_smooth <- function(x, DepthSummary = NULL){
#   # Input from twin please.
#   #Allow passing in either a two elemet list of Eachsize and DepthSummary, or passing in as two variables.
# 
#   x2 <- parse_jac_input2(x, DepthSummary)
#   EachSize = x2[[1]]
#   DepthSummary = x2[[2]]
#   
#   withGamFit <- EachSize %>% group_by(project) %>% nest() %>%
#     mutate(mod = map(data, safe_double_gam),
#            modOnly = map(mod, ~.[[1]]),
#            pred = map2(modOnly, data, safely(predict), se.fit = TRUE),
#            predOnly = map(pred, ~.[[1]]),
#            data01 = map2(data, predOnly,
#                          ~bind_cols(.x, link = .y$fit, lse = .y$se.fit))) %>%
#     select(project, data01) %>%
#     unnest(data01) %>%
#     mutate(link_lower = link - lse,
#            link_upper = link + lse,
#            nnp_smooth = exp(link),
#            nnp_lower = exp(link_lower),
#            nnp_upper = exp(link_upper),
#            np_smooth = nnp_smooth * binsize,
#            tp_smooth = np_smooth * vol,
#            flux_smooth = np_smooth * (C_f_global * lb ^ ag_global)
#     )
#   
#   TotalStuff <- withGamFit %>% group_by(project, profile, time, depth) %>%
#     summarize(smooth_TotParticles = sum(tp_smooth),
#            smooth_nparticles = sum(np_smooth),
#            smooth_nnparticles = sum(nnp_smooth),
#            smooth_flux_fit = sum(flux_smooth)
#            )
#   
#   DepthSummary_B <- left_join(DepthSummary, TotalStuff, by = c("project", "profile", "time", "depth"))
#   
#   return(list(ES = withGamFit, DS = DepthSummary_B))
#   
# }

## Functions for smoothing data, using size and depth informaiton together
## originally coded in SmoothsAndFluxRevisited
## I loose at DRY coding, but also if I try to fix this, it will break everything and I'll spend a week fixing things.
my_double_gam <- function(df){
  gam(TotalParticles ~s(log(lb), log(depth), by = factor(profile)), offset = log(vol * binsize), family = nb(), data = df)
}

safe_double_gam <- safely(my_double_gam)


double_gam_smooth <- function(x, DepthSummary = NULL){
  # Input from twin please.
  #Allow passing in either a two elemet list of Eachsize and DepthSummary, or passing in as two variables.

  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  
  withGamFit <- EachSize %>% group_by(project) %>% nest() %>%
    mutate(mod = map(data, safe_double_gam),
           modOnly = map(mod, ~.[[1]]),
           pred = map2(modOnly, data, safely(predict), se.fit = TRUE),
           predOnly = map(pred, ~.[[1]]),
           data01 = map2(data, pred,
                         ~bind_cols(.x, link = .y$result$fit, lse = .y$result$se.fit))) %>%
    select(project, data01) %>%
    unnest(data01) %>%
    mutate(link_lower = link - lse,
           link_upper = link + lse,
           nnp_smooth = exp(link),
           nnp_lower = exp(link_lower),
           nnp_upper = exp(link_upper),
           np_smooth = nnp_smooth * binsize,
           tp_smooth = np_smooth * vol,
           flux_smooth = np_smooth * (C_f_global * lb ^ ag_global)
    )
  
  TotalStuff <- withGamFit %>% group_by(project, profile, time, depth) %>%
    summarize(smooth_TotParticles = sum(tp_smooth),
           smooth_nparticles = sum(np_smooth),
           smooth_nnparticles = sum(nnp_smooth),
           smooth_flux_fit = sum(flux_smooth)
           )
  
  DepthSummary_B <- left_join(DepthSummary, TotalStuff, by = c("project", "profile", "time", "depth"))
  
  return(list(ES = withGamFit, DS = DepthSummary_B))
  
}

## Sometimes I want to run things on just one profile, especially when building functions that I will later
## run over every profile. This filters out just one profile (aka cast)
filter_profile <- function(x, DepthSummary = NULL, profile = "stn_043"){
  prof2 <- profile
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]] %>% ungroup()
  DepthSummary = x2[[2]] %>% ungroup()
  
  EachSize <- EachSize %>% filter(profile == prof2) %>% select(-c(project, profile, time))
  DepthSummary <- DepthSummary %>% filter(profile == prof2) %>% select(-c(project, profile, time))
  
  return(list(ES = EachSize, DS = DepthSummary))
}

## This calls the eularian prism model and calculates a bunch of relevant parameters
diagnose_disaggregation_one_profile <- function(x, DepthSummary = NULL){
  ## Preamble
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  lb_vec <- sort(unique(EachSize$lb))
  m_vec =  Cm * lb_vec ^ alpha;
  w_vec = Cw * lb_vec ^ gamma;
  llb_01 <- little_lb <- lb_vec[1] - (lb_vec[2] - lb_vec[1])/2
  
  specData <- EachSize %>% select(depth, lb, np_smooth, nnp_smooth, flux_smooth) %>%
    nest(spec_meta = c(lb, np_smooth, nnp_smooth, flux_smooth))
  
  preparedData <- DepthSummary %>% left_join(specData, by = c("depth")) %>%
    arrange(depth) %>%
    mutate(spec_only = map(spec_meta, ~pull(., np_smooth)),
           spec_prev = lag(spec_only),
           flux_prev = lag(smooth_flux_fit),
           DF = smooth_flux_fit - flux_prev,
           #DFP = 1 - DF/flux_prev,  # I was using this for a while.
           DFP = smooth_flux_fit/flux_prev,
           depth_prev = lag(depth),
           DZ = depth - depth_prev,
    )
  
  minDepth = min(preparedData$depth)
  
  saveFirstDepth = preparedData %>% filter(depth == minDepth) %>% select(depth, spec_meta) %>% unnest(spec_meta) 
  
  modelRun <- preparedData %>%
    .[-1,] %>%
    # fix flux leak here
    mutate(use_this_DFP = map2_dbl(spec_prev, DFP, optFun, lbv = lb_vec, mv = m_vec, wv = w_vec, llb = llb_01, alpha = alpha_global, gamma = gamma_global)) %>%
    mutate(spec_pred = map2(spec_prev, use_this_DFP, remin_smooth_shuffle, lbv = lb_vec, mv = m_vec, wv = w_vec, llb = llb_01, alpha = alpha_global, gamma = gamma_global))
  
  #modelRunFixLine1 <- bind_rows(preparedData[1,], modelRun)

  modelConcise <- modelRun %>%
    mutate(spec_meta = map2(spec_meta, spec_prev, ~tibble(.x, np_prev = .y))) %>%
    mutate(spec_meta = map2(spec_meta, spec_pred, ~tibble(.x, np_pred = .y))) %>%
    select(depth, depth_prev, DZ, DF, DFP, use_this_DFP, spec_meta)
  
  modelUnnest <- modelConcise %>%
    #select(depth, spec_meta) %>%
    unnest(spec_meta) %>%
    ungroup()
  
  modelUnnestWithFirstDepth <- bind_rows(saveFirstDepth, modelUnnest)
  
  modelPostCalc <- modelUnnestWithFirstDepth %>%
    mutate(
      flux_prev = np_prev * (C_f_global * lb ^ ag_global),
      flux_pred = np_pred * (C_f_global * lb ^ ag_global)
    )
  
  Tot <- modelPostCalc %>% 
    group_by(depth) %>%
    summarize(depth_prev = first(depth_prev), DZ = first(DZ) ,DF = first(DF), DFP = first(DFP), use_this_DFP =  first(use_this_DFP),
              Flux = sum(flux_smooth), 
              Flux_Prev = sum(flux_prev),
              Flux_Pred = sum(flux_pred))
  
  Small <- modelPostCalc %>% 
    filter(lb <= ParticleSizeCutoff) %>%
    group_by(depth) %>%
    summarize(DF = first(DF), DFP = first(DFP), 
              Flux = sum(flux_smooth), 
              Flux_Prev = sum(flux_prev),
              Flux_Pred = sum(flux_pred))
  Big <- modelPostCalc %>% 
    filter(lb > ParticleSizeCutoff) %>%
    group_by(depth) %>%
    summarize(DF = first(DF), DFP = first(DFP), 
              Flux = sum(flux_smooth), 
              Flux_Prev = sum(flux_prev),
              Flux_Pred = sum(flux_pred))
  
  All <- Tot %>%
    left_join(Small, by = "depth", suffix = c("", "_Small")) %>%
    left_join(Big, by = "depth", suffix = c("", "_Big")) %>%
    mutate(osps = Flux_Small - Flux_Pred_Small,
           obpb = Flux_Pred_Big - Flux_Big, 

           ospsDZ = osps/DZ
           )
  
    DepthSummary_B <- DepthSummary %>%
      left_join(All, by = "depth") %>% rename(Flux_Smooth = Flux)
    
    modelReduced <- modelPostCalc %>% select(
      depth, lb, flux_prev, flux_pred
    )
    
    EachSize_B <- EachSize %>%
      left_join(modelReduced, by = c("depth", "lb"))
  # Code does stuff here
  return(list(ES = EachSize_B, DS = DepthSummary_B))
}

## Sometimes the eularian prism code breaks, in which case, I want to keep running, so I run it "safely"
## with this purr `safely` function
diagnose_disaggregation_one_profile_safe <- safely(diagnose_disaggregation_one_profile)

## Run diagnose_disaggregation_one_profile_safe over every profile and keep track of data
diagnose_disaggregation<- function(x, DepthSummary = NULL){
  ## Preamble
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  ESN <- EachSize %>% 
    group_by(project, profile, time) %>%
    nest() %>%
    rename(ES  = data)
  
  DSN <- DepthSummary %>%
    group_by(project, profile, time) %>%
    nest() %>%
    rename(DS = data)
  
  metaNest <- left_join(ESN, DSN, by = c("project", "profile", "time"))
  
  metaNest <- metaNest %>% mutate(ESDS = map2(ES, DS, ~list(.x, .y)))
  
  metaNest <- metaNest %>%
    mutate(ESDS_Mod_Safe = map(ESDS, diagnose_disaggregation_one_profile_safe))
  
  metaNest01 <- metaNest %>% mutate(ESDS_Mod = map(ESDS_Mod_Safe, ~.[[1]]),
                                    ESDS_Err = map(ESDS_Mod_Safe, ~.[[2]]),
                                    ES01 = map(ESDS_Mod, ~.[[1]]),
                                    DS01 = map(ESDS_Mod, ~.[[2]]))
  

  EachSize01 <- metaNest01 %>% select(project, profile, time, ES01) %>% unnest(ES01)
  DepthSummary01 <- metaNest01 %>% select(project, profile, time, DS01) %>% unnest(DS01)

  return(list(ES = EachSize01, DS = DepthSummary01))
}

## Functions about fitting flux

fit_flux_es <- function(C_f, ag, ES){
ES2 <- ES %>% mutate(flux2 = C_f * nparticles * lb ^ ag)
ES2
}

fit_flux_ds <- function(C_f, ag, ES, DS){
  DS2 <- ES %>%
    group_by(project, profile, time, depth) %>%
    summarize(tot_flux2 = sum(flux2))
  
  DS2 <- left_join(DS, DS2, by = c("project", "profile", "time", "depth"))
}

fit_flux <- function(C_f, ag, ES, DS){
  ES2 <- fit_flux_es(C_f, ag, ES)
  DS2 <- fit_flux_ds(C_f, ag, ES2, DS)
  return(list(ES = ES2, DS = DS2))
}

RMSE <- function(mod){
  RSS <- c(crossprod(mod$residuals))
  MSE = RSS/ length(mod$residuals)
  RMSE = sqrt(MSE)
}

flux_check <- function(DS){
  #diff = log10(DS$tn_flux) - log10(DS$tot_flux2)
  diff = log(DS$tn_flux) - log(DS$tot_flux2)
  squares = diff ^2
  rss = sum(squares)
  mse = rss/length(rss)
  rmse = sqrt(mse)
  rmse
}

fit_check_flux <- function(C_f, ag, ES, DS){
  FF <- fit_flux(C_f, ag, ES, DS)
  FC <- flux_check(FF$DS)
  FC
}

fc_wrap <- function(x, ES, DS){
  C_f <- x[1]
  ag <- x[2]
  FC <- fit_check_flux(C_f, ag, ES, DS)
  FC
}

## Recode time variable

# add_time_data <- function(x, DepthSummary = NULL){
#   require(chron)
#   require(lubridate)
#   x2 <- parse_jac_input2(x, DepthSummary)
#   EachSize = x2[[1]]
#   DepthSummary = x2[[2]]
# 
#   timeDf <- tibble(time = unique(DepthSummary$time)) %>% 
#     mutate(tod <- times(strftime(time,"%H:%M:%S")))
#   
#   # Recode into blocks of time
#   
#   # create breaks
#   breaks <- hour(hm("21:00", "5:00", "9:00", "18:00", "20:59"))
#   # labels for the breaks
#   labels <- c("Night", "Morning", "Afternoon", "Evening")
# 
#   timeDf <- timeDf %>%
#     mutate(timeBlock = cut(x=hour(time), breaks = breaks, labels = labels, include.lowest=TRUE)
# )
#   
#   # hours from noon
#   
#   hour(timeDf$tod)
#   
#   
#   
#   return(list(ES = EachSize, DS = DepthSummary))
# }