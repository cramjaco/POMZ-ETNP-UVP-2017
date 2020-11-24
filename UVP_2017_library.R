### Load stuff

library(transformr)
library(tidyverse)
library(mgcv)
library(broom)
library(furrr)
pass <- function(x, ...) {x}

### Function to load in initial data and prepare a data frame.

bring_in_p2 <- function(){
  
  # bring in metadata that specifies relevant files
  uvpMeta <- read_tsv("data/uvpdata/export_detailed_20190304_23_14_Export_metadata_summary.tsv") %>% rename(time = `yyyy-mm-dd hh:mm`) %>% arrange(time)
  # just take p2 data
  uvpMetaP2 <- uvpMeta %>% filter(Site == "016")
  
  # bring in particle data
  uvp_data_path <- "data/uvpdata"
  particleData <- uvpMetaP2 %>% pull(`Particle filename`) %>%
    map(~ read_tsv(file.path(uvp_data_path, .), locale = locale(encoding = "latin1"))) %>%
    reduce(rbind)
  
  # some initial processing
  particleNumbers <- particleData %>%
    select(profile = Profile, time = `yyyy-mm-dd hh:mm`, depth = `Depth [m]`, vol = `Sampled volume[L]`, `LPM (102-128 µm)[#/L]`:`LPM (>26 mm)[#/L]`) %>%
    gather(key = "sizeclass", value = "nparticles", `LPM (102-128 µm)[#/L]`:`LPM (>26 mm)[#/L]`)
  
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

### 

make_twin_df_list <- function(EachSize){
  
  DepthSummary <- EachSize %>% group_by(profile, time, depth) %>%
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
  
  alpha = 0.52 # Alldredge
  gamma = 0.26 # Alldredge & Gotschalk
  C_f_fit = 10.5 # Normalize_UVP_Flux.Rmd, nonlinear for now
  ag_fit = alpha + gamma # Zerod out for now; I'd like to clean this all up soon.
  
  EachSize2 <- EachSize %>% 
    mutate(
      biovolume = nparticles * lb ^ alpha,
      speed = lb ^ gamma,
      flux = biovolume * speed,
      flux_fit = nparticles * C_f_fit * lb ^ ag_fit
    )
  DepthSummary2 <- EachSize2 %>%
    group_by(profile, time, depth) %>%
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

## Seperate parameters for small < 53 micron and big > 53 micron particles
calc_small_and_big <- function(x, DepthSummary = NULL){
  #Allow passing in either a two elemet list of Eachsize and DepthSummary, or passing in as two variables.
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  small <- EachSize %>%
    filter(lb < 0.53)
  
  big <- EachSize %>%
    filter(lb >= 0.53)
  
  small2 <- small %>%
    group_by(profile, time, depth) %>%
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
    group_by(profile, time, depth) %>%
    summarise(
      big_TotParticles = sum(TotalParticles),
      big_nparticles = sum(nparticles),
      big_nnparticles = sum(n_nparticles),
      big_biovolume = sum(biovolume),
      big_flux = sum(flux),
      big_flux_fit = sum(flux_fit),
      big_speed = big_flux/big_biovolume
    )
  
  DepthSummary2 <- left_join(DepthSummary, small2, by = c("profile", "time", "depth")) %>%
    left_join(big2, by = c("profile", "time", "depth"))
  
  list(ES = EachSize, DS = DepthSummary2)
  
}


### Calculate particle size distribution, intercept and slope
## Takes EachSize, a data frame of particle size specific stuff, and Depth Summary, which is depth specific stuff
## Returns a list with the above.
calc_psd <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  fit_model = function(df) glm(TotalParticles ~ log(lb), offset = log(binsize * vol), data = df, family = "poisson")
  
  psdCalc01 <- EachSize %>% 
    group_by(profile, time, depth) %>%
    nest() %>%
    mutate(model = map(data, fit_model)) %>%
    mutate(tidied = map(model, tidy)) %>%
    select(-data, -model) %>%
    unnest(tidied) %>%
    select(profile:estimate) %>%
    spread(key = "term", value = "estimate") %>%
    rename(icp = `(Intercept)`, psd = `log(lb)`)
  
  DepthSummary2 <- left_join(DepthSummary, psdCalc01, by = c("profile", "time", "depth"))
  
  list(ES = EachSize, DS = DepthSummary2)
  
}

calc_small_psd <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  fit_model = function(df) glm(TotalParticles ~ log(lb), offset = log(binsize * vol), data = df, family = "poisson")
  
  psdCalc01 <- EachSize %>% 
    filter(lb < 0.53) %>%
    group_by(profile, time, depth) %>%
    nest() %>%
    mutate(model = map(data, fit_model)) %>%
    mutate(tidied = map(model, tidy)) %>%
    select(-data, -model) %>%
    unnest(tidied) %>%
    select(profile:estimate) %>%
    spread(key = "term", value = "estimate") %>%
    rename(small_icp = `(Intercept)`, small_psd = `log(lb)`)
  
  DepthSummary2 <- left_join(DepthSummary, psdCalc01, by = c("profile", "time", "depth"))
  
  list(ES = EachSize, DS = DepthSummary2)
}

calc_big_psd <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  fit_model = function(df) glm(TotalParticles ~ log(lb), offset = log(binsize * vol), data = df, family = "poisson")
  
  psdCalc01 <- EachSize %>% 
    filter(lb >= 0.53) %>%
    group_by(profile, time, depth) %>%
    nest() %>%
    mutate(model = map(data, fit_model)) %>%
    mutate(tidied = map(model, tidy)) %>%
    select(-data, -model) %>%
    unnest(tidied) %>%
    select(profile:estimate) %>%
    spread(key = "term", value = "estimate") %>%
    rename(big_icp = `(Intercept)`, big_psd = `log(lb)`)
  
  DepthSummary2 <- left_join(DepthSummary, psdCalc01, by = c("profile", "time", "depth"))
  
  list(ES = EachSize, DS = DepthSummary2)
}

calc_psd_gam <- function(x, DepthSummary = NULL){
  
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

# predict total particles for each size class from the gam
pred_tp_gam <- function(x, DepthSummary = NULL){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  EachSize2 <- DepthSummary %>% select(profile, time, depth, icp_gam, psd_gam) %>%
    right_join(EachSize, by = c("profile", "time", "depth")) %>%
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

tp_quantiles <-  function(x, DepthSummary = NULL,  niter = 10){
  
  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  
  t2Calc01 <- EachSize %>% select(profile, time, depth, binsize, lb, vol, GamPredictTP) %>%
    mutate(vol2 = vol) %>%
    group_by(profile, time, depth) %>%
    nest(data = c(binsize, lb, vol, GamPredictTP))
  
  #t2Calc02 <- t2Calc01 %>% mutate(quantiles = future_map(data, quantilater, niter = niter))
  # the above doesn't work, trying an alternative
  future_map(t2Calc01[["data"]][], quantilater, niter = niter) -> moo
  t2Calc02 <- t2Calc01
  t2Calc02$quantiles = moo
  
  t2Calc03 <- t2Calc02 %>% select(-data) %>% unnest(quantiles)
  
  DepthSummary2 <- DepthSummary %>% left_join(t2Calc03, by = c("profile", "time", "depth"))
  
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
  seq(from = 2000, to = 5000, by = 200)
) %>% unique

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
    group_by(profile, time, DepthBin, lb, ub, binsize) %>%
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
    group_by(profile, time, DepthBin) %>%
    summarize() %>%
    ungroup() %>%
    mutate(depth = as.numeric(as.character(DepthBin))) %>%
    select(-DepthBin)
    pass
  
  return(list(ES = EachSize4, DS = DepthSummary2))
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
#     group_by(profile, time, DepthBin, lb, ub, binsize) %>%
#     summarize(vol = sum(vol), TotalParticles = sum(TotalParticles)) %>%
#     ungroup()
#   
# }

## 24 November 2020

my_double_gam <- function(df){
  gam(TotalParticles ~s(log(lb), log(depth)), offset = log(vol * binsize), family = nb(), data = df)
}

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