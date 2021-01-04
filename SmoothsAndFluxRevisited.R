

# Setup
source("UVP_2017_library.R")
source("ModelStuff.R")
options(readr.default_locale=readr::locale(tz="Mexico/General"))

dataP2 <- bring_in_p2()
dataP16S100 <- bring_in_p16_s100()

dataBoth <- combine_projects(dataP2, dataP16S100) %>% filter(profile != "stn_041") # Station 041 has only one depth, so remvoed

  ## We will always have a data with information about each size, and a summary of each depth

twin01 <- make_twin_df_list(dataBoth)
ES01 <- twin01[[1]]
DS01 <- twin01[[2]]

binnedSFR <- bin_depths(twin01)
bSe <- binnedSFR$ES
bSd <- binnedSFR$DS

## Temp Library


my_double_gam_v2 <- function(df){
  gam(TotalParticles ~s(log(lb), log(depth), by = factor(profile)), offset = log(vol * binsize), family = nb(), data = df)
}

safe_double_gam_v2 <- safely(my_double_gam_v2)


double_gam_smooth_v2 <- function(x, DepthSummary = NULL){
  # Input from twin please.
  #Allow passing in either a two elemet list of Eachsize and DepthSummary, or passing in as two variables.

  x2 <- parse_jac_input2(x, DepthSummary)
  EachSize = x2[[1]]
  DepthSummary = x2[[2]]
  
  
  withGamFit <- EachSize %>% group_by(project) %>% nest() %>%
    mutate(mod = map(data, safe_double_gam_v2),
           modOnly = map(mod, ~.[[1]]),
           pred = map2(modOnly, data, safely(predict), se.fit = TRUE),
           predOnly = map(pred, .[[1]]),
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
           #flux_smooth = np_smooth * (C_f_global * lb ^ ag_global)
    )
  
  TotalStuff <- withGamFit %>% group_by(project, profile, time, depth) %>%
    summarize(smooth_TotParticles = sum(tp_smooth),
           smooth_nparticles = sum(np_smooth),
           smooth_nnparticles = sum(nnp_smooth),
           #smooth_flux_fit = sum(flux_smooth)
           )
  
  DepthSummary_B <- left_join(DepthSummary, TotalStuff, by = c("project", "profile", "time", "depth"))
  
  return(list(ES = withGamFit, DS = DepthSummary_B))
  
}

## Testing
test <- binnedSFR %>% double_gam_smooth_v2()

allETNP <-  binnedSFR$ES %>% filter(project == "ETNP")
# 
# 
testGam <- gam(TotalParticles ~s(log(lb), log(depth), by = factor(profile)), offset = log(vol * binsize), family = nb(), data = allETNP)

# test2Gam <- gam(TotalParticles ~s(log(lb), log(depth), by = factor(profile)), offset = log(vol * binsize), family = nb(),
#                 data = binnedSFR$ES %>% filter(project == "P16"))
# 
# test3Gam <- gam(TotalParticles ~s(log(lb), log(depth), by = factor(profile)), offset = log(vol * binsize), family = nb(),
#                 data = binnedSFR$ES )

bSe2 <- test$ES
bSd2 <- test$DS

# horizontalGamPlot <- dataGamHorizontal %>% ggplot(aes(x = resp_fit, y = depth, col = log(lb), group = lb)) +
# scale_y_reverse() + geom_point() + scale_x_log10(limits = c(10^-8, NA)) + scale_color_viridis_c() + geom_path() +
# geom_vline(xintercept = 1) + geom_vline(xintercept = 5)  + geom_errorbar(aes(xmin = resp_lower, xmax = resp_upper), width = 10, alpha = 0.5)+ theme_bw()

bSe2 %>%
  filter(project == "ETNP") %>%
  #filter(profile == "stn_032") %>%
  ggplot(aes(y = depth, x = nnp_smooth, xmin = nnp_lower, xmax = nnp_upper, col = log(lb), group = lb)) +
  scale_y_reverse(limits = c(500, 0)) + scale_x_log10() +
  geom_point() + geom_errorbar(width = 10) + geom_path() +
  scale_color_viridis_c() +
  facet_wrap(~profile)
# hmm, these actuall ylook a bit different

bSe2 %>%
  filter(project == "ETNP") %>%
  #filter(profile == "stn_032") %>%
  ggplot(aes(y = depth, x = TotalParticles, xmin = nnp_lower, xmax = nnp_upper, col = log(lb), group = lb)) +
  scale_y_reverse(limits = c(500, 0)) + scale_x_log10(limits = c(1e-3, 1e5)) +
  geom_point()  + geom_path() +
  scale_color_viridis_c() +
  facet_wrap(~profile)

# these look different-ish
# some real large particles are somehow larger than slightly smaller bins?

# for actual fitting, I should perhaps sum the etnp data and then fit a smooth on that
# and then fit to the trap flux?
# rationalle being that the flux averages over the week
# but for analyiss I shoudl show each week, I think.

trapFlux <- read_csv("dataOut/fluxMS_distilled.csv")
# There are some low traps in the surface that are anomalously low. How about, if depth is less than 250, and flux is less than 30, remove it
trapFlux2 <- trapFlux %>%
  filter(!(Depth < 250 & C_flux_umol < 30))

# Problem last time I think is that I assumed flux was a power law and tried to fit to that. But actually it increases at depth
# There's lots of good code in Normalize_UVP_Flux.Rmd
# here's the plan to determine C_f and ag in f = C_f * lb ^ ag
# pull from the smooth, estemates of numbers of particles where each trap is
# I'm not sure whether to sum the particles and the volumes sampled, or whether to just average nparticles,
# or to throw all of the profiles in together as seperate data points
# Just wrote combine profiles, which will get us an essentially average profile.
# So we combine, then calculate a smooth (not main smoothing function), then predict out the values where the traps are.
# Then we run a solver on that data frame to minimize distance from traps, using the strategy from Normalize_UVP_Flux.Rmd

forFlux <- binnedSFR %>%  sum_profiles()
forFlux$ES <- forFlux$ES %>% filter(project == "ETNP")
forFlux$DS <- forFlux$DS %>% filter(project == "ETNP")

gamESE <- gam(TotalParticles ~s(log(lb), log(depth)), offset = log(vol * binsize), family = nb(), data = forFlux$ES)

lb_vvv <- sort(unique(forFlux$ES$lb))
binsize_vvv <-sort(unique(forFlux$ES$binsize)) 
lbbs <- bind_cols(lb = lb_vvv, binsize = binsize_vvv)
  
tfD = sort(unique(trapFlux2$Depth))

Expanded <- expand_grid(depth = tfD, lb = lb_vvv) %>% left_join(lbbs, by = "lb")

Predicted <- predict(gamESE, Expanded, type = "link")

EP <- bind_cols(Expanded, link = Predicted)

EP <- EP %>%
  mutate(n_nparticles = exp(link),
         nparticles = n_nparticles * binsize) %>%
  mutate(profile = "multiple", project = "ETNP", time = "togetawatch")

# DP <- EP %>% group_by(depth) %>% summarize() %>% 
#   mutate(profile = "multiple", project = "ETNP", time = "togetawatch")

DP <- trapFlux2 %>% select(depth = Depth, tn_flux = C_flux_umol) %>% 
  mutate(profile = "multiple", project = "ETNP", time = "togetawatch")

PETE <- list(ES = EP, DS = DP)



fit_check_flux(C_f_global, ag_global, ES = PETE$ES, DS = PETE$DS)

opt <- optim(c(C_f_global, ag_global), fc_wrap, ES = PETE$ES, DS = PETE$DS)
opt
# gives non-sensical values

oFF <- fit_flux(opt$par[1], opt$par[2], ES = PETE$ES, DS = PETE$DS) 

oFF$DS %>% pivot_longer(c(tn_flux, tot_flux2)) %>%
  ggplot(aes(x = value, y = depth, col = name)) + geom_point() + scale_y_reverse()


# the tiger is out! Alpha = 2

### Data for figure

BianchiMids <- ((BianchiBins + lag(BianchiBins)) / 2)[2:length(BianchiBins)]

ExpandedB <- expand_grid(depth = BianchiMids, lb = lb_vvv) %>% left_join(lbbs, by = "lb")

PredictedB <- predict(gamESE, ExpandedB, type = "link")

EPB <- bind_cols(ExpandedB, link = PredictedB)

EPB <- EPB %>%
  mutate(n_nparticles = exp(link),
         nparticles = n_nparticles * binsize) %>%
  mutate(profile = "multiple", project = "ETNP", time = "togetawatch") %>%
  mutate(flux = nparticles * (opt$par[1] * lb ^ opt$par[2]))

EDB <- EPB %>% group_by(depth) %>% summarize(Flux = sum(flux))
  

write_csv(EPB, "dataOut/CombinedProfileFluxEst_ED.csv")
write_csv(EDB, "dataOut/CombinedProfileFluxEst_DS.csv")
write_csv(oFF$DS, "dataOut/ObservedVsExpectedFlux.csv")

# DP <- EP %>% group_by(depth) %>% summarize() %>% 
#   mutate(profile = "multiple", project = "ETNP", time = "togetawatch")

