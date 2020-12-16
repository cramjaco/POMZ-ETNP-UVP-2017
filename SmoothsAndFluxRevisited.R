

# Setup

options(readr.default_locale=readr::locale(tz="Mexico/General"))
source("UVP_2017_library.R")
source("ModelStuff.R")
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

# allETNP <-  binnedSFR$ES %>% filter(project == "ETNP")
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

# for actual fitting, I should perhaps sum the etnp data and then fit a smooth on that
# and then fit to the trap flux?
# rationalle being that the flux averages over the week
# but for analyiss I shoudl show each week, I think.
