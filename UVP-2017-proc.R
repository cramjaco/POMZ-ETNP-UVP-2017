## Bring in data

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

## Set up multicore operation

no_cores <- availableCores() - 1
#plan(multiprocess, workers = no_cores)

## Process unbinned data

pt0 = proc.time()
unbinned <- twin01 %>% 
  calc_particle_parameters %>%
  calc_small_and_big() %>%
  calc_psd %>%
  calc_psd_gam_multiprofile() %>%
  pred_tp_gam %>%
  calc_small_psd() %>%
  calc_big_psd() %>%
  #double_gam_smooth() %>%
  #tp_quantiles(niter = 1000) %>%
  pass
pt1 = proc.time()
pt1 - pt0


# Binning by Daniele's bins
## Start with twinned list.
## Danielle's bins from matlab are.
## unique([0:20:100 100:25:200 200:50:1000 1000:100:2000 2000:200:5000]);


binned01 <- bin_depths(twin01)



pt0 = proc.time()
binned02 <- binned01 %>% 
  calc_particle_parameters %>%
  calc_small_and_big() %>%
  calc_psd %>%
  calc_psd_gam_multiprofile() %>%
  pred_tp_gam %>%
  calc_small_psd() %>%
  calc_big_psd() %>%
  double_gam_smooth() %>%
  diagnose_disaggregation() %>%
  #tp_quantiles(niter = 1000) %>%
  pass
pt1 = proc.time()
pt1 - pt0


# Write out data

write_csv(binned02[["ES"]], "dataOut/binned_EachSize.csv")
write_csv(binned02[["DS"]], "dataOut/binned_DepthSummary.csv")
write_csv(unbinned[["ES"]], "dataOut/unbinned_EachSize.csv")
write_csv(unbinned[["DS"]], "dataOut/unbinned_DepthSummary.csv")

save.image("setupProc.RData")


