# make psd gam multi profile

gamTestData <- binned01 %>% 
  calc_particle_parameters %>%
  calc_small_and_big() %>%
  calc_psd()

gtOut <- gamTestData %>%
  calc_psd_gam_multiprofile()

gtOut1 <- gtOut$out

gtMod <- gtOut$psd_gam_mod
