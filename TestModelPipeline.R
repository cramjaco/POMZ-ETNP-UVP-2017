source("UVP-2017-proc.R")
source("ModelStuff.R")

testModPipe <- binned02 %>% 
  filter_profile(profile = "stn_043")

tmpes <- testModPipe$ES
tmpds <- testModPipe$DS
