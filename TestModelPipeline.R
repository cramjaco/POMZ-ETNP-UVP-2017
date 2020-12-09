source("UVP-2017-proc.R")
source("ModelStuff.R")

test_go <- function(){
testModPipe <<- binned02 %>% 
  filter_profile(profile = "stn_043") %>%
  diagnose_disaggregation_one_profile()
}

tmpes <- testModPipe$ES
tmpds <- testModPipe$DS
