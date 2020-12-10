source("UVP-2017-proc.R")
source("ModelStuff.R")

test_go <- function(){
testModPipe <<- binned02 %>% 
  diagnose_disaggregation()
}

test_go_in <- function(){
testModPipe_in <<- binned02 %>%
  filter_profile(profile = "p16n_100") %>%
  diagnose_disaggregation_one_profile()
}

tmpes <- testModPipe$ES
tmpds <- testModPipe$DS
