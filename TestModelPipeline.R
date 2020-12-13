source("UVP-2017-proc.R")
source("ModelStuff.R")

test_go <- function(){
testModPipe <<- binned02 %>% 
  diagnose_disaggregation()

tmpes <- testModPipe$ES
tmpds <- testModPipe$DS
}

test_go_in <- function(){
testModPipe_in <<- binned02 %>%
  filter_profile(profile = "stn_040") %>%
  diagnose_disaggregation_one_profile()


}



# to do -- nix station 041, from all of the analysis -- like very early in the pipeline (it only has one depth)
# pull the non error term from metaNest in diagnose_disaggregation
# unnest everything somehow.