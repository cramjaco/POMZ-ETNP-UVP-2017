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
