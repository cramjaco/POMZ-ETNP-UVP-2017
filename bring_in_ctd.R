library(readxl)
library(lubridate)
#library(reshape2)
ctddata <- read.csv("data/skq201617s-combinedctd.csv")
castdata <- read.csv("data/skq201617s-castinfo.csv")

### Lets bring in CTD data here

# note that P1 is cast 18-30 and P2 is cast 31-43
# This tells us if a cast is from stations P1 or P2, which are the only ones that I care about right now.
pstationator <- function(cast){
  if(cast >=  18 & cast <= 30){
    "P1"
  } else if (cast >= 31){
    "P2"
  } else {
    NA
  }
}

pstation <- sapply(castdata$cast, pstationator)
cast2 <- data.frame(cast = castdata$cast, time = as.POSIXct(castdata$startTime), latitude = castdata$latitude, longitude = castdata$longitude, pstation = pstation)

ctd.m = reshape2::melt(subset(ctddata, select = -X), id = c('cast', 'depth'))

# Also bring in par data

claraData <- read_excel("data/SKQ201617S_CTD_Profile_CF.xlsx", sheet = "PAR")
claraPAR <- claraData %>% select(Date = `mon/dd/yyyy`, Hr = `hh:mm`, PAR, PARPct, depth = `Pressure [db]`) %>%
  mutate(Hr = hour(Hr) + minute(Hr)/60) %>%
  filter(!(Hr > 4 & Hr < 3)) # There was one night cast clara forgot to remove
