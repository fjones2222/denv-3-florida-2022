## Name: Isty Rysava
## Date: 31/07/2023
## Code: A script to estimates dates for unobserved cases withing each identified link

rm(list=ls())
"C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3"

source("code/R/dengue_SI.R")
library(tidyverse)

## Data
cases <- read_csv(file="output/cases_assigned_PrunePi2023-07-18_maxP.csv")

## Get most likely SI
SIs <- rgamma(1000, shape=gamma_shape, scale=gamma_scale)
summary(SIs)
si <- mean(SIs)

## Get the temporal distance
links <- filter(cases, infector>0)
ss <- (as.POSIXlt(cases$Onset.Date)$year-122)*365 + (as.POSIXlt(cases$Onset.Date)$yday) # get exact dates
unobserved <- c()

for(i in 1:nrow(links)){
  id <- links$ID[i]
  infector <- links$infector[i]
  if(ss[id]-ss[infector]>2*si){
    dates <- seq(as.Date(cases$Onset.Date[infector]), as.Date(cases$Onset.Date[id]), by = si)
    fin.dates <- dates[2:(length(dates)-1)]
    unobserved <- rbind(unobserved, data.frame(fin.dates, id, infector))
  }
}

unobserved <- data.frame(unobserved)
colnames(unobserved) <- c("dates_unobserved", "ID", "source")
write.csv(unobserved, "output/unobserved_dates.csv", row.names=FALSE)

