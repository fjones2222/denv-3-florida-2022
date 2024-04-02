## Name: Isty Rysava
## Date: 05/06/2023
## Code: A script to check distances for matched cases

rm(list=ls())
"C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3"

### Data
source("code/R/load_data.R")
# cases <- read.csv(file="output/cases_assigned_PrunePi2023-07-18_maxP.csv")
cases <- read.csv(file="output/cases_assigned_PrunePi2023-07-19_maxPall.csv")
distmatrix <-read.csv("data/gendist_mat.csv", header=T); dim(distmatrix)

### Calculate distances
ss =  (as.POSIXlt(cases$Onset.Date)$year-122)*365 + (as.POSIXlt(cases$Onset.Date)$yday) # get exact dates
cases$time_d <- cases$gen_d <- NA
for(i in 1:nrow(cases)){
  if(cases$infector[i]!=0){
    sourcei <- cases$infector[i]
    cases$time_d[i] <- ss[i] - ss[sourcei]
    cases$gen_d[i] <- distmatrix[i, sourcei]
  }
}
# write.csv(cases, file="output/cases_assigned_PrunePi2023-07-19_maxPall.csv", row.names = F)

### Plot all and case time and gen distances
pdf("figs/dist_hists.pdf")
par(mfrow=c(2, 1), mar=c(5, 5, 0.5, 5))
hist(as.vector(as.matrix(distmatrix)), breaks = seq(0, 50, 1), main = "", xlab= "Gen dist between all cases")
hist(cases$gen_d, breaks = seq(0, 50, 1), main = "", xlab= "Gen dist between linked cases")
hist(as.vector(D_all$dates), breaks = seq(0, 190, 1), main = "", xlab= "Temporal dist between all cases")
hist(cases$time_d, breaks = seq(0, 190, 1), main = "", xlab= "Temporal dist between linked cases")
dev.off()

