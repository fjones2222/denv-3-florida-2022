## Name: Isty Rysava
## Date: 15/06/2023
## Code: A script to run sensitivity analysis on the transmission tree algorithm

rm(list=ls())
#"C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3"

source("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/code/R/functions.R")
source("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/code/R/dengue_SI.R")
source("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/code/R/Trees.R") 
library("scales")

## Import data and prepare for analysis
## 1) dates
cases = read.csv("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/data/cases_sorted.csv", header=T); dim(cases)
cases$ID = 1:nrow(cases)
cases$source <- ifelse(cases$travel.status=="Locally acquired", 0, 1)

## 2) gen matrix
distmatrix = read.csv("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/data/gendist_mat.csv", header=T); dim(distmatrix)

### Get parametrs for epi distributions
SIshape = gamma_shape
SIrate = gamma_rate
SIscale = gamma_scale

mu_year_per_site <- 6.98 * 10^-4 # substitutions per nucleotide per year, 95% highest probability density range 5.39x10-4 – 8.62x10-4
n_sites <- 10170
mutrate = (mu_year_per_site * n_sites / 365) # mutation rate per day and sequence
mut_low <- (5.39 * 10^-4 * n_sites / 365) 
mut_high <- (8.62 * 10^-4 * n_sites / 365) 
mutrates <- c(mut_low, mutrate, mut_high)

pis <- c(0.05, 0.1, 0.15) # reporting probability
quantile = 0.95 # decide on pruning quantile

## 3) Generate distance functions
f_temporal <- fpaircase(type = "temporal", gamma_shape = gamma_shape,
                        gamma_scale = gamma_scale)
f_genetic <- fpaircase(type = "genetic", poisson_rate = mutrate,
                       gamma_shape = gamma_shape, 
                       gamma_scale = gamma_scale)

####----------------------------------------------------------------------------------------------------------
### RUN ALGORITHMS TO GET BOOTSTRAP PROGENITORS
params <- expand.grid(mutrate=mutrates, pi=pis)

for (idx in 1:nrow(params)){
  mutrate <- params[idx,1]
  pi <- params[idx,2]
  
  ## prep cut-offs
  cuts <- c(temporal = get_quantiles(f=f_temporal, p=quantile, pi = pi), 
            genetic = get_quantiles(f=f_genetic, p=quantile, pi = pi))
  
  ## initialize
  n=1000 
  sources = probs = ll = distprobs = SIprobs = matrix(nrow=nrow(cases), ncol=n)
  
  #INFORMED GUESSES - remove cows/livestock
  for (i in 1:n){
    ss =  (as.POSIXlt(cases$Onset.Date)$year-122)*365 + (as.POSIXlt(cases$Onset.Date)$yday) # get exact dates
    SImax = max(f_temporal(5:cuts['temporal'], pi = pi)) # get max P out of all possible SIs
    distmax = max(f_genetic(0:cuts['genetic'], pi = pi)) # get max P out of all possible distances
    
    guesses = unlist(lapply(1:nrow(cases), guess2, possdates=ss, possIDs=cases$ID, pi=pi, prune=quantile,
                            cuts = cuts, SImax = SImax, distmax = distmax,
                            knownsource=as.list(cases$source), distmatrix=distmatrix))
    
    sources[,i] = guesses[seq(1, length(guesses), 5)]
    probs[,i] = guesses[seq(2, length(guesses), 5)]
    ll[,i] = guesses [seq(3, length(guesses), 5)]
    SIprobs[,i] = guesses [seq(4, length(guesses), 5)]
    distprobs[,i] = guesses [seq(5, length(guesses), 5)]
  }
  
  write.csv(sources, file=paste0("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/output/sensitivity/sources_pi", pi, "mut", round(mutrate, digits=3), ".csv"), row.names=FALSE)
  write.csv(probs, file=paste0("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/output/sensitivity/probs_pi", pi, "mut", round(mutrate, digits=3), ".csv"), row.names=FALSE)
  write.csv(ll, file=paste0("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/output/sensitivity/ll_pi", pi, "mut", round(mutrate, digits=3), ".csv"), row.names=FALSE)
  write.csv(SIprobs, file=paste0("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/output/sensitivity/SIprobs_pi", pi, "mut", round(mutrate, digits=3), ".csv"), row.names=FALSE)
  write.csv(distprobs, file=paste0("C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3/output/sensitivity/distprobs_pi", pi, "mut", round(mutrate, digits=3), ".csv"), row.names=FALSE)
  
  print(idx)
}
  
