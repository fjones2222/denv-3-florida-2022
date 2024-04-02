## Name: Isty Rysava
## Date: 05/06/2023
## Code: A script to run transmission tree algorithm

rm(list=ls())
"C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3"

source("code/R/functions.R")
source("code/R/Trees.R") 
library("scales")
library('epicontacts')
library('igraph')
library('tidyverse')

## Init
mu_year_per_site <- 6.98 * 10^-4 # 4 * 10^-4 # substitutions per nucleotide per year
n_sites <- 10170
mutrate = (mu_year_per_site * n_sites / 365) # mutation rate per day and sequence
mut_low <- (5.39 * 10^-4 * n_sites / 365) 
mut_high <- (8.62 * 10^-4 * n_sites / 365) 
mutrates <- c(mut_low, mutrate, mut_high)
pis <- c(0.05, 0.1, 0.15) # reporting probability

####----------------------------------------------------------------------------------------------------------
#### Calculate tree metrics
params <- expand.grid(mutrate=mutrates, pi=pis)
not_identified <- vector("list", nrow(params))
Res <- vector("list", nrow(params))

for (idx in 1:nrow(params)){
  mutrate <- params[idx,1]
  pi <- params[idx,2]
  
  # Read in data
  trees <- read.csv(paste0("output/sensitivity/sources_pi", pi, "mut", round(mutrate, digits=3), ".csv"))
  LL <- read.csv(paste0("output/sensitivity/ll_pi", pi, "mut", round(mutrate, digits=3), ".csv"))
  #distprobs <- read.csv(paste0("output/sensitivity/distprobs_pi", pi, "mut", round(mutrate, digits=3), ".csv"))
  probs <-  read.csv(paste0("output/sensitivity/probs_pi", pi, "mut", round(mutrate, digits=3), ".csv"))
  
  cases = read.csv("data/cases_sorted.csv", header=T)
  cases$ID = 1:nrow(cases)
  cases$source <- ifelse(cases$travel.status=="Locally acquired", 0, 1)
  
  # Examine where tree building algorithm failed:
  nas = which(is.na(apply(trees, 1, sum))) 
  NAS = rep(NA, length(nas))
  for(i in 1:length(nas)){
    NAS[i] <- length(which(!is.na(trees[nas[i],])))
  }
  not_identified[[idx]] <- setdiff(cases$ID[nas[which(NAS == 0)]], which(cases$source==1)) # mostly all travel acquired, except for these 5

  # Find most likely progenitors for each case (and their bootstrap support)
  all <- nrow(cases)
  runs <- ncol(trees)
  infector <- bootstrap <- progll <- prob_tree <- prob_source <- numeric(all)
  for(i in 1:all){ # For every case:
    if(length(which(is.na(trees[i,])))!=runs){   
      infectors = tabulate(t(trees[i,]))
      infector[i] <- which(infectors == max(infectors))[1]
      bootstrap[i] <- max(infectors)
      prob_tree[i] <- max(infectors)/runs
      prob_source[i] <- max(probs[i, which(trees[i,]==infector[i])])
    }
  }
  nc = hist(infector, breaks = -1:max(cases$ID, na.rm=TRUE), plot=F)$counts[-1]
  cases$Re <- nc[match(cases$ID, 1:max(cases$ID))]
  cases$infector <- infector
  cases$prob_tree <- prob_tree
  cases$prob_source <- prob_source
  
  # Find most likely tree
  ML=numeric(ncol(LL))
  for(i in 1:ncol(LL)){ML[i]=sum(as.numeric(LL[,i]), na.rm=TRUE)}
  MLtree=which(ML==max(ML))
  MLtree=MLtree[round(runif(1, 1, length(MLtree)))]
  cases$MLtree_infector <- trees[,MLtree]
  nc2 = hist(trees[,MLtree], breaks = -1:max(cases$ID, na.rm=TRUE), plot=F)$counts[-1]
  cases$Re_MLtree <- nc2[match(cases$ID, 1:max(cases$ID))]
  write.csv(cases, file=paste0("output/sensitivity/cases_assigned__pi", pi, "mut", round(mutrate, digits=3), ".csv"), row.names=FALSE)
  
  # Compute Re for each cluster
  epic <- make_epicontacts(
    linelist = cases %>%
      mutate(Onset.Date=as.Date(Onset.Date,format="%Y-%m-%d")),
    contacts = cases %>%
      filter(infector!=0) %>%
      dplyr::select(ID,infector),
    id = "ID",
    from = "infector",
    to = "ID",
    directed = TRUE
  )
  clusters <- as.igraph(epic) %>% clusters()
  Res[[idx]] <- compute_R(clusters$csize, pi)[[1]] 
  
  print(idx)
}

sensitivity_analist <- list(params=params, not_identified=not_identified, Res=Res)
saveRDS(sensitivity_analist, "output/sensitivity/sensitivity_analist.Rdata")

####----------------------------------------------------------------------------------------------------------
#### Plot Re estimates
pdf("Figs/Re_pi_estimates.pdf", width=6, height=4)
plot(1:3, c(sensitivity_analist$Re[[2]]["central"], sensitivity_analist$Re[[5]]["central"], sensitivity_analist$Re[[8]]["central"]), 
     ylim = c(0, 1), pch = 19, xlab = "Reporting probability", ylab = "Estimated Re", axes = FALSE)
axis(side = 1, at = 1:3,
     labels = c("5%", "10%", "15%"), cex = 0.75)
axis(side = 2)
for(i in 1:3){
  j <- c(2, 5, 8)[i]
  segments(i, sensitivity_analist$Re[[j]]["low"], i, sensitivity_analist$Re[[j]]["up"])
}
dev.off()


