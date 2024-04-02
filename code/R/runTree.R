## Name: Isty Rysava
## Date: 05/06/2023
## Code: A script to run transmission tree algorithm

rm(list=ls())
"C:/Users/tui9/OneDrive - CDC/GitHub/florida_denv_3"

source("code/R/functions.R")
source("code/R/dengue_SI.R")
source("code/R/Trees.R") 
library("scales")

## Import data and prepare for analysis
## 1) dates
cases = read.csv("data/cases_sorted.csv", header=T); dim(cases)
cases$ID = 1:nrow(cases)
cases$source <- ifelse(cases$travel.status=="Locally acquired", 0, 1)
# cases$source <- c(1, rep(0, (nrow(cases)-1)))

## 2) gen matrix
distmatrix = read.csv("data/gendist_mat.csv", header=T); dim(distmatrix)

### Get parameters for epi distributions
SIshape = gamma_shape
SIrate = gamma_rate
SIscale = gamma_scale

mu_year_per_site <- 6.98 * 10^-4 # 4 * 10^-4 # substitutions per nucleotide per year
n_sites <- 10170
mutrate = (mu_year_per_site * n_sites / 365) # mutation rate per day and sequence

pi = 0.01 # decide on reporting probability
quantile = 0.95 # decide on pruning quantile

## 3) Generate distance functions
f_temporal <- fpaircase(type = "temporal", gamma_shape = gamma_shape,
                        gamma_scale = gamma_scale)
f_genetic <- fpaircase(type = "genetic", poisson_rate = mutrate,
                       gamma_shape = gamma_shape, 
                       gamma_scale = gamma_scale)

## 4) prep cut-offs
cuts <- c(temporal = get_quantiles(f=f_temporal, p=quantile, pi = pi), 
          genetic = get_quantiles(f=f_genetic, p=quantile, pi = pi))

####----------------------------------------------------------------------------------------------------------
### RUN ALGORITHMS TO GET BOOTSTRAP PROGENITORS
n=1000 # 1000
sources = probs = ll = distprobs = SIprobs = matrix(nrow=nrow(cases), ncol=n)

system.time(  
  
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
    print(i)
  }
)

write.csv(sources, file="output/sourcesPrunePi01_maxPall_test.csv", row.names=FALSE)
write.csv(probs, file="output/probsPrunePi01_maxPall_test.csv", row.names=FALSE)
write.csv(ll, file="output/llPrunePi01_maxPall_test.csv", row.names=FALSE)
write.csv(SIprobs, file="output/SIprobsPrunePi01_maxPall_test.csv", row.names=FALSE)
write.csv(distprobs, file="output/distprobsPrunePi01_maxPall_test.csv", row.names=FALSE)

####----------------------------------------------------------------------------------------------------------
### Calculate tree metrics:
trees = data.frame(read.csv("output/sourcesPrunePi_maxPall.csv"))
LL = read.csv("output/llPrunePi_maxPall.csv")
distprobs = read.csv("output/distprobsPrunePi_maxPall.csv")
probs = read.csv("output/probsPrunePi_maxPall.csv")

# Examine where tree building algorithm failed
nas = which(is.na(apply(trees, 1, sum))) 
NAS = rep(NA, length(nas))
for(i in 1:length(nas)){
  NAS[i] <- length(which(!is.na(trees[nas[i],])))
}
cases$ID[nas[which(NAS == 0)]] # The IDs that have no ancestor 
setdiff(cases$ID[nas[which(NAS == 0)]], which(cases$source==1)) # travel cases
setdiff(which(cases$source==0), cases$ID[nas[which(NAS == 0)]]) 

setdiff(cases$ID[nas[which(NAS == 0)]], which(cases$travel.status=="Locally acquired")) 
setdiff(cases$ID[nas[which(NAS == 0)]], which(cases$travel.status!="Locally acquired")) 

# Check trees
par(mfrow=c(2,1))
hist(as.numeric(trees[33,]), breaks = 0:max(cases$ID))$counts  # explore potential progenitors of each case
hist(as.numeric(LL[33,]), breaks = seq(-70, 0, 0.01))  # and the log likelihood of each case

# Find most likely progenitors for each case 
# Look at average Re based on 1000s of tree reconstructions
all <- nrow(cases)
runs <- ncol(trees)
infector <- bootstrap <- progll <- prob_source <- numeric(all)
for(i in 1:all){ # For every case:
  if(length(which(is.na(trees[i,])))!=runs){   # If not all results are NA, assign most likely progenitor
    infectors = tabulate(t(trees[i,]))
    infector[i] <- which(infectors == max(infectors))[1]
    bootstrap[i] <- max(infectors)
    prob_source[i] <- max(probs[i, which(trees[i,]==infector[i])])
    print(i)
  }
}
length(which(bootstrap!=0))/all   # % of assigned progenitors 
nc = hist(infector, breaks = -1:max(cases$ID, na.rm=TRUE), plot=F)$counts[-1]
cases$Re <- nc[match(cases$ID, 1:max(cases$ID))]
cases$infector <- infector
cases$prob_source <- prob_source

### Find most likely tree
ML=numeric(ncol(LL))
for(i in 1:ncol(LL)){ML[i]=sum(as.numeric(LL[,i]), na.rm=TRUE)}
range(ML)
(MLtree=which(ML==max(ML, na.rm=T))) 

# assign Re based on ML tree
cases$MLtree_infector <- trees[,MLtree]
nc2 = hist(trees[,MLtree], breaks = -1:max(cases$ID, na.rm=TRUE), plot=F)$counts[-1]
cases$Re_MLtree <- nc2[match(cases$ID, 1:max(cases$ID))]
write.csv(cases, file=paste("output/cases_assigned_PrunePi01_", Sys.Date(), "_maxPall.csv", sep=""), row.names=FALSE)
















