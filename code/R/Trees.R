### ALGORITHMS AND PARAMETERS FOR TREE CONSTRUCTION ###

guess=function(case,
               possdates,
               possIDs,
               distmatrix,
               knownsource=rep(0, length(case)),
               exclude=TRUE){

# find progenitors for cases - each case corresponds to a row in the case table
# default:    no known source - can assign known progenitor where possible

# Arguments:  the case (row in the case table),
#             dates of the cases,
#             IDs of possible progenitors
#             genetic distance matrix between all cases (no NAs to ensure assignments....)
  
  # print(case)

  if (sum(unlist(knownsource[case]))==0){  # if no known source, look for progenitor
		sourcei = which(possdates < possdates[case]) # or any possible progenitor
		} else { sourcei = match(unlist(knownsource[case]), NA)}  # or known possible sources

  # Probabilities of potential sources
  SIs = possdates[case]-possdates[sourcei]
  # SIprobs = dgamma(SIs, shape=SIshape, rate=SIrate)  # Serial interval probabilities of possible sources
  SIprobs = ifelse(SIs <= 4,  
                   NA,  # accounts for the truncation at 4 days
                   dgamma(SIs, shape=SIshape, rate=SIrate))
  
  sourceIDs = possIDs[sourcei]		# Identities of possible sources

  if(!is.na(distmatrix[case,1])){  # only examine for case where distances have been calculated
    sourcedists = distmatrix[case, sourcei] # determine distances of possible progenitors
    if(length(sourcedists)==0){
      distprobs <- 0  # if no source set P to 0 
    }else{
      # distprobs <- dpois(round(as.numeric(sourcedists), digits=0), lambda=mutrate)
      distprobs <- dnbinom(round(as.numeric(sourcedists), digits=0), size=SIshape, prob=(SIscale*mutrate)/(SIscale*mutrate+1))
    }

    sourceprobs = SIprobs*distprobs
    
    # Assign progenitors:
    unitprobs = sourceprobs/sum(sourceprobs, na.rm=T)  # Scale probablities to one
    unitprobs[is.na(unitprobs)] = 0		# BUT GET RID OF NAs!
    cumunitprobs = cumsum(unitprobs)		# and make them cumulative
    ravr = runif(1,0,1)					# pick a random variate (RV) for assigning the progenitor
    sourceID = which(cumunitprobs>ravr)[1]  # which probability matches the RV
    SOURCE = sourceIDs[sourceID]			# Designate that as the source
    sourceprob = unitprobs[sourceID]		# keep track of the OVERALL SCALED probability
    SIprob = SIprobs[sourceID]  	# keep track of the SI probability
    distprob = distprobs[sourceID]    # keep track of the DISTANCE probability
    sourceLL = log(sourceprobs[sourceID]) # assign log likelihood of progenitor

		} else { SOURCE = sourceprob = sourceLL=NA} # assign NAs if distance not calculate-able
	c(SOURCE, sourceprob, sourceLL, SIprob, distprob)
}

## Alternative version adapted from Cori et al., 2018 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6312344/) 
## Includes reporting probability and pruning
guess2=function(case,
               possdates,
               possIDs,
               distmatrix,
               pi,
               prune,
               cuts,
               SImax,
               distmax,
               knownsource=rep(0, length(case)),
               exclude=TRUE){
  
  if (sum(unlist(knownsource[case]))==0){  # if no known source, look for progenitor
    sourcei = which(possdates < possdates[case]) # or any possible progenitor
  } else {sourcei = match(unlist(knownsource[case]), NA)}  # or known possible sources
  
  # Probabilities of potential sources
  SIs = possdates[case]-possdates[sourcei]
  SIs[which(SIs>cuts['temporal'] | SIs<=4)] <- NA # apply temporal cut-offs
  SIprobs = f_temporal(SIs, pi = pi)
  sourceIDs = possIDs[sourcei]		# Identities of possible sources
  
  if(!is.na(distmatrix[case,1])){  # only examine for case where distances have been calculated
    sourcedists = distmatrix[case, sourcei] # determine distances of possible progenitors
    if(length(sourcedists)==0){
      distprobs <- 0  # if no source set P to 0 
    }else{
      distprobs = f_genetic(round(as.numeric(sourcedists), digits=0), pi = pi)
      # get max P out of all possible distances
    }
    distprobs[which(sourcedists>cuts['genetic'])] <- NA # apply genetic cut-offs
    
    # Get probability of being compatible with a linkage
    sourceprobs <- unitprobs <- SIprobs/SImax * distprobs/distmax
    
    # Assign progenitors
    unitprobs[is.na(unitprobs)] = 0		# ged rid of NAs
    ravr = runif(1,0,1)					# pick a random variate (RV) for assigning the progenitor
    sourceID.temp = which(unitprobs>ravr) # which probability matches the RV
    if(length(sourceID.temp)==0){
      SOURCE = NA			# Designate that as the source
      sourceprob = NA		# keep track of the OVERALL SCALED probability
      SIprob = NA  	# keep track of the SI probability
      distprob = NA    # keep track of the DISTANCE probability
      sourceLL = NA # assign log likelihood of progenitor
    }else{
      sourceID = which(unitprobs==max(unitprobs[sourceID.temp]))[1]
      SOURCE = sourceIDs[sourceID]			# Designate that as the source
      sourceprob = unitprobs[sourceID]		# keep track of the OVERALL SCALED probability
      SIprob = SIprobs[sourceID]  	# keep track of the SI probability
      distprob = distprobs[sourceID]    # keep track of the DISTANCE probability
      sourceLL = log(sourceprobs[sourceID]) # assign log likelihood of progenitor
    }
  } else { SOURCE = sourceprob = sourceLL=NA} # assign NAs if distance not calculate-able
  c(SOURCE, sourceprob, sourceLL, SIprob, distprob)
}

## Test version combining both approaches above + cumulative probs
guess3=function(case,
                possdates,
                possIDs,
                distmatrix,
                pi,
                prune,
                cuts,
                SImax,
                distmax,
                knownsource=rep(0, length(case)),
                exclude=TRUE){
  
  if (sum(unlist(knownsource[case]))==0){  # if no known source, look for progenitor
    sourcei = which(possdates < possdates[case]) # or any possible progenitor
  } else {sourcei = match(unlist(knownsource[case]), NA)}  # or known possible sources
  
  # Probabilities of potential sources
  SIs = possdates[case]-possdates[sourcei]
  SIs[which(SIs>cuts['temporal'] | SIs<=4)] <- NA # apply temporal cut-offs
  SIprobs = f_temporal(SIs, pi = pi)
  sourceIDs = possIDs[sourcei]		# Identities of possible sources
  
  if(!is.na(distmatrix[case,1])){  # only examine for case where distances have been calculated
    sourcedists = distmatrix[case, sourcei] # determine distances of possible progenitors
    if(length(sourcedists)==0){
      distprobs <- 0  # if no source set P to 0 
    }else{
      distprobs = f_genetic(round(as.numeric(sourcedists), digits=0), pi = pi)
      # get max P out of all possible distances
    }
    distprobs[which(sourcedists>cuts['genetic'])] <- NA # apply genetic cut-offs
    
    # Get probability of being compatible with a linkage
    unitprobs <- SIprobs/SImax * distprobs/distmax
    
    # Assign progenitors
    rtotal <- sum(unitprobs, na.rm=T)
    nonNA <- sourceIDs[which(!is.na(unitprobs))]
    sourceprobs <- unitprobs[!is.na(unitprobs)]
    cr <- cumsum(sourceprobs)
    ravr <-  runif(1,0,1)*rtotal
    sourceID <- nonNA[which(cr>ravr)[1]]
    if(is.na(sourceID)==T){
      SOURCE = NA			# Designate that as the source
      sourceprob = NA		# keep track of the OVERALL SCALED probability
      SIprob = NA  	# keep track of the SI probability
      distprob = NA    # keep track of the DISTANCE probability
      sourceLL = NA # assign log likelihood of progenitor
    }else{
      SOURCE = sourceIDs[sourceID]			# Designate that as the source
      sourceprob = unitprobs[sourceID]		# keep track of the overall probability
      SIprob = SIprobs[sourceID]  	# keep track of the SI probability
      distprob = distprobs[sourceID]    # keep track of the DISTANCE probability
      sourceLL = log(unitprobs[sourceID]) # assign log likelihood of progenitor
    }
  } else { SOURCE = sourceprob = sourceLL=NA} # assign NAs if distance not calculate-able
  c(SOURCE, sourceprob, sourceLL, SIprob, distprob)
}

## Calculate Re
calcRe = function(start, end, interval, Re_vec, dates){
  #   Function to calculate Re timeseries given subset of data (say from within specific gridcells)
  #   Arguments: 
  #   start and end of timeseries and intervals in the timeseries
  #   list of individual Re
  #   dates of cases
  
  times = seq(start, end, interval)
  Re_mean <- ReUCI <- ReLCI <- Re90CI <- Re10CI <- ReN <- rep(NA, length=length(times)) # Set up vectors
  
  for(i in 1:length(times)){ # Loop thru timeseries to calculate metrics
    Res = Re_vec[which(dates >= (times[i]-interval) & dates < times[i])] # Select cases from correct times
    if(length(Res)>0){
      Re_mean[i] = mean(Res)
      ReUCI[i] = quantile(Res, 0.975)
      ReLCI[i] = quantile(Res, 0.025)
      Re90CI[i] = quantile(Res, 0.9)
      Re10CI[i] = quantile(Res, 0.1)
      ReN[i] = length(Res)
    }
  }
  Re = data.frame(mean = Re_mean, UCI = ReUCI, LCI = ReLCI, CI_90 = Re90CI, CI_10=Re10CI, N = ReN)
  Re
}

## Compute Re (from Cori et. al, 2018)
compute_R  <- function(cl_size, rho) {
  profile <- profile_likelihood(y_obs = cl_size, 
                                rho = rho, 0.01, 20)
  R_estimate <- theta_max_likelihood(profile$theta,
                                     profile$Likelihood, 
                                     0.95)
  R <- c(central = R_estimate$theta_max_likelihood, 
         low = R_estimate$lower_theta, 
         up = R_estimate$upper_theta)
  
  import <- import(y_obs = cl_size, 
                   rho = rho, 
                   profile, 1e3, 1e3, 0.95)
  unobs <- c(central = import$theta_max_likelihood,
             low = import$lower_theta,
             up = import$upper_theta)
  
  return(list(R, unobs))
}

