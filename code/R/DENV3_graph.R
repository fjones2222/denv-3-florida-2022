## Graph-based model to estimate DENV-3 clusters
## Code: Adapted from Cori et al., 2018 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6312344/)
## Date: 05/22/2023

rm(list=ls())

source("code/R/functions.R")
source("code/R/dengue_SI.R")
source("code/R/load_data.R")

## Combine data
plot(D_all, nclass = 60) # play around  with nclass
range(D_all$dates) # 0-182, with travel 0-222
range(D_all$dna) # 0-44.41978, with travel 0-45.14413

## Get model parameters
# serial interval distribution parameters
# mutation rate
mu_year_per_site <- 6.98 * 10^-4  # substitutions per nucleotide per year
n_sites <- genome_size
mu_day_whole <- (mu_year_per_site * n_sites / 365) # mutation rate per day and sequence

# reporting rate
pi <- 0.1 #  try 0.1, 0.15?

# quantiles
q <- c(.50, .75, .90, .95, .95^(1/3), .99, .999, .9995, .9999)
cutoff_choice <- q[4] ## 95% quantile; play  around

# remove lower bound for temporal data
# n <- as.numeric(summary(as.vector(D_all$dates))[2]/2)
inds=which(D_all$dates<5)
D_all$dates[inds] <- 0

## distance functions for each of the 2 types
f_temporal <- fpaircase(type = "temporal", gamma_shape = gamma_shape,
                        gamma_scale = gamma_scale)
par(mfrow=c(3,1))
plot(f_temporal, xlim = c(0,365), pi=1)
plot(f_temporal, xlim = c(0,365), pi=0.5)
plot(f_temporal, xlim = c(0,365), pi=0.1)

f_genetic <- fpaircase(type = "genetic", poisson_rate = mu_day_whole,
                       gamma_shape = gamma_shape, 
                       gamma_scale = gamma_scale)

## Plotting these
par(mfrow=c(2,1))
plot(f_temporal, xlim = c(0,365))
plot(f_genetic, xlim = c(0,5))

## Plot with cut offs
cols <- "red"
plot_overlay <- function(dist, f, q, pi, xlab, breaks, resol = 1,
                         q_color = cols, hist_bordercolor = "grey",
                         hist_color = "lightgrey"){
  
  qtl <- get_quantiles(f, q, pi = pi)
  hist(dist, col = hist_color, border = hist_bordercolor, 
       main = "", xlab = xlab, 
       breaks = breaks)
  par(new = TRUE)
  x <- seq(min(breaks), max(breaks), resol)
  y <- f(x, pi = pi)
  plot(x, y, type = "l", axes = FALSE, main = "", xlab = "", ylab = "")
  
  y <- f(x, pi = 1)
  lines(x, y, col="darkblue")
  
  ## add vertical lines corresponding to quantiles
  abline(v = qtl, col = q_color, lwd = 1)
  # abline(v = 5, col = q_color, lwd = 1)
}

### use the function above to create our plot: 
pdf("Figs/Dist_functions.pdf", width=9, height=7)
par(mfrow=c(2, 1), mar=c(5, 5, 0.5, 5))
## temporal
# plot_overlay(dist = as.vector(D_all$dates), 
#              f = f_temporal, 
#              q = q[4], 
#              pi = pi, 
#              xlab = "Pairwise distance in time (days)", 
#              breaks = seq(0,750, 10),
#              resol = 1)
qtl <- get_quantiles(f_temporal, q[4], pi = pi)
hist(as.vector(D_all$dates), col = "lightgray", border = "gray", 
     main = "", xlab = "Pairwise distance in time (days)", 
     breaks = seq(0,750, 10))
par(new = TRUE)
x <- seq(min(seq(0,750, 10)), max(seq(0,750, 10)), 1)
y1 <- f_temporal(x, pi = 1)
y01 <- f_temporal(x, pi = pi)
plot(x, y01, type = "l", axes = FALSE, main = "", xlab = "", ylab = "", ylim=c(0, max(y1)))
lines(x, y1, col="darkblue")
## add vertical lines corresponding to quantiles
abline(v = qtl, col = "red", lwd = 1)
abline(v = 5, col = "red", lwd = 1)

## genetic
# plot_overlay(dist = as.vector(D_all$dna), 
#              f = f_genetic, 
#              q = q[4], 
#              pi = pi, 
#              xlab = "Pairwise genetic distance", 
#              breaks = seq(0,50,1),
#              resol = 1)

qtl <- get_quantiles(f_genetic, q[4], pi = pi)
hist(as.vector(D_all$dna), col = "lightgray", border = "gray", 
     main = "", xlab = "Pairwise genetic distance", 
     breaks = seq(0,50, 1))
par(new = TRUE)
x <- seq(min(seq(0,50, 1)), max(seq(0,50, 1)), 1)
y1 <- f_genetic(x, pi = 1)
y01 <- f_genetic(x, pi = pi)
plot(x, y01, type = "l", axes = FALSE, main = "", xlab = "", ylab = "", ylim=c(0, max(y1)))
lines(x, y1, col="darkblue")
## add vertical lines corresponding to quantiles
abline(v = qtl, col = "red", lwd = 1)
dev.off()

#########################################################################################################################
## function to get results for a certain cutoff and reporting rate
get_res <- function(D_all, q, pi, f_temporal, f_genetic,
                    type = c("all", "temporal", "genetic")) {
  
  type <- match.arg(type)
  
  ## get the cutoffs
  cuts <- c(temporal = get_quantiles(f_temporal, q, pi = pi), 
            genetic = get_quantiles(f_genetic, q, pi = pi))
  
  if (type == "all") {
    ## use vimes
    out <- vimes(D_all, cutoff = cuts,
                 graph.opt = vimes.graph.opt(col.pal = funky))
  } else if (type == "temporal") {
    out <- vimes(vimes_data(dates = D_all$dates), cutoff = cuts["temporal"],
                 graph.opt = vimes.graph.opt(col.pal = funky))
  } else if (type == "genetic") {
    out <- vimes(vimes_data(dna = D_all$dna), cutoff = cuts["genetic"],
                 graph.opt = vimes.graph.opt(col.pal = funky))
  }
  
  return(out)
}

## prep  combinations of p (q) and pi 
q <- c(0.95)
pi <- c(0.05, 0.1, 0.15)
combi <- expand.grid(p = q,
                     pi = pi)
quantile_pretty <- pi*100
quantile_pretty <- paste0(quantile_pretty, "%")

#res <- vector(9L, mode = "list")
res <- vector(3L, mode = "list")

## run analysis
for (i in 1:nrow(combi)) {
  res[[i]] <- get_res(D_all, q=combi[i, 1],
                      pi=combi[i, 2], f_temporal,
                      f_genetic, type="all")
}
  
## visualise the output
par(mfrow = c(1, 3), mar=c(1,1,3,1))
for (i in 1:length(res)) {
  plot(res[[i]]$graph, vertex.label = "",
       main = paste("surveillance:", quantile_pretty[i]))
}

#########################################################################################################################
## Estimate the underlying reproduction number and number of imported cases
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

i <- which(combi$p %in% 0.95) # or different p
clust_size <- lapply(res, function(i) i$clusters$size) 
rho <- combi[,2]

R_estimate_and_imports <- lapply(1:length(clust_size),
                                 function(i) compute_R(clust_size[[i]], rho[i]))

R_estimates <- sapply(1:length(R_estimate_and_imports), 
                      function(i) R_estimate_and_imports[[i]][[1]])

N_unobs_estimates <- sapply(1:length(R_estimate_and_imports), 
                            function(i) R_estimate_and_imports[[i]][[2]])

N_tot_estimates <- N_unobs_estimates + 
  matrix(rep(lengths(clust_size), 3), 
         nrow = 3, byrow = TRUE)

## plot R estimates
pdf("Figs/Re_pi_estimates.pdf", width=6, height=4)
plot(R_estimates["central",], ylim = c(0, max(c(1.2, max(R_estimates)))),
     pch = 19, xlab = "Reporting probability", ylab = "Estimated Re", axes = FALSE)
axis(side = 1, at = 1:length(R_estimate_and_imports),
     labels = quantile_pretty, cex = 0.75)
axis(side = 2, cex.axis = 0.8)
for (i in 1:length(R_estimate_and_imports)) {
  segments(i, R_estimates["low",i],i, R_estimates["up",i])
}
abline(h = 1, col = "red", lty = 2) # adding R=1
dev.off()
# Re values"
# central 0.99 0.97 0.92
# low     0.92 0.85 0.78
# up      1.08 1.08 1.06

### transmission chain
library(epicontacts)
library(igraph)

contact <- res[[4]]$graph %>% as_data_frame() %>%
            mutate(infector=vimes_analysis$seq_id[as.integer(from)],
                   case_id = vimes_analysis$seq_id[as.integer(to)]) 
            # %>% select(from, to)
contact <- contact[,-c(1:2)]

epic <- make_epicontacts(
  linelist = vimes_analysis %>% mutate(case_id=seq_id),
  contacts = contact,
  id = "case_id",
  from = "infector",
  to = "case_id",
  directed = TRUE
 )


plot(
  epic,
  x_axis = "Onset Date",
  arrow_size = 0.5,
  node_size = 13,
  label = FALSE,
  height = 700,
  width = 700
)

