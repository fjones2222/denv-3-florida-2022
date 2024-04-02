## Prepare dengue serial interval (for graph-based analysis of DENV-3 clusters)
## Code: Simulates distributions of extrinsic and intrinsic incubation period of dengue, and mosquito mortality
## Date: 05/23/2023

rm(list=ls())

## Library
library(fitdistrplus)

## (1) Extrinsic Incubation period
# parameters taken from Chan & Johansson, 2012 (Table 2)
tau=4.9
b0=2.9
bt=-0.08
EIPfc <- function(T) exp(exp(b0 + bt*T) + 1/(2*tau))
temps <- seq(20,30, length=10)
# plot(temps, EIPfc(temps), ylim=c(0, 60))
# EIP mean
# sd 1/(2*tau)

EIP <- c()
for(i in 1:length(temps)){
  EIP.temp <- rlnorm(100, log(EIPfc(temps[i])), 1/(2*tau))
  EIP <- c(EIP, EIP.temp)
}

## (2) Intrinsic Incubation period
# parameters taken from Chan & Johansson, 2012 (Table 3)
beta0.g <- 1.78
v.g <- 16 # gamma shape
lambda.g <- v.g / exp(beta0.g) # gamma rate
IIP <- rgamma(1000, shape=v.g, rate=lambda.g)

## (3) Mosquito longevity
# parameters taken from Johansson et al., 2014
sd <- 2
ML <- c()
for(i in 1:length(temps)){
mju <- 0.3967 - 0.03912*temps[i] + 2.442e-03*temps[i]^2 - 7.479e-05*temps[i]^3 + 9.298e-07*temps[i]^4
l <- 1/(mju*temps[i])
# ML.temp <- rnbinom(100, mu=l, size=sd)
ML.temp <- rexp(100, rate=l)
ML <- c(ML, ML.temp)
}

################################################################################
### Combine and fit to gamma distribution to estimate rate and shape
SI.dat <- EIP + IIP + 0.5*ML # half-time of mosquito life
fit <- fitdist(SI.dat, distr = "gamma", method = "mle")
# summary(fit); plot(fit)
gamma_shape <- fit$estimate[1]
gamma_rate <- fit$estimate[2]
gamma_scale <- 1/gamma_rate
# test1=rgamma(1000, shape=gamma_shape, scale=gamma_scale)
# test2=rgamma(1000, shape=gamma_shape, rate=gamma_rate)
gamma_mean <- mean(rgamma(1000, shape=gamma_shape, scale=gamma_scale))
gamma_sd <- sd(rgamma(1000, shape=gamma_shape, scale=gamma_scale))

rm(fit, SI.dat, IIP, EIP, tau, b0, bt, EIPfc, temps, beta0.g, v.g, lambda.g, i, EIP.temp, sd, mju, l, ML.temp, ML)
