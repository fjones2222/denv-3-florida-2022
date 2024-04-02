pacman::p_load(
  tidyverse,
  readxl,
  lubridate,
  tidycensus,
  sf,
  cowplot,
  RColorBrewer,
  gtsummary
)

library(vimes)
library(branchr)


# give.n <- function(x){
#   return(c(y = mean(x), label = length(x)))
# }


give.n <- function(x){
  return(list(y = -3, label = glue::glue("n={length(x)}")))
}
