### Process Results stoch conc ei
library(tidyverse)
library(tidybayes)
library(posterior)
library(fs)
source("src/utility_functions.R")

args <- commandArgs(trailingOnly=TRUE)


if (length(args) == 0) {
  init_count = 100.0
  seed = 4
} else {
  sim <- as.integer(args[1])
  seed <- as.integer(args[2])
  
}


priors_only = sim == 0




# priors only -------------------------------------------------------------
  priorname <- paste0("stc_params_constrt_pop", init_count, ".0_seed", seed, ".csv")
  prior_gq_samples_all <- read_csv(here::here("data",
                                              priorname)) %>%
    pivot_longer(-c(iteration, chain)) %>%
    select( name, value)
  
  
  priors <- make_fixed_posterior_samples(prior_gq_samples_all)
  
  
  prior_timevarying_quantiles <- make_timevarying_posterior_quantiles(prior_gq_samples_all)
  
  prior_timevarying_name <- paste0("stc_timevarying_quantiles_constrt_pop", init_count, ".0_seed", seed, ".csv")
  write_csv(prior_timevarying_quantiles, here::here("data", prior_timevarying_name))
