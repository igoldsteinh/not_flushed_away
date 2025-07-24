### Process Results det ei normal for newshedding 
library(tidyverse)
library(tidybayes)
library(posterior)
library(fs)
source("src/utility_functions.R")
key_frame <- read_csv(here::here("data", "newshedding_data", "sim_key.csv"))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  sim = 115
  seed = 1
} else {
  sim <- as.integer(args[1])
  seed <- as.integer(args[2])
  
}

datafilename = key_frame %>% 
    filter(sim_num == sim) %>% 
    pull(obsdata_filename) %>% 
    unique() %>%
    as.character()
in_order_val = FALSE
  

priors_only = FALSE
if(priors_only == TRUE) {
  
  priorname <- paste0("prior_generated_quantities_scenario", sim, "_seed", seed, ".csv")
  prior_gq_samples_all <- read_csv(here::here("results",
                                              "det_ei_normal",
                                              priorname)) %>%
    pivot_longer(-c(iteration, chain)) %>%
    select( name, value)
  
  
  priors <- make_fixed_posterior_samples(prior_gq_samples_all)
  
  
  prior_timevarying_quantiles <- make_timevarying_posterior_quantiles(prior_gq_samples_all)
  
  prior_samp_name <- paste0("prior_samples_scenario", sim , "_seed", seed, ".csv")
  prior_timevarying_name <- paste0("prior_timevaryingquantiles_scenario", sim, "_seed", seed, ".csv")
  write_csv(priors, here::here("results", "det_ei_normal", prior_samp_name))
  write_csv(prior_timevarying_quantiles, here::here("results", "det_ei_normal", prior_timevarying_name))
  
  # make prior predictive intervals 
  prior_pred_address <- paste0("results/det_ei_normal/prior_predictive_scenario",
                               sim,
                               "_seed",
                               seed,
                               ".csv")
  prior_pred <- read_csv(prior_pred_address)
  
    simdata <- read_csv(here::here(datafilename)) %>% 
      mutate(total_conc = 1)
    

  prior_pred_intervals <- make_post_pred_intervals(prior_pred, simdata)
  
  prior_pred_interval_address <- paste0("results/det_ei_normal/prior_predictive_intervals_scenario",
                                        sim,
                                        "_seed",
                                        seed,
                                        ".csv")
  
  write_csv(prior_pred_intervals, prior_pred_interval_address)
  
  
  quit()
}


# posterior ---------------------------------------------------------------


# calculate MCMC diagnostics after burnin
gq_address <- paste0("results/det_ei_normal/generated_quantities/generated_quantities_scenario", 
                     sim, 
                     "_seed", 
                     seed,
                     ".csv")

posterior_samples <- read_csv(gq_address) %>%
  rename(.iteration = iteration,
         .chain = chain) %>%
  as_draws()



subset_samples <- subset_draws(posterior_samples, chain = c(1,2,3,4))

mcmc_summary <- summarise_draws(subset_samples)

mcmc_summary_address <- paste0("results/det_ei_normal/mcmc_summaries/mcmc_summary_scenario", 
                               sim, 
                               "_seed",
                               seed,
                               ".csv")
write_csv(mcmc_summary, mcmc_summary_address)

# lp trace plot -----------------------------------------------------------
lp_df <-read_csv(here::here("results", "det_ei_normal", "generated_quantities", paste0("posterior_df_scenario",
                                                                                       sim, 
                                                                                       "_seed",
                                                                                       seed,
                                                                                       ".csv"))) 

trace_plot <- lp_df %>%
  ggplot(aes(x = iteration, y = lp, color = as.factor(chain))) + 
  geom_line() +
  theme_bw() + 
  ggtitle("Stoch Conc Ei Normal LP Trace")

# ggsave(here::here("results", "det_ei_normal", "mcmc_summaries", paste0("trace_scenario", sim, "_seed", seed, ".png" )), trace_plot, width = 5, height = 5)

# create long format fixed samples and time-varying quantiles -----------------
posterior_gq_samples_all <- subset_samples  %>%
  pivot_longer(-c(.iteration, .chain)) %>%
  dplyr::select( name, value)


posterior_fixed_samples <- make_fixed_posterior_samples(posterior_gq_samples_all)

fixed_samples_address <- paste0("results/det_ei_normal/generated_quantities/posterior_fixed_samples_scenario",
                                sim, 
                                "_seed",
                                seed,
                                ".csv")

write_csv(posterior_fixed_samples, fixed_samples_address)


posterior_timevarying_quantiles <- make_timevarying_posterior_quantiles(posterior_gq_samples_all)


timevarying_quantiles_address <- paste0("results/det_ei_normal/generated_quantities/posterior_timevarying_quantiles_scenario",
                                        sim,
                                        "_seed",
                                        seed,
                                        ".csv")

write_csv(posterior_timevarying_quantiles, timevarying_quantiles_address)

rm(posterior_timevarying_quantiles)

rm(posterior_gq_samples_all)


# create posterior predictive quantiles -----------------------------------
# preserve if needed, but comment out for now due to possiblity of whacky numerical errors (not our fault)
post_pred_address <- paste0("results/det_ei_normal/posterior_predictive/posterior_predictive_scenario",
                            sim,
                            "_seed",
                            seed,
                            ".csv")
post_pred <- read_csv(post_pred_address)

  simdata <- read_csv(here::here(datafilename)) %>%
    mutate(total_conc = 1)


post_pred_intervals <- make_post_pred_intervals(post_pred, simdata, in_order_val)

post_pred_interval_address <- paste0("results/det_ei_normal/posterior_predictive/posterior_predictive_intervals_scenario",
                                     sim,
                                     "_seed",
                                     seed,
                                     ".csv")

write_csv(post_pred_intervals, post_pred_interval_address)



