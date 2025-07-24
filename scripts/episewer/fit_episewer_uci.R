# fit EpiSewer
library(EpiSewer)
library(data.table)
library(tidyverse)
library(patchwork)
# load data ---------------------------------------------------------------
real_data <- read_csv(here::here("data", "uci_ww_data.csv"))
measurements_frame <- real_data %>% 
  mutate(rep_id = new_time) %>%
  mutate(concentration = `N2 Calculated Cp/mL`) %>%
  filter(concentration > 0) %>%
  dplyr::select(date,place, rep_id, concentration) 

# sample settings ---------------------------------------------------------
ww_sampling <- model_sampling(
  outliers = outliers_none(),
  sample_effects = sample_effects_none()
)

# sewer settings ----------------------------------------------------------
ww_sewage <- model_sewage(
  flows = flows_assume(1),
  residence_dist = residence_dist_assume(residence_dist = c(1))
)

# shedding settings ----------------------------------------------------------
ww_shedding <- model_shedding(
  shedding_dist = shedding_dist_assume(get_discrete_gamma(gamma_shape = 2.183783, gamma_scale = 1.836103), 
                                       shedding_reference = "symptom_onset"),
  incubation_dist = incubation_dist_assume(get_discrete_gamma(gamma_shape = 1, gamma_scale = 4)),
  load_variation = load_variation_estimate(),
  
)
# latent infection settings ---------------------------------------------------
ww_infections <- model_infections(
  generation_dist = generation_dist_assume(get_discrete_gamma_shifted(gamma_mean = 9.7, 
                                                                      gamma_sd = sqrt(2*(9.7/2)^2))/sum(get_discrete_gamma_shifted(gamma_mean = 9.7, gamma_sd = sqrt(2*(9.7/2)^2)))
  ),  
  R = R_estimate_rw(intercept_prior_mu = 0.5,
                    intercept_prior_sigma = 0.1),
  seeding = seeding_estimate_rw(),
  infection_noise = infection_noise_estimate()
)

# measurement settings G1----------------------------------------------------
frame_G1 <- measurements_frame %>% 
  filter(place == "G1")
ww_measurements_G1 <- model_measurements(
  concentrations = concentrations_observe(measurements = frame_G1,
                                          replicate_col = "rep_id"),
  noise = noise_estimate_constant_var(),
  LOD = LOD_none()
)

frame_G2 <- measurements_frame %>% 
  filter(place == "G2")
ww_measurements_G2 <- model_measurements(
  concentrations = concentrations_observe(measurements = frame_G2,
                                          replicate_col = "rep_id"),
  noise = noise_estimate_constant_var(),
  LOD = LOD_none()
)

frame_E <- measurements_frame %>% 
  filter(place == "E")
ww_measurements_E <- model_measurements(
  concentrations = concentrations_observe(measurements = frame_E,
                                          replicate_col = "rep_id"),
  noise = noise_estimate_constant_var(),
  LOD = LOD_none()
)


# stan settings -----------------------------------------------------------
ww_fit_opts <- set_fit_opts(
  model = model_stan_opts(package = "EpiSewer"),
  sampler = sampler_stan_mcmc(
    iter_warmup = 500,
    iter_sampling = 500,
    chains = 4,
    parallel_chains = 4, # run all chains in parallel
    seed = 42
  )
)

ww_results_opts <- set_results_opts(
  fitted = TRUE,
  summary_intervals = c(0.5, 0.8, 0.95),
  samples_ndraws = 50
)

# fit model ---------------------------------------------------------------
options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
test_G1 <- EpiSewer(
  measurements = ww_measurements_G1,
  sampling = ww_sampling,
  sewage = ww_sewage,
  shedding = ww_shedding,
  infections = ww_infections,
  fit_opts = ww_fit_opts,
  results_opts = ww_results_opts
)
write_rds(test_G1, here::here("scripts", "episewer", "episewer_G1.rds"))
plot_R(test_G1)


# fit model ---------------------------------------------------------------
options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
test_G2 <- EpiSewer(
  measurements = ww_measurements_G2,
  sampling = ww_sampling,
  sewage = ww_sewage,
  shedding = ww_shedding,
  infections = ww_infections,
  fit_opts = ww_fit_opts,
  results_opts = ww_results_opts
)
write_rds(test_G2, here::here("scripts", "episewer", "episewer_G2.rds"))
plot_R(test_G2)

# fit model ---------------------------------------------------------------
options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
test_E <- EpiSewer(
  measurements = ww_measurements_E,
  sampling = ww_sampling,
  sewage = ww_sewage,
  shedding = ww_shedding,
  infections = ww_infections,
  fit_opts = ww_fit_opts,
  results_opts = ww_results_opts
)
write_rds(test_E, here::here("scripts", "episewer", "episewer_E.rds"))
plot_R(test_E)