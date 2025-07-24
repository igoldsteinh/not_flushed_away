# fit EpiSewer
library(EpiSewer)
library(data.table)
library(tidyverse)
library(scales)

# load data ---------------------------------------------------------------
key_frame <- read_csv(here::here("data", "newshedding_data", "sim_key.csv"))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  sim = 102
  seed_val = 1
} else {
  sim <- as.integer(args[1])
  seed_val <- as.integer(args[2])
  
}

datafilename = key_frame %>% 
  filter(sim_num == sim) %>% 
  pull(obsdata_filename) %>% 
  unique() 
sim_data <- read_csv(here::here(datafilename)) %>%
  filter(seed == seed_val) %>%
  dplyr::select(time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
  pivot_longer(cols = -time, names_to = "sample", values_to = "concentration") %>%
  dplyr::select(time, concentration) %>%
  mutate(concentration = exp(concentration))
fake_dates <- data.frame(date = seq.Date(from = as.Date("2020-01-01"), to = as.Date("2020-12-31"), by = "day"),
                         time = seq(1, 366, by = 1)) %>%
  mutate(weekday = "Monday")
measurements_frame <- sim_data %>% 
  left_join(fake_dates, by = "time") %>% 
  dplyr::select(date, time, concentration, weekday) %>%
  rename(rep_id = time)

full_data_filename <- key_frame %>% 
  filter(sim_num == sim) %>% 
  pull(fulldata_filename) %>% 
  unique() %>%
  as.character()

full_simdata <- read_csv(full_data_filename) %>% rename("true_rt" = "rt")

# sewer assumptions -------------------------------------------------------
# ww_assumptions = sewer_assumptions(
#   generation_dist = get_discrete_gamma_shifted(gamma_mean = 11.01098, gamma_sd = 4.80051),
#   shedding_dist = get_discrete_gamma(gamma_shape = 6.440142, gamma_scale = 2.258266),
#   shedding_reference = "infection",
#   incubation_dist = NULL,
#   min_cases = 1,
#   load_per_case = NULL,
#   flows = flows_assume(1),
#   residence_dist = c(1),
#   # R = R_estimate_rw(intercept_prior_mu = 1.9,
#   #                   sd_prior_mu = 0.2,
#   #                   sd_changepoint_dist = 7,
#   #                   sd_changepoint_sd = 0),
#   seeding = seeding_estimate_rw(),
#   infection_noise = infection_noise_estimate(),
#   replicate_col = "rep_id"
# )

# measurement settings ----------------------------------------------------
ww_measurements <- model_measurements(
  concentrations = concentrations_observe(measurements = measurements_frame,
                                          distribution = "log-normal",
                                          replicate_col = "rep_id"),
  noise = noise_estimate_constant_var(),
  LOD = LOD_none()
)

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
  shedding_dist = shedding_dist_assume(get_discrete_gamma(gamma_shape = 6.440142, gamma_scale = 2.258266), 
                                       shedding_reference = "symptom_onset"),
  incubation_dist = incubation_dist_assume(get_discrete_gamma(gamma_shape = 1, gamma_scale = 4)),
  load_variation = load_variation_estimate(),
  
)

# latent infection settings ---------------------------------------------------
ww_infections <- model_infections(
  generation_dist = generation_dist_assume(get_discrete_gamma_shifted(gamma_mean = 11.01098, gamma_sd = 4.80051)),
  seeding = seeding_estimate_rw(),
  R = R_estimate_rw(intercept_prior_mu = 1.93999,
                    intercept_prior_sigma = 0.190114,
                    sd_prior_sigma = 0.2),
  infection_noise = infection_noise_estimate()
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
test_result <- EpiSewer(
  measurements = ww_measurements,
  sampling = ww_sampling,
  sewage = ww_sewage,
  shedding = ww_shedding,
  infections = ww_infections,
  fit_opts = ww_fit_opts,
  results_opts = ww_results_opts
)
true_rt <- full_simdata %>% 
  filter(seed == seed_val) %>% 
  dplyr::select(seed, new_time, true_rt)
rt_res <- test_result[["summary"]][["R"]] %>%
  left_join(fake_dates, by = "date") %>%
  mutate(seed = seed_val) %>%
  left_join(true_rt, by = c("seed" = "seed", "time" = "new_time"))
test_result
my_theme <- list(
  scale_fill_brewer(name = "CI Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw(),
  theme())

episewer_102_seed1_long <- rt_res %>% 
  group_by(time) %>%
  pivot_longer(cols = c("lower_0.95", "lower_0.8", "lower_0.5"), names_to = ".width", values_to = ".lower") %>% 
  pivot_longer(cols = c("upper_0.95", "upper_0.8", "upper_0.5"), names_to = ".width2", values_to = ".upper") %>%
  mutate(.width = ifelse(.width == "lower_0.5" & .width2 == "upper_0.5", 0.5,
                         ifelse(.width == "lower_0.8" & .width2 == "upper_0.8", 0.8, 
                                ifelse(.width == "lower_0.95" & .width2 == "upper_0.95", 0.95, NA)))) %>%
  filter(.width == 0.5 | .width == 0.8 | .width == 0.95) 
episewer_rt <- episewer_102_seed1_long %>% 
  filter(time >= 1) %>%
  ggplot(aes(x = time, y = median, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = time, y = true_rt), color = "orange") +
  theme_bw() +
  ggtitle("") + 
  my_theme + 
  xlab("Time") +
  ylab("Rt") +
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank()) + 
  ylim(c(0, 4)) + 
  ggtitle("Episewer")
# write_csv(rt_res, here::here("results", "episewer", paste0("episewer_sim_", sim, "_seed_", seed_val, ".csv")))
