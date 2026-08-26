# fit off the shelf case models to case data
library(EpiNow2)
library(tidyverse)
library(sdprisk)
source("src/utility_functions.R")
# G data ------------------------------------------------------------------

G_data <- read_csv(here::here("data", "G_case_data.csv"))
data_length <- dim(G_data)[1]

mean_time = (4 + 7)/7
GI_var = 2*(mean_time/2)^2

date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length) -1, by = "days")

epinow2_data <-  data.frame(
  confirm = G_data$cases,
  date = date
)
# delay distribution should be the latent period distribution
gamma_delay <- EpiNow2::delay_opts(
  EpiNow2::Gamma(1, rate = 7/4, max = 42)
)
# gen time 
mean_time = (4 + 7)/7
GI_var = 2*(mean_time/2)^2

gentime <- EpiNow2::gt_opts(
  EpiNow2::Gamma(mean = mean_time, sd = sqrt(GI_var),max = 20))
  
ascertainment <- EpiNow2::obs_opts(
    scale = EpiNow2::Normal(mean = 0.25, sd = 0.2),
    week_effect = FALSE
  )
  
init_rt_prior <- EpiNow2::rt_opts(
  prior = EpiNow2::LogNormal(mean = 1, sd = 0.1)
)
ls_mn_sd = c(6 * 7, 7)
alpha_mn_sd = c(0.2, 0.1) 

gp_rt_prior <- EpiNow2::gp_opts(
  # ls = EpiNow2::LogNormal(mean = ls_mn_sd[1], sd = ls_mn_sd[2]),
  # alpha = EpiNow2::LogNormal(mean = alpha_mn_sd[1], sd = alpha_mn_sd[2])
)
samples = 12000
warmup = 3000
thin = 2
chains = 4
cores = 4
adapt_delta = 0.99
max_treedepth = 12

stan_opts <- EpiNow2::stan_opts(
  samples = samples, warmup = warmup, thin = thin,
  chains = chains, 
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
)
en2_args <- list(
  epinow2_data,
  generation_time = gentime,
  truncation = EpiNow2::trunc_opts(EpiNow2::Fixed(0)),
  delays = gamma_delay,
  obs = ascertainment,
  rt = rt_opts(prior = LogNormal(mean = 1, sd = 0.1), rw = 1),
  gp = NULL,
  stan = stan_opts,
  forecast = NULL,
  CrIs = c(0.5, 0.8, 0.95)
  
)

en2_fit <- do.call(epinow, en2_args)
summaries <- en2_fit[["estimates"]][["summarised"]]

date_crosswalk <- data.frame(date = date, real_date = G_data$date)
R_summaries <- summaries %>% filter(variable == "R") %>%
  left_join(date_crosswalk, by = "date")
write_csv(R_summaries, here::here("results", "epinow2", "epinow2_G_rt1_rw.csv"))


# E data ------------------------------------------------------------------

E_data <- read_csv(here::here("data", "E_case_data.csv")) %>%
  filter(date > as.Date("2022-02-01"))
data_length <- dim(E_data)[1]

mean_time = (4 + 7)/7
GI_var = 2*(mean_time/2)^2

date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length) -1, by = "days")

epinow2_data <-  data.frame(
  confirm = E_data$cases,
  date = date
)
# delay distribution should be the latent period distribution
gamma_delay <- EpiNow2::delay_opts(
  EpiNow2::Gamma(1, rate = 7/4, max = 42)
)
# gen time 
mean_time = (4 + 7)/7
GI_var = 2*(mean_time/2)^2

gentime <- EpiNow2::gt_opts(
  EpiNow2::Gamma(mean = mean_time, sd = sqrt(GI_var),max = 20))

ascertainment <- EpiNow2::obs_opts(
  scale = EpiNow2::Normal(mean = 0.25, sd = 0.2),
  week_effect = FALSE
)

init_rt_prior <- EpiNow2::rt_opts(
  prior = EpiNow2::LogNormal(mean = 1, sd = 0.1)
)
ls_mn_sd = c(6 * 7, 7)
alpha_mn_sd = c(0.2, 0.1) 

gp_rt_prior <- EpiNow2::gp_opts(
  # ls = EpiNow2::LogNormal(mean = ls_mn_sd[1], sd = ls_mn_sd[2]),
  # alpha = EpiNow2::LogNormal(mean = alpha_mn_sd[1], sd = alpha_mn_sd[2])
)
samples = 12000
warmup = 3000
thin = 2
chains = 4
cores = 4
adapt_delta = 0.99
max_treedepth = 12

stan_opts <- EpiNow2::stan_opts(
  samples = samples, warmup = warmup, thin = thin,
  chains = chains, 
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
)
en2_args <- list(
  epinow2_data,
  generation_time = gentime,
  truncation = EpiNow2::trunc_opts(EpiNow2::Fixed(0)),
  delays = gamma_delay,
  obs = ascertainment,
  rt = rt_opts(prior = LogNormal(mean = 1, sd = 0.1), rw = 1),
  gp = NULL,
  stan = stan_opts,
  forecast = NULL,
  CrIs = c(0.5, 0.8, 0.95)
  
)

en2_fit <- do.call(epinow, en2_args)
summaries <- en2_fit[["estimates"]][["summarised"]]

date_crosswalk <- data.frame(date = date, real_date = E_data$date)
R_summaries <- summaries %>% filter(variable == "R") %>%
  left_join(date_crosswalk, by = "date")
write_csv(R_summaries, here::here("results", "epinow2", "epinow2_E_rt1_rw.csv"))




