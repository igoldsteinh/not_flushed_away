# fit off the shelf case models to case data
library(epidemia)
library(tidyverse)
library(sdprisk)
source("src/utility_functions.R")
# G data ---------------------------------------------------------------
G_data <- read_csv(here::here("data", "G_case_data.csv"))

# Run Epidemia -----------------------------------------------
data_length <- dim(G_data)[1]
dur_latent_days = 4
dur_infec_days = 7
rates <- c(7 / dur_latent_days, 7 / dur_infec_days) # pick the generation time

epidemia_weights <- epidemia_hypoexp(data_length, rates) %>% `/`(., sum(.))
delay_weights <- epidemia_gamma(data_length, 1, 7 / dur_latent_days) %>% `/`(., sum(.)) # pick the delay distribution

date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length), by = "days")
policy = c(NA, rep(0,7), rep(1,9))
epidemia_data <- data.frame(
  city = "Irvine",
  cases = c(NA, G_data$cases),
  date = date,
  day = weekdays(date),
  policy = policy
)

rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.15),
  prior_intercept = normal(log(0.5), 0.1),
  link = "log"
)

obs <- epiobs(
  formula = cases ~ 1 + policy,
  prior_intercept = rstanarm::normal(location = qlogis(0.25), scale = 0.2),
  prior = rstanarm::normal(location = 0, scale = 0.5),
  link = "logit",
  i2o = delay_weights[1:data_length]
)
args <- list(
  rt = rt,
  inf = epiinf(gen = epidemia_weights[1:data_length]),
  obs = obs,
  data = epidemia_data,
  iter = 6000,
  thin = 6,
  seed = 225
)

args$inf <- epiinf(
  gen = epidemia_weights[1:data_length],
  latent = TRUE,
  prior_aux = normal(10, 2)
)
estimnormal_posterior <- do.call(epim, args)

start_date <- min(G_data$time)
max_date <- max(G_data$time)

estimnormal_posterior_rt <-
  posterior_rt(estimnormal_posterior)[["draws"]] %>%
  data.frame() %>%
  `colnames<-`(start_date:(max_date + 1)) %>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws,
               names_to = "epidemia_time",
               values_to = "value",
               names_transform = list(epidemia_time = as.integer)
  ) %>%
  mutate(variable = "rt") %>%
  dplyr::select(variable, epidemia_time, value) %>%
  group_by(variable, epidemia_time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  filter(epidemia_time != 1) %>%
  mutate(time = epidemia_time - 1) %>%
  mutate(
    method = "estim_normal",
    name = "Rt"
  ) %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)

write_csv(estimnormal_posterior_rt, here::here("results", "case_models", "epidemia_mc_cases_posterior_rtv2.csv"))

# Run Epidemia centered at rt = 1-----------------------------------------------
data_length <- dim(G_data)[1]
dur_latent_days = 4
dur_infec_days = 7
rates <- c(7 / dur_latent_days, 7 / dur_infec_days) # pick the generation time

epidemia_weights <- epidemia_hypoexp(data_length, rates) %>% `/`(., sum(.))
delay_weights <- epidemia_gamma(data_length, 1, 7 / dur_latent_days) %>% `/`(., sum(.)) # pick the delay distribution

date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length), by = "days")
policy = c(NA, rep(0,7), rep(1,9))
epidemia_data <- data.frame(
  city = "Irvine",
  cases = c(NA, G_data$cases),
  date = date,
  day = weekdays(date),
  policy = policy
)

rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.15),
  prior_intercept = normal(log(1), 0.1),
  link = "log"
)


obs <- epiobs(
  formula = cases ~ 1 + policy,
  prior_intercept = rstanarm::normal(location = qlogis(0.25), scale = 0.2),
  prior = rstanarm::normal(location = 0, scale = 0.5),
  link = "logit",
  i2o = delay_weights[1:data_length]
)
args <- list(
  rt = rt,
  inf = epiinf(gen = epidemia_weights[1:data_length]),
  obs = obs,
  data = epidemia_data,
  iter = 6000,
  thin = 4,
  seed = 225
)

args$inf <- epiinf(
  gen = epidemia_weights[1:data_length],
  latent = TRUE,
  prior_aux = normal(10, 2)
)
estimnormal_posterior <- do.call(epim, args)

start_date <- min(G_data$time)
max_date <- max(G_data$time)

estimnormal_posterior_rt <-
  posterior_rt(estimnormal_posterior)[["draws"]] %>%
  data.frame() %>%
  `colnames<-`(start_date:(max_date + 1)) %>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws,
               names_to = "epidemia_time",
               values_to = "value",
               names_transform = list(epidemia_time = as.integer)
  ) %>%
  mutate(variable = "rt") %>%
  dplyr::select(variable, epidemia_time, value) %>%
  group_by(variable, epidemia_time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  filter(epidemia_time != 1) %>%
  mutate(time = epidemia_time - 1) %>%
  mutate(
    method = "estim_normal",
    name = "Rt"
  ) %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)

write_csv(estimnormal_posterior_rt, here::here("results", "case_models", "epidemia_mc_cases_rt=1_posterior_rtv2.csv"))
# E data ---------------------------------------------------------------

E_data <- read_csv(here::here("data", "E_case_data.csv"))

# Run Epidemia -----------------------------------------------
data_length <- dim(E_data)[1]
dur_latent_days = 4
dur_infec_days = 7
rates <- c(7 / dur_latent_days, 7 / dur_infec_days) # pick the generation time

epidemia_weights <- epidemia_hypoexp(data_length, rates) %>% `/`(., sum(.))
delay_weights <- epidemia_gamma(data_length, 1, 7 / dur_latent_days) %>% `/`(., sum(.)) # pick the delay distribution

date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length), by = "days")
policy = c(NA, rep(0,7), rep(1,9))

epidemia_data <- data.frame(
  city = "Irvine",
  cases = c(NA, E_data$cases),
  date = date,
  day = weekdays(date), 
  policy = policy
)

rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.15),
  prior_intercept = normal(log(0.5), 0.1),
  link = "log"
)

obs <- epiobs(
  formula = cases ~ 1 + policy,
  prior_intercept = rstanarm::normal(location = qlogis(0.25), scale = 0.2),
  prior = rstanarm::normal(location = 0, scale = 0.5),
  link = "logit",
  i2o = delay_weights[1:data_length]
)
args <- list(
  rt = rt,
  inf = epiinf(gen = epidemia_weights[1:data_length]),
  obs = obs,
  data = epidemia_data,
  iter = 4000,
  thin = 4,
  seed = 225
)

args$inf <- epiinf(
  gen = epidemia_weights[1:data_length],
  latent = TRUE,
  prior_aux = normal(10, 2)
)
estimnormal_posterior <- do.call(epim, args)

start_date <- min(E_data$time)
max_date <- max(E_data$time)

E_estimnormal_posterior_rt <-
  posterior_rt(estimnormal_posterior)[["draws"]] %>%
  data.frame() %>%
  `colnames<-`(start_date:(max_date + 1)) %>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws,
               names_to = "epidemia_time",
               values_to = "value",
               names_transform = list(epidemia_time = as.integer)
  ) %>%
  mutate(variable = "rt") %>%
  dplyr::select(variable, epidemia_time, value) %>%
  group_by(variable, epidemia_time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  filter(epidemia_time != 1) %>%
  mutate(time = epidemia_time - 1) %>%
  mutate(
    method = "estim_normal",
    name = "Rt"
  ) %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)

write_csv(E_estimnormal_posterior_rt, here::here("results", 
                                                    "case_models", 
                                                    "epidemia_E_cases_posterior_rtv2.csv"))

# Run Epidemia rt centered at 1-----------------------------------------------
data_length <- dim(E_data)[1]
dur_latent_days = 4
dur_infec_days = 7
rates <- c(7 / dur_latent_days, 7 / dur_infec_days) # pick the generation time

epidemia_weights <- epidemia_hypoexp(data_length, rates) %>% `/`(., sum(.))
delay_weights <- epidemia_gamma(data_length, 1, 7 / dur_latent_days) %>% `/`(., sum(.)) # pick the delay distribution

date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length), by = "days")
policy = c(NA, rep(0,7), rep(1,9))

epidemia_data <- data.frame(
  city = "Irvine",
  cases = c(NA, E_data$cases),
  date = date,
  day = weekdays(date),
  policy = policy
)

rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.15),
  prior_intercept = normal(log(1), 0.1),
  link = "log"
)

obs <- epiobs(
  formula = cases ~ 1 + policy,
  prior_intercept = rstanarm::normal(location = qlogis(0.25), scale = 0.2),
  prior = rstanarm::normal(location = 0, scale = 0.5),
  link = "logit",
  i2o = delay_weights[1:data_length]
)
args <- list(
  rt = rt,
  inf = epiinf(gen = epidemia_weights[1:data_length]),
  obs = obs,
  data = epidemia_data,
  iter = 8000,
  thin = 8,
  seed = 124
)

args$inf <- epiinf(
  gen = epidemia_weights[1:data_length],
  latent = TRUE,
  prior_aux = normal(10, 2)
)
estimnormal_posterior <- do.call(epim, args)

start_date <- min(E_data$time)
max_date <- max(E_data$time)

E_estimnormal_posterior_rt <-
  posterior_rt(estimnormal_posterior)[["draws"]] %>%
  data.frame() %>%
  `colnames<-`(start_date:(max_date + 1)) %>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws,
               names_to = "epidemia_time",
               values_to = "value",
               names_transform = list(epidemia_time = as.integer)
  ) %>%
  mutate(variable = "rt") %>%
  dplyr::select(variable, epidemia_time, value) %>%
  group_by(variable, epidemia_time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  filter(epidemia_time != 1) %>%
  mutate(time = epidemia_time - 1) %>%
  mutate(
    method = "estim_normal",
    name = "Rt"
  ) %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)

write_csv(E_estimnormal_posterior_rt, here::here("results", 
                                                    "case_models", 
                                                    "epidemia_E_cases_rt=1_posterior_rtv2.csv"))

