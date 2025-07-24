# simulate observed gene counts from our individual seiiiiiiirrr engine
# using normal emisssions distribution
# new shedding 
# NOTE
# the sanest way to simulate shedding is to start everyone in the E compartment
# this will mean having to filter so that we're starting when shedding is reasonable
library(tidyverse)
library(glmnet)
library(patchwork)
library(GGally)
library(scales)
library(cowplot)
library(gridExtra)

source("src/simulate_stochastic_seiiiiiiirrr_new.R")
# make a changing vector for r0
# pick a vector for R0 (right now, it's not changing, but it could!)
r03 = 2 
r0_init = r03[1]
r0_vec = r03
unique_r0 = r03
# we start at 0.9 at the start of r03, so that is the first time point
# choose how many times R0 changes
change_points = 0

# pick initial number of individuals in the compartments
# choose the populations size
N = 1000
# set the initial number of individuals in the I compartments
I1_init = 0
I2_init = 0
I3_init = 0
I4_init = 0
I5_init = 0
I6_init = 0
I7_init = 0
# set the initial number of individuals in the E comparmtnet
E_init = 20
R1_init = 0
R2_init = 0
# pick the average time spent infected but not infectious 
# 4 days
gamma = 1/4
# pick the average time spent infectious
# 7 days
nu = 7/7
total_nu = 1/7
# avg time spent recovered but shedding RNA
# 18 days total
eta = 2/18
total_eta = 1/18
#transform r0 vector into beta, pick the very first beta
beta_init = r0_init * total_nu
# create the vector of betas (should be same length as change points)
beta_vec = unique_r0 * total_nu 
# how many simulations?
num_sims = 100
#list to hold the simulated epidemics
epidemics <- vector(mode='list', length=num_sims)
num = 1
i = 1

t_stop = 15 * 7
while (i <= num_sims) {
  set.seed(num)
  # this is exhausting itself too quickly
  # first step is to double check r0 internal calculation
  # maybe add an Rt calculation
  # 
  potential_epidemic = sim_SEIIIIIIIRRR_nonconst_new(N, 
                                                 E_init, 
                                                 I1_init, 
                                                 I2_init, 
                                                 I3_init, 
                                                 I4_init,
                                                 I5_init,
                                                 I6_init,
                                                 I7_init,
                                                 R1_init,
                                                 R2_init,
                                                 beta_vec, 
                                                 change_points, 
                                                 gamma, 
                                                 nu1 = nu,
                                                 nu2 = nu,
                                                 nu3 = nu,
                                                 nu4 = nu,
                                                 nu5 = nu, 
                                                 nu6 = nu,
                                                 nu7 = nu,
                                                 eta1 = eta, 
                                                 eta2 = eta,
                                                 t_stop)
  epidemics[[i]] = potential_epidemic
  
    i = i + 1
  
  
  num = num + 1
  
  
}

#the simulation engine simulates individuals 
# so for 2000 people, it tracks all infection and recovery times for 2000 people
# we can calculate the population level counts from the individual simulation
#individual level data
individ_data = map(epidemics, pluck, 1)
#population level data 
state_data = map(epidemics, pluck, 2) %>% 
            map(~.x %>% create_daily_data_seiiiiiiirrr_new() %>% ungroup()) 

# write_rds(individ_data, here::here("data",  paste0("seirr_normalnewinfpop", init_count,"_individ_data.rds")))
# write_rds(state_data, here::here("data", paste0("seirr_normalnewinfpop", init_count, "_truecurve.rds")))

# plot the data -----------------------------------------------------------

plot_rt = state_data[[1]] %>% 
          filter(time <= 7 * 14) %>%
          ggplot(aes(x = time, y = rt)) +
          geom_point() + 
          theme_bw()

plot_states <- state_data[[1]] %>% 
               filter(time <= 7 * 14) %>%
               dplyr::select(time,  I1) %>% 
               pivot_longer(-time) %>%
               ggplot(aes(x = time, y = value, color = name)) + 
               geom_point()

plot_rt
plot_states
# creating true gene count data ------------------------------------------------------
# write_csv(individ_data, here::here("data", "sim_data", "scenario1_individ_data.csv"))
# individ_data = read_csv(here::here("data", "sim_data", "scenario1_individ_data.csv"))

# turn individual level data into something readable
# label is an individual
# each row is an individual and the times they changed to a new state
set.seed(1234)
wide_format = map(individ_data, ~.x %>% 
                    rowwise() %>%
                    mutate(individ_weight = rnorm(1, 0, 1.09)))

# if a value in wide_format is NA replace it with t_stop
wide_format = map(wide_format, ~.x %>% mutate_all(~replace_na(., t_stop)))
# create list of times at which we want to calculate the population level RNA concentration
times = map(state_data, ~seq(0, max(.x$integer_day), by = 1))


# create true RNA conc data ----------------------------------------------
# create 9 weights that sum to one
dumb_weights <- c(0.05, 0.4, 0.8, 0.4, 0.2, 0.1, 0.05, 0.025, 0.0125)
norm_weights <- dumb_weights / sum(dumb_weights)
comp_weights = log(norm_weights * 10e6, base = 10)
# create the true concentration of RNA before observation noise
true_conc_data <- vector(mode='list', length=num_sims)
for (i in (1:num_sims)) {
  #choose the times to calculate the population level RNA concentration before environmental/measurement noise
  current_times = times[[i]]
  
  true_conc_data[[i]] = map(current_times, ~new_shedding(time = .x, wide_format[[i]], comp_weights  ,N = N)) %>%
    bind_rows() %>%
    mutate(sim = i)
  
}


# plot the concentration
plot_conc = true_conc_data[[1]] %>% 
            ggplot(aes(x = time, y = total_conc)) +
            geom_point() + 
            theme_bw() + 
            xlab("Time") +
            ylab("Conc.")
# simulate data sets from the true values (Scenario 1)------------------------
obs_start = 1
obs_end = obs_start + (14 * 7) - 1
seeds = 1:num_sims
tau = 0.5
obs_true_data = map(true_conc_data, ~.x %>% filter(time >= obs_start & time <= obs_end))

# now simulate the actually observed RNA concentration data
obs_data = map2(seeds, obs_true_data, ~simulate_gene_data_normal(.y, 
                              seed = .x, 
                              rho = 0.0009, 
                              sd = tau))
# take every other day 
final_obs_data = map(obs_data, ~.x %>%
  mutate(new_time = time - obs_start + 1,
         everyother = new_time %% 2 == 1)  %>% 
    mutate(week = floor((new_time - 1)/7)) %>%
  filter(everyother == TRUE))


full_obs_data <- map2(obs_data, state_data, ~.x %>%
  left_join(.y, by = c("time" = "integer_day")) %>%
    mutate(new_time = time - obs_start + 1) %>% 
    mutate(week = floor((new_time - 1)/7)))

obs_true_data <- map2(obs_true_data, state_data, ~.x %>% 
  left_join(.y, by = c("time" = "integer_day")))


full_obs_data <- full_obs_data %>% 
                 bind_rows(.id = "seed") 


final_obs_data <- final_obs_data %>% 
                  bind_rows(.id = "seed")



# plot data ---------------------------------------------------------------
data_plot <- final_obs_data %>%
             filter(seed == 1) %>%
             dplyr::select(time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
             pivot_longer(-time) %>%
             ggplot(aes(x = time, y = value)) + 
             geom_point() + 
             theme_bw() + 
  xlab("Time") +
  ylab("Log Conc.")

# save files --------------------------------------------------------------
# write_rds(obs_true_data, here::here("data", paste0("seirr_normalnewinfpop", init_count, "_truegenecounts.rds")))
write_csv(full_obs_data, here::here("data", "newshedding_data", paste0("pop1000_fulldata.csv")))
write_csv(final_obs_data, here::here("data", "newshedding_data", paste0("pop1000_data.csv")))

sim_key <- read_csv(here::here("data", "newshedding_data", "sim_key.csv"))

key_frame <- data.frame(sim_num = 102,
                      sim_name = "newshedding_pop1000",
                      sim_filename = "simulate_sei7_newshedding_pop1k.R",
                      obsdata_filename = "data/newshedding_data/pop1000_data.csv",
                      fulldata_filename = "data/newshedding_data/pop1000_fulldata.csv",
                      N = N,
                      E_init = E_init,
                      gamma = gamma,
                      nu = nu,
                      eta = eta,
                      r0 = as.character(r0_vec),
                      change_points = as.character(change_points))


sim_key <- bind_rows(sim_key, key_frame)

write_csv(sim_key, here::here("data", "newshedding_data", "sim_key.csv"))
# create overdisp data ----------------------------------------------------
# overdispcase_data = simulate_case_data(obs_E2I_data[[1]],
#                                            seed = 1011,
#                                            rho = 0.5,
#                                            phi = 57.55)
# # aggregate to the week
# overdispfinal_case_data = overdispcase_data %>%
#                         mutate(new_time = time - obs_start + 1) %>%
#                         mutate(week = floor((new_time - 1)/7)) %>%
#                         group_by(week) %>%
#                         summarise(total_cases = sum(cases),
#                                   total_E2I = sum(E2I_transitions),
#                                   num_days = n()) %>%
#                         mutate(new_week = week + 1)
# 
# 
# 
# write_csv(overdispfinal_case_data, here::here("data", paste0("seirr_normalnewinfpop", init_count, "total_pop", N, "_overdisp.csv")))

# initial conditions ------------------------------------------------------
# state_data = read_rds(here::here("data", "seirr_normal_truecurve.rds"))
# 
# initial_states = map(state_data, ~.x %>% mutate(diff = (obs_start-1) - integer_day) %>% 
#                                          filter(diff > 0) %>%
#                                          filter(diff == min(diff))) %>% bind_rows(.id = "seed")
# 
# testing = state_data %>% bind_rows(.id = "seed") %>% filter(seed == 3)
# write_csv(initial_states, here::here("data", "seirr_normal_initstates.csv"))


# visualize simulations ---------------------------------------------------

# visualize results
# all credit to Damon Bayer for plot functions 
# my_theme <- list(
#   scale_fill_brewer(name = "Credible Interval Width",
#                     labels = ~percent(as.numeric(.))),
#   guides(fill = guide_legend(reverse = TRUE)),
#   theme_minimal_grid(),
#   theme(legend.position = "bottom"))
# 
# sim = 1
# make_rt_plot <- function(seed_val) {
#   full_obs_data %>%
#     filter(seed == seed_val) %>%
#     ggplot(aes(time, Rt)) +
#     geom_point(color = "coral1") + 
#     scale_y_continuous("Rt", label = comma) +
#     scale_x_continuous(name = "Time") +
#     ggtitle(str_c("EI Rt Scenario ", sim, " Seed ", seed_val)) +
#     my_theme
# }
# 
# ggsave2(filename = here::here("data",  paste0("stoch_conc_ei_truert_scenario", sim, ".pdf")),
#         plot = full_obs_data %>%
#           distinct(seed) %>%
#           arrange(seed) %>%
#           pull(seed) %>%
#           map(make_rt_plot) %>%
#           marrangeGrob(ncol = 1, nrow = 1),
#         width = 12,
#         height = 8)
# 


