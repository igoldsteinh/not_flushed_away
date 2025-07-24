# summarise results det_ei_normal
library(tidyverse)
library(tidybayes)
library(posterior)
library(fs)
library(gridExtra)
library(ggplot2)
library(scales)
library(cowplot)
source("src/utility_functions.R")
key_frame <- read_csv(here::here("data", "newshedding_data", "sim_key.csv"))
args <- commandArgs(trailingOnly=TRUE)
ln = FALSE
first_rt = 0.9
if (length(args) == 0) {
  sim = 105
  seed = 1
} else {
  sim <- as.integer(args[1])
  seed <- as.integer(args[2])
  
}



# read in data and results ---------------------------------------------------------

timevarying_suffix <- paste0("posterior_timevarying_quantiles_scenario",
                             sim)
timevarying_list <- list.files(path = path("results", "det_ei_normal", "generated_quantities"),  pattern = timevarying_suffix)

mcmc_suffix <- paste0("mcmc_summary_scenario", 
                      sim)
mcmc_list <- list.files(path = path("results", "det_ei_normal", "mcmc_summaries"),  pattern = mcmc_suffix)

full_stan_diag <- map(mcmc_list, ~read_csv(here::here("results", "det_ei_normal",  "mcmc_summaries", .x)) %>% 
                        mutate(address = .x)) %>%
  bind_rows() %>%
  mutate(seed = as.numeric(stringr::str_extract(address, stringr::regex("(\\d+)(?!.*\\d)"))),
         scenario = stringr::str_match(address, "scenario(\\w+)_seed")[,2]) %>%  
  filter(scenario == sim) %>%
  group_by(seed) %>% 
  filter(variable != "R2[1]" & variable != "C[1]") %>%
  summarise(min_rhat = min(rhat),
            max_rhat = max(rhat),
            min_ess_bulk = min(ess_bulk),
            max_ess_bulk = max(ess_bulk),
            min_ess_tail = min(ess_tail),
            max_ess_tail = max(ess_tail))

write_csv(full_stan_diag, here::here("results", "det_ei_normal", paste0("det_ei_normal_scenario", sim,  "_allseeds_stan_diag.csv")))
# create final rt frame ---------------------------------------------------
full_data_filename <- key_frame %>% 
  filter(sim_num == sim) %>% 
  pull(fulldata_filename) %>% 
  unique() %>%
  as.character()

full_simdata <- read_csv(full_data_filename) %>% rename("true_rt" = "rt")


# the data set is simulated for a certain period of time
# but then based on how we choose to space apart observations (every two days, every seven etc)
# there is a max observed time in the data set, we should not judge the model beyond the fitted data (for now)
# this time should be the same across models (it will be the same for the case models even though they're slightly different)
fitted_simdata_address <- key_frame %>%
  filter(sim_num == sim) %>%
  pull(obsdata_filename) %>%
  unique() %>%
  as.character()

fitted_simdata <- read_csv(fitted_simdata_address)
max_time <- fitted_simdata %>% group_by(seed) %>% summarise(max_time = max(new_time))

timevarying_quantiles <- map(timevarying_list, ~read_csv(here::here("results", "det_ei_normal", "generated_quantities", .x)) %>% 
                                 mutate(address = .x,
                                        seed = as.numeric(stringr::str_extract(address, stringr::regex("(\\d+)(?!.*\\d)"))),
                                        scenario = stringr::str_match(address, "scenario(\\w+)_seed")[,2]))
  
  rt_quantiles <- timevarying_quantiles %>%
    map(~.x %>% filter(name == "rt_t_values") %>%
          rename(week = time) %>% 
          dplyr::select(seed, week, scenario, value, .lower, .upper, .width,.point, .interval)) %>%
    bind_rows() %>% 
    filter(scenario == sim) %>%
    right_join(full_simdata, by = c("week", "seed")) %>%
    left_join(max_time, by = "seed") %>%
    filter(time <= max_time,
           week >= 0)
  
write_csv(rt_quantiles, here::here("results", "det_ei_normal", paste0("det_ei_normal_scenario", sim,  "_allseeds_rt_quantiles.csv")))
# write_csv(I_quantiles, here::here("results", "det_ei_normal", paste0("eirr_scenario", sim,  "_allseeds_prevI_quantiles.csv")))

# visualize results
# all credit to Damon Bayer for plot functions 
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

make_rt_plot <- function(seed_val) {
  rt_quantiles %>%
    filter(seed == seed_val) %>%
    ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    geom_point(aes(time, true_rt), color = "coral1") + 
    scale_y_continuous("Rt", label = comma) +
    scale_x_continuous(name = "Time") +
    ggtitle(str_c("EI Posterior Rt Scenario ", sim, " Seed ", seed_val)) +
    my_theme
}
ggsave2(filename = here::here("results", "det_ei_normal", paste0("det_ei_normal_rt_plots_scenario", sim, ".pdf")),
        plot = rt_quantiles %>%
          distinct(seed) %>%
          arrange(seed) %>%
          pull(seed) %>%
          map(make_rt_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)


