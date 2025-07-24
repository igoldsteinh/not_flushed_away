# summarise results stoch_conc_ei_normal
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
if (length(args) == 0) {
  sim = 115
  seed = 1
} else {
  sim <- as.integer(args[1])
  seed <- as.integer(args[2])
  
}

# read in data and results ---------------------------------------------------------

timevarying_suffix <- paste0("episewer_sim_",
                             sim)
timevarying_list <- list.files(path = path("results", "episewer"),  pattern = timevarying_suffix)


# diagnostics -------------------------------------------------------------

timevarying_suffix_diag <- paste0("episewer_diagnostics_sim_",
                             sim)
timevarying_list_diag <- list.files(path = path("results", "episewer"),  pattern = timevarying_suffix_diag)


# create final rt frame ---------------------------------------------------

rt_quantiles <- map(timevarying_list, ~read_csv(here::here("results", "episewer",.x)) %>% 
                                 mutate(address = .x)) %>%
  bind_rows()
  

write_csv(rt_quantiles, here::here("results", "episewer", paste0("episewer", sim,  "_allseeds_rt_quantiles.csv")))


# create final diagnostics frame ---------------------------------------------------

diagnostics <- map(timevarying_list_diag, ~read_csv(here::here("results", "episewer",.x)) %>% 
                      mutate(address = .x)) %>%
  bind_rows()

problems <- diagnostics %>%
  group_by(seed, sim) %>%
  summarise(
    num_diverge = sum(num_divergent),
    num_treedepth = sum(max_treedepth),
    min_bfmi = min(bfmi)
  ) %>%
  filter(min_bfmi < 0.3 | num_diverge > 0 | num_treedepth > 0)

write_csv(diagnostics, here::here("results", "episewer", paste0("episewer", sim,  "_allseeds_diagnostics.csv")))
