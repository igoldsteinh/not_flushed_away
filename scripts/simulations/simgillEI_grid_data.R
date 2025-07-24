# sim gillEI data grid
# idea is we're going to make many single data sets at different 
# initial pop sizes for stochastic Gillespie EI
# to try to get a better sense on when our model starts failing



# libraries ---------------------------------------------------------------

library(tidyverse)
library(glmnet)
library(patchwork)
library(GGally)
library(scales)
library(cowplot)
library(gridExtra)

source("src/simulate_stochastic_models.R")

# set the grid size variable and model params ----------------------------------------------

init_count = 10
num_epidemics = 100
I_init = init_count
E_init = init_count
gamma = 1/4
nu = 1/7
rho = 200
t_df = 2.99
t_sd = 0.5


## use same rt curve from previous scenario
sim = 1
seed_val = 4
full_simdata_address <- paste0("data/scenario", sim, "_pop1000_full_genecount_obsdata.csv")
full_simdata <- read_csv(here::here(full_simdata_address), show_col_types = FALSE)

fullseed_data <- full_simdata %>% filter(seed == seed_val) %>% 
  fill(Rt, .direction = "down") 

alpha_dat <- fullseed_data %>% 
  mutate(alpha = Rt * nu) %>% 
  dplyr::select(time, alpha)


alpha_init = alpha_dat$alpha[1]
alpha_vec = alpha_dat$alpha
change_points = alpha_dat$time
epidemics <- vector(mode='list', length=num_epidemics)
num = 1
i = 1
t_stop = 7*15
while (i <= num_epidemics) {
  set.seed(num)
  potential_epidemic = gillespie_ei_nonconst(E_init, 
                                             I_init, 
                                             alpha_init,
                                             alpha_vec,
                                             change_points,
                                             gamma, 
                                             nu, 
                                             t_stop)
  # we want at least 14 weeks worth of data
  if (max(potential_epidemic$t) > (7*14)) {
    epidemics[[i]] = potential_epidemic
    i = i + 1
  }
  
  num = num + 1
  
  
}
obs_start = 1
obs_end = obs_start + (14 * 7) - 1
seeds = 1:num_epidemics

daily_counts = map(epidemics, create_daily_data, 7*14) %>%
  map(~.x %>% filter(time >= obs_start & time <= obs_end))

obs_data = map(daily_counts, ~.x %>% 
                 mutate(log_genes_mean = log(I) + log(rho)) %>%
                 rowwise() %>%
                 mutate(
                   log_gene_copies1 = log_genes_mean  + (t_sd * rt(1,t_df)),
                   log_gene_copies2 = log_genes_mean  + (t_sd * rt(1,t_df)),
                   log_gene_copies3 = log_genes_mean  + (t_sd * rt(1,t_df)),
                   log_gene_copies4 = log_genes_mean  + (t_sd * rt(1,t_df)),
                   log_gene_copies5 = log_genes_mean  + t_sd * rt(1,t_df),
                   log_gene_copies6 = log_genes_mean  + t_sd * rt(1,t_df),
                   log_gene_copies7 = log_genes_mean  + t_sd * rt(1,t_df),
                   log_gene_copies8 = log_genes_mean  + t_sd * rt(1,t_df),
                   log_gene_copies9 = log_genes_mean  + t_sd * rt(1,t_df),
                   log_gene_copies10 = log_genes_mean + t_sd * rt(1,t_df)))
# take every other day 
final_obs_data = map(obs_data, ~.x %>%
                       mutate(new_time = time - obs_start + 1,
                              everyother = new_time %% 2 == 1)  %>% 
                       mutate(week = floor((new_time - 1)/7)) %>%
                       filter(everyother == TRUE))

state_data <- map(daily_counts, ~.x %>% left_join(alpha_dat, by = "time") %>%
                    mutate(Rt = alpha/nu)%>%
                    dplyr::select(integer_day, alpha, Rt)) 

full_obs_data <- map2(obs_data, state_data, ~.x %>%
                        left_join(.y, by = c("time" = "integer_day")) %>%
                        mutate(new_time = time - obs_start + 1) %>% 
                        mutate(week = floor((new_time - 1)/7)))



full_obs_data <- full_obs_data %>% 
  bind_rows(.id = "seed")

final_obs_data <- final_obs_data %>% 
  bind_rows(.id = "seed")

write_csv(full_obs_data, here::here("data", paste0("gillEIpop", init_count, "_full_genecount_obsdata.csv")))
write_csv(final_obs_data, here::here("data", paste0("gillEIpop", init_count, "_fitted_genecount_obsdata.csv")))


# make an id table --------------------------------------------------------
sim_name <- c("gillEIpop5", 
              "gillEIpop10", 
              "gillEIpop20",
              "gillEIpop30", 
              "gillEIpop40",
              "gillEIpop50",
              "gillEIpop60",
              "gillEIpop70",
              "gillEIpop80",
              "gillEIpop90",
              "gillEIpop100",
              "gillEIpop150",
              "gillEIpop200",
              "gillEIpop250",
              "gillEIpop300", 
              "gillEIpop5_truert",
              "gillEIpop10_truert",
              "gillEIpop20_truert",
              "gillEIpop30_truert",
              "gillEIpop40_truert",
              "gillEIpop50_truert",
              "gillEIpop60_truert",
              "gillEIpop70_truert",
              "gillEIpop80_truert",
              "gillEIpop90_truert",
              "gillEIpop100_truert")

scenario_name <- c("gillEIpop5", 
                   "gillEIpop10", 
                   "gillEIpop20",
                   "gillEIpop30", 
                   "gillEIpop40",
                   "gillEIpop50",
                   "gillEIpop60",
                   "gillEIpop70",
                   "gillEIpop80",
                   "gillEIpop90",
                   "gillEIpop100",
                   "gillEIpop150",
                   "gillEIpop200",
                   "gillEIpop250",
                   "gillEIpop300", 
                   "gillEIpop5",
                   "gillEIpop10",
                   "gillEIpop20",
                   "gillEIpop30",
                   "gillEIpop40",
                   "gillEIpop50",
                   "gillEIpop60",
                   "gillEIpop70",
                   "gillEIpop80",
                   "gillEIpop90",
                   "gillEIpop100")

sim_num <- c(105, 110, 120, 130, 140, 150, 160, 170, 180, 190, 1100,
             1150, 1200, 1250, 1300, 205, 210, 220, 230, 240, 250, 260,
             270, 280, 290, 2100)

sim_id <- data.frame(sim_num, sim_name, scenario_name)

write_csv(sim_id, here::here("data", "gillEI_sim_dict.csv"))
