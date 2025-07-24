# sample prior posterior plot for paper 


# params ------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
library(tidybayes)
suppressPackageStartupMessages(library(posterior))
library(fs)
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
source(here::here("src", "utility_functions.R"))
key_frame <- read_csv(here::here("data", "newshedding_data", "sim_key.csv"))
sim = 102
seed_val = 1
prior_seed = 1
lump_val = 7
true_init_rt = 1.958
# True Frame

fixed_names <- c("gamma",
                 "nu",
                 "rho_conc",
                 "rt_init",
                 "sigma_rt",
                 "tau",
                 "E_init",
                 "I_init")

true_vals <- c(1/4,
               1/7,
               NA,
               true_init_rt,
               NA,
               0.5,
               20.0,
               0.0)

true_frame <- data.frame(name = fixed_names, true_value = true_vals)

# samples -----------------------------------------------------------------


stc_fixed_posterior_samples <- read_csv(here::here("results",
                                                   "stoch_ei_normal",
                                                   "generated_quantities",
                                                   paste0("posterior_fixed_samples_scenario", sim, "_seed", seed_val, ".csv")), show_col_types = FALSE) %>%
  left_join(true_frame, by = "name") %>%
  mutate(type = "posterior") %>%
  filter(name %in% fixed_names)

stc_prior_samp_name <- paste0("prior_samples_scenario", sim , "_seed", prior_seed, ".csv")
stc_fixed_prior_samp <- read_csv(here::here("results", 
                                            "stoch_ei_normal", 
                                            stc_prior_samp_name)) %>%
  mutate(type = "Prior") %>% 
  left_join(true_frame, by = "name")



# plot --------------------------------------------------------------------



priors_and_posteriors <- rbind(stc_fixed_posterior_samples %>% dplyr::select(name, value, type, true_value), stc_fixed_prior_samp)

levels_list <- c("E_init", "I_init", "gamma", "nu", "rt_init", "sigma_rt", "rho_conc", "tau")
labels_list <- c("E(0)", "I(0)", "Gamma", "Nu", "Initial Rt", "Sigma", "Rho", "Tau")

priors_and_posteriors$name <- factor(priors_and_posteriors$name, levels = levels_list, labels = labels_list)
priors_and_posteriors$type <- factor(priors_and_posteriors$type, levels = c("posterior", "Prior"), labels = c("Posterior", "Prior"))

param_plot <- priors_and_posteriors %>%
  filter(!(name == "Tau" & value > 10)) %>%
  rename(Type = type) %>%
  ggplot(aes(value, Type, fill = Type)) +
  stat_halfeye(normalize = "xy")  +
  geom_vline(aes(xintercept = true_value), linetype = "dotted", size = 1) + 
  facet_wrap(. ~ name, scales = "free_x") +
  theme_bw() +
  theme(text = element_text(size = 18),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ggtitle("Stc Fixed Parameter Prior and Posteriors (Fixed Scenario)") +
  ylab("Type") + 
  xlab("Value")

# ggsave(here::here("figures", "paper_example_fixedpp_plot.pdf"), param_plot, height = 12, width = 12)
priors_and_posteriors <- priors_and_posteriors %>% 
                         mutate(facet_type = type)

priors_and_posteriors$facet_type <- factor(priors_and_posteriors$facet_type, levels = c("Prior", "Posterior"))
stc_tau_plot <- priors_and_posteriors %>%
  rename(Type = type) %>%
  filter(name == "Tau") %>%
  filter(!(name == "Tau" & value > 10)) %>%
  ggplot(aes(value, Type, fill = Type)) +
  stat_halfeye(normalize = "xy")  +
  geom_vline(aes(xintercept = true_value), linetype = "dotted", size = 1) + 
  facet_wrap(. ~ facet_type, scales = "free_x") +
  theme_bw() +
  theme(text = element_text(size = 18),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ggtitle("Stc Tau Prior and Posterior (Fixed Scenario)") + 
  xlab("Value")

stc_tau_plot

# ggsave(here::here("figures", "paper_example_taupp_plot.pdf"), stc_tau_plot, height = 4, width = 8)