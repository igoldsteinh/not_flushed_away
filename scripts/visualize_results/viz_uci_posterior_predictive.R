# visualize uci posterior predictives

# data --------------------------------------------------------------------

#| output: false
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

my_theme <- list(
  scale_fill_brewer(name = "Quantile",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme())

lump_val = 7
seed_val = 1
prior_seed = 1
dorms <- c("G1", "G2", "E")

uci_data <- read_csv(here::here("data", "uci_ww_data.csv"), 
                     show_col_types = FALSE)

count <- uci_data %>% 
  group_by(place) %>%
  summarise(num_days = n_distinct(date))


# create time crosswalks for visualizing compartment counts
G1_data <- uci_data %>% filter(place == "G1")
G1_param_change_times <- seq(lump_val, max(G1_data$new_time), by = lump_val)
G1_solve_times <- sort(union(unique(G1_data$new_time), G1_param_change_times))
G1_num_timepoints <- length(G1_solve_times)
G1_min_date <- min(G1_data$date) - ddays(1)
G1_time_crosswalk <- data.frame(time = 0:(G1_num_timepoints-1), 
                                 new_time = G1_solve_times) %>%
  mutate(date = G1_min_date + ddays(new_time),
         lump = floor((new_time)/lump_val))

G2_data <- uci_data %>% filter(place == "G2")

G2_param_change_times <- seq(lump_val, max(G2_data$new_time), by = lump_val)
G2_solve_times <- sort(union(unique(G2_data$new_time), G2_param_change_times))
G2_num_timepoints <- length(G2_solve_times)
G2_min_date <- min(G2_data$date) - ddays(1)
G2_time_crosswalk <- data.frame(time = 0:(G2_num_timepoints-1), 
                                 new_time = G2_solve_times) %>%
  mutate(date = G2_min_date + ddays(new_time),
         lump = floor((new_time)/lump_val))

E_data <- uci_data %>% filter(place == "E")

E_param_change_times <- seq(lump_val, max(E_data$new_time), by = lump_val)
E_solve_times <- sort(union(unique(E_data$new_time), E_param_change_times))
E_num_timepoints <- length(E_solve_times)
E_min_date <- min(E_data$date) - ddays(1)
E_time_crosswalk <- data.frame(time = 0:(E_num_timepoints-1), 
                                    new_time = E_solve_times) %>%
  mutate(date = E_min_date + ddays(new_time),
         lump = floor((new_time)/lump_val))




# G1 ---------------------------------------------------------------------

sim = "G1"
G1_real_data <- G1_data %>%
  dplyr::select(new_time, log_conc)

# stoch
G1_stc_post_pred_intervals <- read_csv(here::here("results",
                                                   "stoch_ei_normal",
                                                   "posterior_predictive",
                                                   paste0("posterior_predictive_intervals_scenario", sim, "_seed", seed_val, ".csv")), show_col_types = FALSE)


G1_stc_post_pred_intervals <- G1_stc_post_pred_intervals %>%
  left_join(G1_time_crosswalk, by = c("new_time")) %>%
  left_join(G1_real_data, by = "new_time")

G1_stc_posterior_predictive_plot <- G1_stc_post_pred_intervals %>%
  ggplot() +
  geom_ribbon(aes(x = date, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
  geom_line(aes(x = date, y = value)) +
  geom_point(mapping = aes(x = date, y = log_conc), color = "orange") +
  scale_fill_brewer(name = "Quantile") +
  theme_bw() + 
  ggtitle("Stc G1") +
  xlab("") +
  ylab("")

# det
G1_post_pred_intervals <- read_csv(here::here("results",
                                               "det_ei_normal",
                                               "posterior_predictive",
                                               paste0("posterior_predictive_intervals_scenario", sim, "_seed", seed_val, ".csv")), show_col_types = FALSE)


G1_post_pred_intervals <- G1_post_pred_intervals %>%
  left_join(G1_time_crosswalk, by = c("new_time"))


G1_post_pred_intervals <- G1_post_pred_intervals %>%
  left_join(G1_real_data, by = "new_time")
G1_posterior_predictive_plot <- G1_post_pred_intervals %>%
  ggplot() +
  geom_ribbon(aes(x = date, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
  geom_line(aes(x = date, y = value)) +
  geom_point(mapping = aes(x = date, y = log_conc), color = "orange") +
  scale_fill_brewer(name = "Quantile") +
  theme_bw() + 
  ggtitle("Det G1") +
  xlab("") +
  ylab("Log. Conc.")


G1_plot <- G1_posterior_predictive_plot + G1_stc_posterior_predictive_plot  + 
  plot_annotation(title = "G1 Posterior Predictive") & theme(text = element_text(size = 18))


# G2 ---------------------------------------------------------------------

sim = "G2"
G2_real_data <- G2_data %>%
  dplyr::select(new_time, log_conc)

# stoch
G2_stc_post_pred_intervals <- read_csv(here::here("results",
                                                   "stoch_ei_normal",
                                                   "posterior_predictive",
                                                   paste0("posterior_predictive_intervals_scenario", sim, "_seed", seed_val, ".csv")), show_col_types = FALSE)


G2_stc_post_pred_intervals <- G2_stc_post_pred_intervals %>%
  left_join(G2_time_crosswalk, by = c("new_time")) %>%
  left_join(G2_real_data, by = "new_time")

G2_stc_posterior_predictive_plot <- G2_stc_post_pred_intervals %>%
  ggplot() +
  geom_ribbon(aes(x = date, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
  geom_line(aes(x = date, y = value)) +
  geom_point(mapping = aes(x = date, y = log_conc), color = "orange") +
  scale_fill_brewer(name = "Quantile") +
  theme_bw() + 
  ggtitle("Stc G2") +
  ylab("") +
  xlab("Date")

# det
G2_post_pred_intervals <- read_csv(here::here("results",
                                               "det_ei_normal",
                                               "posterior_predictive",
                                               paste0("posterior_predictive_intervals_scenario", sim, "_seed", seed_val, ".csv")), show_col_types = FALSE)


G2_post_pred_intervals <- G2_post_pred_intervals %>%
  left_join(G2_time_crosswalk, by = c("new_time")) %>%
  left_join(G2_real_data, by = "new_time")

G2_posterior_predictive_plot <- G2_post_pred_intervals %>%
  ggplot() +
  geom_ribbon(aes(x = date, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
  geom_line(aes(x = date, y = value)) +
  geom_point(mapping = aes(x = date, y = log_conc), color = "orange") +
  scale_fill_brewer(name = "Quantile") +
  theme_bw() + 
  ggtitle("Det G2") +
  ylab("Log Conc.") +
  xlab("Date") 


G2_plot <- G2_posterior_predictive_plot + G2_stc_posterior_predictive_plot + 
  plot_annotation(title = "G2 Posterior Predictive") & theme(text = element_text(size = 18))


# E ------------------------------------------------------------------

sim = "E"
E_real_data <- E_data %>%
  dplyr::select(new_time, log_conc)

# stoch
E_stc_post_pred_intervals <- read_csv(here::here("results",
                                                      "stoch_ei_normal",
                                                      "posterior_predictive",
                                                      paste0("posterior_predictive_intervals_scenario", sim, "_seed", seed_val, ".csv")), show_col_types = FALSE)


E_stc_post_pred_intervals <- E_stc_post_pred_intervals %>%
  left_join(E_time_crosswalk, by = c("new_time")) %>%
  left_join(E_real_data, by = "new_time")

E_stc_posterior_predictive_plot <- E_stc_post_pred_intervals %>%
  ggplot() +
  geom_ribbon(aes(x = date, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
  geom_line(aes(x = date, y = value)) +
  geom_point(mapping = aes(x = date, y = log_conc), color = "orange") +
  scale_fill_brewer(name = "Quantile") +
  theme_bw() + 
  ggtitle("Stc E") +
  xlab("") +
  ylab("")

# det
E_post_pred_intervals <- read_csv(here::here("results",
                                                  "det_ei_normal",
                                                  "posterior_predictive",
                                                  paste0("posterior_predictive_intervals_scenario", sim, "_seed", seed_val, ".csv")), show_col_types = FALSE)


E_post_pred_intervals <- E_post_pred_intervals %>%
  left_join(E_time_crosswalk, by = c("new_time")) %>%
  left_join(E_real_data, by = "new_time")

E_posterior_predictive_plot <- E_post_pred_intervals %>%
  ggplot() +
  geom_ribbon(aes(x = date, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
  geom_line(aes(x = date, y = value)) +
  geom_point(mapping = aes(x = date, y = log_conc), color = "orange") +
  scale_fill_brewer(name = "Quantile") +
  theme_bw() + 
  ggtitle("Det E") +
  xlab("") +
  ylab("Log Conc.")

E_plot <- E_posterior_predictive_plot + E_stc_posterior_predictive_plot  +
                plot_annotation(title = "E Posterior Predictive") & theme(text = element_text(size = 18))

full_plot <- E_plot / G1_plot / G2_plot + plot_layout(guides = "collect") + plot_annotation(title = "UCI Posterior Predictive")

# ggsave(here::here("figures", "uci_posterior_predictive.pdf"), full_plot, height = 11, width = 11)



