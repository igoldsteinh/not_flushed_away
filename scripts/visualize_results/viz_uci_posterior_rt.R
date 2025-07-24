# visualize real data from case and wastewater models (separately)
suppressPackageStartupMessages(library(tidyverse))
library(tidybayes)
suppressPackageStartupMessages(library(posterior))
library(fs)
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

seed_val = 1
lump_val = 7
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw(),
  theme())


# visualize uci data ------------------------------------------------------
# ww data
uci_data <- read_csv(here::here("data", "uci_ww_data.csv"), 
                     show_col_types = FALSE)
dorms <- c("G1", "G2", "E")
dorm_labels <- c("E", "G1", "G2")
count <- uci_data %>% 
  group_by(place) %>%
  summarise(num_days = n_distinct(date))
uci_plot_data <- uci_data %>% 
  rename(`Sub-Comm` = place) %>%
  filter(`Sub-Comm` %in% dorms)

uci_plot_data$`Sub-Comm` <- factor(uci_plot_data$`Sub-Comm`, levels = c("E", "G1", "G2"), labels = dorm_labels)
plot_data <- uci_plot_data %>% 
  group_by(`Sub-Comm`, date) %>%
  mutate(mean_log = mean(log_conc)) %>%
  ggplot(aes(x = date, y = log_conc)) +
  geom_point(size = 2) + 
  geom_line(aes(x = date, y = mean_log), linewidth = 1) + 
  ylab("Log Conc") + 
  xlab("Date") + 
  ylim(c(0,10)) +
  theme_bw() + 
  facet_wrap(vars(`Sub-Comm`))

plot_data

# case data
G_cases = read_csv(here::here("data", "G_case_data.csv")) %>%
  mutate(Community = "G")

E_cases = read_csv(here::here("data", "E_case_data.csv")) %>%
  mutate(Community = "E")

weekly_cases = bind_rows(G_cases, E_cases)

comm_labels <- c( "E", "G")
weekly_cases$Community <- factor(weekly_cases$Community, levels = c("E", "G"), labels = comm_labels)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_uci_cases <- weekly_cases %>% 
  ggplot(aes(x = date, y = cases)) +
  geom_point(size = 2) + 
  geom_line(linewidth = 1) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) +
  xlab("Date") + 
  ylab("Cases") + 
  facet_wrap(vars(Community))
plot_uci_cases     


uci_data_plot <- plot_data / plot_uci_cases + 
  plot_annotation(title = "UCI Case and Wastewater Data Feb 2022 - May 2022") &
  theme(text = element_text(size = 18))

# ggsave(here::here("figures", "uci_data_plot.pdf"), uci_data_plot, height = 10, width = 10)


# rt at 0.5 ww models -----------------------------------------------------
# create time crosswalks for visualizing compartment counts
uci_data <- read_csv(here::here("data", "uci_ww_data.csv"), 
                     show_col_types = FALSE)

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
G1_stc_timevarying_quantiles <- read_csv(here::here("results", 
                                                     "stoch_ei_normal", 
                                                     "generated_quantities", 
                                                     paste0("posterior_timevarying_quantiles_scenario", 
                                                            "G1", "_seed", seed_val, ".csv")),
                                          show_col_types = FALSE) 


G1_stc_rt_quantiles <- G1_stc_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(G1_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

# visualize results
G1_stc_plot_posterior_rt <- G1_stc_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Stc-ww G1") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = "none")


# redo for det
G1_det_timevarying_quantiles <- read_csv(here::here("results", 
                                                     "det_ei_normal", 
                                                     "generated_quantities", 
                                                     paste0("posterior_timevarying_quantiles_scenario", 
                                                            "G1", "_seed", seed_val, ".csv")),
                                          show_col_types = FALSE) 


G1_det_rt_quantiles <- G1_det_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(G1_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

# visualize results
G1_det_plot_posterior_rt <- G1_det_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Det-ww G1") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("") + 
  theme(legend.position = "none")



# G2 ---------------------------------------------------------------------
G2_stc_timevarying_quantiles <- read_csv(here::here("results", 
                                                     "stoch_ei_normal", 
                                                     "generated_quantities", 
                                                     paste0("posterior_timevarying_quantiles_scenario", 
                                                            "G2", "_seed", seed_val, ".csv")),
                                          show_col_types = FALSE) 


G2_stc_rt_quantiles <- G2_stc_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(G2_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

# visualize results
G2_stc_plot_posterior_rt <- G2_stc_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Stc-ww G2") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("Date") + 
  theme(legend.position = "none")



# redo for det
G2_det_timevarying_quantiles <- read_csv(here::here("results", 
                                                     "det_ei_normal", 
                                                     "generated_quantities", 
                                                     paste0("posterior_timevarying_quantiles_scenario", 
                                                            "G2", "_seed", seed_val, ".csv")),
                                          show_col_types = FALSE) 


G2_det_rt_quantiles <- G2_det_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(G2_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

# visualize results
G2_det_plot_posterior_rt <- G2_det_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Det-ww G2") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("Date") + 
  theme(legend.position = "none")


# E ------------------------------------------------------------------

E_stc_timevarying_quantiles <- read_csv(here::here("results", 
                                                        "stoch_ei_normal", 
                                                        "generated_quantities", 
                                                        paste0("posterior_timevarying_quantiles_scenario", 
                                                               "E", "_seed", seed_val, ".csv")),
                                             show_col_types = FALSE) 


E_stc_rt_quantiles <- E_stc_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(E_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

extra_rows <- G2_stc_rt_quantiles %>% 
             filter(date == min(date))

extra_rows$value <- NA
extra_rows$.lower <- NA
extra_rows$.upper <- NA

E_stc_rt_quantiles <- bind_rows(extra_rows, E_stc_rt_quantiles)
# visualize results
E_stc_plot_posterior_rt <- E_stc_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Stc-ww E") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = "none")



# redo for det
E_det_timevarying_quantiles <- read_csv(here::here("results", 
                                                        "det_ei_normal", 
                                                        "generated_quantities", 
                                                        paste0("posterior_timevarying_quantiles_scenario", 
                                                               "E", "_seed", seed_val, ".csv")),
                                             show_col_types = FALSE) 


E_det_rt_quantiles <- E_det_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(E_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

E_det_rt_quantiles <- bind_rows(E_det_rt_quantiles, extra_rows)
# visualize results
E_det_plot_posterior_rt <- E_det_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Det-ww E") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("") + 
  theme(legend.position = "none")


E_stoch_vs_det <- E_det_plot_posterior_rt + E_stc_plot_posterior_rt + plot_layout(guides = "collect")


# epidemia case fits ------------------------------------------------------
# G  --------------------------------------------------------------
G_case_posteriort = read_csv(here::here("results", 
                                         "case_models", 
                                         "epidemia_G_cases_posterior_rtv2.csv"))

G_case_posteriort =  G_case_posteriort %>%
  mutate(lump = time-1) %>% 
  dplyr::select(lump, value, .lower, .upper, .width) %>%
  right_join(G1_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

G_case_rt_plot <- G_case_posteriort %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Epidemia-cases G") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = "none")


G_case_rt_plot_bottom <- G_case_posteriort %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Epidemia-cases G") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("Date") + 
  theme(legend.position = "none")


# E ------------------------------------------------------------------

E_case_posteriort = read_csv(here::here("results", 
                                         "case_models", 
                                         "epidemia_E_cases_posterior_rtv2.csv"))

E_case_posteriort =  E_case_posteriort %>%
  mutate(lump = time-1) %>% 
  dplyr::select(lump, value, .lower, .upper, .width) %>%
  right_join(E_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

E_case_posteriort = bind_rows(E_case_posteriort, extra_rows)
E_case_rt_plot <- E_case_posteriort %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Epidemia-cases E") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = c(0.55, 0.75), legend.background = element_blank())


# vis episewer posteriors -------------------------------------------------

episewer_E <- read_rds(here::here("results", "episewer", "episewer_E.rds"))
episewer_G2 <- read_rds(here::here("results", "episewer", "episewer_G2.rds"))
episewer_G1 <- read_rds(here::here("results", "episewer", "episewer_G1.rds"))

episewer_rtE <- episewer_E[["summary"]][["R"]] 
episewer_rtG2 <- episewer_G2[["summary"]][["R"]] 
episewer_rtG1 <- episewer_G1[["summary"]][["R"]] 

episewer_rtE_long <- episewer_rtE %>% 
  group_by(date) %>%
  pivot_longer(cols = c("lower_0.95", "lower_0.8", "lower_0.5"), names_to = ".width", values_to = ".lower") %>% 
  pivot_longer(cols = c("upper_0.95", "upper_0.8", "upper_0.5"), names_to = ".width2", values_to = ".upper") %>%
  mutate(.width = ifelse(.width == "lower_0.5" & .width2 == "upper_0.5", 0.5,
                         ifelse(.width == "lower_0.8" & .width2 == "upper_0.8", 0.8, 
                                ifelse(.width == "lower_0.95" & .width2 == "upper_0.95", 0.95, NA)))) %>%
  filter(.width == 0.5 | .width == 0.8 | .width == 0.95) 


episewer_rtG2_long <- episewer_rtG2 %>% 
  group_by(date) %>%
  pivot_longer(cols = c("lower_0.95", "lower_0.8", "lower_0.5"), names_to = ".width", values_to = ".lower") %>% 
  pivot_longer(cols = c("upper_0.95", "upper_0.8", "upper_0.5"), names_to = ".width2", values_to = ".upper") %>%
  mutate(.width = ifelse(.width == "lower_0.5" & .width2 == "upper_0.5", 0.5,
                         ifelse(.width == "lower_0.8" & .width2 == "upper_0.8", 0.8, 
                                ifelse(.width == "lower_0.95" & .width2 == "upper_0.95", 0.95, NA)))) %>%
  filter(.width == 0.5 | .width == 0.8 | .width == 0.95) 

episewer_rtG1_long <- episewer_rtG1 %>% 
  group_by(date) %>%
  pivot_longer(cols = c("lower_0.95", "lower_0.8", "lower_0.5"), names_to = ".width", values_to = ".lower") %>% 
  pivot_longer(cols = c("upper_0.95", "upper_0.8", "upper_0.5"), names_to = ".width2", values_to = ".upper") %>%
  mutate(.width = ifelse(.width == "lower_0.5" & .width2 == "upper_0.5", 0.5,
                         ifelse(.width == "lower_0.8" & .width2 == "upper_0.8", 0.8, 
                                ifelse(.width == "lower_0.95" & .width2 == "upper_0.95", 0.95, NA)))) %>%
  filter(.width == 0.5 | .width == 0.8 | .width == 0.95) 

# graphs ------------------------------------------------------------------
episewer_E_posterior_rt <- episewer_rtE_long %>% 
  bind_rows(extra_rows) %>%
  filter(min(E_stc_rt_quantiles$date) <= date) %>%
  ggplot(aes(x = date, y = median,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Episewer E") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("Date") + 
  theme(legend.position = "none")

episewer_G2_posterior_rt <- episewer_rtG2_long %>% 
  filter(min(G2_stc_rt_quantiles$date) <= date) %>%
  ggplot(aes(x = date, y = median,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Episewer G2") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("Date") + 
  theme(legend.position = "none")
episewer_G1_posterior_rt <- episewer_rtG1_long %>% 
  filter(min(G1_stc_rt_quantiles$date) <= date) %>%
  ggplot(aes(x = date, y = median,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Episewer G1") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("Date") + 
  theme(legend.position = "none")
E_plot <- (E_det_plot_posterior_rt +  E_stc_plot_posterior_rt + episewer_E_posterior_rt +  E_case_rt_plot) + 
  plot_layout(ncol = 4, nrow = 1)&
  theme(text = element_text(size = 20)) + 
  plot_layout(nrow = 3)
G1_plot <-  (G1_det_plot_posterior_rt + G1_stc_plot_posterior_rt + episewer_G1_posterior_rt+ G_case_rt_plot) + 
  plot_layout(ncol = 4, nrow = 1)&
  theme(text = element_text(size = 20)) + 
  plot_layout(nrow = 3)
G2_plot <- (G2_det_plot_posterior_rt + G2_stc_plot_posterior_rt + episewer_G2_posterior_rt +  G_case_rt_plot_bottom) + 
  plot_layout(ncol = 4, nrow = 1)&
  theme(text = element_text(size = 20)) + 
  plot_layout(nrow = 3)
episewer_posterior_grid <- plot_grid(
  E_plot,
  G1_plot,
  G2_plot,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1, 1)) 

# ggsave(here::here("figures", "uci_all_posterior_fig_episewer.pdf"), episewer_posterior_grid, height = 15, width = 15)


# rt=1 --------------------------------------------------------------------

# G1 ---------------------------------------------------------------------
G1_stc_timevarying_quantiles <- read_csv(here::here("results", 
                                                    "stoch_ei_normal", 
                                                    "generated_quantities", 
                                                    paste0("posterior_timevarying_quantiles_scenario", 
                                                           "G1_rt1", "_seed", seed_val, ".csv")),
                                         show_col_types = FALSE) 


G1_stc_rt_quantiles <- G1_stc_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(G1_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

# visualize results
G1_stc_plot_posterior_rt <- G1_stc_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Stc-ww G1") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = "none")


# redo for det
G1_det_timevarying_quantiles <- read_csv(here::here("results", 
                                                    "det_ei_normal", 
                                                    "generated_quantities", 
                                                    paste0("posterior_timevarying_quantiles_scenario", 
                                                           "G1_rt1", "_seed", seed_val, ".csv")),
                                         show_col_types = FALSE) 


G1_det_rt_quantiles <- G1_det_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(G1_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

# visualize results
G1_det_plot_posterior_rt <- G1_det_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Det-ww G1") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("") + 
  theme(legend.position = "none")



# G2 ---------------------------------------------------------------------
G2_stc_timevarying_quantiles <- read_csv(here::here("results", 
                                                    "stoch_ei_normal", 
                                                    "generated_quantities", 
                                                    paste0("posterior_timevarying_quantiles_scenario", 
                                                           "G2_rt1", "_seed", seed_val, ".csv")),
                                         show_col_types = FALSE) 


G2_stc_rt_quantiles <- G2_stc_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(G2_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

# visualize results
G2_stc_plot_posterior_rt <- G2_stc_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Stc-ww G2") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("Date") + 
  theme(legend.position = "none")



# redo for det
G2_det_timevarying_quantiles <- read_csv(here::here("results", 
                                                    "det_ei_normal", 
                                                    "generated_quantities", 
                                                    paste0("posterior_timevarying_quantiles_scenario", 
                                                           "G2_rt1", "_seed", seed_val, ".csv")),
                                         show_col_types = FALSE) 


G2_det_rt_quantiles <- G2_det_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(G2_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

# visualize results
G2_det_plot_posterior_rt <- G2_det_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Det-ww G2") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("Date") + 
  theme(legend.position = "none")


# E ------------------------------------------------------------------

E_stc_timevarying_quantiles <- read_csv(here::here("results", 
                                                   "stoch_ei_normal", 
                                                   "generated_quantities", 
                                                   paste0("posterior_timevarying_quantiles_scenario", 
                                                          "E_rt1", "_seed", seed_val, ".csv")),
                                        show_col_types = FALSE) 


E_stc_rt_quantiles <- E_stc_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(E_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

extra_rows <- G2_stc_rt_quantiles %>% 
  filter(date == min(date))

extra_rows$value <- NA
extra_rows$.lower <- NA
extra_rows$.upper <- NA

E_stc_rt_quantiles <- bind_rows(extra_rows, E_stc_rt_quantiles)
# visualize results
E_stc_plot_posterior_rt <- E_stc_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Stc-ww E") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = "none")



# redo for det
E_det_timevarying_quantiles <- read_csv(here::here("results", 
                                                   "det_ei_normal", 
                                                   "generated_quantities", 
                                                   paste0("posterior_timevarying_quantiles_scenario", 
                                                          "E_rt1", "_seed", seed_val, ".csv")),
                                        show_col_types = FALSE) 


E_det_rt_quantiles <- E_det_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(E_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

E_det_rt_quantiles <- bind_rows(E_det_rt_quantiles, extra_rows)
# visualize results
E_det_plot_posterior_rt <- E_det_rt_quantiles %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Det-ww E") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("Rt") + 
  xlab("") + 
  theme(legend.position = "none")


E_stoch_vs_det <- E_det_plot_posterior_rt + E_stc_plot_posterior_rt + plot_layout(guides = "collect")

# epidemia case fits ------------------------------------------------------
# G  --------------------------------------------------------------
G_case_posteriort = read_csv(here::here("results", 
                                        "case_models", 
                                        "epidemia_G_cases_rt=1_posterior_rtv2.csv"))

G_case_posteriort =  G_case_posteriort %>%
  mutate(lump = time-1) %>% 
  dplyr::select(lump, value, .lower, .upper, .width) %>%
  right_join(G1_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

G_case_rt_plot <- G_case_posteriort %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Epidemia-cases G") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = "none")


G_case_rt_plot_bottom <- G_case_posteriort %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Epidemia-cases G") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("Date") + 
  theme(legend.position = "none")


# E ------------------------------------------------------------------

E_case_posteriort = read_csv(here::here("results", 
                                        "case_models", 
                                        "epidemia_E_cases_rt=1_posterior_rtv2.csv"))

E_case_posteriort =  E_case_posteriort %>%
  mutate(lump = time-1) %>% 
  dplyr::select(lump, value, .lower, .upper, .width) %>%
  right_join(E_time_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, new_time, date, value, .lower, .upper, .width) %>%
  distinct()

E_case_posteriort = bind_rows(E_case_posteriort, extra_rows)
E_case_rt_plot <- E_case_posteriort %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = date, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Epidemia-cases E") + 
  my_theme + 
  ylim(c(0,3.5)) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = c(0.55, 0.75), legend.background = element_blank())
# every posterior together ------------------------------------------------

all_posterior_fig <-  (E_det_plot_posterior_rt +  E_stc_plot_posterior_rt + E_case_rt_plot)/
  (G1_det_plot_posterior_rt + G1_stc_plot_posterior_rt + G_case_rt_plot)/
  (G2_det_plot_posterior_rt + G2_stc_plot_posterior_rt + G_case_rt_plot_bottom) + plot_annotation(
    title = 'Estimated Rt from Wastewater and Cases Initial Rt ~ 1'
  ) &
  theme(text = element_text(size = 18))

# ggsave(here::here("figures", "uci_all_posterior_rt=1_fig.pdf"), all_posterior_fig, height = 12, width = 12)

