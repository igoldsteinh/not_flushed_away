#visualize individual fits for paper
library(tidyverse)
library(patchwork)
library(scales)
library(tidybayes)
# stoch 
stoch_102 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario102_allseeds_rt_quantiles.csv")) 
stoch_115 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario115_allseeds_rt_quantiles.csv")) 
det_102 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario102_allseeds_rt_quantiles.csv"))
det_115 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario115_allseeds_rt_quantiles.csv")) 
my_theme <- list(
  scale_fill_brewer(name = "CI Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw(),
  theme())
episewer_115 <- read_csv(here::here("results", "episewer", "episewer115_allseeds_rt_quantiles.csv")) 
episewer_102_seed1 <- read_csv(here::here("results", "episewer", "episewer_sim_102_seed_1.csv"))

# steep visualization -----------------------------------------------------

# data
key_frame <- read_csv(here::here("data", "newshedding_data", "sim_key.csv"))
datafilename = key_frame %>% 
  filter(sim_num == 115) %>% 
  pull(obsdata_filename) %>% 
  unique() %>%
  as.character()

steep_seed = 1
fitted_simdata <- read_csv(here::here(datafilename), show_col_types = FALSE)


seed_data <- fitted_simdata %>% 
  filter(seed == steep_seed)

data_plot <- seed_data %>% 
  dplyr::select(new_time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>% 
  pivot_longer(-new_time) %>% 
  group_by(new_time) %>%
  mutate(mean_val = mean(value)) %>%
  ggplot(aes(x = new_time, y = value)) + 
  geom_point() + 
  geom_line(aes(x = new_time, y = mean_val)) +
  theme_bw() +
  ggtitle(paste0("Wastewater Data")) + 
  ylab("Log Conc.") +
  xlab("") + 
  theme(text = element_text(size = 20),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())


stc_steep_quantiles <- stoch_115 %>%
                       filter(seed == steep_seed)

det_steep_quantiles <- det_115 %>%
                       filter(seed == steep_seed)


stc_steep_posterior_rt <- stc_steep_quantiles %>% 
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  theme_bw() +
  ggtitle("Stc Posterior Rt") + 
  my_theme + 
  xlab("") +
  ylab("Rt") +
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) + 
  ylim(c(0, 4))





det_steep_posterior_rt <- det_steep_quantiles %>% 
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  theme_bw() +
  ggtitle("Det Posterior Rt") + 
  my_theme + 
  xlab("") +
  ylab("Rt") +
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) + 
  ylim(c(0, 4))

# add a plot for episewer
episewer_115_seed1_long <- episewer_115 %>% 
  filter(seed == steep_seed) %>%
  group_by(time) %>%
  pivot_longer(cols = c("lower_0.95", "lower_0.8", "lower_0.5"), names_to = ".width", values_to = ".lower") %>% 
  pivot_longer(cols = c("upper_0.95", "upper_0.8", "upper_0.5"), names_to = ".width2", values_to = ".upper") %>%
  mutate(.width = ifelse(.width == "lower_0.5" & .width2 == "upper_0.5", 0.5,
                         ifelse(.width == "lower_0.8" & .width2 == "upper_0.8", 0.8, 
                                ifelse(.width == "lower_0.95" & .width2 == "upper_0.95", 0.95, NA)))) %>%
  filter(.width == 0.5 | .width == 0.8 | .width == 0.95) 
episewer_115_rt_plot <- episewer_115_seed1_long %>% 
  filter(time >= 1) %>%
  ggplot(aes(x = time, y = median, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = time, y = true_rt), color = "orange") +
  theme_bw() +
  ggtitle("") + 
  my_theme + 
  xlab("") +
  ylab("Rt") +
  theme(text = element_text(size = 20),
        legend.position = c(0.8,0.8),
        legend.background = element_blank()) + 
  ylim(c(0, 4)) + 
  ggtitle("Episewer")

steep_plot <- (data_plot + det_steep_posterior_rt + stc_steep_posterior_rt + episewer_115_rt_plot) + plot_annotation(title = "Steep Rt", 
                                                                                            theme = theme(plot.title = element_text(size = 20))) + 
  plot_layout(ncol = 4, nrow = 1)


stc_steep_posterior_rt <- stc_steep_quantiles %>% 
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  theme_bw() +
  ggtitle("Stc Posterior Rt") + 
  my_theme + 
  xlab("") +
  ylab("Rt") +
  theme(text = element_text(size = 20),
        legend.position = c(0.8,0.8),
        legend.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) + 
  ylim(c(0, 4))




# newinf ------------------------------------------------------------------

# data
datafilename = key_frame %>% 
  filter(sim_num == 102) %>% 
  pull(obsdata_filename) %>% 
  unique() %>%
  as.character()

fixed_seed = 1
fitted_simdata <- read_csv(here::here(datafilename), show_col_types = FALSE)



seed_data <- fitted_simdata %>% 
  filter(seed == fixed_seed)

newinfdata_plot <- seed_data %>% 
  dplyr::select(new_time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>% 
  pivot_longer(-new_time) %>% 
  group_by(new_time) %>%
  mutate(mean_val = mean(value)) %>%
  ggplot(aes(x = new_time, y = value)) + 
  geom_point() + 
  geom_line(aes(x = new_time, y = mean_val)) +
  theme_bw() +
  ggtitle(paste0("")) + 
  ylab("Log Conc.") +
  xlab("Time") + 
  theme(text = element_text(size = 20))


stc_newinfquantiles <- stoch_102 %>%
  filter(seed == fixed_seed)

det_newinfquantiles <- det_102 %>%
  filter(seed == fixed_seed)


stc_newinfposterior_rt <- stc_newinfquantiles %>% 
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  theme_bw() +
  ggtitle("") + 
  my_theme + 
  xlab("Time") +
  ylab("Rt") +
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank()) + 
  ylim(c(0, 4))





det_newinfposterior_rt <- det_newinfquantiles %>% 
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  theme_bw() +
  ggtitle("") + 
  my_theme + 
  xlab("Time") +
  ylab("Rt") +
  theme(text = element_text(size = 20),
        legend.position = "none") + 
  ylim(c(0, 4))

episewer_102_seed1_long <- episewer_102_seed1 %>% 
  group_by(time) %>%
  pivot_longer(cols = c("lower_0.95", "lower_0.8", "lower_0.5"), names_to = ".width", values_to = ".lower") %>% 
  pivot_longer(cols = c("upper_0.95", "upper_0.8", "upper_0.5"), names_to = ".width2", values_to = ".upper") %>%
  mutate(.width = ifelse(.width == "lower_0.5" & .width2 == "upper_0.5", 0.5,
                         ifelse(.width == "lower_0.8" & .width2 == "upper_0.8", 0.8, 
                                ifelse(.width == "lower_0.95" & .width2 == "upper_0.95", 0.95, NA)))) %>%
  filter(.width == 0.5 | .width == 0.8 | .width == 0.95) 
episewer_rt_102_plot <- episewer_102_seed1_long %>% 
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
  ylim(c(0, 4)) 


newinfplot <- (newinfdata_plot + det_newinfposterior_rt + stc_newinfposterior_rt + episewer_rt_102_plot) + plot_annotation(title = "Fixed Rt", 
                                                                                                          theme = theme(plot.title = element_text(size = 20))) +
  plot_layout(ncol = 4, nrow = 1)

full_plot <- steep_plot / newinfplot + plot_annotation(title = "Simulated Data and Posterior Rt Estimates",
                                                                      theme = theme(plot.title = element_text(size = 20)))


ggsave(here::here("figures", "example_rt_plot.pdf"), full_plot, height = 8, width = 15)
