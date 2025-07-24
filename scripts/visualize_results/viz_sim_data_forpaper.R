# visualize simulated data for paper
library(tidyverse)
library(patchwork)
library(GGally)
library(gridExtra)
library(cowplot)
library(scales)
library(latex2exp)

seed_val = 25
sim = "seirr_normalshallowpop10"
scenario_sim = "seirr_normalshallowpop10"
fitted_simdata_address <-"data/newshedding_data/pop1000_data.csv"
full_simdata_address <- "data/newshedding_data/pop1000_fulldata.csv"
fitted_simdata <- read_csv(here::here(fitted_simdata_address), show_col_types = FALSE)
full_simdata <- read_csv(here::here(full_simdata_address), show_col_types = FALSE) 


seed_data <- fitted_simdata %>% 
  filter(seed == seed_val)
lump_val = 7
fullseed_data <- full_simdata %>% filter(seed == seed_val) %>%
  mutate(lump = floor((time - 1)/lump_val))

data_plot <- seed_data %>% 
  dplyr::select(new_time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>% 
  pivot_longer(-new_time) %>% 
  group_by(new_time) %>%
  mutate(mean_val = mean(value)) %>%
  ggplot(aes(x = new_time, y = value)) + 
  geom_point(size = 4) + 
  geom_line(aes(x = new_time, y = mean_val), linewidth = 2) +
  theme_bw() +
  ggtitle(paste0("Wastewater Data")) + 
  ylab("Log Concentration") +
  xlab("Time") 

true_curve <- fullseed_data %>%
  mutate(I = I1 + I2 + I3 + I4 + I5 + I6 + I7) %>%
  dplyr::select(new_time, S, E, I, R1, R2) %>%
  pivot_longer(cols = -c(new_time)) %>%
  rename(Compartment = name)

level_list <- c("S", "E", "I", "R1", "R2")

true_curve$Compartment <- factor(true_curve$Compartment, levels=level_list)

epi_curve <- true_curve %>%
  ggplot(aes(x = new_time, y = value, color = Compartment)) + 
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("Simulated Epidemic") +
  geom_point(size = 4) + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.7),
        legend.background = element_blank()) 


true_rt_plot <- fullseed_data %>% 
  dplyr::select(new_time, rt) %>%
  ggplot(aes(x = new_time, y = rt)) +
  theme_bw() +
  geom_point(color = "orange", size = 4) +
  ggtitle(paste0("True Rt")) + 
  xlab("Time")


combined_plot <- true_rt_plot + epi_curve + data_plot + ggtitle("Simulated Data") &
  theme(text = element_text(size = 30))


combined_plot
# ggsave(here::here("figures", "sim_ww_paper.pdf"), combined_plot, height = 8, width = 24)
