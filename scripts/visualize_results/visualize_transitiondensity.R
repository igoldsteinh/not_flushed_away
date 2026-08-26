#compare ei MJP vs ei LN
library(tidyverse)
library(patchwork)
library(tidybayes)


# read in sim values ------------------------------------------------------
mjp_10 <- read_csv(here::here("data", "compare_data", "ei_mjp10.csv")) %>%
  mutate(Model = "MJP",
         Count = 10)
LN_10 <- read_csv(here::here("data", "compare_data", "ei_SDE10.csv")) %>%
  mutate(Model = "LN",
         Count = 10)
mjp_5 <- read_csv(here::here("data", "compare_data", "ei_mjp5.csv")) %>%
  mutate(Model = "MJP",
         Count = 5)
LN_5 <- read_csv(here::here("data", "compare_data", "ei_SDE5.csv")) %>%
  mutate(Model = "LN",
         Count = 5)
mjp_20 <- read_csv(here::here("data", "compare_data", "ei_mjp20.csv")) %>%
  mutate(Model = "MJP",
         Count = 20)
LN_20 <- read_csv(here::here("data", "compare_data", "ei_SDE20.csv")) %>%
  mutate(Model = "LN",
         Count = 20)


# visualize the results ---------------------------------------------------
E_plot <- mjp_10 %>%
  bind_rows(LN_10) %>%
  bind_rows(mjp_5) %>%
  bind_rows(LN_5) %>%
  bind_rows(mjp_20) %>%
  bind_rows(LN_20) %>%
  ggplot(aes(x = as.factor(Count), y = E, fill = Model)) + 
  geom_boxplot() + 
  theme_minimal() +
  xlab("Initial Count")

I_plot <- mjp_10 %>%
  bind_rows(LN_10) %>%
  bind_rows(mjp_5) %>%
  bind_rows(LN_5) %>%
  bind_rows(mjp_20) %>%
  bind_rows(LN_20) %>%
  ggplot(aes(x = as.factor(Count), y = I, fill = Model)) + 
  geom_boxplot() + 
  theme_minimal() +
  xlab("Initial Count")

combined_plot <- E_plot + I_plot + 
  plot_layout(guides = 'collect') & theme(text = element_text(size = 18))
# ggsave(here::here("figures", "LN_vs_Gill_E_dissertation.pdf"), dissertation_fig, height = 8, width = 12)


# 2d plots ----------------------------------------------------------------

twod_plot <- mjp_10 %>%
  bind_rows(LN_10) %>%
  bind_rows(mjp_5) %>%
  bind_rows(LN_5) %>%
  bind_rows(mjp_20) %>%
  bind_rows(LN_20) %>%
  ggplot(aes(x = E, y = I)) + 
  geom_density_2d(aes(colour = Model), contour_var = "ndensity") +
  facet_wrap(vars(Count), scales ="free") +
  ggtitle("Transition Densities: MJP vs LN") +
  theme(text = element_text(size = 20))
# ggsave(here::here("figures", "LN_vs_MJP_density_plot.pdf"), height = 8, width = 12)