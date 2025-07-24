# visualizing trajectories of gillespie vs our log-normal model 
# simulation files are simgillEI_grid_data.R
# and stc_constant_grid_data.jl and process_stc_sims.R
library(tidyverse)
library(patchwork)
library(tidybayes)
# gillespie simulations -------------------------------------------
init_counts <- c(10,20,30)
other_gill_sims <- map(init_counts, ~read_csv(here::here("data", paste0("gillEIpop", .x, "_constrt_states.csv"))) %>% 
                         mutate(init_count = .x) %>%
                         rename(time = integer_day))

gill_sims <- bind_rows(other_gill_sims)
gill_medians <- gill_sims %>%
                 pivot_longer(-c(time, id, init_count)) %>% 
                 group_by(init_count, name, time) %>% 
                 dplyr::select(-id) %>% 
                 median_qi(.width = c(0.5, 0.8, 0.95)) %>%
                 rename(new_time = time) %>% 
                 mutate(model = "Gill")

gillE_plot <- gill_medians %>% filter(name == "E") %>%
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper)) +
  scale_fill_brewer(name = "Credible Interval Width") +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("Gillespie E compartment")


# log normal simulations --------------------------------------------------
init_counts = c(10,20,30)
seed_val = 4

prior_timevarying_quantiles <- map(init_counts, ~read_csv(here::here("data", 
                                                   paste0("stc_timevarying_quantiles_constrt_pop", 
                                                          .x, ".0_seed", seed_val, ".csv")),
                                   show_col_types = FALSE))

traj_names <- c("E", "I")

time_crosswalk <- data.frame(time = 0:96, new_time = 1:97)

ln_medians <- map2(prior_timevarying_quantiles, init_counts, ~.x %>%
  left_join(time_crosswalk, by = "time") %>%
  filter(name %in% traj_names) %>%
  mutate(init_count = .y) %>%
  dplyr::select(init_count, name, new_time, value, .lower, .upper, .width, .point, .interval) %>% 
  mutate(model = "LN"))


E_plot <- ln_medians %>%
  bind_rows() %>% 
  filter(name == "E") %>% 
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper)) +
  scale_fill_brewer(name = "Interval Width") +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("Log-Normal E Curve") 

E_plot

I_plot <- ln_medians %>%
  bind_rows() %>% 
  filter(name == "E") %>% 
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper)) +
  scale_fill_brewer(name = "Interval Width") +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("Log-Normal E Curve") 

I_plot


# combine em --------------------------------------------------------------
combined_E <- bind_rows(ln_medians, gill_medians) %>% 
            filter(.width == 0.95) %>%
            filter(name == "E") %>%
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, color = model, fill = model), alpha = 0.5) +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("Log-Normal vs Gill E Curve") 

combined_I <- bind_rows(ln_medians, gill_medians) %>% 
  filter(.width == 0.95) %>%
  filter(name == "I") %>%
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, color = model, fill = model), alpha = 0.5) +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("Log-Normal vs Gill I Curve") 


time_cut = 30

combined_E_trunc <- bind_rows(ln_medians, gill_medians) %>% 
  filter(.width == 0.95) %>%
  filter(name == "E") %>%
  filter(new_time <= time_cut) %>% 
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, color = model, fill = model), alpha = 0.5) +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("Log-Normal vs Gill E Curve (Trunc)") 

combined_I_trunc <- bind_rows(ln_medians, gill_medians) %>% 
  filter(.width == 0.95) %>%
  filter(name == "I") %>%
  filter(new_time <= time_cut) %>%
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, color = model, fill = model), alpha = 0.5) +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("Log-Normal vs Gill I Curve (Trunc)") 
# dissertation  --------------------------------------------------------

combined_I_trunc <- bind_rows(ln_medians, gill_medians) %>% 
  filter(.width == 0.95) %>%
  filter(name == "I") %>%
  mutate(model = ifelse(model == "Gill", "MJP", model)) %>%
  rename(Model = model) %>%
  filter(new_time <= time_cut) %>%
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, color = Model, fill = Model), alpha = 0.5) +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("") +
  ylab("I") +
  xlab("Time") + 
  theme(text = element_text(size = 18),
        legend.position = c(0.8, 0.6), 
        legend.background = element_blank())

combined_E_trunc <- bind_rows(ln_medians, gill_medians) %>% 
  filter(.width == 0.95) %>%
  filter(name == "E") %>%
  mutate(model = ifelse(model == "Gill", "MJP", model)) %>%
  rename(Model = model) %>%
  filter(new_time <= time_cut) %>%
  ggplot() +
  geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, color = Model, fill = Model), alpha = 0.5) +
  theme_bw() +
  facet_wrap(vars(init_count), scales = "free_y") +
  ggtitle("") +
  ylab("E") +
  xlab("Time") + 
  theme(text = element_text(size = 18),
        legend.position = "none", 
        legend.background = element_blank())

dissertation_fig <- combined_E_trunc / combined_I_trunc + plot_annotation(title = "Marginal Quantiles of E and I counts for LN vs MJP") & theme(text = element_text(size = 18))

# ggsave(here::here("figures", "LN_vs_Gill_E_dissertation.pdf"), dissertation_fig, height = 8, width = 12)
