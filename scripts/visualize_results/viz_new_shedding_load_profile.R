suppressPackageStartupMessages(library(tidyverse))
library(tidybayes)
suppressPackageStartupMessages(library(posterior))
library(fs)
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(kableExtra))

source(here::here("src", "utility_functions.R"))

my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme())


dumb_weights <- c(0.05, 0.4, 0.8, 0.4, 0.2, 0.1, 0.05, 0.025, 0.0125)
norm_weights <- dumb_weights / sum(dumb_weights)
comp_weights = log(norm_weights * 10e6, base = 10)
comp_names = c("I1", "I2", "I3", "I4", "I5", "I6", "I7", "R1", "R2")
weight_frame <- data.frame(comp_names, comp_weights)

# plot the weights
weight_plot <- ggplot(weight_frame, aes(x = comp_names, y = comp_weights)) +
  geom_point(stat = "identity") +
  labs(title = "Shedding in each compartment",
       x = "Compartment",
       y = "Shedding (Log base 10)") +
  theme_minimal()

ggsave(here::here("figures", "compartment_shedding_weights.pdf"), weight_plot, width = 6, height = 4)



# example shedding --------------------------------------------------------

individ_conc_data <- read_csv(here::here("data", "newshedding_data", "example_shedding.csv"))

individ_shedding <- individ_conc_data %>% 
  bind_rows() %>%
  group_by(sim) %>%
  filter(sim %in% c(1:9)) %>%
  dplyr::select(sim, time, 
                total_conc, 
                I1 = total_conc_I1,
                I2 = total_conc_I2,
                I3 = total_conc_I3,
                I4 = total_conc_I4,
                I5 = total_conc_I5,
                I6 = total_conc_I6,
                I7 = total_conc_I7) %>%
  pivot_longer(-c(time,sim)) %>%
  filter(name != "total_conc") %>%
  filter(value > 0) %>%
  rename(Compartment = name) %>%
  ggplot(aes(x = time, y = value, color = Compartment)) +
  geom_point() +
  geom_line() + 
  scale_x_continuous(breaks = seq(0, 14, by = 1)) +
  facet_wrap(vars(sim), scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Time since start of Infectiousness") +
  ylab("Shedding") + 
  theme(text = element_text(size = 18)) + 
  ggtitle("Simulated shedding for 9 individuals")

ggsave(here::here("figures", "example_shedding.pdf"), individ_shedding, width = 12, height = 8)
