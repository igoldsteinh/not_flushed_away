# freq mets with new shedding
library(tidyverse)
library(patchwork)
source(here::here("src", "utility_functions.R"))
# first check the stan diags ----------------------------------------------
# stoch
stan102 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario102_allseeds_stan_diag.csv"))
stan115 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario115_allseeds_stan_diag.csv"))
stan106 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario106_allseeds_stan_diag.csv"))
stan127 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario127_allseeds_stan_diag.csv"))
fails102 <- stan102 %>%
  filter(max_rhat > 1.05 | min_ess_bulk < 100 | min_ess_tail < 100) %>% pull(seed)
fails115 <- stan115 %>%
  filter(max_rhat > 1.05 | min_ess_bulk < 100 | min_ess_tail < 100) %>% pull(seed)
fails106 <- stan106 %>%
  filter(max_rhat > 1.05 | min_ess_bulk < 100 | min_ess_tail < 100) %>% pull(seed)
fails127 <- stan127 %>%
  filter(max_rhat > 1.05 | min_ess_bulk < 100 | min_ess_tail < 100) %>% pull(seed)
# det
det_stan102 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario102_allseeds_stan_diag.csv"))
det_stan115 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario115_allseeds_stan_diag.csv"))
det_stan106 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario106_allseeds_stan_diag.csv"))
det_stan127 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario127_allseeds_stan_diag.csv"))
det_fails102 <- det_stan102 %>%
  filter(max_rhat > 1.05 | min_ess_bulk < 100 | min_ess_tail < 100) %>% pull(seed)
det_fails115 <- det_stan115 %>%
  filter(max_rhat > 1.05 | min_ess_bulk < 100 | min_ess_tail < 100) %>% pull(seed)
det_fails106 <- det_stan106 %>%
  filter(max_rhat > 1.05 | min_ess_bulk < 100 | min_ess_tail < 100) %>% pull(seed)
det_fails127 <- det_stan127 %>%
  filter(max_rhat > 1.05 | min_ess_bulk < 100 | min_ess_tail < 100) %>% pull(seed)


# read in the rt quantiles ------------------------------------------------
width = 0.8
# stoch 
stoch_102 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario102_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(seed %in% fails102))
stoch_115 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario115_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(seed %in% fails115)) %>%
  filter(seed != 84)
stoch_106 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario106_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(seed %in% fails106))
stoch_127 <- read_csv(here::here("results", "stoch_ei_normal", "stoch_ei_normal_scenario127_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(seed %in% fails127))
# det
det_102 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario102_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(seed %in% det_fails102))
det_115 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario115_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(seed %in% det_fails115))
det_106 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario106_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(seed %in% det_fails106))
det_127 <- read_csv(here::here("results", "det_ei_normal", "det_ei_normal_scenario127_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(seed %in% det_fails127))
# calculate metrics -------------------------------------------------------
#stoch
stoch_102_metrics = NULL
stoch_102_seeds = unique(stoch_102$seed)
for (i in 1:length(stoch_102_seeds)) {
  seed_val = stoch_102_seeds[i]
  sub_frame = stoch_102 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = seed_val,
                                                                                            model = "Stoch Fixed")
  stoch_102_metrics = bind_rows(stoch_102_metrics, metrics)
}
stoch_115_metrics = NULL
stoch_115_seeds = unique(stoch_115$seed)
for (i in 1:length(stoch_115_seeds)) {
  seed_val = stoch_115_seeds[i]
  sub_frame = stoch_115 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = seed_val,
                                                                                            model = "Stoch Steep")
  stoch_115_metrics = bind_rows(stoch_115_metrics, metrics)
}
stoch_106_metrics = NULL
stoch_106_seeds = unique(stoch_106$seed)
for (i in 1:length(stoch_106_seeds)) {
  seed_val = stoch_106_seeds[i]
  sub_frame = stoch_106 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = seed_val,
                                                                                            model = "Stoch Total500")
  stoch_106_metrics = bind_rows(stoch_106_metrics, metrics)
}
stoch_127_metrics = NULL
stoch_127_seeds = unique(stoch_127$seed)
for (i in 1:length(stoch_127_seeds)) {
  seed_val = stoch_127_seeds[i]
  sub_frame = stoch_127 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = seed_val,
                                                                                            model = "Stoch Total2000")
  stoch_127_metrics = bind_rows(stoch_127_metrics, metrics)
}
#det
det_102_metrics = NULL
det_102_seeds = unique(det_102$seed)
for (i in 1:length(det_102_seeds)) {
  seed_val = det_102_seeds[i]
  sub_frame = det_102 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = seed_val,
                                                                                            model = "Det Fixed")
  det_102_metrics = bind_rows(det_102_metrics, metrics)
}
det_115_metrics = NULL
det_115_seeds = unique(det_115$seed)
for (i in 1:length(det_115_seeds)) {
  seed_val = det_115_seeds[i]
  sub_frame = det_115 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = seed_val,
                                                                                            model = "Det Steep")
  det_115_metrics = bind_rows(det_115_metrics, metrics)
}
det_106_metrics = NULL
det_106_seeds = unique(det_106$seed)
for (i in 1:length(det_106_seeds)) {
  seed_val = det_106_seeds[i]
  sub_frame = det_106 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = seed_val,
                                                                                            model = "Det Total500")
  det_106_metrics = bind_rows(det_106_metrics, metrics)
}
det_127_metrics = NULL
det_127_seeds = unique(det_127$seed)
for (i in 1:length(det_127_seeds)) {
  seed_val = det_127_seeds[i]
  sub_frame = det_127 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = seed_val,
                                                                                            model = "Det Total2000")
  det_127_metrics = bind_rows(det_127_metrics, metrics)
}

# episewer comparison -----------------------------------------------------
episewer_102 <- read_csv(here::here("results", "episewer", "episewer102_allseeds_rt_quantiles.csv")) %>%
  filter(seed != 71)
episewer_115 <- read_csv(here::here("results", "episewer", "episewer115_allseeds_rt_quantiles.csv")) 



episewer_102_metrics <- NULL
episewer_102_seeds = unique(episewer_102$seed)
for (i in 1:length(episewer_102_seeds)) {
  seed_val = episewer_102_seeds[i]
  sub_frame = episewer_102 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = median, upper = upper_0.8, lower = lower_0.8) %>% mutate(seed = seed_val,
                                                                                                   model = "Episewer")
  episewer_102_metrics = bind_rows(episewer_102_metrics, metrics)
}
episewer_102_metrics <- episewer_102_metrics %>% mutate(model = "Episewer Fixed")
episewer_115_metrics <- NULL
episewer_115_seeds = unique(episewer_115$seed)
for (i in 1:length(episewer_115_seeds)) {
  seed_val = episewer_115_seeds[i]
  sub_frame = episewer_115 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = median, upper = upper_0.8, lower = lower_0.8) %>% mutate(seed = seed_val,
                                                                                                   model = "Episewer Steep")
  episewer_115_metrics = bind_rows(episewer_115_metrics, metrics)
}

# merge metrics -----------------------------------------------------------
all_metrics <- bind_rows(stoch_102_metrics,
                         stoch_115_metrics,
                         stoch_106_metrics,
                         stoch_127_metrics,
                         det_102_metrics,
                         det_115_metrics,
                         det_106_metrics,
                         det_127_metrics,
                         episewer_102_metrics,
                         episewer_115_metrics) %>% 
  pivot_longer(cols = -c(seed, model), names_to = "metric")
all_metrics_summary <- all_metrics %>% 
  group_by(model, metric) %>%
  summarise(mean_value = mean(value),
            median_value = median(value))
level_list <- c("mean_dev", "mean_env", "MCIW", "MASV", "true_MASV")
label_list <- c("Deviation", "Envelope", "MCIW", "MASV", "True MASV")
all_metrics$metric <- factor(all_metrics$metric, levels=level_list, labels=label_list)
model_level_list <- c(
                      "Det Steep",
                      "Stoch Steep",
                      "Episewer Steep",
                      "Det Fixed",
                      "Stoch Fixed",
                      "Episewer Fixed",
                      "Det Total500",
                      "Stoch Total500",
                      "Det Total2000",
                      "Stoch Total2000")
all_metrics$model <- factor(all_metrics$model, levels = model_level_list)
all_metrics$useful_value <- 0
all_metrics$useful_value[all_metrics$metric == "Envelope"] <- width
all_metrics$useful_value[all_metrics$metric == "MCIW"] <- NA
all_metrics$useful_value[all_metrics$metric == "MASV"] <- NA

# base plot ---------------------------------------------------------------
base_list <- c(
               "Det Steep",
               "Stoch Steep",
               "Episewer Steep",
               "Det Fixed",
               "Stoch Fixed",
               "Episewer Fixed")
all_metric_dev_plot <- all_metrics %>% 
  filter(model %in% base_list) %>% 
  filter(metric == "Deviation") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Deviation") + 
  ylab("Deviation") + 
  xlab("")
all_metric_env_plot <- all_metrics %>% 
  filter(model %in% base_list) %>% 
  filter(metric == "Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Envelope") + 
  ylab("Envelope") + 
  xlab("")
all_metric_mciw_plot <- all_metrics %>% 
  filter(model %in% base_list) %>% 
  filter(metric == "MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("MCIW") + 
  ylab("MCIW") + 
  xlab("Model")
true_masv <- all_metrics %>% 
  filter(model %in% base_list) %>%
  filter(metric == "True MASV") %>%
  mutate(model = ifelse(model == "Det Steep" | model == "Stoch Steep" | model == "Episewer Steep", "True Steep",
                               "True Fixed"),
         metric = "MASV")
masv_plot_data <- all_metrics %>% 
  filter(model %in% base_list) %>% 
  filter(metric == "MASV") %>%
  bind_rows(true_masv)
masv_label_list = c(
                    "Det Steep",
                    "Stoch Steep",
                    "Episewer Steep",
                    "True Steep",
                    "Det Fixed",
                    "Stoch Fixed",
                    "Episewer Fixed",
                    "True Fixed")

masv_plot_data$model <- factor(masv_plot_data$model, levels = masv_label_list, labels = masv_label_list)
all_metric_masv_plot <- masv_plot_data %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("MASV") + 
  ylab("MASV") + 
  xlab("Model")
new_all_metric_plot <- all_metric_dev_plot + all_metric_env_plot + all_metric_mciw_plot + all_metric_masv_plot + plot_annotation(
  title = 'Frequentist Metrics for Stc vs Det vs Episewer Models'
) &
  theme(text = element_text(size = 18))
ggsave(here::here("figures", "freq_metrics_baseline_newshedding2.pdf"), new_all_metric_plot, height = 11, width = 11)
# varying stoch fig -------------------------------------------------------
stoch_list <- c( 
                 "Det Total500",
                 "Stoch Total500",
                 "Det Fixed",
                 "Stoch Fixed",
                 "Det Total2000",
                 "Stoch Total2000")
all_metrics$model <- factor(all_metrics$model, levels = stoch_list)

stoch_metric_dev_plot <- all_metrics %>% 
  filter(model %in% stoch_list) %>% 
  filter(metric == "Deviation") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Deviation") + 
  ylab("Deviation") + 
  xlab("")
stoch_metric_env_plot <- all_metrics %>% 
  filter(model %in% stoch_list) %>% 
  filter(metric == "Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Envelope") + 
  ylab("Envelope") + 
  xlab("")
stoch_metric_mciw_plot <- all_metrics %>% 
  filter(model %in% stoch_list) %>% 
  filter(metric == "MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("MCIW") + 
  ylab("MCIW") + 
  xlab("Model")
stoch_true_masv <- all_metrics %>% 
  filter(model %in% stoch_list) %>%
  filter(metric == "True MASV") %>%
  mutate(model = ifelse(model == "Det Fixed" | model == "Stoch Fixed", 
                        "True Fixed",
                        ifelse(model == "Det Total500" | model == "Stoch Total500", 
                                      "True Total500",
                                             "True Total2000")),
         metric = "MASV")
stoch_masv_plot_data <- all_metrics %>% 
  filter(model %in% stoch_list) %>% 
  filter(metric == "MASV") %>%
  bind_rows(stoch_true_masv)
stoch_masv_label_list = c(
                          "Det Total500",
                          "Stoch Total500",
                          "True Total500",
                          "Det Fixed", 
                          "Stoch Fixed", 
                          "True Fixed",
                          "Det Total2000",
                          "Stoch Total2000",
                          "True Total2000")

stoch_masv_plot_data$model <- factor(stoch_masv_plot_data$model, levels = stoch_masv_label_list, labels = stoch_masv_label_list)
stoch_metric_masv_plot <- stoch_masv_plot_data %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("MASV") + 
  ylab("MASV") + 
  xlab("Model")
new_stoch_metric_plot <- stoch_metric_dev_plot + stoch_metric_env_plot + stoch_metric_mciw_plot + stoch_metric_masv_plot + plot_annotation(
  title = 'Frequentist Metrics for Stoch vs Det Models (Varying Stochasticity)'
) &
  theme(text = element_text(size = 18))
ggsave(here::here("figures", "freq_metrics_stoch_newshedding.pdf"), new_stoch_metric_plot, height = 11, width = 11)


# episewer comparison -----------------------------------------------------
episewer_102 <- read_csv(here::here("results", "episewer", "episewer102_allseeds_rt_quantiles_mod.csv")) 
episewer_115 <- read_csv(here::here("results", "episewer", "episewer115_allseeds_rt_quantiles.csv")) 

episewer_102_metrics <- NULL
episewer_102_seeds = unique(episewer_102$seed)
for (i in 1:length(episewer_102_seeds)) {
  seed_val = episewer_102_seeds[i]
  sub_frame = episewer_102 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = median, upper = upper_0.8, lower = lower_0.8) %>% mutate(seed = seed_val,
                                                                                            model = "Episewer")
  episewer_102_metrics = bind_rows(episewer_102_metrics, metrics)
}
episewer_102_metrics <- episewer_102_metrics %>% mutate(model = "Episewer Fixed")
episewer_115_metrics <- NULL
episewer_115_seeds = unique(episewer_115$seed)
for (i in 1:length(episewer_115_seeds)) {
  seed_val = episewer_115_seeds[i]
  sub_frame = episewer_115 %>% filter(seed == seed_val) 
  metrics = rt_metrics(sub_frame, value = median, upper = upper_0.8, lower = lower_0.8) %>% mutate(seed = seed_val,
                                                                                                   model = "Episewer Steep")
  episewer_115_metrics = bind_rows(episewer_115_metrics, metrics)
}

# merge metrics -----------------------------------------------------------
episewer_metrics <- bind_rows(stoch_102_metrics,
                         episewer_102_metrics, 
                         stoch_115_metrics,
                         episewer_115_metrics) %>% 
  pivot_longer(cols = -c(seed, model), names_to = "metric")
episewer_metrics_summary <- episewer_metrics %>% 
  group_by(model, metric) %>%
  summarise(mean_value = mean(value),
            median_value = median(value))
level_list <- c("mean_dev", "mean_env", "MCIW", "MASV", "true_MASV")
label_list <- c("Deviation", "Envelope", "MCIW", "MASV", "True MASV")
episewer_metrics$metric <- factor(episewer_metrics$metric, levels=level_list, labels=label_list)
model_level_list <- c(
  "Stoch Fixed",
  "Episewer Fixed", 
  "Stoch Steep", 
  "Episewer Steep")
episewer_metrics$model <- factor(episewer_metrics$model, levels = model_level_list)
episewer_metrics$useful_value <- 0
episewer_metrics$useful_value[episewer_metrics$metric == "Envelope"] <- width
episewer_metrics$useful_value[episewer_metrics$metric == "MCIW"] <- NA
episewer_metrics$useful_value[episewer_metrics$metric == "MASV"] <- NA

# base plot ---------------------------------------------------------------
base_list <- c(
  "Stoch Fixed",
  "Episewer Fixed",
  "Stoch Steep",
  "Episewer Steep")
all_metric_dev_plot <- episewer_metrics %>% 
  filter(model %in% base_list) %>% 
  filter(metric == "Deviation") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Deviation") + 
  ylab("Deviation") + 
  xlab("")
all_metric_env_plot <- episewer_metrics %>% 
  filter(model %in% base_list) %>% 
  filter(metric == "Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Envelope") + 
  ylab("Envelope") + 
  xlab("")
all_metric_mciw_plot <- episewer_metrics %>% 
  filter(model %in% base_list) %>% 
  filter(metric == "MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("MCIW") + 
  ylab("MCIW") + 
  xlab("Model")
true_masv <- episewer_metrics %>% 
  filter(model %in% base_list) %>%
  filter(metric == "True MASV") %>%
  mutate(model = ifelse(model == "Det Steep" | model == "Stoch Steep", "True Steep",
                        ifelse(model == "Det Steep" | model == "Stoch Steep", "True Steep",
                               "True Fixed")),
         metric = "MASV") %>%
  filter(!is.na(value))
masv_plot_data <- episewer_metrics %>% 
  filter(model %in% base_list) %>% 
  filter(metric == "MASV") %>%
  bind_rows(true_masv)
masv_label_list = c(
  "Stoch Fixed",
  "Episewer Fixed", 
  "True Fixed",
  "Stoch Steep", 
  "Episewer Steep",
  "True Steep"
)

masv_plot_data$model <- factor(masv_plot_data$model, levels = masv_label_list, labels = masv_label_list)
all_metric_masv_plot <- masv_plot_data %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("MASV") + 
  ylab("MASV") + 
  xlab("Model")
episewer_all_metric_plot <- all_metric_dev_plot + all_metric_env_plot + all_metric_mciw_plot + all_metric_masv_plot + plot_annotation(
  title = 'Frequentist Metrics for Stc vs Episewer Models'
) &
  theme(text = element_text(size = 18))
# ggsave(here::here("figures", "freq_metrics_episewer.pdf"), episewer_all_metric_plot, height = 11, width = 11)

