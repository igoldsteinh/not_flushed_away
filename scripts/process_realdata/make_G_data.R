# file for combining G1 and G2
library(tidyverse)


# read in G1 and G2 -------------------------------------------------------
G1_data <- read_csv(here::here("data", "G1_data.csv"))

G2_data <- read_csv(here::here("data", "G2_data.csv"))

G_data <- G1_data %>%
  bind_rows(G2_data)%>%
  mutate(conc = exp(log_conc)) %>%
  group_by(date, new_time) %>%
  summarise(
    mean_conc = mean(conc),
    log_mean_conc = log(mean_conc),
    n_places = n_distinct(place)
  ) 

write_csv(G_data, here::here("data", "Gcombined_data.csv"))


# shorten G data ----------------------------------------------------------
end_date <- as.Date("2022-05-01")

G_data_trunc <- G_data %>% 
  filter(date <= end_date)
write_csv(G_data_trunc, here::here("data", "Gcombined_data_trunc.csv"))
