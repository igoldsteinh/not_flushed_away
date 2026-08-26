# file for truncating data to May 1st
library(tidyverse)
# E ww --------------------------------------------------------------------
end_date <- as.Date("2022-05-01")
E_data <- read_csv(here::here("data", "E_data.csv"))

E_data_trunc <- E_data %>% 
  filter(date <= end_date)

write_csv(E_data_trunc, here::here("data", "E_data_trunc.csv"))


# E cases -----------------------------------------------------------------
E_case_data <- read_csv(here::here("data", "E_case_data.csv"))

# correcting dates 
E_case_data$date[E_case_data$lump == 1] <- "2022-02-07"
E_case_data$date[E_case_data$lump == 3] <- "2022-02-21"
E_case_data$date[E_case_data$lump == 5] <- "2022-03-07"
E_case_data$date[E_case_data$lump == 6] <- "2022-03-14"
E_case_data$date[E_case_data$lump == 9] <- "2022-04-04"
E_case_data$date[E_case_data$lump== 12] <- "2022-04-25"
E_case_data_trunc <- E_case_data %>%
  filter(date <= end_date)

write_csv(E_case_data_trunc, here::here("data", "E_case_data_trunc.csv"))
write_csv(E_case_data, here::here("data", "E_case_data_corrected.csv"))

# G ww  -------------------------------------------------------------------
G_data <- read_csv(here::here("data", "Gcombined_data.csv"))

G_data_trunc <- G_data %>% 
  filter(date <= end_date)

write_csv(G_data_trunc, here::here("data", "Gcombined_data_trunc.csv"))


# G case data truncated ---------------------------------------------------
G_case_data <- read_csv(here::here("data", "G_case_data.csv"))
# rectifying labeling error in the dates
# the weeks are correct
# the corresponding date is too early and should be changed for visualization 
# this does not impact inference at all
G_case_data$date[G_case_data$date == "2022-02-17"] = "2022-02-21"
G_case_data$date[G_case_data$lump == 5] = "2022-03-07"
G_case_date_trunc <- G_case_data %>% 
  filter(date <= end_date)

write_csv(G_case_date_trunc, here::here("data", "G_case_data_trunc.csv"))
write_csv(G_case_data, here::here("data", "G_case_data_corrected.csv"))
