# for sim 102 seed 32
# the original run had 12 divergent transitions
# I reran it on a different seed (33) to get a fit with fewer divergent transitions
# this file removes the output of the original run to the output of the new run 
library(tidyverse)
episewer_102 <- read_csv(here::here("results", "episewer", "episewer102_allseeds_rt_quantiles.csv")) 
new_sim32 <- read_csv(here::here("results", "episewer", "episewerlocal_sim_102_seed_32.csv"))

episewer_102 <- episewer_102 %>% 
  filter(!(sim == 102 & seed == 32)) %>% 
  bind_rows(new_sim32) 
# save 
write_csv(episewer_102, here::here("results", "episewer", "episewer102_allseeds_rt_quantiles_mod.csv"))