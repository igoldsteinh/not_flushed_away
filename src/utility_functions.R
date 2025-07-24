# Utility functions
library(tidyverse)
library(lubridate)
library(tidybayes)
set.seed(1234)



# make fixed posterior samples --------------------------------------------

make_fixed_posterior_samples <- function(posterior_gq) {
  posterior_gq_samples <- posterior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]", negate = T))
  
  return(posterior_gq_samples)
}


# make time varying posterior quantiles -----------------------------------
make_timevarying_posterior_quantiles <- function(posterior_gq) {
  timevarying_posterior_quantiles <-
    posterior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]")) %>%
    mutate(time = name %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.numeric(),
           name = name %>%
             str_extract("^.+(?=\\[)") %>%
             str_remove("data_")) %>%
    group_by(name, time) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    left_join(.,tibble(time = 0:max(.$time)))
  
  return(timevarying_posterior_quantiles)
  
}



# make posterior predictive intervals -------------------------------------
# posterior_predictive = trueish_postpred
# sim_data = simdata
# cases = FALSE
# ten_sim = ten_sim_val
# three_mean = three_mean_val
make_post_pred_intervals <- function(posterior_predictive, sim_data, in_order = FALSE, num_samps = 10 ){
  if (in_order == FALSE) {
    obs_time <- sim_data %>%
      filter(total_conc > 0) %>%
      mutate(obs_index = row_number()) %>%
      dplyr::select(obs_index, new_time)
    
    posterior_predictive_samples <- posterior_predictive %>%
      pivot_longer(-c(iteration, chain)) %>%
      mutate(obs_index = name %>%
               str_extract("(?<=\\[)\\d+(?=\\])") %>%
               as.numeric(),
             name = name %>%
               str_extract("^.+(?=\\[)") %>%
               str_remove("data_")) %>%
      bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                  summarize(value = sum(value),
                            .groups = "drop"))  %>%
      left_join(obs_time, by = "obs_index")
    
  } else {
    obs_time <- sim_data %>%
      filter(total_conc > 0) %>%
      mutate(obs_index = row_number()) %>%
      dplyr::select(obs_index, new_time)
    
    posterior_predictive_samples <- posterior_predictive %>%
      pivot_longer(-c(iteration, chain)) %>%
      mutate(samp_index = name %>%
               str_extract("(?<=\\[)\\d+(?=\\])") %>%
               as.numeric(),
             name = name %>%
               str_extract("^.+(?=\\[)") %>%
               str_remove("data_")) %>%
      bind_rows(., group_by(., chain, iteration, samp_index, name) %>%
                  summarize(value = sum(value),
                            .groups = "drop"))  %>%
      mutate(obs_index = ceiling(samp_index/num_samps)) %>%
      left_join(obs_time, by = "obs_index")
    
  }
  
  posterior_predictive_intervals <- posterior_predictive_samples %>%
    select(new_time, name, value) %>%
    group_by(new_time, name) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    select(new_time, name, value, starts_with(".")) %>%
    distinct()
  return(posterior_predictive_intervals)
}


# make posterior predictive plot ------------------------------------------
make_post_pred_plot <- function(posterior_predictive_intervals, 
                                sim_data, 
                                cases = FALSE,
                                ten_sim = FALSE,
                                three_mean = FALSE) {
  
  if (cases == FALSE & ten_sim == FALSE & three_mean == FALSE) {
    true_data <- sim_data %>%
      dplyr::select(new_time, 
                    log_gene_copies1, 
                    log_gene_copies2, 
                    log_gene_copies3) %>%
      rename("log_copies1" = "log_gene_copies1",
             "log_copies2" = "log_gene_copies2",
             "log_copies3" = "log_gene_copies3") %>%
      pivot_longer(cols = - new_time) %>%
      rename("true_value" = "value") %>%
      filter(true_value > 0)
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_time"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_time, y = value)) + 
      geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      # scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue3"), name="fill") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
  } else if (cases == FALSE & ten_sim == TRUE & three_mean == FALSE) {
    true_data <- sim_data %>%
      dplyr::select(new_time, 
                    log_gene_copies1, 
                    log_gene_copies2, 
                    log_gene_copies3,
                    log_gene_copies4,
                    log_gene_copies5,
                    log_gene_copies6,
                    log_gene_copies7,
                    log_gene_copies8,
                    log_gene_copies9,
                    log_gene_copies10,
      ) %>%
      rename("log_copies1" = "log_gene_copies1",
             "log_copies2" = "log_gene_copies2",
             "log_copies3" = "log_gene_copies3",
             "log_copies4" = "log_gene_copies4",
             "log_copies5" = "log_gene_copies5",
             "log_copies6" = "log_gene_copies6",
             "log_copies7" = "log_gene_copies7",
             "log_copies8" = "log_gene_copies8",
             "log_copies9" = "log_gene_copies9",
             "log_copies10" = "log_gene_copies10") %>%
      pivot_longer(cols = - new_time) %>%
      rename("true_value" = "value") %>%
      filter(true_value > 0)
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_time"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_time, y = value)) + 
      geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      # scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue3"), name="fill") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  } else if (cases == FALSE & ten_sim == FALSE & three_mean == TRUE) {
    # not sure what is going to happen here, wait until we have a posterior to work with before finishing
    true_data <- sim_data %>%
      dplyr::select(new_time, 
                    log_mean_copiesthree
      ) %>%
      rename("log_mean_copies" = "log_mean_copiesthree") %>%
      pivot_longer(cols = - new_time) %>%
      rename("true_value" = "value") %>%
      filter(true_value > 0)
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_time"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_time, y = value)) + 
      geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      # scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue3"), name="fill") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  }
  else {
    true_data <- sim_data %>%
      dplyr::select(new_week, 
                    total_cases) %>%
      rename("true_value" = "total_cases")
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_week"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_week, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_week, y = value)) + 
      geom_point(mapping = aes(x = new_week, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      # scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue3"), name="fill") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  }
  
  return(posterior_predictive_plot)
  
}




# make prior posterior plot -----------------------------------------------

make_fixed_param_plot <- function(posterior_gq_samples, prior_gq_samples, sim_data = TRUE){
  
  if (sim_data == TRUE) {
    priors_and_posteriors <- rbind(posterior_gq_samples %>% dplyr::select(name, value, type, true_value), prior_gq_samples)
    
    
    param_plot <- priors_and_posteriors %>%
      ggplot(aes(value, type, fill = type)) +
      stat_halfeye(normalize = "xy")  +
      geom_vline(aes(xintercept = true_value), linetype = "dotted", size = 1) + 
      facet_wrap(. ~ name, scales = "free_x") +
      theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle("Fixed parameter prior and posteriors")
    
  } else {
    priors_and_posteriors <- rbind(posterior_gq_samples %>% dplyr::select(name, value, type), prior_gq_samples)
    
    
    param_plot <- priors_and_posteriors %>%
      ggplot(aes(value, type, fill = type)) +
      stat_halfeye(normalize = "xy")  +
      facet_wrap(. ~ name, scales = "free_x") +
      theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle("Fixed parameter prior and posteriors")
    
  }
  
  return(param_plot)
  
}

# rt_metrics --------------------------------------------------------------
# used for calculating frequentist characteristics of inference for rt
rt_metrics<- function(data, value, upper, lower) {
  metric_one <- data %>%
    mutate(dev = abs({{ value }} - true_rt),
           CIW = abs({{ upper }} - {{ lower }}),
           envelope = true_rt >= {{ lower }} & true_rt <=  {{ upper }}) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(mean_dev = mean(dev),
              MCIW = mean(CIW),
              mean_env = mean(envelope))
  
  metrics_two <- data %>%
    mutate(prev_val = lag({{ value }}),
           prev_rt = lag(true_rt),
           sv = abs({{ value }} - prev_val),
           rt_sv = abs(true_rt - prev_rt)) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(MASV = mean(sv),
              true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}

# Discretize  Distributions ------------------------------------------
# Epidemia style discretization of gamma
epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, y)
  pmf[1] <- pgamma(1.5, alpha, rate = beta)
  for (i in 2:y) {
    pmf[i] <- pgamma(i + .5, alpha, rate = beta) - pgamma(i - .5, alpha, rate = beta)
  }
  
  pmf
}

zero_epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, (y + 1))
  pmf[1] <- pgamma(0.5, alpha, rate = beta)
  for (i in 2:(y + 1)) {
    pmf[i] <- pgamma(i - 1 + .5, alpha, rate = beta) - pgamma(i - 1 - .5, alpha, rate = beta)
  }
  
  pmf
}


epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, y)
  pmf[1] <- phypoexp(1.5, rates)
  for (i in 2:y) {
    pmf[i] <- phypoexp(i + .5, rates) - phypoexp(i - .5, rates)
  }
  
  pmf
}


epidemia_lognormal <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- plnorm(1.5, meanlog = params[1], sdlog = params[2])
  for (i in 2:y) {
    pmf[i] <- plnorm(i + .5, meanlog = params[1], sdlog = params[2]) -
      plnorm(i - .5, meanlog = params[1], sdlog = params[2])
  }
  
  pmf
}

epidemia_weibull <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- pweibull(1.5, shape = params[1], scale = params[2])
  for (i in 2:y) {
    pmf[i] <- pweibull(i + .5, shape = params[1], scale = params[2]) -
      pweibull(i - .5, shape = params[1], scale = params[2])
  }
  
  pmf
}



zero_epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, (y + 1))
  pmf[1] <- phypoexp(0.5, rates)
  for (i in 2:y + 1) {
    pmf[i] <- phypoexp(i - 1 + .5, rates) - phypoexp(i - 1 - .5, rates)
  }
  
  pmf
}

