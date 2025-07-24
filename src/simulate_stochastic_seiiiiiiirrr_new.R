# simulate stochastic seiiiiiiirr 
# individual level
library(fields)
library(tidyverse)

# stochastic seiiiirr with non-constant beta--------------------------------------------------------
# beta_vec[1] is the initial beta
# change_points[1] is 0
sim_SEIIIIIIIRRR_nonconst_new <- function(N, 
                                  E_init, 
                                  I1_init, 
                                  I2_init, 
                                  I3_init, 
                                  I4_init,
                                  I5_init,
                                  I6_init,
                                  I7_init,
                                  R1_init,
                                  R2_init,
                                  beta_vec, 
                                  change_points, 
                                  gamma, 
                                  nu1,
                                  nu2,
                                  nu3,
                                  nu4,
                                  nu5, 
                                  nu6,
                                  nu7,
                                  eta1, 
                                  eta2,
                                  t_stop) {
  # initial number of comp counts 
  E <- E_init
  I1 <- I1_init
  I2 <- I2_init
  I3 <- I3_init
  I4 <- I4_init
  I5 <- I5_init
  I6 <- I6_init
  I7 <- I7_init
  R1 <- R1_init
  R2 <- R2_init
  S <- N - E_init -
           I1_init - 
           I2_init - 
           I3_init - 
           I4_init -
           I5_init - 
           I6_init - 
           I7_init - 
           R1_init - 
           R2_init
  # set up the rate frame
  S_vec <- rep("S", S)
  E_vec <- rep("E", E)
  I1_vec <- rep("I1", I1)
  I2_vec <- rep("I2", I2)
  I3_vec <- rep("I3", I3)
  I4_vec <- rep("I4", I4)
  I5_vec <- rep("I5", I5)
  I6_vec <- rep("I6", I6)
  I7_vec <- rep("I7", I7)
  R1_vec <- rep("R1", R1)
  R2_vec <- rep("R2", R2)
  pop_vec <- c(S_vec, E_vec, I1_vec, I2_vec, 
               I3_vec, I4_vec, I5_vec, I6_vec, I7_vec, R1_vec, R2_vec)
  rate_vec = rep(0, N)
  id_vec = seq(1, N, by = 1)
  rate_frame = data.frame(id = id_vec, state = pop_vec, rate = rate_vec)
  
  # set up the individual frame output
  individ_frame = data.frame(id = id_vec) %>%
    mutate(infec_time = NA,
           infectious_time = NA,
           i2_time = NA,
           i3_time = NA,
           i4_time = NA,
           i5_time = NA,
           i6_time = NA,
           i7_time = NA,
           recover_time = NA,
           r2_time = NA,
           stopshed_time = NA)
  # set the initial times to be 0
  individ_frame$infec_time[rate_frame$id[rate_frame$state == "E"]] <- 0
  individ_frame$infectious_time[rate_frame$id[rate_frame$state == "I1"]] <- 0
  individ_frame$i2_time[rate_frame$id[rate_frame$state == "I2"]] <- 0
  individ_frame$i3_time[rate_frame$id[rate_frame$state == "I3"]] <- 0
  individ_frame$i4_time[rate_frame$id[rate_frame$state == "I4"]] <- 0
  individ_frame$i5_time[rate_frame$id[rate_frame$state == "I5"]] <- 0
  individ_frame$i6_time[rate_frame$id[rate_frame$state == "I6"]] <- 0
  individ_frame$i7_time[rate_frame$id[rate_frame$state == "I7"]] <- 0
  individ_frame$recover_time[rate_frame$id[rate_frame$state == "R1"]] <- 0
  individ_frame$r2_time[rate_frame$id[rate_frame$state == "R2"]] <- 0
  
  # set up initial state frame
  beta_idx = 1
  beta_t = beta_vec[beta_idx]
  next_cp = if (length(beta_vec) > 1 ) {
    change_points[beta_idx + 1]
  } else {
    Inf
  }
  r0 = (beta_t/nu1) + (beta_t/nu2) + (beta_t/nu3) + (beta_t/nu4) + (beta_t/nu7) + (beta_t/nu5) + (beta_t/nu6)
  rt = (S/N) * r0 
  t <- 0
  initial_state = data.frame(time = t,
                             S = S,
                             E = E,
                             I1 = I1, 
                             I2 = I2, 
                             I3 = I3, 
                             I4 = I4, 
                             I5 = I5, 
                             I6 = I6, 
                             I7 = I7, 
                             R1 = R1, 
                             R2 = R2, 
                             beta_t = beta_t, 
                             r0 = r0, 
                             rt = rt)
  state_frame = initial_state
  # update the rate_frame
  rate_frame$rate <- (rate_frame$state == "S") * ((beta_t/N) * (I1 + I2 + I3 + I4 + I5 + I6 + I7)) +
    (rate_frame$state == "E") * gamma + 
    (rate_frame$state == "I1") * nu1 + 
    (rate_frame$state == "I2") * nu2 + 
    (rate_frame$state == "I3") * nu3 + 
    (rate_frame$state == "I4") * nu4 + 
    (rate_frame$state == "I5") * nu5 +
    (rate_frame$state == "I6") * nu6 +
    (rate_frame$state == "I7") * nu7 + 
    (rate_frame$state == "R1") * eta1 + 
    (rate_frame$state == "R2") * eta2 
  
  
  # start simulating
  while ((E + I1 + I2 + I3 + I4 + I5 + I6 + I7) > 0 & t < t_stop) {
    # time to next event is min of all possible events
    total_rate = sum(rate_frame$rate)
    next_event = rexp(1, rate = total_rate)
    
    # check to see if the next event is a beta change
    if (t + next_event >= next_cp){
      # update time
      t = next_cp
      # update beta and next_cp
      beta_idx = beta_idx + 1
      beta_t = beta_vec[beta_idx]
      if (beta_idx < length(beta_vec)) {
        next_cp = change_points[beta_idx + 1]
      } else {
        next_cp = Inf  
      }
    } else {
      # choose which event happens propoportional to the rates
      which_id = sample(rate_frame$id, 1, prob = rate_frame$rate)
      # update time 
      t = t + next_event
      if (rate_frame$state[which_id] == "S") {
        rate_frame$state[which_id] = "E"
        E = E + 1
        S = S - 1
        individ_frame$infec_time[which_id] = t
      } else if (rate_frame$state[which_id] == "E") {
        rate_frame$state[which_id] = "I1"
        I1 = I1 + 1
        E = E - 1
        individ_frame$infectious_time[which_id] = t
      } else if (rate_frame$state[which_id] == "I1") {
        rate_frame$state[which_id] = "I2"
        I2 = I2 + 1
        I1 = I1 - 1
        individ_frame$i2_time[which_id] = t
      } else if (rate_frame$state[which_id] == "I2") {
        rate_frame$state[which_id] = "I3"
        I3 = I3 + 1
        I2 = I2 - 1
        individ_frame$i3_time[which_id] = t
      } else if (rate_frame$state[which_id] == "I3") {
        rate_frame$state[which_id] = "I4"
        I4 = I4 + 1
        I3 = I3 - 1
        individ_frame$i4_time[which_id] = t
      } else if (rate_frame$state[which_id] == "I4") {
        rate_frame$state[which_id] = "I5"
        I5 = I5 + 1
        I4 = I4 - 1
        individ_frame$i5_time[which_id] = t
      } else if (rate_frame$state[which_id] == "I5") {
        rate_frame$state[which_id] = "I6"
        I6 = I6 + 1
        I5 = I5 - 1
        individ_frame$i6_time[which_id] = t
      } else if (rate_frame$state[which_id] == "I6") {
        rate_frame$state[which_id] = "I7"
        I7 = I7 + 1
        I6 = I6 - 1
        individ_frame$i7_time[which_id] = t
      } else if (rate_frame$state[which_id] == "I7") {
        rate_frame$state[which_id] = "R1"
        R1 = R1 + 1
        I7 = I7 - 1
        individ_frame$recover_time[which_id] = t
      } else if (rate_frame$state[which_id] == "R1") {
        rate_frame$state[which_id] = "R2"
        R2 = R2 + 1
        R1 = R1 - 1
        individ_frame$r2_time[which_id] = t
      } else {
        rate_frame$state[which_id] = "R3"
        R2 = R2 - 1
        individ_frame$stopshed_time[which_id] = t
      }
    }
    # update rates
    rate_frame$rate <- (rate_frame$state == "S") * ((beta_t/N) * (I1 + I2 + I3 + I4 + I5 + I6 + I7)) +
      (rate_frame$state == "E") * gamma + 
      (rate_frame$state == "I1") * nu1 + 
      (rate_frame$state == "I2") * nu2 + 
      (rate_frame$state == "I3") * nu3 + 
      (rate_frame$state == "I4") * nu4 + 
      (rate_frame$state == "I5") * nu5 +
      (rate_frame$state == "I6") * nu6 +
      (rate_frame$state == "I7") * nu7 + 
      (rate_frame$state == "R1") * eta1 + 
      (rate_frame$state == "R2") * eta2 
    # update params
    r0 = (beta_t/nu1) + (beta_t/nu2) + (beta_t/nu3) + (beta_t/nu4) + (beta_t/nu7) + (beta_t/nu5) + (beta_t/nu6)
    rt = (S/N) * r0 
    new_state = data.frame(time = t,
                               S = S,
                               E = E,
                               I1 = I1, 
                               I2 = I2, 
                               I3 = I3, 
                               I4 = I4, 
                               I5 = I5, 
                               I6 = I6, 
                               I7 = I7, 
                               R1 = R1, 
                               R2 = R2, 
                               beta_t = beta_t, 
                               r0 = r0, 
                               rt = rt)
    state_frame = bind_rows(state_frame, new_state)
  }
  res <- list(individ_frame, state_frame)
  res
}

# function for creating day level state data seiiiirr ------------------------------
create_daily_data_seiiiiiiirrr_new <- function(state_frame) {
  day_data <- state_frame %>% 
    mutate(integer_day = ceiling(time),
           time_diff = integer_day - time) %>%
    group_by(integer_day) %>% 
    filter(time_diff == min(time_diff)) 
  
  return(day_data)
}

# calc_individ_counts_nonoise ---------------------------------------------
# function for calculating individual gene counts at time t 
# with no individual variation
x <- c(-3, -1, 1, 3, 5, 7, 9, 13, 17, 21, 25, 29)
y <- c(5, 6.9, 6.7, 6.5, 6.2, 5.9, 5.5, 4.75, 3.9, 3.3, 2.1, 1.3)
x_adj <- x + 3

#translate from log base 10 scale to real scale
exp_y <- 10^(y)
#tps cant do this for some reason, so instead do the spline on log base 10
data <- data.frame(x_adj, exp_y, y)
# tp_spline <- Tps(data$x_adj, data$y)
set.seed(1234)
tp_spline <- Tps(data$x_adj, data$y)

calc_individ_counts_nonoise <- function(t) {
  mean_log10_prediction<- predict(tp_spline, t)
  
  # hard cutoff at zero
  if (mean_log10_prediction < 0){
    prediction <- 0
  } else {
    prediction <- 10^(mean_log10_prediction) # should be 1.09
  }
  
  return(prediction)
}


# calculate total gene concentration with no noise and seiiiirr -----------
calc_total_gene_concs_nonoise_seiiiiiiirrr <- function(epi_curve, time, N) {
  # first collect all active individuals
  # epi_curve = wide_format
  # time = 5
  active_individuals <- epi_curve %>% 
    ungroup() %>% 
    mutate(i1_individuals = infectious_time <= time & i2_time > time,
           i2_individuals = i2_time <= time & i3_time > time,
           i3_individuals = i3_time <= time & i4_time > time,
           i4_individuals = i4_time <= time & i5_time > time,
           i5_individuals = i5_time <= time & i6_time > time,
           i6_individuals = i6_time <= time & i7_time > time,
           i7_individuals = i7_time <= time & recover_time > time,
           r1_individuals = recover_time <= time & r2_time > time,
           r2_individuals = r2_time <= time & stopshed_time > time) %>%
    filter(i1_individuals == TRUE | i2_individuals == TRUE | 
             i3_individuals == TRUE | i4_individuals == TRUE |
             i5_individuals == TRUE | i6_individuals == TRUE |
             i7_individuals == TRUE |
             r1_individuals == TRUE | r2_individuals == TRUE)
  
  if (dim(active_individuals)[1] > 0) {
    active_individuals <- active_individuals %>%
      mutate(time_since_infectious = time - infectious_time) %>% 
      rowwise() %>% 
      mutate(gene_counts = calc_individ_counts_nonoise(time_since_infectious))
    
    # report the total gene counts, and the total number in each compartment
    res <- data.frame("time" = time,
                      "total_conc" = sum(active_individuals$gene_counts)/N, 
                      "total_conc_I1" = sum(active_individuals$gene_counts[active_individuals$i1_individuals == TRUE])/N,
                      "total_conc_I2" = sum(active_individuals$gene_counts[active_individuals$i2_individuals == TRUE])/N,
                      "total_conc_I3" = sum(active_individuals$gene_counts[active_individuals$i3_individuals == TRUE])/N,
                      "total_conc_I4" = sum(active_individuals$gene_counts[active_individuals$i4_individuals == TRUE])/N,
                      "total_conc_I5" = sum(active_individuals$gene_counts[active_individuals$i5_individuals == TRUE])/N,
                      "total_conc_I6" = sum(active_individuals$gene_counts[active_individuals$i6_individuals == TRUE])/N,
                      "total_conc_I7" = sum(active_individuals$gene_counts[active_individuals$i7_individuals == TRUE])/N,
                      "total_conc_R1" = sum(active_individuals$gene_counts[active_individuals$r1_individuals == TRUE])/N,
                      "total_conc_R2" = sum(active_individuals$gene_counts[active_individuals$r2_individuals == TRUE])/N,
                      "num_I1" = sum(active_individuals$i1_individuals),
                      "num_I2" = sum(active_individuals$i2_individuals),
                      "num_I3" = sum(active_individuals$i3_individuals),
                      "num_I4" = sum(active_individuals$i4_individuals),
                      "num_I5" = sum(active_individuals$i5_individuals),
                      "num_I6" = sum(active_individuals$i6_individuals),
                      "num_I7" = sum(active_individuals$i7_individuals),
                      "num_R1" = sum(active_individuals$r1_individuals),
                      "num_R2" = sum(active_individuals$r2_individuals))
    
    return(res)
  } else {
    res <- data.frame("time" = time,
                      "total_conc" = 0, 
                      "total_conc_I1" = 0,
                      "total_conc_I2" = 0 ,
                      "total_conc_I3" = 0, 
                      "total_conc_I4" = 0, 
                      "total_conc_I5" = 0, 
                      "total_conc_I6" = 0, 
                      "total_conc_I7" = 0, 
                      "total_conc_R1" = 0,
                      "total_conc_R2" = 0,
                      "num_I1" = 0,
                      "num_I2" = 0,
                      "num_I3" = 0,
                      "num_I4" = 0,
                      "num_I5" = 0,
                      "num_I6" = 0,
                      "num_I7" = 0,
                      "num_R1" = 0,
                      "num_R2" = 0)
    
    return(res)
    
  }
}


# create new shedding function based on infection stages and indiv --------
new_shedding <- function(time, individ_data, comp_weights, N) {
    # individ_data = wide_format[[1]]
    # time = 5
    individ_data <- individ_data %>% 
      ungroup() %>% 
      mutate(i1_individuals = infectious_time <= time & i2_time > time,
             i2_individuals = i2_time <= time & i3_time > time,
             i3_individuals = i3_time <= time & i4_time > time,
             i4_individuals = i4_time <= time & i5_time > time,
             i5_individuals = i5_time <= time & i6_time > time,
             i6_individuals = i6_time <= time & i7_time > time,
             i7_individuals = i7_time <= time & recover_time > time,
             r1_individuals = recover_time <= time & r2_time > time,
             r2_individuals = r2_time <= time & stopshed_time > time) %>%
      mutate(shedding = ifelse(i1_individuals == TRUE, 10^(comp_weights[1] + individ_weight),
                              ifelse(i2_individuals == TRUE, 10^(comp_weights[2] + individ_weight),
                              ifelse(i3_individuals == TRUE, 10^(comp_weights[3] + individ_weight),
                              ifelse(i4_individuals == TRUE, 10^(comp_weights[4] + individ_weight),
                              ifelse(i5_individuals == TRUE, 10^(comp_weights[5] + individ_weight),
                              ifelse(i6_individuals == TRUE, 10^(comp_weights[6] + individ_weight),
                              ifelse(i7_individuals == TRUE, 10^(comp_weights[7] + individ_weight),
                              ifelse(r1_individuals == TRUE, 10^(comp_weights[8] + individ_weight),
                              ifelse(r2_individuals == TRUE, 10^(comp_weights[9] + individ_weight), 0))))))))))
    
      # report the total gene counts, and the total number in each compartment
      res <- data.frame("time" = time,
                        "total_conc" = sum(individ_data$shedding)/N, 
                        "total_conc_I1" = sum(individ_data$shedding[individ_data$i1_individuals == TRUE])/N,
                        "total_conc_I2" = sum(individ_data$shedding[individ_data$i2_individuals == TRUE])/N,
                        "total_conc_I3" = sum(individ_data$shedding[individ_data$i3_individuals == TRUE])/N,
                        "total_conc_I4" = sum(individ_data$shedding[individ_data$i4_individuals == TRUE])/N,
                        "total_conc_I5" = sum(individ_data$shedding[individ_data$i5_individuals == TRUE])/N,
                        "total_conc_I6" = sum(individ_data$shedding[individ_data$i6_individuals == TRUE])/N,
                        "total_conc_I7" = sum(individ_data$shedding[individ_data$i7_individuals == TRUE])/N,
                        "total_conc_R1" = sum(individ_data$shedding[individ_data$r1_individuals == TRUE])/N,
                        "total_conc_R2" = sum(individ_data$shedding[individ_data$r2_individuals == TRUE])/N,
                        "num_I1" = sum(individ_data$i1_individuals),
                        "num_I2" = sum(individ_data$i2_individuals),
                        "num_I3" = sum(individ_data$i3_individuals),
                        "num_I4" = sum(individ_data$i4_individuals),
                        "num_I5" = sum(individ_data$i5_individuals),
                        "num_I6" = sum(individ_data$i6_individuals),
                        "num_I7" = sum(individ_data$i7_individuals),
                        "num_R1" = sum(individ_data$r1_individuals),
                        "num_R2" = sum(individ_data$r2_individuals))
      
      return(res)
}


# simulate observed concentrations from true ------------------------------

simulate_gene_data_normal <- function(true_gene_counts, seed, rho, sd){
  set.seed(seed)
  sim_data <- true_gene_counts %>% 
    mutate(log_genes_mean = log(total_conc) + log(rho)) %>%
    rowwise() %>%
    mutate(
      log_gene_copies1 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies2 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies3 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies4 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies5 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies6 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies7 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies8 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies9 = log_genes_mean  + (sd * rnorm(1)),
      log_gene_copies10 = log_genes_mean + (sd * rnorm(1))) %>%
    mutate(log_mean_copiesten = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3) +
                                       exp(log_gene_copies4) + exp(log_gene_copies5) + exp(log_gene_copies6) +
                                       exp(log_gene_copies7) + exp(log_gene_copies8) + exp(log_gene_copies9) +
                                       exp(log_gene_copies10))/10),
           log_mean_copiesthree = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3))/3))
  
  return(sim_data)
  
}

