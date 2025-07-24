# calculate mean gen time and sd of gen time 
# also find a gamma distribution to match the SLD
library(sdprisk)


# mean and sd of gen time -------------------------------------------------
set.seed(1)
latent_period = rexp(10000, rate = 1/4)
inf_period = rgamma(10000, shape = 7, rate = 1)
gen_time = latent_period + inf_period
mean(gen_time)
sd(gen_time)
# find SLD approximate ----------------------------------------------------
dumb_weights <- c(0, 0.05, 0.4, 0.8, 0.4, 0.2, 0.1, 0.05, 0.025, 0.0125)
norm_weights <- dumb_weights / sum(dumb_weights)
comp_weights = log(norm_weights * 10e6, base = 10)
times <- c(0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 11.5, 20.5)

compare_gamma_pdf <- function(candidate_params, true_weights, true_times) {
  candidate_vals <- dgamma(true_times, shape = candidate_params[1], rate = candidate_params[2])
  
  loss <- sum((candidate_vals - true_weights)^2)
  return(loss)
}

start_shape = 3
start_scale = 1
start_params <- c(start_shape, start_scale)
optim_params <- optim(par = start_params,
                        fn = compare_gamma_pdf, 
                        true_weights = norm_weights, 
                        true_times = times)

optim_shape = optim_params$par[1]
optim_scale = optim_params$par[2]
optim_vals = dgamma(times, shape = optim_shape, rate = optim_scale)
compare_frame = data.table(times, norm_weights, optim_vals)
compare_frame %>%
  pivot_longer(cols = -times, names_to = "type", values_to = "value") %>%
  ggplot(aes(x = times, y = value, color = type)) + 
  geom_point()
sld_mean = optim_shape / optim_scale
sld_sd = sqrt(optim_shape/optim_scale^2)
