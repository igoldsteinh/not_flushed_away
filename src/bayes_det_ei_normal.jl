# Turing model for EIRR-ww model (flexible resolution)
using LinearAlgebra
using ForwardDiff
using NaNMath
@model function bayes_det_ei_normal!(outs_tmp,data_log_copies, obstimes, param_change_times, grid_size, index)
  # Calculate number of observed datapoints timepoints
  l_copies = length(obstimes)
  l_param_change_times = length(param_change_times)
  # Priors
  rt_params_non_centered ~ MvNormal(zeros(l_param_change_times + 2), Diagonal(ones(l_param_change_times + 2))) # +2, 1 for var, 1 for init
  I_init_non_centered ~ Normal()
  E_init_non_centered ~ Normal()
  gamma_non_centered ~ Normal() # rate to I
  nu_non_centered ~ Normal() # rate to R1
  rho_conc_non_centered ~ Normal() # scale concentrations
  tau_non_centered ~ Normal() # standard deviation for log scale data
  # Transformations
  gamma = exp(gamma_non_centered * gamma_sd + gamma_mean)
  nu = exp(nu_non_centered * nu_sd + nu_mean)
  rho_conc = exp(rho_conc_non_centered * rho_conc_sd + rho_conc_mean)
  tau = exp(tau_non_centered * tau_sd + tau_mean)
  sigma_rt_non_centered = rt_params_non_centered[1]
  sigma_rt = exp(sigma_rt_non_centered * sigma_rt_sd + sigma_rt_mean)
  rt_init_non_centered = rt_params_non_centered[2]
  rt_init = exp(rt_init_non_centered * rt_init_sd + rt_init_mean)
  alpha_init = rt_init * nu
  log_rt_steps_non_centered = rt_params_non_centered[3:end]
  I_init = exp(I_init_non_centered * I_init_sd + I_init_mean)
  E_init = exp(E_init_non_centered * E_init_sd + E_init_mean)
  u0 = [E_init, I_init] 
  # Time-varying parameters
  alpha_t_values_no_init = exp.(log(rt_init) .+ cumsum(vec(log_rt_steps_non_centered) * sigma_rt)) * nu
  alpha_t_values_with_init = vcat(alpha_init, alpha_t_values_no_init)
  # solve ODE
  sol_reg_scale_array = newnew_ei_closed_solution!(outs_tmp, 1:maximum(obstimes), param_change_times, grid_size, 0.0, alpha_t_values_with_init, u0, gamma, nu)
  log_genes_mean = NaNMath.log.(sol_reg_scale_array[3, 2:end]) .+ log(rho_conc) # first entry is the initial conditions, we want 2:end
  # Likelihood
  for i in 1:l_copies
    data_log_copies[i] ~ Normal(log_genes_mean[round(Int64,index[i])], tau) 
  end
  # Generated quantities
  rt_t_values_with_init = alpha_t_values_with_init/ nu
  return (
    gamma = gamma,
    nu = nu,
    rho_conc = rho_conc,
    rt_init = rt_init,
    sigma_rt = sigma_rt,
    tau = tau,
    alpha_t_values = alpha_t_values_with_init,
    rt_t_values = rt_t_values_with_init,
    E_init,
    I_init,
    E = sol_reg_scale_array[2, :],
    I = sol_reg_scale_array[3, :],
    log_genes_mean = log_genes_mean
  )
end
