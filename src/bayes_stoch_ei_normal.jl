# Turing model for stochastic EI-ww model
using LinearAlgebra
using ForwardDiff
using LogExpFunctions
# function for creating compartment counts
# using our log normal approximation 
function solve_log_compartments(__varinfo__,
                                __context__,
                                log_comp_counts_non_centered, 
                                l_timepoints,               
                                solve_times, 
                                param_change_times, 
                                l_param_change_times,
                                u0, 
                                alpha_t_values_with_init, 
                                gamma, 
                                nu)
  log_comp_counts = zeros(Base.promote_eltype(log_comp_counts_non_centered[1]), 2, l_timepoints)
  alpha = alpha_t_values_with_init[1]
  alpha_index = 1
  next_change = param_change_times[1]
  change_index = 1
  # assume the first the initial condition is time 0
  for i in 1:l_timepoints 
    # calculate the time which has ellapsed from the current state to the next state
    # assume we initialize at time 0
    if i == 1
      t = solve_times[i]
    else 
      t = solve_times[i] - solve_times[i-1]
    end
    # solve for the conditional means at time solve_times[i] given the compartment counts at time solve_times[i-1]
    # before solving closed_ei_moments check if exp(sqrt(4*alpha*gamma + power(gamma - nu,2))*t) is infinite, and if so reject
    # this equivalent to checking if the sum of the moments is over 8 billion, but catches it at an earlier stage before it can turn into a NAN ouput
    if isinf(exp(sqrt(4*alpha*gamma + power(gamma - nu,2))*t)) == true || isinf(6*exp((gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) + nu)*t)) == true
        Turing.@addlogprob! -Inf
        return nothing
    end
    moment_array = closed_ei_moments(t,u0, alpha, gamma, nu)
    # reject if I compartment mean is less than 1e-10 or the sum of the means is larger than 8 billion 
    if moment_array[2] < 1e-10 || sum(moment_array[1:2]) > 8e9 
        Turing.@addlogprob! -Inf
        return nothing 
    end
    # calculate the covariance matrix of the log compartments
    # using the delta method approximation
     log_var_matrix = reshape(vcat(moment_array[3]/(moment_array[1])^2, #varE/muE^2
     moment_array[5]/((moment_array[1])*(moment_array[2])), #covEI/muE * muI
     moment_array[5]/((moment_array[1])*(moment_array[2])), #covEI/muE * muI
     moment_array[4]/(moment_array[2])^2), (2,2)) #varI/muI^2
      # cholesky decomposition of covariance matrix of log compartments with noise on the diagonal to help pos def
     final_matrix = cholesky(Hermitian(log_var_matrix + Diagonal(fill(1e-12, 2))))
     #create compartment counts at time solve_times[i]
     log_comp_counts[1:2, i]= log.(moment_array[1:2]) + final_matrix.L * log_comp_counts_non_centered[(1+(i-1)*2):i*2]
     # update u0, the initial conditions for time solve_times[i+1]
     # since the compartments are conditioned on, there is no variance in the initial conditions
     u0 = vcat(exp.(log_comp_counts[1:2,i]), 0,0,0)
     # check if this was a time when we need to update alpha
     # if it is, update alpha and update the next time we change 
     if solve_times[i] == next_change
      # if we are at the last time point we are solving for, do nothing
      if i == l_timepoints
        # do nothing
      else
        # update alpha
        alpha_index = alpha_index + 1
        alpha = alpha_t_values_with_init[alpha_index]
        # update the next time we change alpha
        # if we're at the last change, set the next change to be a time that is larger than the last time point
        if change_index == l_param_change_times
          next_change = maximum(solve_times) + 1
        else 
          # update the next change time 
          change_index = change_index + 1
          next_change = param_change_times[change_index]
        end
      end 
    end
    # reject if sum of first three u0 entries is more than 8 billion
    if sum(u0[1:2]) > 8e9
        Turing.@addlogprob! -Inf
        return  
    end
  end 
  return(log_comp_counts)
end
# the full model
@model function bayes_stoch_ei_normal(data_log_copies, obstimes, solve_times, index, param_change_times)
  # Calculate number of observed datapoints timepoints
  l_copies = length(obstimes)
  l_timepoints = length(solve_times)
  l_param_change_times = length(param_change_times)
  log_comp_counts_non_centered ~ MvNormal(zeros((l_timepoints * 2)), Diagonal(ones((l_timepoints * 2)))) 
  rt_params_non_centered ~  MvNormal(zeros(l_param_change_times + 2), Diagonal(ones(l_param_change_times + 2))) # +2, 1 for var, 1 for init
  gamma_non_centered ~ Normal() # rate to I
  nu_non_centered ~ Normal() # rate to R
  I_init_non_centered ~ Normal()
  E_init_non_centered ~ Normal()
  rho_conc_non_centered ~ Normal() # scaling factor for gene concs
  tau_non_centered ~ Normal() # noise parameter
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
  u0 = [E_init, I_init, 0,0,0] 
  alpha_t_values_no_init = exp.(log(rt_init) .+ cumsum(vec(log_rt_steps_non_centered) * sigma_rt)) * nu
  alpha_t_values_with_init = vcat(alpha_init, alpha_t_values_no_init)
  # Solve for the mean and variance
  # for each time point, we need to model the conditional density of the compartments given the previous compartments
  log_comp_counts = solve_log_compartments(__varinfo__,
                                            __context__,
                                          log_comp_counts_non_centered, 
                                          l_timepoints,
                                          solve_times,
                                          param_change_times,
                                          l_param_change_times,
                                          u0,
                                          alpha_t_values_with_init,
                                          gamma,
                                          nu)  
  if log_comp_counts == nothing
    Turing.@addlogprob! -Inf
    return
  end
  # create mean of log genes  
  log_genes_mean = log_comp_counts[2, :] .+ log(rho_conc)
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
    I_init,
    E_init,
    log_genes_mean = vcat(log(I_init) + log(rho_conc), log_genes_mean),
    E = vcat(E_init,exp.(log_comp_counts[1, :])),
    I = vcat(I_init, exp.(log_comp_counts[2, :]))
  )
end
