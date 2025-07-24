using DrWatson
using Revise
using JLD2
using CSV
using DataFrames
using Turing
using LogExpFunctions
using Random
using ForwardDiff
using Optim
using LineSearches
using Logging
using PreallocationTools
using not_flushed_away


sim =
if length(ARGS) == 0
  75
else
  parse(Int64, ARGS[1])
end

seed = 
if length(ARGS) == 0
  3
else 
  parse(Int64, ARGS[2])
end 

Logging.disable_logging(Logging.Warn)

## Control Parameters
n_chains = 4
priors_only = false

Logging.disable_logging(Logging.Warn)

mkpath(resultsdir("det_ei_normal"))
mkpath(resultsdir("det_ei_normal", "posterior_samples"))
## Control Parameters
n_samples = 250
n_chains = 4

# Load Data and priors
include(projectdir("src/load_dp_det_ei_normal.jl"))

obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)

# pick the change times 
change_grid = 7.0
  if maximum(obstimes) % change_grid == 0
    param_change_max = maximum(obstimes) - change_grid
  else 
    param_change_max = maximum(obstimes)
  end 
param_change_times = collect(change_grid:change_grid:param_change_max)

full_time_series = collect(minimum(obstimes):grid_size:maximum(obstimes))
outs_tmp = dualcache(zeros(3,length(full_time_series)), 10)

index = zeros(length(obstimes))
for i in 1:length(index)
    time = obstimes[i]
    index[i] = indexin(round(Int64,time), full_time_series)[1]
end 

## Define closed form solution
include(projectdir("src/newnew_closed_soln_ei.jl"))


## Load Model
  include(projectdir("src/bayes_det_ei_normal.jl"))



my_model = bayes_det_ei_normal!(
    outs_tmp, 
    data_log_copies,
    obstimes, 
    param_change_times,
    grid_size,
    index)

missing_log_copies = repeat([missing], length(data_log_copies))

my_model_forecast_missing = bayes_det_ei_normal!(
  outs_tmp, 
  missing_log_copies,
  obstimes, 
  param_change_times,
  grid_size,
  index)


if priors_only
    prior_samples = load(resultsdir("det_ei_normal", string("prior_samples_scenario", sim, "_seed", seed, ".jld2")))["prior_samples"]
    
    indices_to_keep = .!isnothing.(generated_quantities(my_model, prior_samples));
    
    prior_samples_randn = ChainsCustomIndex(prior_samples, indices_to_keep);
    
    
    Random.seed!(seed)
    
    
    prior_predictive_randn = predict(my_model_forecast_missing, prior_samples_randn)
    CSV.write(resultsdir("det_ei_normal", string("prior_predictive_scenario", sim, "_seed", seed, ".csv")), DataFrame(prior_predictive_randn))
    
    Random.seed!(seed)
    prior_gq_randn = get_gq_chains(my_model, prior_samples_randn);
    CSV.write(resultsdir("det_ei_normal", string("prior_generated_quantities_scenario", sim, "_seed", seed, ".csv")), DataFrame(prior_gq_randn))
    
    
        exit()
end

posterior_samples = load(resultsdir("det_ei_normal", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")))["posterior_samples"]

indices_to_keep = .!isnothing.(generated_quantities(my_model, posterior_samples));

posterior_samples_randn = ChainsCustomIndex(posterior_samples, indices_to_keep);

Random.seed!(seed)
predictive_randn = predict(my_model_forecast_missing, posterior_samples_randn)
CSV.write(resultsdir("det_ei_normal", "posterior_predictive", string("posterior_predictive", "_scenario", sim, "_seed", seed,  ".csv")), DataFrame(predictive_randn))

Random.seed!(seed)
gq_randn = get_gq_chains(my_model, posterior_samples_randn);
CSV.write(resultsdir("det_ei_normal", "generated_quantities", string("generated_quantities", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(gq_randn))

posterior_df = DataFrame(posterior_samples)
CSV.write(resultsdir("det_ei_normal", "generated_quantities", string("posterior_df", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(posterior_samples))
