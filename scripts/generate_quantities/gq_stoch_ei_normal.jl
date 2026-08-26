# fitting the stochastic ei normal model 
using Logging
using DrWatson
using JLD2
using CSV
using DataFrames
using Random
using Turing
using not_flushed_away
# using GenericLinearAlgebra
sim =
if length(ARGS) == 0
  "E_rt0.75"
else
  parse(Int64, ARGS[1])
end

seed = 
if length(ARGS) == 0
  97
else 
  parse(Int64, ARGS[2])
end 

Logging.disable_logging(Logging.Warn)

## Control Parameters
n_chains = 4
priors_only = false

print(sim)
print(seed)

# Load Data and priors
include(projectdir("src/load_dp_stoch_ei_normal.jl"))

# create times
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

unique_times = unique(obstimes)
index = zeros(length(obstimes))
solve_times = sort(union(unique_times, param_change_times))
l_timepoints = length(solve_times)

for i in 1:length(index)
    time = obstimes[i]
    index[i] = indexin(round(Int64,time), solve_times)[1]
end 


# load model 
include(projectdir("src/closed_form_ei_moments.jl"))
include(projectdir("src/bayes_stoch_ei_normal.jl"))
my_model = bayes_stoch_ei_normal(data_log_copies, obstimes, solve_times, index, param_change_times)


missing_log_copies = repeat([missing], length(data_log_copies))

my_model_forecast_missing = bayes_stoch_ei_normal(missing_log_copies, obstimes, solve_times, index, param_change_times)


if priors_only
    prior_samples = load(resultsdir("stoch_ei_normal", string("prior_samples_scenario", sim, "_seed", seed, ".jld2")))["prior_samples"]
    
    indices_to_keep = .!isnothing.(generated_quantities(my_model, prior_samples));
    
    prior_samples_randn = ChainsCustomIndex(prior_samples, indices_to_keep);
    
    
    Random.seed!(seed)
    
    
    prior_predictive_randn = predict(my_model_forecast_missing, prior_samples_randn)
    CSV.write(resultsdir("stoch_ei_normal", string("prior_predictive_scenario", sim, "_seed", seed, ".csv")), DataFrame(prior_predictive_randn))
    
    Random.seed!(seed)
    prior_gq_randn = get_gq_chains(my_model, prior_samples_randn);
    CSV.write(resultsdir("stoch_ei_normal", string("prior_generated_quantities_scenario", sim, "_seed", seed, ".csv")), DataFrame(prior_gq_randn))
    
    
        exit()
end

posterior_samples = load(resultsdir("stoch_ei_normal", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")))["posterior_samples"]

indices_to_keep = .!isnothing.(generated_quantities(my_model, posterior_samples));

posterior_samples_randn = ChainsCustomIndex(posterior_samples, indices_to_keep);

Random.seed!(seed)
predictive_randn = predict(my_model_forecast_missing, posterior_samples_randn)
CSV.write(resultsdir("stoch_ei_normal", "posterior_predictive", string("posterior_predictive", "_scenario", sim, "_seed", seed,  ".csv")), DataFrame(predictive_randn))

Random.seed!(seed)
gq_randn = get_gq_chains(my_model, posterior_samples_randn);
CSV.write(resultsdir("stoch_ei_normal", "generated_quantities", string("generated_quantities", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(gq_randn))

posterior_df = DataFrame(posterior_samples)
CSV.write(resultsdir("stoch_ei_normal", "generated_quantities", string("posterior_df", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(posterior_samples))
