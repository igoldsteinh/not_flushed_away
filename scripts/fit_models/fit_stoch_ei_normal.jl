# fit the ei model with normal emissions distribution
using Logging
using DrWatson
using JLD2
using CSV
using DataFrames
using Random
using Turing
using not_flushed_away
sim =
if length(ARGS) == 0
  "Gcombined_trunc"
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
# generate prior samples 
if priors_only
  Random.seed!(seed)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
  wsave(resultsdir("stoch_ei_normal", string("prior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict prior_samples)

  exit()
end
# create initial values using MAP plus noise
Random.seed!(seed)
MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]
Random.seed!(seed)
MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]
init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
# lets try using a sample from the deterministic posterior as the starting point
# if we use zeros, then we are initializing where the deterministic posterior is at
if sim == "G1" 
  det_posterior_samples = load(resultsdir("det_ei_normal", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")))["posterior_samples"]
  det_data_frame = DataFrame(det_posterior_samples)
  median_vals = vcat(zeros(length(solve_times) * 2), [median(det_data_frame[:,i]) for i in 3:28])
  init = repeat([median_vals], n_chains) 
end 
# fit the model 
Turing.setadbackend(:forwarddiff)
Random.seed!(seed)
n_samples = 250
posterior_samples = sample(my_model, NUTS(-1, 0.9), MCMCThreads(), n_samples, n_chains, discard_initial = n_samples, init_params = init)
wsave(resultsdir("stoch_ei_normal", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict posterior_samples)
