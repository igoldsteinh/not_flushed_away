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
  "E_rt0.75"
else
  parse(Int64, ARGS[1])
end
seed = 
if length(ARGS) == 0
  1
else 
  parse(Int64, ARGS[2])
end 
Logging.disable_logging(Logging.Warn)
## Control Parameters
n_samples = 250
n_chains = 4
priors_only = false
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
# sample prior
if priors_only
  Random.seed!(seed)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
  wsave(resultsdir("det_ei_normal",string("prior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict prior_samples)
  exit()
end
Random.seed!(seed)
# initialize at MAP plus noise
MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]
Random.seed!(seed)
MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]
init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
# sample posterior
Random.seed!(seed)
posterior_samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_samples, init_params = init)
wsave(resultsdir("det_ei_normal", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict posterior_samples)


