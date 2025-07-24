# simulating epidemic curves from the stochastic concentration model
using Logging
using DrWatson
using JLD2
using CSV
using DataFrames
using Random
using Turing
using DifferentialEquations
using stoch_wastewater
# using GenericLinearAlgebra
init_count =
if length(ARGS) == 0
  30.0
else
  parse(Int64, ARGS[1])
end

seed = 
if length(ARGS) == 0
  4
else 
  parse(Int64, ARGS[2])
end 

Logging.disable_logging(Logging.Warn)

## Control Parameters
n_chains = 4
priors_only = false
tau = 0.5
rt = 1.1

sim = "gillEIconst"
# Load Data and priors
include(projectdir("src/load_dp_stoch_conc_ei.jl"))


obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)

# pick the change times 
# if maximum(obstimes) % 7 == 0
#   param_change_max = maximum(obstimes) - 7
# else 
#   param_change_max = maximum(obstimes)
# end 
# param_change_times = collect(7:7.0:param_change_max)
# full_time_series = collect(minimum(obstimes):grid_size:maximum(obstimes))

# pick the change times as the times used to simulate gillespie data
rt_times = CSV.read("data/gillEIconst_scenario1_pop1000_full_genecount_obsdata.csv", DataFrame)

rt_times = subset(rt_times, :seed => ByRow(x -> x == seed))
param_change_times = rt_times[:,:time][1:end-1]
# lets use 1.2 as the rt 
# repeat 1.2 for the length of param_change_times
# repeat rt for the length of param_change_times

alpha_t_values_no_init = fill(rt, length(param_change_times)) * (1/7)
alpha_init = rt * (1/7)

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
include(projectdir("src/bayes_stoch_ei_simulator.jl"))

E_init = init_count
I_init = init_count
my_model = bayes_stoch_ei_simulator(data_log_copies, obstimes, solve_times, index, param_change_times, alpha_init, alpha_t_values_no_init, E_init, I_init, tau) 

# generate prior samples 
Random.seed!(seed)
prior_samples = sample(my_model, Prior(), MCMCThreads(), 400, 1)
 

indices_to_keep = .!isnothing.(generated_quantities(my_model, prior_samples));
    
prior_samples_randn = ChainsCustomIndex(prior_samples, indices_to_keep);


missing_log_copies = repeat([missing], length(data_log_copies))

my_model_forecast_missing = bayes_stoch_ei_simulator(missing_log_copies, obstimes, solve_times, index, param_change_times, alpha_init, alpha_t_values_no_init, E_init, I_init, tau)

Random.seed!(seed)


prior_predictive_randn = predict(my_model_forecast_missing, prior_samples_randn)
CSV.write(projectdir("data", string("ln_pp_pop", init_count, "_tau", tau, "_rt", rt, "_seed", seed, ".csv")), DataFrame(prior_predictive_randn))

Random.seed!(seed)
prior_gq_randn = get_gq_chains(my_model, prior_samples_randn);
CSV.write(projectdir("data", string("ln_params_pop", init_count, "_tau", tau, "_rt", rt, "_seed", seed, ".csv")), DataFrame(prior_gq_randn))

