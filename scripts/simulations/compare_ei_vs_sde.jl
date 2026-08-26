# this file is for comparing the transition densities of the EI MJP vs the SDE used in the model
using DrWatson
using Distributions
using Random
using DataFrames
using StatsBase
using LinearAlgebra
using CSV
include(srcdir("sim_ei_mjp.jl"))
include(srcdir("sim_ei_sde.jl"))
include(srcdir("closed_form_ei_moments.jl"))
# let's start with 10 in each
init_E = 10
init_I = 10
gamma = 1/4
nu = 1/7
rt = 1.3
alpha = rt * nu
t_end = 3.0
### mjp
mjp_frame10 = Array{Float64}(undef, 10000, 3) 
for s in 1:10000
    Random.seed!(s)
    sim_frame = sim_ei_mjp(init_E, init_I, alpha, nu, gamma, t_end)
    last_row = sim_frame[sim_frame.time .<= t_end,:][end,:]
    mjp_frame10[s,1] = s
    mjp_frame10[s,2] = last_row.E
    mjp_frame10[s,3] = last_row.I
end 
column_names = [:iteration, :E, :I]
# Convert the array to a DataFrame
mjp_df10 = DataFrame(mjp_frame10, column_names)
### sde 
sde_frame10 = Array{Float64}(undef, 10000, 3) 
t_diff = t_end - 0.0
for s in 1:10000
    Random.seed!(s)
    comp_counts = sim_ei_sde(init_E, init_I, alpha, nu, gamma, t_diff)
    sde_frame10[s,1] = s
    sde_frame10[s,2] = comp_counts[1]
    sde_frame10[s,3] = comp_counts[2]
end 
sde_df10 = DataFrame(sde_frame10, column_names)
CSV.write(datadir("compare_data", "ei_mjp10.csv"), mjp_df10)
CSV.write(datadir("compare_data", "ei_sde10.csv"), sde_df10)

# let's start with 5 in each
init_E = 5
init_I = 5
gamma = 1/4
nu = 1/7
rt = 1.3
alpha = rt * nu
t_end = 3.0
### mjp
mjp_frame5 = Array{Float64}(undef, 10000, 3) 
for s in 1:10000
    Random.seed!(s)
    sim_frame = sim_ei_mjp(init_E, init_I, alpha, nu, gamma, t_end)
    last_row = sim_frame[sim_frame.time .<= t_end,:][end,:]
    mjp_frame5[s,1] = s
    mjp_frame5[s,2] = last_row.E
    mjp_frame5[s,3] = last_row.I
end 
column_names = [:iteration, :E, :I]
# Convert the array to a DataFrame
mjp_df5 = DataFrame(mjp_frame5, column_names)
### sde 
sde_frame5 = Array{Float64}(undef, 10000, 3) 
t_diff = t_end - 0.0
for s in 1:10000
    Random.seed!(s)
    comp_counts = sim_ei_sde(init_E, init_I, alpha, nu, gamma, t_diff)
    sde_frame5[s,1] = s
    sde_frame5[s,2] = comp_counts[1]
    sde_frame5[s,3] = comp_counts[2]
end 
sde_df5 = DataFrame(sde_frame5, column_names)
CSV.write(datadir("compare_data", "ei_mjp5.csv"), mjp_df5)
CSV.write(datadir("compare_data", "ei_sde5.csv"), sde_df5)

# let's start with 20 in each
init_E = 20
init_I = 20
gamma = 1/4
nu = 1/7
rt = 1.3
alpha = rt * nu
t_end = 3.0
### mjp
mjp_frame20 = Array{Float64}(undef, 10000, 3) 
for s in 1:10000
    Random.seed!(s)
    sim_frame = sim_ei_mjp(init_E, init_I, alpha, nu, gamma, t_end)
    last_row = sim_frame[sim_frame.time .<= t_end,:][end,:]
    mjp_frame20[s,1] = s
    mjp_frame20[s,2] = last_row.E
    mjp_frame20[s,3] = last_row.I
end 
column_names = [:iteration, :E, :I]
# Convert the array to a DataFrame
mjp_df20 = DataFrame(mjp_frame20, column_names)
### sde 
sde_frame20 = Array{Float64}(undef, 10000, 3) 
t_diff = t_end - 0.0
for s in 1:10000
    Random.seed!(s)
    comp_counts = sim_ei_sde(init_E, init_I, alpha, nu, gamma, t_diff)
    sde_frame20[s,1] = s
    sde_frame20[s,2] = comp_counts[1]
    sde_frame20[s,3] = comp_counts[2]
end 
sde_df20 = DataFrame(sde_frame20, column_names)
CSV.write(datadir("compare_data", "ei_mjp20.csv"), mjp_df20)
CSV.write(datadir("compare_data", "ei_sde20.csv"), sde_df20)
