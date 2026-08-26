function sim_ei_sde(E_init::Int64, I_init::Int64, alpha::Float64, nu::Float64, gamma::Float64, t_diff::Float64)
    # Setting up initial states and frames
    u0 = [E_init, I_init, 0,0,0]

    moment_array = closed_ei_moments(t_diff,u0, alpha, gamma, nu)
    # calculate the covariance matrix of the log compartments
    # using the delta method approximation
     log_var_matrix = reshape(vcat(moment_array[3]/(moment_array[1])^2, #varE/muE^2
     moment_array[5]/((moment_array[1])*(moment_array[2])), #covEI/muE * muI
     moment_array[5]/((moment_array[1])*(moment_array[2])), #covEI/muE * muI
     moment_array[4]/(moment_array[2])^2), (2,2)) #varI/muI^2
      # cholesky decomposition of covariance matrix of log compartments with noise on the diagonal to help pos def
     final_matrix = cholesky(Hermitian(log_var_matrix + Diagonal(fill(1e-12, 2))))
     random_normals = randn(2)
     #create compartment counts at time solve_times[i]
     log_comp_counts = log.(moment_array[1:2]) + final_matrix.L * random_normals
     comp_counts = exp.(log_comp_counts)

    return comp_counts
end


