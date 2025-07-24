## Load Data and priors for stoch conc ei normal model
if sim == 102
  sim_key = CSV.read("data/newshedding_data/sim_key.csv", DataFrame)
  # filter sim_key by sim 
  sim_vals = filter(row -> row[:sim_num] == sim, sim_key)
  all_dat = CSV.read(sim_vals.obsdata_filename[1], DataFrame)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.2)
  const rt_init_sd = 0.1
  const rt_init_mean = log(1.9)
  const tau_sd = 1
  const tau_mean = log(1)
  const E_init_sd = 0.05
  const E_init_mean = log(20)
  const I_init_sd = 0.05
  const I_init_mean = log(1)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
  long_dat = DataFrames.stack(dat, [:log_gene_copies1, 
  :log_gene_copies2, 
  :log_gene_copies3])
  data_log_copies = long_dat[:, :value]
  grid_size = 1.0
end
if sim == 106
  sim_key = CSV.read("data/newshedding_data/sim_key.csv", DataFrame)
  # filter sim_key by sim 
  sim_vals = filter(row -> row[:sim_num] == sim, sim_key)
  all_dat = CSV.read(sim_vals.obsdata_filename[1], DataFrame)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.2)
  const rt_init_sd = 0.1
  const rt_init_mean = log(1.9)
  const tau_sd = 1
  const tau_mean = log(1)
  const E_init_sd = 0.05
  const E_init_mean = log(10)
  const I_init_sd = 0.05
  const I_init_mean = log(1)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
  long_dat = DataFrames.stack(dat, [:log_gene_copies1, 
  :log_gene_copies2, 
  :log_gene_copies3])
  data_log_copies = long_dat[:, :value]
  grid_size = 1.0
end
if sim == 115
  sim_key = CSV.read("data/newshedding_data/sim_key.csv", DataFrame)
  # filter sim_key by sim 
  sim_vals = filter(row -> row[:sim_num] == sim, sim_key)
  all_dat = CSV.read(sim_vals.obsdata_filename[1], DataFrame)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.2)
  const rt_init_sd = 0.1
  const rt_init_mean = log(1.0)
  const tau_sd = 1
  const tau_mean = log(1)
  const E_init_sd = 0.05
  const E_init_mean = log(20)
  const I_init_sd = 0.05
  const I_init_mean = log(1)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
  long_dat = DataFrames.stack(dat, [:log_gene_copies1, 
  :log_gene_copies2, 
  :log_gene_copies3])
  data_log_copies = long_dat[:, :value]
  grid_size = 1.0
end
if sim == 127
  sim_key = CSV.read("data/newshedding_data/sim_key.csv", DataFrame)
  # filter sim_key by sim 
  sim_vals = filter(row -> row[:sim_num] == sim, sim_key)
  all_dat = CSV.read(sim_vals.obsdata_filename[1], DataFrame)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.2)
  const rt_init_sd = 0.1
  const rt_init_mean = log(1.9)
  const tau_sd = 1
  const tau_mean = log(1)
  const E_init_sd = 0.05
  const E_init_mean = log(40)
  const I_init_sd = 0.05
  const I_init_mean = log(1)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
  long_dat = DataFrames.stack(dat, [:log_gene_copies1, 
  :log_gene_copies2, 
  :log_gene_copies3])
  data_log_copies = long_dat[:, :value]
  grid_size = 1.0
end
if sim == "G1"
  all_dat = CSV.read(projectdir("data", "uci_ww_data.csv"), DataFrame)
  long_dat = filter(row -> row.place == sim, all_dat)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const eta_sd = 0.2
  const eta_mean = log(1/18)
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.15)
  const rt_init_sd = 0.1
  const rt_init_mean = log(0.5)
  const tau_sd = 1
  const tau_mean = 0.0
  data_log_copies = long_dat[:, :log_conc]
  grid_size = 1.0
  const E_init_sd = 0.05
  const E_init_mean = log(6)
  const I_init_sd = 0.05
  const I_init_mean = log(10)
end 
if sim == "G2"
  all_dat = CSV.read(projectdir("data", "uci_ww_data.csv"), DataFrame)
  long_dat = filter(row -> row.place == sim, all_dat)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const eta_sd = 0.2
  const eta_mean = log(1/18)
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.15)
  const rt_init_sd = 0.1
  const rt_init_mean = log(0.5)
  const tau_sd = 1
  const tau_mean = 0.0
  data_log_copies = long_dat[:, :log_conc]
  grid_size = 1.0
  const E_init_sd = 0.5
  const E_init_mean = log(7)
  const I_init_sd = 0.5
  const I_init_mean = log(12)
end 
if sim == "E"
  all_dat = CSV.read(projectdir("data", "uci_ww_data.csv"), DataFrame)
  long_dat = filter(row -> row.place == sim, all_dat)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const eta_sd = 0.2
  const eta_mean = log(1/18)
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.15)
  const rt_init_sd = 0.1
  const rt_init_mean = log(0.5)
  const tau_sd = 1
  const tau_mean = 0.0
  data_log_copies = long_dat[:, :log_conc]
  grid_size = 1.0
  const E_init_sd = 0.05
  const E_init_mean = log(5)
  const I_init_sd = 0.05
  const I_init_mean = log(9)
end 
if sim == "G1_rt1"
  all_dat = CSV.read(projectdir("data", "uci_ww_data.csv"), DataFrame)
  long_dat = filter(row -> row.place == "G1", all_dat)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const eta_sd = 0.2
  const eta_mean = log(1/18)
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.15)
  const rt_init_sd = 0.1
  const rt_init_mean = log(1)
  const tau_sd = 1
  const tau_mean = 0.0
  data_log_copies = long_dat[:, :log_conc]
  grid_size = 1.0
  const E_init_sd = 0.05
  const E_init_mean = log(6)
  const I_init_sd = 0.05
  const I_init_mean = log(10)
end 
if sim == "G2_rt1"
  all_dat = CSV.read(projectdir("data", "uci_ww_data.csv"), DataFrame)
  long_dat = filter(row -> row.place == "G2", all_dat)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const eta_sd = 0.2
  const eta_mean = log(1/18)
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.15)
  const rt_init_sd = 0.1
  const rt_init_mean = log(1)
  const tau_sd = 1
  const tau_mean = 0.0
  data_log_copies = long_dat[:, :log_conc]
  grid_size = 1.0
  const E_init_sd = 0.5
  const E_init_mean = log(7)
  const I_init_sd = 0.5
  const I_init_mean = log(12)
end 
if sim == "E_rt1"
  all_dat = CSV.read(projectdir("data", "uci_ww_data.csv"), DataFrame)
  # filter all_dat for place=="G1"
  long_dat = filter(row -> row.place == "E", all_dat)
  ## Define Priors
  const gamma_sd = 0.2
  const gamma_mean = log(1/4) 
  const nu_sd = 0.2
  const nu_mean = log(1/7) 
  const eta_sd = 0.2
  const eta_mean = log(1/18)
  const rho_conc_sd =  1.0
  const rho_conc_mean = log(1)
  const sigma_rt_sd = 0.1
  const sigma_rt_mean = log(0.15)
  const rt_init_sd = 0.1
  const rt_init_mean = log(1)
  const tau_sd = 1
  const tau_mean = 0.0
  data_log_copies = long_dat[:, :log_conc]
  grid_size = 1.0
  const E_init_sd = 0.05
  const E_init_mean = log(5)
  const I_init_sd = 0.05
  const I_init_mean = log(9)
end 
