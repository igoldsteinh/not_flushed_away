# not_flushed_away
This repository has all code needed to recreate the analyses conducted in the paper The Signal is Not Flushed Away: Inferring the Effective Reproduction Number From Wastewater Data in Small Populations.
Models were fit in `Julia`, while simulation of synthetic data, and visualization of results was done in `R`. 
All results files needed to reproduce the figures are in the repo, individual simulation results are excluded for the sake of storage.

## Navigation
```
├── data                          <- Processed real and simulated data
│   └── newshedding_data                  <- Simulated data
│
├── figures                       <- Paper figures
│
├── results                       <- Model outputs, organized by model, and then by output type
│   │                                Example for the stochastic EI-ww model is shown
│   │                                Structure is the same for main models
│   └── stoch_ei_normal           <- Summaries of simulation results
│        ├── generated_quantities <- Correctly scaled posteriors, quantiles and mcmc samples
│        ├── posterior_predictive <- Posterior predictive mcmc samples and quantiles
│        └── posterior_samples    <- Raw Julia posterior samples
│
├── scripts                       <- Paper code 
│   ├── fit_models                <- Fit models to real and simulated data
│   ├── generated_quantities      <- Turn raw Julia MCMC output into more useable csv files
│   ├── process_results           <- Final processing of mcmc/summarising simulation results
│   ├── simulated_data            <- Simulate data
│   └── visualize_results         <- Turn summaries of model results into paper figures
│   
├── src                           <- Models, priors, simulation engines, utility functions
│   
├── vignettes                     <- Example code for fitting models and processing results
└──     
```

## Setting up the `Julia` environment. 
The results from this project were generated using `Julia 1.9.3`. 
If you want to use a more recent version of `Julia`, delete the [Manifest.toml](https://github.com/igoldsteinh/not_flushed_away/blob/main/Manifest.toml) file.
We recommend installing `Julia` via [`juliaup`](https://github.com/JuliaLang/juliaup).
Once you have `Julia` installed, from the terminal, navigate to the project root directory then type `julia`. 
Your terminal will look like:
```
julia>
```
Now type `]`. Your terminal should now look like:
```
(@v1.9.3) pkg>
```
Then use the following commands
```
activate .
```
and 
```
instantiate
```
More information on `Julia` environments is available in the [Environments documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project).

## Quarto and Julia
The [vignettes folder](https://github.com/igoldsteinh/not_flushed_away/tree/main/vignettes) has two Quarto vignettes which condense the model fitting workflow into one [`Julia` vignette](https://github.com/igoldsteinh/not_flushed_away/blob/main/vignettes/fit_ei_ww.qmd) that demonstrates how to fit the stochastic EI-ww model to the UCI wastewater data via MCMC and one [`R` vignette](https://github.com/igoldsteinh/not_flushed_away/blob/main/vignettes/process_ei_ww.qmd) that uses the results of the `Julia` vignette to visualize the saved MCMC results. 

To execute the vignettes, we recommend using the IDE [VS Code](https://code.visualstudio.com) with the [`Julia`](https://code.visualstudio.com/docs/languages/julia) and [Quarto](https://quarto.org/docs/tools/vscode.html) extensions. 
Additional information on compiling `Julia` Quarto files is available [here](https://quarto.org/docs/computations/julia.html). 
You can also run each chunk of `Julia` code by copy pasting it into the REPL (the interactive Julia environment that opens when you type `julia` in the terminal).

## Model fitting workflow
The original workflow for the main models involves multiple files. 
As an example, to generate results from the the stochastic EI-ww model, use [fit_stoch_ei_normal.jl](https://github.com/igoldsteinh/not_flushed_away/blob/main/scripts/fit_models/fit_stoch_ei_normal.jl) to fit the model, then [gq_stoch_ei_normal.jl](https://github.com/igoldsteinh/not_flushed_away/blob/main/scripts/generate_quantities/gq_stoch_ei_normal.jl) to re-scale the posterior and generate posterior predictive values, 
finally [process_results_stoch_ei_normal_newshedding.R](https://github.com/igoldsteinh/not_flushed_away/blob/main/scripts/process_results/process_results_stoch_ei_normal_newshedding.R) creates tidy versions of the posterior and posterior predictive summaries.
When summarising results from multiple simulations, [summarise_results_stoch_ei_normal_newshedding.R](https://github.com/igoldsteinh/not_flushed_away/blob/main/scripts/process_results/summarise_results_stoch_ei_normal_newshedding.R) creates summary outputs. 
Similarly named files exist for the deterministic EI-ww model. 

## Simulation name key
When executing scripts, the `sim` parameter controls what simulation is being used, the `seed` parameter controls the seed and also the specific data set used. 
For simulations, we used values of `seed` from 1 to 100. 
Here is a key translating the values of `sim`:
* `sim=102` = `Fixed`
* `sim=115` = `Steep`
* `sim=106` = `Total500`
* `sim=127` = `Total2000`
More detailed information on simulations are stored in the [simulation key file](https://github.com/igoldsteinh/not_flushed_away/blob/main/data/newshedding_data/sim_key.csv).

## Model name key
The model names used in the code are not the same as those used in the paper. 
Here is a key:
* `stoch_ei_normal` = `stochastic EI-ww`
* `det_ei_normal` = `deterministic EI-ww`
* `epidemia` = `Epidemia-cases`
* `episewer` = `Episewer`

## Fitting Episewer
To install Episewer, please follow the [package instructions](https://adrian-lison.github.io/EpiSewer/).
To fit the model use [fit_episewer.R](https://github.com/igoldsteinh/not_flushed_away/blob/main/scripts/episewer/fit_episewer.R) for simulated data and [fit_episewer_uci.R](https://github.com/igoldsteinh/not_flushed_away/blob/main/scripts/episewer/fit_episewer_uci.R) for real data. 

## Fitting Epidemia
The instructions for installing Epidemia are [here](https://imperialcollegelondon.github.io/epidemia/articles/install.html). 
Unfortunately, Epidemia fails to install on recent versions of R.
The file used to run it for this paper is [fit_uci_casesv2.R](https://github.com/igoldsteinh/not_flushed_away/blob/main/scripts/fit_models/fit_uci_casesv2.R).
