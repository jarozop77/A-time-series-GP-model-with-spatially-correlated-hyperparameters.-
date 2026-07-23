# Shared simulated data

`simulated_data.rds` is a stored R data frame with 960 rows and 13 variables, comprising 16 districts and 60 monthly observations per district. Variables are `distcode` (district code), `t` (time index), `sim.cases` (simulated case count), `date`, `year`, `month`, `distname`, `total_population`, `inc` (incidence proportion), `inc_10000`, `log_inc`, `sin12`, and `cos12` (annual harmonic covariates). Both implementations model `log_inc`, use `sin12` and `cos12`, train on `t <= 53`, and reserve seven times for testing.

`district_boundaries.rds` is a 16-feature `sf` object used to construct district adjacency. `initial_values.rds` and `initial_covariance_matrices.rds` contain GP-MLE starting/prior quantities required only by the tailor-made sampler.

The data were stored in the original project as `data/dataset_toyexample.RDS`. No simulation-generation script was found, and the generating parameter values and random seed cannot be recovered reliably from the inspected files. Therefore no speculative generator has been created. The stored dataset is used unchanged and is shared rather than duplicated. Both `tailor_made_algorithm/run_model.R` and `standard_implementation/run_model.R` read it directly.
