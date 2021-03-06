% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/musppr_mapping_simulation_functions.R
\name{run_bayesboot}
\alias{run_bayesboot}
\title{Run scans of Bayesian bootstrap samples to estimate a positional confidence interval for a QTL}
\usage{
run_bayesboot(
  qtl_table,
  sim_data,
  genoprobs,
  K,
  map,
  num_samples = 100,
  prob = seq(0.8, 0.99, by = 0.01)
)
}
\description{
This function ...
}
\examples{
run_bayesboot()
}
