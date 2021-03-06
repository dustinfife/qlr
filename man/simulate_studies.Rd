% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qlr_functions.R
\name{simulate_studies}
\alias{simulate_studies}
\title{Simulate a chosen number of studies based on a fixed correlation matrix}
\usage{
simulate_studies(
  rho_matrix,
  number_of_studies = 10,
  prob_any_missing = 0.5,
  prop_missing = 0.6,
  imputations = 5,
  return.list = FALSE
)
}
\arguments{
\item{rho_matrix}{The correlation matrix (obtained from \link{random_cor_cov}).}

\item{number_of_studies}{The total number of studies to simulate}

\item{prob_any_missing}{The probability that any given study will have at least one missing correlation}

\item{prop_missing}{The proportion of variables from each study that will have missing data}
}
\value{

}
\description{
Simulate a chosen number of studies based on a fixed correlation matrix
}
