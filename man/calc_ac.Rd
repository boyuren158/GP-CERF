% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_ac.R
\name{calc_ac}
\alias{calc_ac}
\title{Calculate Covariate Balance}
\usage{
calc_ac(w, X, weights)
}
\arguments{
\item{w}{A vector exposure values across all subjects.}

\item{X}{A matrix of covariate values. Subjects in rows and covariates in columns.}

\item{weights}{A vector of weights assigned to all subjects based on the trained GP.}
}
\value{
A vector of correlations between w and each column of X.
}
\description{
Calculate weighted correlation between a list of covariates and an exposure.
Weights are defined by the covariance function of the GP.
}
\examples{

set.seed(429)

# generate data
data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

# generate random weights
weights <- runif(nrow(data))
weights <- weights/sum(weights)

# covariate matrix
design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])

cb <- GPCERF:::calc_ac(w = data$treat, X = design_mt, weights=weights)

}
\keyword{internal}
