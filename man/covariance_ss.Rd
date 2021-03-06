% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_covariance.R
\name{covariance_ss}
\alias{covariance_ss}
\title{Perform an estimation of the covariance with smoothing splines.}
\usage{
covariance_ss(curves, grid, nbasis = 5, center = TRUE, nodiag = TRUE)
}
\arguments{
\item{curves}{list, where each element represents a curve. Each curve have to
be defined as a list with two entries:
\itemize{
\item \strong{$t} Sampling points
\item \strong{$x} Observed points.
}}

\item{grid}{Vector, sampling points at which estimate the covariance.}

\item{nbasis}{Integer (default = 5), number of basis to use for the splines.}

\item{center}{Boolean (default = TRUE), center the data?}

\item{nodiag}{Boolean (default = TRUE), should the diagonal be removed from
the fitting of the model?}
}
\value{
Matrix representing the covariance surface.
}
\description{
This function performs the estimation of the covariance of a set of curves
using smoothing splines.
}
\references{
Cai, T., Yuan, M. (2010) - Nonparametric covariance function
estimation for functional and longitudinal data. University of Pennsylvania
and Georgia institute of technology
}
