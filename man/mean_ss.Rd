% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_mean.R
\name{mean_ss}
\alias{mean_ss}
\title{Perform an estimation of the mean with smoothing splines.}
\usage{
mean_ss(curves, grid)
}
\arguments{
\item{curves}{List, where each element represents a curve. Each curve have to
be defined as a list with two entries:
\itemize{
\item \strong{$t} Sampling points
\item \strong{$x} Observed points.
}}

\item{grid}{Vector, sampling points at which estimate the mean.}
}
\value{
Vector representing the mean curve.
}
\description{
This function performs the estimation of the mean of a set of curves using
smoothing splines.
}
\references{
Cai T., Yuan M. (2011) - Optimal Estimation of the mean function
based on discretely sampled functional data: phase transition, The Annals of
Statistics
}
