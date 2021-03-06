% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_bandwidth.R
\name{estimate_bandwidths_mean}
\alias{estimate_bandwidths_mean}
\title{Perform an estimation of the bandwidth for the estimation of the mean.}
\usage{
estimate_bandwidths_mean(
  curves,
  grid_param = c(0.25, 0.5, 0.75),
  grid_bandwidth = NULL,
  delta_f = NULL,
  n_obs_min = 2,
  kernel_name = "epanechnikov"
)
}
\arguments{
\item{curves}{List, where each element represents a curve. Each curve have to
be defined as a list with two entries:
\itemize{
\item \strong{$t} Sampling points
\item \strong{$x} Observed points.
}}

\item{grid_param}{Vector (default = c(0.25, 0.5, 0.75)), the sampling points
at which we estimate the parameters.}

\item{grid_bandwidth}{Vector (default = NULL), grid of bandwidths.}

\item{delta_f}{Function (default = NULL), function to determine the delta.}

\item{n_obs_min}{Integer (default = 2), minimum number of points in the
neighborhood to keep the curve in the estimation.}

\item{kernel_name}{String (default = 'epanechnikov'), the kernel used for the
estimation:
\itemize{
\item epanechnikov
\item uniform
\item biweight
}}
}
\value{
Dataframe, with elements:
\itemize{
\item \strong{sigma} Estimation of the standard deviation of the noise
\item \strong{variance} Estimation of the variance of the process
\item \strong{hursts} Estimation of \eqn{H_0}
\item \strong{constants} Estimation of \eqn{L_0}
\item \strong{bandwidths} Estimation of the bandwidth
}
}
\description{
This function performs an estimation of the bandwidth to be used in the
Nadaraya-Watson estimator.
}
\references{
Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
estimation of irregular mean and covariance functions.
}
