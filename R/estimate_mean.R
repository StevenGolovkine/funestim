################################################################################
#                       Functions for mean estimation                          #
################################################################################

#' Perform an estimation of the mean with local linear smoothers.
#' 
#' This function performs the estimation of the mean of a set of curves using 
#' local linear smoothers where the bandwidth is estimated using the 
#' methodology from Golovkine et al. (2021).
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#' be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points.
#'  } 
#' @param grid Vector (default = seq(0, 1, length.out = 101)), sampling points
#' at which estimate the curves.
#' @param grid_param Vector (default = c(0.25, 0.5, 0.75)), sampling points at
#' which we estimate the parameters.
#' @param grid_bandwidth Vector (default = NULL), grid of bandwidths.
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param n_obs_min Integer (default = 2), minimum number of observation for the
#' smoothing.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#' 
#' @return List of with two entries:
#'  \itemize{
#'   \item \strong{$parameter} Estimated parameters.
#'   \item \strong{$mu} Estimated mean.
#'  }
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#' estimation of irregular mean and covariance functions.
#' @export
mean_ll <- function(
    curves,
    grid = seq(0, 1, length.out = 101),
    grid_param = c(0.25, 0.5, 0.75),
    grid_bandwidth = NULL, delta_f = NULL,
    n_obs_min = 2, kernel_name = 'epanechnikov'
){
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  curves_smooth <- smooth_curves_mean(
    curves, grid = grid, grid_param = grid_param, 
    grid_bandwidth = grid_bandwidth, delta_f = delta_f,
    kernel_name = kernel_name, n_obs_min = n_obs_min)
  mu <- curves_smooth$curves_smooth |> 
    (\(x) do.call("rbind", x))() |> 
    colMeans(na.rm = TRUE)
  list("parameter" = curves_smooth$parameter, "mu" = mu)
}

#' Perform an estimation of the mean with smoothing splines.
#' 
#' This function performs the estimation of the mean of a set of curves using
#' smoothing splines.
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#' be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points.
#'  }
#' @param grid Vector, sampling points at which estimate the mean.
#' 
#' @return Vector representing the mean curve.
#'  
#' @references Cai T., Yuan M. (2011) - Optimal Estimation of the mean function 
#' based on discretely sampled functional data: phase transition, The Annals of 
#' Statistics
#' @export
mean_ss <- function(curves, grid){
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  curves_ <- list2cai(curves)
  mod <- stats::smooth.spline(curves_$time, curves_$x)
  stats::predict(mod, grid)$y
}

#' Perform an estimation of the mean with local linear smoothers.
#' 
#' This function performs the estimation of the mean of a set of curves using 
#' local linear smoothers where the bandwidth is estimated using the 
#' methodology from Zhang et Wang (2016).
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#' be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points.
#'  } 
#' @param grid Vector, sampling points at which estimate the mean.
#' 
#' @return Vector representing the mean curve.
#' 
#' @references Zhang X. and Wang J.-L. (2016) - From sparse to dense functional
#' data and beyond, The Annals of Statistics
#' @export
mean_lll <- function(curves, grid) {
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  curves_ <- list2cai(curves)
  L3 <- fdapace::MakeFPCAInputs(
    IDs = curves_$obs, tVec = curves_$time, yVec = curves_$x, deduplicate = TRUE)
  fdapace::GetMeanCurve(
    L3$Ly, L3$Lt,
    list(kernel = 'epan', nRegGrid = length(grid), methodBwMu = 'GCV')
  )$mu
}
# ----
