################################################################################
#        Functions that performs kernel smoothing over a set of curves         #
################################################################################

#' Perform the smoothing of an individual curve.
#'
#' This function performs the smoothing of a curve using the Nadaraya-Watson 
#' estimator given a particular kernel.
#' 
#' @param curve List, with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points.
#'   \item \strong{$x} Observed points.
#'  } 
#' @param grid Vector, sampling points at which the curve is estimated.
#' @param bandwith Vector, estimation of the bandwidth. If a unique element is 
#' provided, we use a  unique bandwidth for the curve. However, if a vector is 
#' given, the bandwidth changes depending on the sampling points. 
#' @param bandwith_times Vector (default = NULL), times at which the bandwidths 
#' have been estimated. Only used if the parameter \code{bandwidth} is a vector.
#' @param kernel String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#' @param n_obs_min Integer (default = 1), minimum number of observation for 
#' the smoothing.
#' @useDynLib funestim
#'
#' @return List, with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points.
#'   \item \strong{$x} Estimated points.
#'  }
#'  
#' @export
estimate_curve <- function(curve, grid, bandwidth, bandwidth_times = NULL,
                           kernel = "epanechnikov", n_obs_min = 1) {
  if (length(bandwidth) == 1) {
    bandwidth <- rep(bandwidth, length(grid))
  } else if ((length(bandwidth) != length(grid)) & !is.null(bandwidth_times)) {
    bandwidth <- stats::approx(
      bandwidth_times, bandwidth, xout = grid,
      yleft = bandwidth[1], yright = bandwidth[length(bandwidth)],
      ties = 'ordered')$y
  } else if (length(bandwidth) == length(grid)) {
    bandwidth <- bandwidth
  } else {
    stop("Issues with the bandwidth parameter.")
  }
  
  if (kernel == "epanechnikov") {
    x_hat <- epaKernelSmoothingCurve(
      grid, curve$t, curve$x, bandwidth, n_obs_min)
  } else if (kernel == "uniform") {
    x_hat <- uniKernelSmoothingCurve(
      grid, curve$t, curve$x, bandwidth, n_obs_min)
  } else if (kernel == "biweight") {
    x_hat <- biweightKernelSmoothingCurve(
      grid, curve$t, curve$x, bandwidth, n_obs_min)
  } else {
    print("Wrong kernel name")
    x_hat <- rep(0, length(grid))
  }
  as.vector(x_hat)
}

#' Perform a non-parametric smoothing of a set of curves for mean estimation.
#'
#' This function performs a non-parametric smoothing of a set of curves using 
#' the Nadaraya-Watson estimator.
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#' be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points.
#'   \item \strong{$x} Observed points.
#'  } 
#' @param grid Vector (default = NULL), sampling points at which estimate the 
#' curves. If NULL, the sampling points for the estimation are the same than the
#' observed ones.
#' @param grid_param Vector (default = 0.5), sampling points at which we 
#' estimate the parameters.
#' @param grid_bandwidth Vector (default = NULL), grid of bandwidths.
#' @param kernel String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#' @param n_obs_min Integer (default = 2), minimum number of observation for 
#' the smoothing.
#'
#' @return A list, which contains two elements. The first one is a list which 
#'  contains the estimated parameters:
#'  \itemize{
#'   \item \strong{sigma} Estimation of the standard deviation of the noise.
#'   \item \strong{variance} Estimation of the variance of the process.
#'   \item \strong{H0} Estimation of \eqn{H_0}.
#'   \item \strong{L0} Estimation of \eqn{L_0}.
#'   \item \strong{b} Estimation of the bandwidth.
#'  }
#'  The second one is another list which contains the estimation of the curves:
#'  \itemize{
#'   \item \strong{$t} Sampling points.
#'   \item \strong{$x} Estimated points.
#'  }
#'  
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
smooth_curves <- function(curves, grid = NULL,
                          grid_param = 0.5, grid_bandwidth = NULL, 
                          kernel = 'epanechnikov', n_obs_min = 2){
  
  if (kernel == 'uniform') type_k = 1
  else if (kernel == 'epanechnikov') type_k = 2
  else if (kernel == 'biweight') type_k = 3
  else type_k = 1
  
  # Estimation of the different parameters
  param_estim <- estimate_bandwidths(curves, t0_list = grid_param, 
                                     grid = grid_bandwidth,
                                     nb_obs_minimal = n_obs_min,
                                     type_k = type_k)

  # Estimation of the curves
  if (is.null(grid)) {
    curves_estim <- curves |> 
      purrr::map(~ estimate_curve(.x, U = .x$t, b = param_estim$b,
                                  t0_list = grid_param, kernel = kernel,
                                  n_obs_min = n_obs_min))
  } else {
    curves_estim <- curves |> 
      purrr::map(~ estimate_curve(.x, U = grid, b = param_estim$b,
                                  t0_list = grid_param, kernel = kernel,
                                  n_obs_min = n_obs_min))
  }
  list("parameter" = param_estim, "smooth" = curves_estim)
}
# ----
