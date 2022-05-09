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
#' @param bandwidth Vector, estimation of the bandwidth. If a unique element is 
#' provided, we use a  unique bandwidth for the curve. However, if a vector is 
#' given, the bandwidth changes depending on the sampling points. 
#' @param bandwidth_times Vector (default = NULL), times at which the bandwidths 
#' have been estimated. Only used if the parameter \code{bandwidth} is a vector.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the 
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
                           kernel_name = "epanechnikov", n_obs_min = 1) {
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
  
  if (kernel_name == "epanechnikov") {
    x_hat <- epaKernelSmoothingCurve(
      grid, curve$t, curve$x, bandwidth, n_obs_min)
  } else if (kernel_name == "uniform") {
    x_hat <- uniKernelSmoothingCurve(
      grid, curve$t, curve$x, bandwidth, n_obs_min)
  } else if (kernel_name == "biweight") {
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
#' @param grid_param Vector (default = c(0.25, 0.5, 0.75)), sampling points at
#' which we estimate the parameters.
#' @param grid_bandwidth Vector (default = NULL), grid of bandwidths.
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param n_obs_min Integer (default = 2), minimum number of observation for 
#' the smoothing.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#'
#' @return A list, which contains two elements. The first one is a list which 
#'  contains the estimated parameters:
#'  \itemize{
#'   \item \strong{sigma} Estimation of the standard deviation of the noise.
#'   \item \strong{variance} Estimation of the variance of the process.
#'   \item \strong{H0} Estimation of \eqn{H_0}.
#'   \item \strong{L0} Estimation of \eqn{L_0}.
#'   \item \strong{bandwidth} Estimation of the bandwidth.
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
smooth_curves_mean <- function(
    curves, grid = NULL, grid_param = c(0.25, 0.5, 0.75),
    grid_bandwidth = NULL, delta_f = NULL,
    kernel_name = 'epanechnikov', n_obs_min = 2
){
  # Estimation of the different parameters
  param_estim <- estimate_bandwidths_mean(
    curves, grid_param = grid_param, grid_bandwidth = grid_bandwidth,
    delta_f = delta_f, n_obs_min = n_obs_min, kernel_name = kernel_name)
  bandwidth_estim <- unlist(param_estim$bandwidth)
  
  # Estimation of the curves
  if (is.null(grid)) {
    curves_estim <- curves |> lapply(function(curve) {
      estimate_curve(curve, grid = curve$t, bandwidth = bandwidth_estim,
                     bandwidth_times = grid_param, kernel_name = kernel_name,
                     n_obs_min = n_obs_min)
    })
  } else {
    curves_estim <- curves |> lapply(function(curve) {
      estimate_curve(curve, grid = grid, bandwidth = bandwidth_estim,
                     bandwidth_times = grid_param, kernel_name = kernel_name,
                     n_obs_min = n_obs_min)
    })
  }
  
  list("parameters" = param_estim, "curves_smooth" = curves_estim)
}

#' Perform a non-parametric smoothing of a set of curves for covariance estimation.
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
#' @param grid Vector (default = seq(0, 1, length.out = 101)), sampling points 
#' at which estimate the curves.
#' @param grid_param Vector (default = c(0.25, 0.5, 0.75)), sampling points at
#' which we estimate the parameters.
#' @param grid_bandwidth Vector (default = NULL), grid of bandwidths.
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param n_obs_min Integer (default = 2), minimum number of observation for 
#' the smoothing.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#'
#' @return A list, which contains three elements. The first one is a list which 
#'  contains the estimated parameters:
#'  \itemize{
#'   \item \strong{sigma} Estimation of the standard deviation of the noise.
#'   \item \strong{variance} Estimation of the variance of the process.
#'   \item \strong{H0} Estimation of \eqn{H_0}.
#'   \item \strong{L0} Estimation of \eqn{L_0}.
#'   \item \strong{bandwidth} Estimation of the bandwidth.
#'  }
#'  The second one is the bandwidths matrix. And the last one is the estimation
#'  of the covariance.
#'  
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
smooth_curves_covariance <- function(
    curves, grid = seq(0, 1, length.out = 101),
    grid_param = c(0.25, 0.5, 0.75), grid_bandwidth = NULL, 
    delta_f = NULL, n_obs_min = 2, kernel_name = 'epanechnikov'
){
  # Inner function to compute the covariance on a particular point (s, t)
  gamma_st <- function(curves, point_s, point_t, bandwidth, n_obs_min = 2){
    curves |>
      lapply(function(curve) {
        estimate_curve(curve, c(point_s, point_t), 
                       bandwidth, bandwidth_times = NULL,
                       kernel_name = kernel_name, n_obs_min = n_obs_min)
      }) |> 
      sapply(function(curve) {
        prod(curve)
      }) |> 
      mean(na.rm = TRUE)
  }
  
  # Estimation of the different parameters
  param_estim <- estimate_bandwidths_covariance(
    curves, grid_param = grid_param, grid_bandwidth = grid_bandwidth,
    delta_f = delta_f, n_obs_min = n_obs_min, kernel_name = kernel_name)

  # Create the bandwidth vector
  bb <- matrix(0, nrow = length(grid_param), ncol = length(grid_param))
  bb[upper.tri(bb, diag = TRUE)] <- unlist(param_estim$bandwidth)
  bb <- bb + t(bb) - diag(diag(bb))
  bb_large <- approx_2D(grid_param, bb, grid)
  
  # Estimation of the curves
  zz <- expand.grid(point_s = grid, point_t = grid)
  zz$bandwidth <- as.vector(bb_large)
  zz_upper <- zz[zz$point_t <= zz$point_s, ]
  zz_upper$cov <- zz_upper |> apply(1, function(row) {
    gamma_st(curves, 
             row['point_s'], row['point_t'], row['bandwidth'], 
             n_obs_min = n_obs_min)
  })
  
  list("parameters" = param_estim, 
       "bandwidths_mat" = bb_large,
       "covariance" = zz_upper$cov)
}
# ----
