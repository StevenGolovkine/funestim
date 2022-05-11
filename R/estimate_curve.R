################################################################################
#        Functions that performs kernel smoothing over a set of curves         #
################################################################################

#' Perform the smoothing of an individual curve.
#'
#' This function performs the smoothing of a curve using the Nadaraya-Watson 
#' estimator given a particular kernel.
#' 
#' @importFrom magrittr %>%
#' @importFrom Rcpp evalCpp
#' 
#' @param curve A list, with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U Vector, sampling points at which estimate the curve. 
#' @param b Vector, estimation of the bandwidth. If one is provided, we use a 
#' unique bandwidth for the curve. However, if a vector is given, the bandwidth 
#' changes depending on the sampling points. 
#' @param t0_list Vector, times at which the bandwidths have been 
#'  estimated. Only used if the parameter \code{b} is a vector.
#' @param kernel Character string, the kernel used for the estimation:
#'  \itemize{
#'   \item epanechnikov (default)
#'   \item uniform
#'   \item beta
#'   \item mBeta 
#'  }
#' @param n_obs_min Integer, minimum number of observation for the smoothing.
#' @useDynLib funestim
#'
#' @return List, with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The estimated points.
#'  }
#'  
#' @export
estimate_curve <- function(curve, U, b, t0_list = NULL,
                           kernel = "epanechnikov", n_obs_min = 1) {
  if (length(b) == 1) {
    bandwidth <- rep(b, length(U))
  } else if ((length(b) != length(U)) & !is.null(t0_list)) {
    bandwidth <- stats::approx(t0_list, b, xout = U,
                               yleft = b[1], yright = b[length(b)],
                               ties = 'ordered')$y
  } else if (length(b) == length(U)) {
    bandwidth <- b
  } else {
    stop("Issues with the bandwidth parameter.")
  }
  
  if (kernel == "epanechnikov") {
    x_hat <- epaKernelSmoothingCurve(U, curve$t, curve$x, bandwidth, n_obs_min)
  } else if (kernel == "uniform") {
    x_hat <- uniKernelSmoothingCurve(U, curve$t, curve$x, bandwidth, n_obs_min)
  } else if (kernel == "biweight") {
    x_hat <- biweightKernelSmoothingCurve(U, curve$t, curve$x, bandwidth, n_obs_min)
  } else {
    print("Wrong kernel name")
    x_hat <- rep(0, length(U))
  }
  list(t = U, x = as.vector(x_hat))
}

#' Perform a non-parametric smoothing of a set of curves for mean estimation.
#'
#' This function performs a non-parametric smoothing of a set of curves using 
#' the Nadaraya-Watson estimator.
#' 
#' @importFrom magrittr %>%
#'
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U Vector, default=NULL. Sampling points at which estimate the curves.
#'  If NULL, the sampling points for the estimation are the same than the 
#'  observed ones.
#' @param t0_list Vector, the sampling points at which we estimate the 
#'  parameters.
#' @param grid Vector, default=NULL. A grid of bandwidths.
#' @param nb_obs_minimal Integer, minimum number of observation for the smoothing.
#' @param kernel Character string, the kernel used for the estimation:
#'  \itemize{
#'   \item epanechnikov (default)
#'   \item uniform
#'   \item biweight
#'  }
#'
#' @return A list, which contains two elements. The first one is a list which 
#'  contains the estimated parameters:
#'  \itemize{
#'   \item \strong{sigma} An estimation of the standard deviation of the noise
#'   \item \strong{variance} An estimation of the variance of the process
#'   \item \strong{H0} An estimation of \eqn{H_0}
#'   \item \strong{L0} An estimation of \eqn{L_0}
#'   \item \strong{b} An estimation of the bandwidth
#'  }
#'  The second one is another list which contains the estimation of the curves:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The estimated points.
#'  }
#'  
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
smooth_curves <- function(data, U = NULL, t0_list = 0.5, grid = NULL, 
                          nb_obs_minimal = 2, kernel = 'epanechnikov',
                          H_true = NULL){
  
  if (kernel == 'uniform')
    type_k = 1
  else if (kernel == 'epanechnikov')
    type_k = 2
  else if (kernel == 'biweight')
    type_k = 3
  else
    type_k = 1
  
  # Estimation of the different parameters
  param_estim <- estimate_bandwidths(data, t0_list = t0_list, grid = grid,
                                     nb_obs_minimal = nb_obs_minimal,
                                     type_k = type_k, H_true = H_true)

  # Estimation of the curves
  if (is.null(U)) {
    curves <- data %>% 
      purrr::map(~ estimate_curve(.x, U = .x$t, b = param_estim$b,
                                  t0_list = t0_list, kernel = kernel,
                                  n_obs_min = nb_obs_minimal))
  } else {
    curves <- data %>% 
      purrr::map(~ estimate_curve(.x, U = U, b = param_estim$b,
                                  t0_list = t0_list, kernel = kernel,
                                  n_obs_min = nb_obs_minimal))
  }
  list("parameter" = param_estim, "smooth" = curves)
}
# ----
