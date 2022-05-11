################################################################################
#                       Functions for mean estimation                          #
################################################################################

#' Perform an estimation of the mean with local linear smoothers.
#' 
#' This function performs the estimation of the mean of a set of curves using 
#' local linear smoothers where the bandwidth is estimated using the 
#' methodology from Golovkine et al. (2021).
#' 
#' @importFrom magrittr %>%
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U Vector, sampling points at which estimate the curves.
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
#' @return A list of with two entries:
#'  \itemize{
#'   \item \strong{$parameter} The estimated parameters.
#'   \item \strong{$mu} The estimated mean.
#'  }
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
mean_ll <- function(data, U = seq(0, 1, length.out = 101), t0_list = 0.5, 
                    grid = NULL, nb_obs_minimal = 2, kernel = 'epanechnikov',
                    H_true = NULL){
  if (!inherits(data, 'list')) data <- checkData(data)
  data_smooth <- smooth_curves(data, U = U, t0_list = t0_list, grid = grid,
                               nb_obs_minimal = nb_obs_minimal, kernel = kernel,
                               H_true = H_true)
  mu <- data_smooth$smooth %>% 
    purrr::map_dfc(~ .x$x) %>% 
    rowMeans(na.rm = TRUE)
  list("parameter" = data_smooth$parameter, "mu" = mu)
}

#' Perform an estimation of the mean with smoothing splines.
#' 
#' This function performs the estimation of the mean of a set of curves using
#' smoothing splines.
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  }
#' @param U Vector, sampling points at which estimate the mean.
#' 
#' @return A vector representing the mean curve.
#'  
#' @references Cai T., Yuan M. (2011) - Optimal Estimation of the mean function 
#' based on discretely sampled functional data: phase transition, The Annals of 
#' Statistics
#' @export
mean_ss <- function(data, U){
  if (!inherits(data, 'list')) data <- checkData(data)
  data_ <- list2cai(data)
  mod <- stats::smooth.spline(data_$time, data_$x)
  stats::predict(mod, U)$y
}

#' Perform an estimation of the mean with local linear smoothers.
#' 
#' This function performs the estimation of the mean of a set of curves using 
#' local linear smoothers where the bandwidth is estimated using the 
#' methodology from Zhang et Wang (2016).
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U Vector, sampling points at which estimate the mean.
#' 
#' @return A vector representing the mean curve.
#' 
#' @references Zhang X. and Wang J.-L. (2016) - From sparse to dense functional
#'  data and beyond, The Annals of Statistics
#' @export
mean_lll <- function(data, U) {
  if (!inherits(data, 'list')) data <- checkData(data)
  data_ <- list2cai(data)
  L3 <- fdapace::MakeFPCAInputs(IDs = data_$obs,
                                tVec = data_$time,
                                yVec = data_$x,
                                deduplicate = TRUE)
  fdapace::GetMeanCurve(L3$Ly, L3$Lt, 
                        list(kernel = 'epan', 
                             nRegGrid = length(U),
                             methodBwMu = 'GCV'))$mu
}
# ----
