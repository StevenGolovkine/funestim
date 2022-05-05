################################################################################
#                Functions for bandwidth parameter estimation                  #
################################################################################

#' Perform an estimation of the bandwidth for the estimation of the mean.
#' 
#' This function performs an estimation of the bandwidth for a univariate kernel
#' regression estimator defined over continuous data.
#'
#' @family estimate bandwidth
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#' be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points.
#'  } 
#' @param point Numeric, sampling point at which we estimate the bandwidth.
#' @param hurst Numeric, estimation of the Hurst coefficient \eqn{H_0}.
#' @param constant Numeric, estimation of the constant \eqn{L_0}.
#' @param sigma Numeric, estimation of the std of the noise \eqn{\sigma}.
#' @param variance Numeric, estimation of the variance of the process.
#' @param grid_bandwidth Vector (default = lseq(0.001, 0.1, length.out = 11)), 
#' grid of bandwidths.
#' @param n_obs_min Integer (default = 2), minimum number of points in the 
#'  neighborhood to keep the curve in the estimation.
#' @param type_kernel Integer (default = 2), Used kernel.
#'  \itemize{
#'   \item \strong{1} uniform kernel
#'   \item \strong{2} epanechnikov kernel
#'   \item \strong{3} biweight kernel
#'  }
#'
#' @return Numeric, estimation of the bandwidth.
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_bandwidth <- function(
    curves, point, 
    hurst = 0.5, constant = 1, 
    sigma = 0, variance = 0,
    grid_bandwidth = lseq(0.001, 0.1, length.out = 11),
    n_obs_min = 2, type_kernel = 2) {
  # Define constants
  cst_kernel <- switch(
    type_kernel,  
    1 / (1 + 2 * hurst), 
    1.5 * (1 / (1 + 2 * hurst) - 1 / (3 + 2 * hurst)),
    1.875 * (1 / (1 + 2 * hurst) - 2 / (3 + 2 * hurst) + 1 / (5 + 2 * hurst))
  )
  q1 <- constant / factorial(floor(hurst)) * sqrt(cst_kernel)
  q2 <- sigma
  q3 <- sqrt(variance)
  
  risk <- rep(NA, length(grid_bandwidth))
  for (idx in 1:length(grid_bandwidth)) {
    current_bandwidth <- grid_bandwidth[idx]
    
    wi <- curves |> 
      purrr::map_dbl(~ neighbors(.x$t, point, current_bandwidth, n_obs_min))
    WN <- sum(wi)
    if (WN == 0) next
    
    temp <- curves |> 
      purrr::map(~ kernel((.x$t - point) / current_bandwidth, type_kernel))
    Wi <- temp |> purrr::map(~ .x / sum(.x))
    Ni <- wi / purrr::map_dbl(Wi, ~ max(.x))
    Ni[Ni == 0] <- NA
    Nmu <- WN / mean(1/Ni, na.rm = TRUE)
    
    risk[idx] <- q1**2 * current_bandwidth**(2 * hurst) +
      q2**2 / Nmu +
      q3**2 / WN
  }
  grid_bandwidth[which.min(risk)]
}

#' Perform an estimation of the bandwidth for the estimation of the mean.
#'
#' This function performs an estimation of the bandwidth to be used in the
#' Nadaraya-Watson estimator.
#'
#' @param curves List, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points.
#'  } 
#' @param grid_param Vector, the sampling points at which we estimate the 
#'  parameters.
#' @param grid_bandwidth Vector (default = NULL), grid of bandwidths.
#' @param n_obs_min Integer (default = 2), minimum number of points in the 
#'  neighborhood to keep the curve in the estimation.
#' @param type_kernel Integer (default = 2), used kernel:
#'  \itemize{
#'   \item \strong{1} uniform kernel
#'   \item \strong{2} epanechnikov kernel
#'   \item \strong{3} biweight kernel
#'  }
#'
#' @return List, with elements:
#'  \itemize{
#'   \item \strong{sigma} Estimation of the standard deviation of the noise
#'   \item \strong{variance} Estimation of the variance of the process
#'   \item \strong{hursts} Estimation of \eqn{H_0}
#'   \item \strong{constants} Estimation of \eqn{L_0}
#'   \item \strong{bandwidths} Estimation of the bandwidth
#'  }
#'
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_bandwidths <- function(
    curves, grid_param = 0.5, grid_bandwidth = NULL,
    n_obs_min = 2, type_kernel = 2) {
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  
  # Estimation of the parameters -- TO MODIFY FOR RECURSIVE ESTIMATION
  M <- curves |> purrr::map_dbl(~ length(.x$t)) |> mean()
  data_presmooth <- presmoothing(curves, grid_param, gamma = 0.5)
  sigma_estim <- estimate_sigma(curves, grid_param, k0_list = 2)
  variance_estim <- estimate_var(data_presmooth)
  hurst_estim <- estimate_H0(data_presmooth)
  constant_estim <- estimate_L0(data_presmooth, hurst_estim, M)
  
  if (is.null(grid_bandwidth)) {
    N <- length(curves)
    Mi <- curves |> purrr::map_dbl(~ length(.x$t))
    aa <- log(1/(N*max(Mi))) / min(2 * hurst_estim + 1) - log(1)
    bb <- log(1/(N*min(Mi))) / max(2 * hurst_estim + 1) + log(5)
    grid_bandwidth <- exp(seq(aa, bb, length.out = 151))
  }
  
  # Estimation of the bandwidth -- CODE IN C++
  bandwidth_estim <- list(grid_param, hurst_estim, constant_estim, 
                          sigma_estim, variance_estim) |>
    purrr::pmap_dbl(function(t0, H0, L0, s, v){
      estimate_bandwidth(curves, point = t0, hurst = H0, constant = L0, 
                         sigma = s, variance = v, 
                         grid_bandwidth = grid_bandwidth, n_obs_min = n_obs_min,
                         type_kernel = type_kernel)
      }
    )

  list(
    "sigma" = sigma_estim,
    "variance" = variance_estim,
    "H0" = hurst_estim,
    "L0" = constant_estim,
    "b" = bandwidth_estim
  )
}

#' Perform an estimation of the bandwidth for the estimation of the covariance.
#' 
#' This function performs an estimation of the bandwidth for a univariate kernel
#' regression estimator defined over continuous data.
#'
#' @family estimate bandwidth
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param s0 Numeric, the sampling point at which we estimate the bandwidth.
#' @param t0 Numeric, the sampling point at which we estimate the bandwidth.
#' @param H0 Vector, an estimation of \eqn{H_0} at \eqn{s_0} and \eqn{t_0}.
#' @param L0 Vector, an estimation of \eqn{L_0} at \eqn{s_0} and \eqn{t_0}.
#' @param moment2 Vector, an estimation of \eqn{EX^2} at \eqn{s_0} and \eqn{t_0}.
#' @param sigma Numeric, an estimation of \eqn{\sigma}.
#' @param variance Numeric, an estimation of \eqn{Var(X_{s_0}X_{t_0)}}.
#' @param grid Vector, default=lseq(0.001, 0.1, length.out = 101). A grid of 
#'  bandwidths.
#' @param nb_obs_minimal Integer, default=2. Minimum number of points in the 
#'  neighborhood to keep the curve in the estimation.
#' @param type_k Integer, default=2. Used kernel.
#'  \itemize{
#'   \item \strong{1} Uniform kernel
#'   \item \strong{2} Epanechnikov kernel
#'   \item \strong{3} Biweight kernel
#'  }
#'
#' @return Numeric, an estimation of the bandwidth.
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_bandwidth_covariance <- function(data, s0, t0,
                                          H0 = c(0.5, 0.5),
                                          L0 = c(1, 1), 
                                          moment2 = c(1, 1),
                                          sigma = 0, 
                                          variance = 0,
                                          grid = NULL,
                                          nb_obs_minimal = 2, 
                                          type_k = 2) {
  if (!inherits(data, 'list')) data <- checkData(data)
  if (is.null(grid)) grid <- lseq(0.01, 0.1, length.out = 41)
  
  # Constant definition
  cst_k <- switch(type_k,
                  1 / (1 + 2 * H0),
                  1.5 * (1 / (1 + 2 * H0) - 1 / (3 + 2 * H0)),
                  1.875 * (1 / (1 + 2 * H0) - 2 / (3 + 2 * H0) + 1 / (5 + 2 * H0)))
  q1_s <- sqrt(2 * moment2[2] * cst_k[1]) * L0[1] / factorial(floor(H0[1]))
  q1_t <- sqrt(2 * moment2[1] * cst_k[2]) * L0[2] / factorial(floor(H0[2]))
  q2_s <- sigma * sqrt(moment2[2])
  q2_t <- sigma * sqrt(moment2[1])
  q3 <- sqrt(variance)
  
  risk <- rep(NA, length(grid))
  for (b in 1:length(grid)) {
    current_b <- grid[b]
    
    wis <- data |> purrr::map_dbl(~ neighbors(.x$t, s0, current_b, nb_obs_minimal))
    wit <- data |> purrr::map_dbl(~ neighbors(.x$t, t0, current_b, nb_obs_minimal))
    wi <- wis * wit
    WN <- sum(wi)
    if (WN == 0) next
    
    temps <- data |> purrr::map(~ kernel((.x$t - s0) / current_b, type_k))
    Wis <- temps |> purrr::map(~ .x / sum(.x))
    Nis <- wi / purrr::map_dbl(Wis, ~ max(.x))
    Nis[Nis == 0] <- NaN
    Ngammas <- WN / mean(1/Nis, na.rm = TRUE)
    
    tempt <- data |> purrr::map(~ kernel((.x$t - t0) / current_b, type_k))
    Wit <- tempt |> purrr::map(~ .x / sum(.x))
    Nit <- wi / purrr::map_dbl(Wit, ~ max(.x))
    Nit[Nit == 0] <- NaN
    Ngammat <- WN / mean(1/Nit, na.rm = TRUE)

    risk[b] <- q1_s**2 * current_b**(2 * H0[1]) +
      q1_t**2 * current_b**(2 * H0[2]) +
      q2_s**2 / Ngammas + q2_t**2 / Ngammat +
      q3**2 / WN
  }
  grid[which.min(risk)]
}
# ----
