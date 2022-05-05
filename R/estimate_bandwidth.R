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
#' @param params List, estimation of the different parameters:
#'  \itemize{
#'   \item \strong{$point} Time point where the smoothing has been done.
#'   \item \strong{$H} Estimated regularity.
#'   \item \strong{$L} Estimated constant.
#'   \item \strong{$var} Estimated variance.
#'  }
#' @param sigma Numeric, estimation of the std of the noise \eqn{\sigma}.
#' @param grid_bandwidth Vector (default = lseq(0.001, 0.1, length.out = 101)), 
#' grid of bandwidths.
#' @param n_obs_min Integer (default = 2), minimum number of points in the 
#'  neighborhood to keep the curve in the estimation.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#'
#' @return Numeric, estimation of the bandwidth.
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_bandwidth_mean <- function(
    curves, params, sigma = 0, 
    grid_bandwidth = lseq(0.001, 0.1, length.out = 101),
    n_obs_min = 2, kernel_name = 'epanechnikov') {
  # Define constants
  cst_kernel <- switch(
    kernel_name,  
    uniform = 1 / (1 + 2 * params$H), 
    epanechnikov = 1.5 * (1 / (1 + 2 * params$H) - 1 / (3 + 2 * params$H)),
    biweight = 1.875 * (
      1 / (1 + 2 * params$H) - 2 / (3 + 2 * params$H) + 1 / (5 + 2 * params$H)
    )
  )
  q1 <- params$L / factorial(floor(params$H)) * sqrt(cst_kernel)
  q2 <- sigma
  q3 <- sqrt(params$var)
  
  risk <- rep(NA, length(grid_bandwidth))
  for (idx in 1:length(grid_bandwidth)) {
    current_bandwidth <- grid_bandwidth[idx]
    
    wi <- curves |> 
      sapply(function(curve) {
        neighbors(curve$t, params$point, current_bandwidth, n_obs_min)
      })
    WN <- sum(wi)
    if (WN == 0) next
    
    temp <- curves |> 
      lapply(function(curve) {
        kernel((curve$t - params$point) / current_bandwidth, type = kernel_name)
      })
    Wi <- temp |> 
      lapply(function(curve) {
      curve / sum(curve)
      })
    Ni <- wi / sapply(Wi, function(curve) max(curve))
    Ni[Ni == 0] <- NA
    Nmu <- WN / mean(1/Ni, na.rm = TRUE)
    
    risk[idx] <- q1**2 * current_bandwidth**(2 * params$H) +
      q2**2 / Nmu +
      q3**2 / WN
  }
  
  list(
    "point" = params$point,
    "H" = params$H,
    "L" = params$L,
    "var" = params$var,
    "bandwidth" = grid_bandwidth[which.min(risk)]
  )
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
#' @param grid_param Vector (default = c(0.25, 0.5, 0.75)), the sampling points 
#' at which we estimate the parameters.
#' @param grid_bandwidth Vector (default = NULL), grid of bandwidths.
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param n_obs_min Integer (default = 2), minimum number of points in the 
#'  neighborhood to keep the curve in the estimation.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
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
estimate_bandwidths_mean <- function(
    curves, grid_param = c(0.25, 0.5, 0.75), grid_bandwidth = NULL,
    delta_f = NULL, n_obs_min = 2, kernel_name = 'epanechnikov') {
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  
  # Estimation of the parameters -- TO MODIFY FOR RECURSIVE ESTIMATION
  m <- curves |> sapply(function(curve) length(curve$t)) |> mean()
  sigma_estim <- estimate_sigma(curves, delta = delta_f(m))
  params_estim <- estimate_parameters(
    curves, grid = grid_param, delta_f = delta_f, kernel = kernel_name, beta = 1)
  
  if (is.null(grid_bandwidth)) {
    N <- length(curves)
    Mi <- curves |> sapply(function(curve) length(curve$t))
    hurst_estim <- sapply(params_estim, function(param) param$H)
    aa <- log(1 / m)
    bb <- log(1/(N * min(Mi))) / max(2 * hurst_estim + 1) + log(5)
    grid_bandwidth <- exp(seq(aa, bb, length.out = 151))
  }
  
  # Estimation of the bandwidth -- CODE IN C++
  params_estim |>
    lapply(function(param) {
      estimate_bandwidth_mean(
        curves, param, sigma = sigma_estim, 
        grid_bandwidth = grid_bandwidth,
        n_obs_min = n_obs_min, kernel_name = kernel_name)
    })
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
