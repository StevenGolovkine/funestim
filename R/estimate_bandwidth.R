################################################################################
#                Functions for bandwidth parameter estimation                  #
################################################################################

################################################################################
# FOR MEAN
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
#' @return List, estimation of the bandwidth.
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
        kernel_f((curve$t - params$point) / current_bandwidth, type = kernel_name)
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
    "mom" = params$mom,
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
#' @return Dataframe, with elements:
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
  
  # Estimation of the parameters
  m <- curves |> sapply(function(curve) length(curve$t)) |> mean()
  sigma_estim <- estimate_sigma(curves, delta = delta_f(m))
  params_estim <- estimate_parameters_mean(
    curves, grid = grid_param, 
    delta_f = delta_f, kernel_name = kernel_name, beta = 1)
  
  if (is.null(grid_bandwidth)) {
    N <- length(curves)
    Mi <- curves |> sapply(function(curve) length(curve$t))
    hurst_estim <- unlist(params_estim$H)
    aa <- log(1 / m)
    bb <- log(1/(N * min(Mi))) / max(2 * hurst_estim + 1) + log(5)
    grid_bandwidth <- exp(seq(aa, bb, length.out = 151))
  }
  
  # Estimation of the bandwidth
  params_estim |>
    apply(1, function(param) {
      estimate_bandwidth_mean(
        curves, param, sigma = sigma_estim, 
        grid_bandwidth = grid_bandwidth,
        n_obs_min = n_obs_min, kernel_name = kernel_name)
    }) |> 
    (\(x) do.call("rbind", x))() |> 
    as.data.frame()
}

################################################################################
# FOR COVARIANCE
################################################################################

#' Perform an estimation of the bandwidth for the estimation of the covariance.
#' 
#' This function performs an estimation of the bandwidth for a univariate kernel
#' regression estimator defined over continuous data.
#'
#' @family estimate bandwidth
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param params List, estimation of the different parameters for the data
#' points pair to estimate:
#'  \itemize{
#'   \item \strong{$point} Time point where the smoothing has been done.
#'   \item \strong{$H} Estimated regularity.
#'   \item \strong{$L} Estimated constant.
#'   \item \strong{$var} Estimated variance.
#'   \item \strong{$mom} Estimated \eqn{E(X^{2}_{t_0})}.
#'   \item \strong{$var_st} Estimated \eqn{E(X^{2}_{t_0})}.
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
#' @return Numeric, an estimation of the bandwidth.
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_bandwidth_covariance <- function(
    curves, params, sigma = 0, 
    grid_bandwidth = lseq(0.001, 0.1, length.out = 101),
    n_obs_min = 2, kernel_name = 'epanechnikov') {
  # Constant definition
  H <- unlist(c(params$H_s, params$H_t))
  L <- unlist(c(params$L_t, params$L_t))
  mom <- unlist(c(params$mom_s, params$mom_t))
  # Define constants
  cst_k <- switch(
    kernel_name,  
    uniform = 1 / (1 + 2 * H), 
    epanechnikov = 1.5 * (1 / (1 + 2 * H) - 1 / (3 + 2 * H)),
    biweight = 1.875 * (1 / (1 + 2 * H) - 2 / (3 + 2 * H) + 1 / (5 + 2 * H))
  )
  q1_s <- sqrt(2 * mom[2] * cst_k[1]) * L[1] / factorial(floor(H[1]))
  q1_t <- sqrt(2 * mom[1] * cst_k[2]) * L[2] / factorial(floor(H[2]))
  q2_s <- sigma * sqrt(mom[2])
  q2_t <- sigma * sqrt(mom[1])
  q3 <- sqrt(params$var_st)
  
  risk <- rep(NA, length(grid_bandwidth))
  for (idx in 1:length(grid_bandwidth)) {
    current_bandwidth <- grid_bandwidth[idx]

    wis <- curves |> 
      sapply(function(curve) {
        neighbors(curve$t, params$point_s, current_bandwidth, n_obs_min)
      })
    wit <- curves |> 
      sapply(function(curve) {
        neighbors(curve$t, params$point_t, current_bandwidth, n_obs_min)
      })
    wi <- wis * wit
    WN <- sum(wi)
    if (WN == 0) next

    temps <- curves |> 
      lapply(function(curve) {
        kernel_f((curve$t - params$point_s) / current_bandwidth, type = kernel_name)
      })
    Wis <- temps |> lapply(function(curve) {
      curve / sum(curve)
    })
    Nis <- wi / sapply(Wis, function(curve) max(curve))
    Nis[Nis == 0] <- NaN
    Ngammas <- WN / mean(1/Nis, na.rm = TRUE)
    
    tempt <- curves |> 
      lapply(function(curve) {
        kernel_f((curve$t - params$point_t) / current_bandwidth, type = kernel_name)
      })
    Wit <- tempt |> lapply(function(curve) {
      curve / sum(curve)
    })
    Nit <- wi / sapply(Wit, function(curve) max(curve))
    Nit[Nit == 0] <- NaN
    Ngammat <- WN / mean(1/Nit, na.rm = TRUE)

    risk[idx] <- q1_s**2 * current_bandwidth**(2 * H[1]) +
      q1_t**2 * current_bandwidth**(2 * H[2]) +
      q2_s**2 / Ngammas + q2_t**2 / Ngammat +
      q3**2 / WN
  }
  
  list(
    "point_s" = params$point_s,
    "H_s" = params$H_s,
    "L_s" = params$L_s,
    "var_s" = params$var_s,
    "mom_s" = params$mom_s,
    "point_t" = params$point_t,
    "H_t" = params$H_t,
    "L_t" = params$L_t,
    "var_t" = params$var_t,
    "mom_t" = params$mom_t,
    "var_st" = params$var_st,
    "bandwidth" = grid_bandwidth[which.min(risk)]
  )
}


#' Perform an estimation of the bandwidth for the estimation of the covariance
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
estimate_bandwidths_covariance <- function(
    curves, grid_param = c(0.25, 0.5, 0.75), grid_bandwidth = NULL,
    delta_f = NULL, n_obs_min = 2, kernel_name = 'epanechnikov') {
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  
  # Estimation of the parameters
  m <- curves |> sapply(function(curve) length(curve$t)) |> mean()
  sigma_estim <- estimate_sigma(curves, delta = delta_f(m))
  params_estim <- estimate_parameters_covariance(
    curves, grid = grid_param, 
    delta_f = delta_f, kernel_name = kernel_name, beta = 1)
  
  if (is.null(grid_bandwidth)) {
    N <- length(curves)
    Mi <- curves |> sapply(function(curve) length(curve$t))
    hurst_estim <- unlist(params_estim$H_s)
    aa <- log(1 / m)
    bb <- log(1/(N * min(Mi))) / max(2 * hurst_estim + 1) + log(5)
    grid_bandwidth <- exp(seq(aa, bb, length.out = 151))
  }
  
  # Estimation of the bandwidth
  params_estim |>
    apply(1, function(param) {
      estimate_bandwidth_covariance(
        curves, param, sigma = sigma_estim, 
        grid_bandwidth = grid_bandwidth,
        n_obs_min = n_obs_min, kernel_name = kernel_name)
    }) |> 
    (\(x) do.call("rbind", x))() |> 
    as.data.frame()
}
