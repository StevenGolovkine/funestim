################################################################################
#                Functions for bandwidth parameter estimation                  #
################################################################################

#' Perform an estimation of the bandwidth for the estimation of the mean.
#' 
#' This function performs an estimation of the bandwidth for a univariate kernel
#' regression estimator defined over continuous data.
#'
#' @importFrom magrittr %>%
#'
#' @family estimate bandwidth
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param t0 Numeric, the sampling point at which we estimate the bandwidth.
#' @param H0 Numeric, an estimation of \eqn{H_0}.
#' @param L0 Numeric, an estimation of \eqn{L_0}.
#' @param sigma Numeric, an estimation of \eqn{\sigma}.
#' @param variance Numeric, an estimation of the variance of the process.
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
estimate_bandwidth <- function(data, t0, H0 = 0.5, L0 = 1, sigma = 0, 
                               variance = 0,
                               grid = lseq(0.001, 0.1, length.out = 101),
                               nb_obs_minimal = 2, type_k = 2) {
  # Define constants
  cst_k <- switch(type_k,  
                  1 / (1 + 2 * H0), 
                  1.5 * (1 / (1 + 2 * H0) - 1 / (3 + 2 * H0)),
                  1.875 * (1 / (1 + 2 * H0) - 2 / (3 + 2 * H0) + 1 / (5 + 2 * H0)))
  q1 <- L0 / factorial(floor(H0)) * sqrt(cst_k)
  q2 <- sigma
  q3 <- sqrt(variance)
  
  risk <- rep(NA, length(grid))
  for (b in 1:length(grid)) {
    current_b <- grid[b]
    
    wi <- data %>% purrr::map_dbl(~ neighbors(.x$t, t0, current_b, nb_obs_minimal))
    WN <- sum(wi)
    if (WN == 0) next
    
    temp <- data %>% purrr::map(~ kernel((.x$t - t0) / current_b, type_k))
    Wi <- temp %>% purrr::map(~ .x / sum(.x))
    Ni <- wi / purrr::map_dbl(Wi, ~ max(.x))
    Ni[Ni == 0] <- NA
    Nmu <- WN / mean(1/Ni, na.rm = TRUE)
    
    risk[b] <- q1**2 * current_b**(2 * H0) +
      q2**2 / Nmu +
      q3**2 / WN
  }
  grid[which.min(risk)]
}

#' Perform an estimation of the bandwidth for the estimation of the mean.
#'
#' This function performs an estimation of the bandwidth to be used in the
#' Nadaraya-Watson estimator.
#' 
#' @importFrom magrittr %>%
#'
#' @param data List, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param t0_list Vector, the sampling points at which we estimate the 
#'  parameters.
#' @param grid Vector, default=NULL. A grid of bandwidths.
#' @param nb_obs_minimal Integer, default=2. Minimum number of points in the 
#'  neighborhood to keep the curve in the estimation.
#' @param type_k Integer, default=2. Used kernel.
#'  \itemize{
#'   \item \strong{1} Uniform kernel
#'   \item \strong{2} Epanechnikov kernel
#'   \item \strong{3} Biweight kernel
#'  }
#'
#' @return List, with elements:
#'  \itemize{
#'   \item \strong{sigma} An estimation of the standard deviation of the noise
#'   \item \strong{variance} An estimation of the variance of the process
#'   \item \strong{H0} An estimation of \eqn{H_0}
#'   \item \strong{L0} An estimation of \eqn{L_0}
#'   \item \strong{b} An estimation of the bandwidth
#'  }
#'
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_bandwidths <- function(data, t0_list = 0.5, grid = NULL,
                                nb_obs_minimal = 2, type_k = 2,
                                H_true = NULL) {
  if (!inherits(data, 'list')) data <- checkData(data)
  
  # Estimation of the parameters
  M <- data %>% purrr::map_dbl(~ length(.x$t)) %>% mean()
  data_presmooth <- presmoothing(data, t0_list, gamma = 0.5)
  sigma_estim <- estimate_sigma(data, t0_list, k0_list = 2)
  variance_estim <- estimate_var(data_presmooth)
  
  if (is.null(H_true)) {
    H0_estim <- estimate_H0(data_presmooth)
  } else {
    H0_estim <- H_true
  }
  L0_estim <- estimate_L0(data_presmooth, H0_estim, M)  
  
  if (is.null(grid)) {
    N <- length(data)
    Mi <- data %>% purrr::map_dbl(~ length(.x$t))
    aa <- log(1/(N*max(Mi))) / min(2 * H0_estim + 1) - log(1)
    bb <- log(1/(N*min(Mi))) / max(2 * H0_estim + 1) + log(5)
    grid <- exp(seq(aa, bb, length.out = 151))
  }
  
  # Estimation of the bandwidth
  b_estim <- list(t0_list, H0_estim, L0_estim, sigma_estim, variance_estim) %>% 
    purrr::pmap_dbl(function(t0, H0, L0, s, v){
      estimate_bandwidth(data, t0 = t0, H0 = H0, L0 = L0, sigma = s, 
                         variance = v, grid = grid, 
                         nb_obs_minimal = nb_obs_minimal, type_k = type_k)
      }
    )
  
  list(
    "sigma" = sigma_estim,
    "variance" = variance_estim,
    "H0" = H0_estim,
    "L0" = L0_estim,
    "b" = b_estim
  )
}

#' Perform an estimation of the bandwidth for the estimation of the covariance.
#' 
#' This function performs an estimation of the bandwidth for a univariate kernel
#' regression estimator defined over continuous data.
#'
#' @importFrom magrittr %>%
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
    
    wis <- data %>% purrr::map_dbl(~ neighbors(.x$t, s0, current_b, nb_obs_minimal))
    wit <- data %>% purrr::map_dbl(~ neighbors(.x$t, t0, current_b, nb_obs_minimal))
    wi <- wis * wit
    WN <- sum(wi)
    if (WN == 0) next
    
    temps <- data %>% purrr::map(~ kernel((.x$t - s0) / current_b, type_k))
    Wis <- temps %>% purrr::map(~ .x / sum(.x))
    Nis <- wi / purrr::map_dbl(Wis, ~ max(.x))
    Nis[Nis == 0] <- NaN
    Ngammas <- WN / mean(1/Nis, na.rm = TRUE)
    
    tempt <- data %>% purrr::map(~ kernel((.x$t - t0) / current_b, type_k))
    Wit <- tempt %>% purrr::map(~ .x / sum(.x))
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
