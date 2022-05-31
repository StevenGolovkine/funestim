################################################################################
#                     Functions for parameters estimation                      #
################################################################################

# Functions for the estimation of the different parameters that are developed in
# S. Golovkine, N. Klutchnikoff and V. Patilea (2021) - Adaptive optimal
# estimation of irregular mean and covariance functions.


# Estimate sigma -- the standard deviation of the noise ----

#' Perform an estimation of the standard deviation of the noise.
#' 
#' This function performs an estimation of the standard deviation of the noise 
#' in the curves. The following formula is used:
#' \deqn{\hat{\sigma^2} = \frac{1}{N}\sum_{n = 1}^{N} 
#'       \frac{1}{2(M_n - 1)}\sum_{l = 2}^{M_n}(Y_{n, (l)} - Y_{n, (l-1)})^2}
#' 
#' @importFrom magrittr %>%
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' @param t0_list Vector, the sampling points at which we estimate \eqn{H_0}. We 
#'  will consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for the 
#'  estimation of \eqn{H_0} when \eqn{\sigma} is unknown.
#' @param k0_list Vector, the number of neighbors of \eqn{t_0} to consider. 
#' Should be set as \eqn{k0 = M\exp(-(\log(\log(M))**2))}. We can set a 
#' different \eqn{k_0}, but in order to use the same for each \eqn{t_0}, just 
#' put a unique numeric.
#'  
#' @return A list, an estimation of sigma at different \eqn{t_0}.
#' @references S. Golovkine, N. Klutchnikoff, V. Patilea (2020) - Learning the
#'  smoothness of noisy curves with application to online curve estimation.
#' @export
estimate_sigma <- function(data, t0_list, k0_list){
  if(!inherits(data, 'list')) data <- checkData(data)
  t0_list %>% 
    purrr::map2_dbl(k0_list, 
      function(t0, k0) {
        idxs <- data %>% 
          purrr::map_dbl(~ min(order(abs(.x$t - t0))[seq_len(8 * k0 - 6)]))
        df_sub <- data %>% 
          purrr::map2(idxs, ~ list(t = .x$t[.y:(.y + 8 * k0 - 7)],
                                   x = .x$x[.y:(.y + 8 * k0 - 7)]))
        estimateSigma(df_sub)
      }
    )
}
# ----

# Estimate the different quantities using pre-smoothing ----
# The quantities are:
#   * H_0 -> the regularity parameter
#   * L_0 -> the Holder constant
#   * Var(X_t) -> the variance of the process at t

#' Perform a pre-smoothing of the data.
#' 
#' This function performs a pre-smoothing of the data by local linear smoother.
#' We use a Gaussian kernel and a naive bandwidth.
#' 
#' @importFrom magrittr %>%
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  } 
#' @param t0_list Vector, the sampling point at which we pre-smooth the data.
#' @param gamma Numeric, default=0.5. Constant \eqn{\gamma} used in the theorem 
#' 1 in the paper. Should be between 0 and 1.
#' @param order Integer, default=0. Regularity of the input data.
#' @param drv Integer, default=0. Order of derivative to be estimated.
#'
#' @return List of array. Contains the smoothed data.
#' 
#' @references S. Golovkine, N. Klutchnikoff and V. Patilea (2021) - Adaptive 
#'  optimal estimation of irregular mean and covariance functions.
#' @export
presmoothing <- function(data, t0_list = 0.5, gamma = 0.5,
                         order = 1, drv = 0){
  if(!inherits(data, 'list')) data <- checkData(data)
  
  m <- data %>% purrr::map_dbl(~ length(.x$t)) %>% mean()
  delta <- exp(-log(m)**gamma)
  t1_list <- t0_list - delta / 2
  t3_list <- t0_list + delta / 2
  
  #b_naive <- (delta / round(m))**(1 / (2 * order + 1))
  b_naive <- log(m) / m
  
  results <- list()
  for(idx_t0 in 1:length(t0_list)){
    df <- array(dim = c(length(data), 11))
    for(i in 1:length(data)){
      # The bandwidth is divided by 3 because of the Gaussian kernel used in
      # KernSmooth::locpoly
      pred <- KernSmooth::locpoly(data[[i]]$t, data[[i]]$x,
                                  bandwidth = b_naive / 3,
                                  drv = drv, degree = 0, gridsize = 11,
                                  range.x = c(t1_list[idx_t0], t3_list[idx_t0]))
      df[i, ] <- pred$y
    }
    results[[idx_t0]] <- list(t = pred$x, x = df)
  }
  return(results)
}

#' Perform an estimation of \eqn{Var(X_{t_0)}}.
#' 
#' This function performs an estimation of \eqn{Var(X_{t_0})} used for the
#' estimation of the bandwidth for the mean and the covariance by a univariate
#' kernel regression estimator.
#' 
#' @importFrom magrittr %>% 
#' 
#' @param data List of array, resulting from the pre-smoothing function.
#' 
#' @return Vector, estimation of the variance at each \eqn{t_0}.
#' 
#' @references S. Golovkine, N. Klutchnikoff and V. Patilea (2021) - Adaptive 
#' optimal estimation of irregular mean and covariance functions.
#' @export
estimate_var <- function(data){
  data %>% purrr::map_dbl(~ var(.x$x[,6], na.rm = TRUE))
}

#' Perform an estimation of \eqn{H_0}.
#' 
#' This function performs an estimation of \eqn{H_0} used for the estimation of 
#' the bandwidth for a univariate kernel regression estimator defined over 
#' continuous domains data. 
#'
#' @importFrom magrittr %>%
#' @family estimate \eqn{H_0}
#' 
#' @param data List of array, resulting from the pre-smoothing function.
#' 
#' @return Vector, an estimation of \eqn{H_0} at each \eqn{t_0}.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_H0 <- function(data){
  data %>% purrr::map_dbl(function(d) {
      a <- mean((d$x[, 6] - d$x[, 1])**2, na.rm = TRUE)
      b <- mean((d$x[, 11] - d$x[, 1])**2, na.rm = TRUE)
      min(max((log(b) - log(a)) / (2 * log(2)), 0.1), 1)
    }
  )
}

#' Perform the estimation of \eqn{L_0}.
#'
#' This function performs an estimation of \eqn{L_0} used for the estimation of
#' the bandwidth for a univariate kernel regression estimator defined over 
#' continuous domains data.
#'
#' @importFrom magrittr %>%
#' @family estimate \eqn{L_0}
#' 
#' @param data List of array, resulting from the pre-smoothing function.
#' @param H0_list Vector, resulting from estimate_H0 function.
#' @param M Numeric, mean number of sampling points per curve.
#'
#' @return Vector, estimation of \eqn{L_0}.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_L0 <- function(data, H0_list, M) {
  H0 <- H0_list %>% purrr::map_dbl(~ .x - 1 / log(M)**1.01)
  V1 <- data %>% 
    purrr::map2(H0, ~ (.x$x[, 6] - .x$x[, 1])**2 / abs(.x$t[6] - .x$t[1])**(2 * .y))
  V2 <- data %>% 
    purrr::map2(H0, ~ (.x$x[, 11] - .x$x[, 6])**2 / abs(.x$t[11] - .x$t[6])**(2 * .y))
  V_mean <- V1 %>% purrr::map2_dfc(V2, ~ (.x + .y) / 2)
  unname(sqrt(colMeans(V_mean, na.rm = TRUE)))
}

#' Perform the estimation of the moments.
#'
#' This function performs an estimation of the moments used for the estimation 
#' of the bandwidth for a univariate kernel regression estimator defined over 
#' continuous domains data.
#'
#' @importFrom magrittr %>%
#' @family estimate moment
#' 
#' @param data List of array, resulting from the pre-smoothing function.
#' @param order Numeric (default=1), the moment to estimate.
#'
#' @return Vector, estimation of the moments.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_moment <- function(data, order = 1) {
  data %>% purrr::map_dbl(~ mean(.x$x[, 6]**order, na.rm = TRUE))
}

#' Perform the estimation of \eqn{Var(X_{s}X_{t)}}.
#'
#' This function performs an estimation of  \eqn{Var(X_{s}X_{t)}} used for the 
#' estimation of the bandwidth for a univariate kernel regression estimator 
#' defined over continuous domains data.
#'
#' @importFrom magrittr %>%
#' @family estimate variance
#' 
#' @param data List of array, resulting from the pre-smoothing function.
#'
#' @return Vector, estimation of \eqn{Var(X_{s}X_{t)}}.
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
variance <- function(data) {
  var_st <- matrix(NA, nrow = length(data), ncol = length(data))
  for(i in 1:length(data)){
    for(j in 1:length(data)){
      var_st[i, j] <- stats::var(data[[i]]$x[, 6] * data[[j]]$x[,6], na.rm = TRUE)
    }
  }
  as.vector(var_st)
}
# ----
