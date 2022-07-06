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
#' @param delta Numeric (default = 0.1), neighborhood for the estimation.
#'  
#' @return A list, an estimation of sigma at different \eqn{t_0}.
#' @references S. Golovkine, N. Klutchnikoff, V. Patilea (2020) - Learning the
#'  smoothness of noisy curves with application to online curve estimation.
#' @export
estimate_sigma <- function(data, delta = 0.1){
  estimateSigma(data, delta = delta)
}
# ----

# Estimate the minimum of the density ----

#' Perform an estimation of the minimum of the sampling points density.
#' 
#' This function performs an estimation of the minimum of the density of the
#' noise.
#' 
#' @importFrom magrittr %>%
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  
#' @return An estimation of the minimum of the density of the sampling points.
#' @export
estimate_density <- function(data){
  T_all <- data %>% purrr::map(~.x$t) %>% unlist() %>% sort()
  min(density(T_all, from = 0.1, to = 0.9)$y)
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
presmoothing <- function(
    data,
    t0_list = seq(.2, .8, l = 20),
    init_b = 1,
    init_L = 1
){
  if (!inherits(data, 'list')) data <- checkData(data)
  
  NW <- function(t, T_mi, Y_mi, h, alpha = 1) {
    K <- function(x, beta){
      ((1 + beta) / (2 * beta))  * (1 - abs(x)^beta) * (abs(x) <= 1)
    }
    tt <- matrix(rep(t, length(T_mi)), ncol = length(T_mi))
    TT <- t(matrix(rep(T_mi, length(t)), ncol = length(t)))
    Kx <- 2 * K(x = (tt - TT) / h, beta = alpha)
    (Kx %*% Y_mi) / matrix(rowSums(Kx))
  }
  
  m <- data %>% purrr::map_dbl(~ length(.x$t)) %>% mean()
  
  delta <- min(log(m)^(-1.1), 0.2)
  t1_list <- t0_list - delta / 2
  t3_list <- t0_list + delta / 2
  
  sigma <- estimate_sigma(data)
  mu0 <- estimate_density(data)
  
  #b_naive <- (delta / round(m))**(1 / (2 * order + 1))
  #b_naive <- log(m) / m
  aa <- (init_b + 1) / 2 * init_b**2 * mu0
  c <- (sigma**(2*init_b) * init_L * aa**init_b)**(1 / (2*init_b + 1))
  psi_m <- (1 / m)**(init_b / (2 * init_b + 1))
  b_naive <- pmax(pmin((c * psi_m / init_L)**(1 / init_b), delta/4), log(m)/m)
  
  inner_loop <- function(i, data, t_list, b_naive, init_b) {
    sapply(data, function(x) {
      NW(t = t_list[, i],
         T_mi = x$t, Y_mi = x$x,
         h = b_naive, alpha = init_b)
    }) %>% t() 
  }
  
  t_list <- rbind(t1_list, t0_list, t3_list)
  purrr::map(1:ncol(t_list), ~list(
    t_list = t_list[,.x],
    x = inner_loop(.x, data = data, t_list = t_list,
                   b_naive = b_naive, init_b = 1)))
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
  data %>% purrr::map_dbl(~ var(.x$x[,2], na.rm = TRUE))
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
      a <- mean((d$x[, 3] - d$x[, 1])**2, na.rm = TRUE)
      b <- mean((d$x[, 2] - d$x[, 1])**2, na.rm = TRUE)
      c <- mean((d$x[, 3] - d$x[, 2])**2, na.rm = TRUE)
      max(min((2 * log(a) - log(b * c)) / log(16), 1), 0.1)
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
    purrr::map2(H0, 
                ~ (.x$x[, 2] - .x$x[, 1])**2 / abs(.x$t[2] - .x$t[1])**(2 * .y))
  V2 <- data %>% 
    purrr::map2(H0, 
                ~ (.x$x[, 3] - .x$x[, 2])**2 / abs(.x$t[3] - .x$t[2])**(2 * .y))
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
  data %>% purrr::map_dbl(~ mean(.x$x[, 2]**order, na.rm = TRUE))
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
  for (i in 1:length(data)) {
    for (j in 1:length(data)) {
      var_st[i, j] <- stats::var(data[[i]]$x[, 2] * data[[j]]$x[, 2], na.rm = TRUE)
    }
  }
  as.vector(var_st)
}
# ----
