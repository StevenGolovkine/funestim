################################################################################
#                   Functions for covariance estimation                        #
################################################################################

#' Perform an estimation of the covariance with local linear smoothers.
#' 
#' This function performs the estimation of the covariance of a set of curves 
#' using local linear smoothers where the bandwidth is estimated using the 
#' methodology from Golovkine et al. (2021).
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter
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
#' @param nb_obs_minimal Integer, minimum number of observation for the smoothing.
#' @param centered Boolean (default=FALSE), are the data centered?
#' 
#' @return A list of with three entries:
#'  \itemize{
#'   \item \strong{$parameter} The estimated parameters.
#'   \item \strong{$bandwidth} The estimated bandwidths.
#'   \item \strong{$cov} The estimated covariance.
#'  }
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
covariance_ll <- function(data, U = seq(0, 1, length.out = 101),
                          t0_list = seq(0.1, 0.9, by = 0.1),
                          nb_obs_minimal = 2,
                          centered = FALSE, H_true = NULL){
  # Inner function to compute the covariance on a particular point (s, t)
  gamma_st <- function(data, s0, t0, b, n_obs_min = 2){
    data %>%
      purrr::map(~ estimate_curve(.x, c(s0, t0), b, n_obs_min = n_obs_min)$x) %>% 
      purrr::map_dbl(~ prod(.x)) %>% 
      mean(na.rm = TRUE)
  }
  # Inner function to compute the mu(s)mu(t) on a particular point (s, t)
  mu_st <- function(data, s0, t0, b, n_obs_min = 2){
    data %>%
      purrr::map(~ estimate_curve(.x, c(s0, t0), b, n_obs_min = n_obs_min)$x) %>% 
      do.call(rbind, .) %>%
      colMeans(na.rm = TRUE) %>%
      prod(na.rm = TRUE)
  }
  
  if (!centered) {
    mean_global <- mean(purrr::map_dbl(data, ~ mean(.x$x)))
    data <- data %>% purrr::map(~ list(t = .x$t, x = .x$x - mean_global))
  }
  
  if (!inherits(data, 'list')) data <- checkData(data)
  #mu_estim <- mean_ll(data, U = U, t0_list = t0_list, H_true = H_true)
  
  # Estimation of the parameters
  Mi <- data %>% purrr::map_dbl(~ length(.x$t))
  data_presmooth <- presmoothing(data, t0_list)
  sigma_estim <- max(estimate_sigma(data, delta = 0.2))
  
  if (is.null(H_true)) {
    H0_estim <- estimate_H0(data_presmooth)
  } else {
    H0_estim <- H_true
  }
  L0_estim <- estimate_L0(data_presmooth, H0_estim, mean(Mi))
  mom2 <- estimate_moment(data_presmooth, 2)
  var_st <- variance(data_presmooth)
  
  # Estimate the bandwidth on t0_list
  zz <- tidyr::expand_grid(s = t0_list, t = t0_list)
  zz <- zz %>% 
    dplyr::mutate(
      is_upper = t <= s,
      H0_s = rep(H0_estim, each = length(t0_list)),
      H0_t = rep(H0_estim, times = length(t0_list)),
      L0_s = rep(L0_estim, each = length(t0_list)),
      L0_t = rep(L0_estim, times = length(t0_list)),
      mom2_s = rep(mom2, each = length(t0_list)),
      mom2_t = rep(mom2, times = length(t0_list)),
      var_st = var_st
    )
  
  zz_nodiag <- zz %>% 
    dplyr::filter(is_upper) %>% 
    dplyr::mutate(b = purrr::pmap_dbl(
      list(s, t, 
           H0_s, H0_t,
           L0_s, L0_t, 
           mom2_s, mom2_t, 
           var_st),
      function(s, t, H0_s, H0_t, L0_s, L0_t, mom2_s, mom2_t, var_st){
        estimate_bandwidth_covariance(
          data, s, t, 
          H0 = c(H0_s, H0_t), 
          L0 = c(L0_s, L0_t),
          moment2 = c(mom2_s, mom2_t), 
          sigma = sigma_estim, 
          variance = var_st, 
          nb_obs_minimal = nb_obs_minimal, 
          type_k = 3
        )
      }
    ))
  
  # Create the bandwidth vector
  bb <- matrix(0, nrow = length(t0_list), ncol = length(t0_list))
  bb[upper.tri(bb, diag = TRUE)] <- zz_nodiag$b
  bb <- bb + t(bb) - diag(diag(bb))
  bb_large <- approx_2D(t0_list, bb, U)
  
  cov_df <- tidyr::expand_grid(s = U, t = U)
  cov_df <- cov_df %>% 
    dplyr::mutate(is_upper = t <= s, b = as.vector(bb_large)) %>%
    dplyr::filter(is_upper) %>% 
    dplyr::mutate(gamma = purrr::pmap_dbl(
      list(s, t, b), gamma_st, data = data, n_obs_min = nb_obs_minimal)) %>%
    dplyr::mutate(mu = purrr::pmap_dbl(
      list(s, t, b), mu_st, data = data, n_obs_min = nb_obs_minimal)) %>%
    dplyr::mutate(cov = gamma - mu)

  #prod_mu <- mu_estim$mu %*% t(mu_estim$mu)
  #prod_mu <- prod_mu[upper.tri(prod_mu, diag = TRUE)]
  
  # Create the final covariance
  res <- matrix(0, nrow = length(U), ncol = length(U))
  res[upper.tri(res, diag = TRUE)] <- cov_df$cov
  for (t in 1:ncol(res)) {
    s <- 1
    current_cov <- res[s, t - s + 1]
    while (s <= (t - s + 1)) {
      if (abs(U[s] - U[t - s + 1]) > bb_large[s, t - s + 1]) {
        current_cov <- res[s, t - s + 1]
      } else {
        res[s, t - s + 1] <- current_cov 
      }
      s <- s + 1
    }
  }
  for (s in 1:nrow(res)) {
    t <- ncol(res)
    current_cov <- res[ncol(res) + s - t, t]
    while (t >= (ncol(res) + s - t)) {
      if (abs(U[ncol(res) + s - t] - U[t]) > bb_large[ncol(res) + s - t, t]) {
        current_cov <- res[ncol(res) + s - t, t]
      } else {
        res[ncol(res) + s - t, t] <- current_cov 
      }
      t <- t - 1
    }
  }

  cov <- res + t(res) - diag(diag(res))
  list("parameter" = zz, "cov" = cov, "bandwidth" = bb_large)
}

#' Perform an estimation of the covariance with smoothing splines.
#' 
#' This function performs the estimation of the covariance of a set of curves
#' using smoothing splines.
#' 
#' @importFrom gss ssanova
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  }
#' @param U Vector, sampling points at which estimate the covariance.
#' @param nbasis Integer (default=5), number of basis to use for the splines.
#' @param centered Boolean (default=FALSE), are the data centered?
#' @param nodiag Boolean (default=TRUE), should the diagonal be removed from
#'  the fitting of the model?
#' 
#' @return A matrix representing the covariance surface.
#' 
#' @references Cai, T., Yuan, M. (2010) - Nonparametric covariance function 
#' estimation for functional and longitudinal data. University of Pennsylvania 
#' and Georgia inistitute of technology
#' @export
covariance_ss <- function(data, U, nbasis = 5, centered = FALSE, nodiag = TRUE){
  predict_ssanova <- utils::getFromNamespace("predict.ssanova", "gss")
  
  if (!inherits(data, 'list')) data <- checkData(data)
  data_ <- list2cai(data)
  time <- data_$time
  x <- data_$x
  subject <- data_$obs
  
  if (!centered) {
    fit <- stats::smooth.spline(time, x)
    x <- x - stats::fitted(fit)
  }
  gg <- NULL
  for (zz in unique(subject)) {
    if (sum(subject == zz) > 1) {
      tt <- time[subject == zz]
      xx <- x[subject == zz]
      g <- expand.grid(t1 = tt, t2 = tt)
      scov <- xx %*% t(xx)
      if (nodiag) scov <- scov + diag(rep(Inf, length(xx)))
      g$z <- matrix(scov, ncol = 1)
      gg <- rbind(gg, g[g$z < Inf, ])
    }
  }
  
  gg <- unique(gg)
  #nobs <- nrow(gg)
  tt <- min(time) + (max(time) - min(time)) * (1:nbasis)/(nbasis + 1)
  g <- expand.grid(t1 = tt, t2 = tt)
  g$z <- 0
  gg <- rbind(g, gg)
  
  fit <- gss::ssanova(z ~ t1 * t2, data = gg, id.basis = 1:(nbasis * nbasis))
  
  new <- expand.grid(t1 = U, t2 = U)
  estim <- predict_ssanova(fit, newdata = new)
  matrix(estim, ncol = length(U), nrow = length(U))
}

#' Perform an estimation of the covariance with local linear smoothers.
#' 
#' This function performs the estimation of the covariance of a set of curves 
#' using local linear smoothers where the bandwidth is given.
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U Vector, sampling points at which estimate the covariance.
#' @param b Numeric (default=0.1), the bandwidth to use.
#' 
#' @return A matrix representing the covariance surface.
#' 
#' @references Zhang X. and Wang J.-L. (2016) - From sparse to dense functional
#'  data and beyond, The Annals of Statistics
#' @export
covariance_lll <- function(data, U, b = 0.1){
  if (!inherits(data, 'list')) data <- checkData(data)
  data_ <- list2cai(data)
  L3 <- fdapace::MakeFPCAInputs(IDs = data_$obs, 
                                tVec = data_$time, 
                                yVec = data_$x,
                                deduplicate = TRUE)
  fdapace::GetCovSurface(L3$Ly, L3$Lt, 
                         list(kernel = 'epan',
                              nRegGrid = length(U),
                              methodMuCovEst = 'smooth',
                              userBwCov = b))$cov
}
# ----
