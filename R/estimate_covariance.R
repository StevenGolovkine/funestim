################################################################################
#                   Functions for covariance estimation                        #
################################################################################

#' Perform an estimation of the covariance with local linear smoothers.
#'
#' This function performs the estimation of the covariance of a set of curves
#' using local linear smoothers where the bandwidth is estimated using the
#' methodology from Golovkine et al. (2021).
#'
#' @param curves List, where each element represents a curve. Each curve have to
#' be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points.
#'  }
#' @param grid Vector (default = seq(0, 1, length.out = 101)), sampling points
#' at which estimate the curves.
#' @param grid_param Vector (default = seq(0.1, 0.9, by = 0.1)), the sampling
#' points at which we estimate the parameters.
#' @param grid_bandwidth Vector (default = NULL), grid of bandwidths.
#' @param center Boolean (default = TRUE), center the data?
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param n_obs_min Integer (default = 2), minimum number of observation for the
#' smoothing.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#'  
#' @return List of with three entries:
#'  \itemize{
#'   \item \strong{$parameters} Estimated parameters.
#'   \item \strong{$bandwidths_mat} Estimated bandwidths matrix.
#'   \item \strong{$covariance} Estimated covariance.
#'  }
#'
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
covariance_ll <- function(
    curves,
    grid = seq(0, 1, length.out = 101),
    grid_param = seq(0.1, 0.9, by = 0.1),
    grid_bandwidth = NULL,
    center = TRUE,
    delta_f = NULL,
    n_obs_min = 2, kernel_name = 'epanechnikov'
    
){
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  if (center) {
    mean_global <- curves |> 
      sapply(function(curve) curve$x) |> 
      unlist() |> 
      mean()
    curves <- curves |> 
      lapply(function(curve) {
        list(t = curve$t, x = curve$x - mean_global)
      })
  }

  # Compute E(X)
  mu_estim <- mean_ll(
    curves, grid = grid, grid_param = grid_param, 
    grid_bandwidth = grid_bandwidth, delta_f = delta_f,
    n_obs_min = n_obs_min, kernel_name = kernel_name)

  # Compute E(X^2)
  cov_estim <- smooth_curves_covariance(
    curves, grid = grid, grid_param = grid_param,
    grid_bandwidth = grid_bandwidth, delta_f = delta_f,
    n_obs_min = n_obs_min, kernel_name = kernel_name)

  # Compute E(X)^2
  prod_mu <- mu_estim$mu %*% t(mu_estim$mu)
  prod_mu <- prod_mu[upper.tri(prod_mu, diag = TRUE)]

  # Create the final covariance E(X^2) - E(X)^2
  results <- matrix(0, nrow = length(grid), ncol = length(grid))
  results[upper.tri(results, diag = TRUE)] <- cov_estim$covariance - prod_mu
  for (t in 1:ncol(results)) {
    s <- 1
    current_cov <- results[s, t - s + 1]
    while (s <= (t - s + 1)) {
      if (abs(grid[s] - grid[t - s + 1]) > 
          cov_estim$bandwidths_mat[s, t - s + 1]) {
          current_cov <- results[s, t - s + 1]
      } else {
        results[s, t - s + 1] <- current_cov
      }
      s <- s + 1
    }
  }
  for (s in 1:nrow(results)) {
    t <- ncol(results)
    current_cov <- results[ncol(results) + s - t, t]
    while (t >= (ncol(results) + s - t)) {
      if (abs(grid[ncol(results) + s - t] - grid[t]) >
          cov_estim$bandwidths_mat[ncol(results) + s - t, t]) {
        current_cov <- results[ncol(results) + s - t, t]
      } else {
        results[ncol(results) + s - t, t] <- current_cov
      }
      t <- t - 1
    }
  }

  covariance <- results + t(results) - diag(diag(results))
  list("parameters" = cov_estim$parameters, 
       "bandwidths_mat" = cov_estim$bandwidths_mat,
       "covariance" = covariance)
}

#' Perform an estimation of the covariance with smoothing splines.
#'
#' This function performs the estimation of the covariance of a set of curves
#' using smoothing splines.
#'
#' @importFrom gss ssanova
#'
#' @param curves list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points.
#'  }
#' @param grid Vector, sampling points at which estimate the covariance.
#' @param nbasis Integer (default = 5), number of basis to use for the splines.
#' @param center Boolean (default = TRUE), center the data?
#' @param nodiag Boolean (default = TRUE), should the diagonal be removed from
#' the fitting of the model?
#'
#' @return Matrix representing the covariance surface.
#'
#' @references Cai, T., Yuan, M. (2010) - Nonparametric covariance function
#' estimation for functional and longitudinal data. University of Pennsylvania
#' and Georgia institute of technology
#' @export
covariance_ss <- function(curves, grid,
                          nbasis = 5, center = TRUE, nodiag = TRUE){
  predict_ssanova <- utils::getFromNamespace("predict.ssanova", "gss")

  if (!inherits(curves, 'list')) curves <- checkData(curves)
  curves_ <- list2cai(curves)
  time <- curves_$time
  x <- curves_$x
  subject <- curves_$obs

  if (center) {
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
  tt <- min(time) + (max(time) - min(time)) * (1:nbasis)/(nbasis + 1)
  g <- expand.grid(t1 = tt, t2 = tt)
  g$z <- 0
  gg <- rbind(g, gg)

  fit <- gss::ssanova(z ~ t1 * t2, data = gg, id.basis = 1:(nbasis * nbasis),
                      weights = NULL, subset = NULL, offset = NULL)

  new <- expand.grid(t1 = grid, t2 = grid)
  estim <- predict_ssanova(fit, newdata = new)
  matrix(estim, ncol = length(grid), nrow = length(grid))
}

#' Perform an estimation of the covariance with local linear smoothers.
#'
#' This function performs the estimation of the covariance of a set of curves
#' using local linear smoothers where the bandwidth is given.
#'
#' @param curves List, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points.
#'  }
#' @param grid Vector, sampling points at which estimate the covariance.
#' @param bandwidth Numeric (default = 0.1), the bandwidth to use.
#'
#' @return Matrix representing the covariance surface.
#'
#' @references Zhang X. and Wang J.-L. (2016) - From sparse to dense functional
#'  data and beyond, The Annals of Statistics
#' @export
covariance_lll <- function(curves, grid, bandwidth = 0.1){
  if (!inherits(curves, 'list')) curves <- checkData(curves)
  curves_ <- list2cai(curves)
  L3 <- fdapace::MakeFPCAInputs(
    IDs = curves_$obs, tVec = curves_$time, yVec = curves_$x,
    deduplicate = TRUE)
  fdapace::GetCovSurface(
    L3$Ly, L3$Lt,
    list(kernel = 'epan', nRegGrid = length(grid),
         methodMuCovEst = 'smooth', userBwCov = bandwidth)
    )$cov
}
# ----
