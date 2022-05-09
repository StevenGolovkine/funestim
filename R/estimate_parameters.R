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
#' @param curves List, where each element represents a curve. Each curve have to
#'  be defined as list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points.
#'   \item \strong{$x} Observed points.
#'  }
#' @param delta Numeric (default = 0.1), neighborhood for the estimation.
#'  
#' @return Numeric, estimation of the std of the noise \eqn{\sigma}
#' 
#' @references S. Golovkine, N. Klutchnikoff, V. Patilea (2020) - Learning the
#'  smoothness of noisy curves with application to online curve estimation.
#' @export
estimate_sigma <- function(curves, delta = 0.1){
  estimateSigma(curves, delta)
}
# ----

# Estimate the different quantities using pre-smoothing ----
# The quantities are:
#   * hurst -> the regularity parameter (Hurst coefficient).
#   * constant -> the Holder constant.
#   * Var(X_t) -> the variance of the process at point t.

#' Perform a pre-smoothing of the data.
#' 
#' This function performs a pre-smoothing of the data using a Nadaraya-Watson
#' estimator. We use an Epanechnikov kernel and a naive bandwidth.
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} Sampling points
#'   \item \strong{$x} Observed points
#'  } 
#' @param point Numeric (default = 0.5), sampling point at which the data is 
#' pre-smoothed.
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param kernel String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#' @param beta Numeric (default = 1), pre-specified regularity of the curves.
#'  The default value is 1, which correspond to at least one time
#'  differentiable curves.
#' @param bandwidth_naive Numeric (default = 0), bandwidth to use for the
#' presmoothing. If set to 0, the bandwidth will be defined as 
#' \deqn{\frac{\delta}{m}^{1 / (2\beta + 1)}}
#' where
#'  \itemize{
#'   \item \eqn{m} is the mean number of sampling points per curve.
#'   \item \eqn{\delta} is the length of the interval where the smoothing is
#'    done.
#'  \item \eqn{\beta} represents the regularity of the curves.
#'  }
#' 
#' @return List, with two entries:
#'  \itemize{
#'   \item \strong{$grid} Grid on which the smoothing has been done.
#'   \item \strong{$x_smooth} The smoothed data.
#'  }
#' 
#' @references S. Golovkine, N. Klutchnikoff and V. Patilea (2021) - Adaptive 
#'  optimal estimation of irregular mean and covariance functions.
#' @export
presmoothing <- function(curves, point = 0.5, delta_f = NULL,
                         kernel = 'epanechnikov',
                         beta = 1, bandwidth_naive = 0){
  m <- curves |> sapply(function(curve) length(curve$t)) |> mean()
  delta <- delta_f(m)
  t_vec <- c(point - delta / 2, point, point + delta / 2)
  
  if (bandwidth_naive == 0)
    bandwidth_naive <- (delta / m)**(1 / (2 * beta + 1))
  
  list(
    grid = t_vec,
    x_smooth = sapply(curves, function(curve) {
      estimate_curve(curve, grid = t_vec, bandwidth = bandwidth_naive)
    })
  )
}

#' Perform an estimation of \eqn{Var(X_{t_0)}}.
#' 
#' This function performs an estimation of \eqn{Var(X_{t_0})} used for the
#' estimation of the bandwidth for the mean and the covariance by a univariate
#' kernel regression estimator.
#' 
#' @param curves_smooth List, resulting from the `presmoothing` function.
#' 
#' @return Numeric, estimation of the variance at \eqn{t_0}.
#' 
#' @references S. Golovkine, N. Klutchnikoff and V. Patilea (2021) - Adaptive 
#' optimal estimation of irregular mean and covariance functions.
#' @export
estimate_var <- function(curves_smooth){
  stats::var(curves_smooth$x_smooth[2,], na.rm = TRUE)
}

#' Perform an estimation of the regularity \eqn{H_0}.
#' 
#' This function performs an estimation of \eqn{H_0} used for the estimation of 
#' the bandwidth for a univariate kernel regression estimator defined over 
#' continuous domains data. 
#' 
#' @family estimate \eqn{H_0}
#' 
#' @param curves_smooth List, resulting from the `presmoothing` function.
#' 
#' @return Numeric, an estimation of \eqn{H_0} at \eqn{t_0}.
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_regularity <- function(curves_smooth){
  current_smooth <- curves_smooth$x_smooth
  a <- mean((current_smooth[3,] - current_smooth[1,])**2, na.rm = TRUE)
  b <- mean((current_smooth[2,] - current_smooth[1,])**2, na.rm = TRUE)
  c <- mean((current_smooth[3,] - current_smooth[2,])**2, na.rm = TRUE)
  max(min((2 * log(a) - log(b * c)) / log(16), 1), 0.1)
}

#' Perform the estimation of the constant \eqn{L_0}.
#'
#' This function performs an estimation of \eqn{L_0} used for the estimation of
#' the bandwidth for a univariate kernel regression estimator defined over 
#' continuous domains data.
#'
#' @family estimate \eqn{L_0}
#' 
#' @param curves_smooth List, resulting from the `presmoothing` function.
#' @param regularity Numeric, estimation of the regularity of the curves and
#' resulting from the the `estimate_regularity` function.
#'
#' @return Numeric, an estimation of \eqn{L_0} at \eqn{t_0}.
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_constant <- function(curves_smooth, regularity) {
  current_grid <- curves_smooth$grid
  current_smooth <- curves_smooth$x_smooth
  a <- mean((current_smooth[3,] - current_smooth[1,])**2, na.rm = TRUE)
  b <- abs(current_grid[3] - current_grid[1])**(2 * regularity)
  sqrt(a / b)
}

#' Perform the estimation of the moments.
#'
#' This function performs an estimation of the moments \eqn{E(X^{\alpha}_{t_0})}
#' used for the estimation of the bandwidth for a univariate kernel regression
#' estimator defined over continuous domains data.
#'
#' @family estimate moment
#' 
#' @param curves_smooth List, resulting from the `presmoothing` function.
#' @param order Numeric (default = 1), the moment to estimate.
#'
#' @return Numeric, an estimation of the moments \eqn{E(X^{\alpha}_{t_0})}
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_moment <- function(curves_smooth, order = 1) {
  mean(curves_smooth$x_smooth[2,]**order, na.rm = TRUE)
}

#' Perform the estimation of the variance \eqn{Var(X_{s}X_{t)}}.
#'
#' This function performs an estimation of the variance \eqn{Var(X_{s}X_{t)}}
#' used for the estimation of the bandwidth for a univariate kernel regression
#' estimator defined over continuous domains data.
#'
#' @family estimate variance
#' 
#' @param curves_smooth_s List, smoothing of the curves at point \eqn{s} and
#' resulting from the `presmoothing` function.
#' @param curves_smooth_t List, smoothing of the curves at point \eqn{t} and
#' resulting from the `presmoothing` function.
#'
#' @return Numeric, estimation of \eqn{Var(X_{s}X_{t)}}.
#' 
#' @references Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
#'  estimation of irregular mean and covariance functions.
#' @export
estimate_crossvar <- function(curves_smooth_s, curves_smooth_t) {
  current_smooth_s <- curves_smooth_s$x_smooth
  current_smooth_t <- curves_smooth_t$x_smooth
  stats::var(current_smooth_s[2,] * current_smooth_t[2,], na.rm = TRUE)
}
# ----

# Recursive estimation of the parameters ----

#' Perform a recursive estimation of the parameters.
#' 
#' This function performs a recursive estimation of the different parameters
#' used for the estimation of the mean and covariance estimation of functional
#' data. The recursion is made by small step onto the estimation of the
#' regularity of the curves. The pre-smoothing of the data is done using a
#' Nadaraya-Watson estimator and the used bandwidth modified using each new 
#' estimation of the regularity.
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  } 
#' @param point Numeric (default = 0.5), sampling point at which the data is 
#' pre-smoothed.
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#' @param beta Numeric (default = 1), pre-specified regularity of the curves to
#' start the recursion. The default value is 1, which correspond to at least one
#'time differentiable curves.
#'
#' @return List, with six entries:
#'  \itemize{
#'   \item \strong{$point} Time point where the smoothing has been done.
#'   \item \strong{$curves} Smoothed curves
#'   \item \strong{$H} Estimated regularity.
#'   \item \strong{$L} Estimated constant.
#'   \item \strong{$var} Estimated variance.
#'   \item \strong{$mom} Estimated \eqn{E(X^{2}_{t_0})}
#'  }
#' 
#' @references S. Golovkine, N. Klutchnikoff and V. Patilea (2021) - Adaptive 
#'  optimal estimation of irregular mean and covariance functions.
#' @export
estimate_parameters_recursion <- function(
    curves, point = 0.5, 
    delta_f = NULL, kernel_name = 'epanechnikov', beta = 1) {
  n_loop <- 0
  H_estim <- 0
  H_prev <- beta
  while (abs(H_prev - H_estim) > 0.1) {
    H_prev <- beta - 0.1 * n_loop
    curves_smooth <- presmoothing(
      curves, point, delta_f, kernel = kernel_name, beta = H_prev)
    H_estim <- estimate_regularity(curves_smooth)
    n_loop <- n_loop + 1
  }
  L_estim <- estimate_constant(curves_smooth, H_estim)
  var_estim <- estimate_var(curves_smooth)
  mom_estim <- estimate_moment(curves_smooth, 2)
  
  list(
    'point' = point,
    'curves' = curves_smooth,
    'H' = H_estim,
    'L' = L_estim,
    'var' = var_estim,
    'mom' = mom_estim
  )
}

#' Perform a recursive estimation of the parameters over a grid of points for
#' the estimation of the mean.
#' 
#' This function performs a recursive estimation of the different parameters
#' used for the estimation of the mean estimation of functional data. The
#' recursion is made by small step onto the estimation of the regularity of the
#' curves. The pre-smoothing of the data is done using a Nadaraya-Watson
#' estimator and the used bandwidth modified using each new estimation of the
#' regularity.
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  } 
#' @param grid Vector (default = c(0.25, 0.5, 0.75)), sampling points at which
#' the data is pre-smoothed.
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#' @param beta Numeric (default = 1), pre-specified regularity of the curves to
#' start the recursion. The default value is 1, which correspond to at least one
#' time differentiable curves.
#' @param compute_crossvar Boolean (default = FALSE), should
#' \eqn{Var(X_{s}X_{t)}} be computed for covariance estimation?
#'
#' @return Dataframe, with columns:
#'  \itemize{
#'   \item \strong{$point} Time point where the smoothing has been done.
#'   \item \strong{$curves} Smoothed curves.
#'   \item \strong{$H} Estimated regularity.
#'   \item \strong{$L} Estimated constant.
#'   \item \strong{$var} Estimated variance.
#'   \item \strong{$mom} Estimated \eqn{E(X^{2}_{t_0})}
#'  }
#' 
#' @references S. Golovkine, N. Klutchnikoff and V. Patilea (2021) - Adaptive 
#'  optimal estimation of irregular mean and covariance functions.
#' @export
estimate_parameters_mean <- function(
    curves, grid = c(0.25, 0.5, 0.75), delta_f = NULL, 
    kernel_name = 'epanechnikov', beta = 1){
  lapply(1:length(grid), function(idx){
    estimate_parameters_recursion(
      curves, point = grid[idx], delta_f = delta_f,
      kernel_name = kernel_name, beta = beta)
  }) |> 
    (\(x) do.call("rbind", x))() |> 
    as.data.frame()
}

#' Perform a recursive estimation of the parameters over a grid of points for
#' the estimation of the covariance.
#' 
#' This function performs a recursive estimation of the different parameters
#' used for the estimation of the covariance estimation of functional data. The
#' recursion is made by small step onto the estimation of the regularity of the #' curves. The pre-smoothing of the data is done using a Nadaraya-Watson
#' estimator and the used bandwidth modified using each new estimation of the
#' regularity.
#' 
#' @param curves List, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  } 
#' @param grid Vector (default = c(0.25, 0.5, 0.75)), sampling points at which
#' the data is pre-smoothed.
#' @param delta_f Function (default = NULL), function to determine the delta.
#' @param kernel_name String (default = 'epanechnikov'), the kernel used for the 
#' estimation:
#'  \itemize{
#'   \item epanechnikov
#'   \item uniform
#'   \item biweight
#'  }
#' @param beta Numeric (default = 1), pre-specified regularity of the curves to
#' start the recursion. The default value is 1, which correspond to at least one
#' time differentiable curves.
#' @param compute_crossvar Boolean (default = FALSE), should
#' \eqn{Var(X_{s}X_{t)}} be computed for covariance estimation?
#'
#' @return Dataframe, with columns:
#'  \itemize{
#'   \item \strong{$point} Time point where the smoothing has been done.
#'   \item \strong{$curves} Smoothed curves.
#'   \item \strong{$H} Estimated regularity.
#'   \item \strong{$L} Estimated constant.
#'   \item \strong{$var} Estimated variance.
#'   \item \strong{$mom} Estimated \eqn{E(X^{2}_{t_0})}
#'   \item \strong{$var_st} \eqn{Var(X_{s}X_{t)}}
#'  }
#' 
#' @references S. Golovkine, N. Klutchnikoff and V. Patilea (2021) - Adaptive 
#'  optimal estimation of irregular mean and covariance functions.
#' @export
estimate_parameters_covariance <- function(
    curves, grid = c(0.25, 0.5, 0.75), delta_f = NULL, 
    kernel_name = 'epanechnikov', beta = 1){
  params_estim <- estimate_parameters_mean(
    curves, grid = grid, delta_f = delta_f, 
    kernel_name = kernel_name, beta = beta)
  zz <- expand.grid(point_s = grid, point_t = grid) |> 
    merge(params_estim, 
          by.x = "point_s", by.y = "point", all.x = TRUE, 
          suffixes = c("", "_s"), sort = FALSE) |> 
    merge(params_estim, 
          by.x = "point_t", by.y = "point", all.x = TRUE, 
          suffixes = c("_s", "_t"), sort = FALSE)
  zz_upper <- zz[zz$point_t <= zz$point_s, ]
  zz_upper$var_st <- zz_upper |> apply(1, function(rows) {
    estimate_crossvar(rows$curves_s, rows$curves_t)
  })
  zz_upper[order(unlist(zz_upper$point_s), unlist(zz_upper$point_t)), ]
}


