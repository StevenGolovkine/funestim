################################################################################
#                     Utility functions for the package                        #
################################################################################

# Functions for the conversion of the data between the different formats. ----
#
# See C. Happ-Kurz (2020) Object-Oriented Software for Functional Data. 
# Journal of Statistical Software, 93(5): 1-38 .
# See Cai and Yuan, Nonparametric Covariance Function Estimation for Functional 
# and Longitudinal Data - 2010 (Technical Report).

#' Convert \code{funData::funData} objects into the right format.
#' 
#' @param data An object of the class \code{funData::funData}.
#' @param norm Boolean, if TRUE, the sampling points are normalized on 
#'  \eqn{[0, 1]}.
#' 
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  
#' @references C. Happ-Kurz (2020) Object-Oriented Software for Functional Data. 
#'  Journal of Statistical Software, 93(5): 1-38 .
#' @export
funData2list <- function(data, norm = TRUE){
  t <- funData::argvals(data)[[1]]
  x <- funData::X(data)
  
  data_list <- list()
  for (i in 1:funData::nObs(data)) {
    if (norm) t <- (t - min(t)) / (max(t) - min(t))
    data_list[[i]] <- list(t = t, x = x[i,])
  }
  
  data_list
}

#' Convert comprehensive lists into \code{funData::funData} objects.
#' 
#' We assume that we \strong{know} that the curves are on the same interval.
#'  
#' @param data_list A list, where each element represents a curve. Each curve is
#'  defined as list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'   
#' @return An object of the class \code{funData::funData}.
#' 
#' @references C. Happ-Kurz (2020) Object-Oriented Software for Functional Data. 
#'  Journal of Statistical Software, 93(5): 1-38 .
#' @export
list2funData <- function(data_list){
  argvals <- data_list[[1]]$t
  obs <- data_list |>
    purrr::map_dfc("x") |>
    as.matrix()
  funData::funData(argvals = argvals, X = t(obs))
}

#' Convert \code{funData::irregFunData} objects into the right format.
#' 
#' @param data An object of the class \code{funData::irregFunData}.
#' @param norm Boolean, if TRUE, the sampling points are normalized on 
#'  \eqn{[0, 1]}
#' 
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item
#'    \strong{$x} The observed points
#'  }
#'  
#' @references C. Happ-Kurz (2020) Object-Oriented Software for Functional Data. 
#'  Journal of Statistical Software, 93(5): 1-38.
#' @export
irregFunData2list <- function(data, norm = TRUE){
  t <- data@argvals
  x <- data@X
  
  data_list <- list()
  if (norm) {
    data_list <- purrr::map2(t, x, 
                             ~ list(t = (.x - min(.x)) / (max(.x) - min(.x)),
                                    x = .y))
  } else {
    data_list <- purrr::map2(t, x, ~ list(t = .x, x = .y))
  }
  
  data_list
}

#' Convert comprehensive lists into \code{funData::irregFunData} objects.
#'  
#' @param data_list A list, where each element represents a curve. Each curve is
#'  defined as list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'   
#' @return  An object of the class \code{funData::irregFunData}
#' 
#' @references C. Happ-Kurz (2020) Object-Oriented Software for Functional Data. 
#'  Journal of Statistical Software, 93(5): 1-38 .
#' @export
list2irregFunData <- function(data_list){
  argvalsList <- data_list |> purrr::map("t")
  obsList <- data_list |> purrr::map("x")
  funData::irregFunData(argvals = argvalsList, X = obsList)
}

#' Convert \code{funData::multiFunData} objects into the right format.
#' 
#' @param data An object of the class \code{funData::multiFunData}.
#' @param norm Boolean, if TRUE, the sampling points are normalized on 
#'  \eqn{[0, 1]}.
#'
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' 
#' @references C. Happ-Kurz (2020) Object-Oriented Software for Functional Data. 
#'  Journal of Statistical Software, 93(5): 1-38 .
#' @export
multiFunData2list <- function(data, norm = TRUE){
  
  data_list <- list()
  cpt <- 1
  for (fun_data in data) {
    if (inherits(fun_data, 'funData')) {
      data_list[[cpt]] <- funData2list(fun_data, norm)
    } else if (inherits(fun_data, 'irregFunData')) {
      data_list[[cpt]] <- irregFunData2list(fun_data, norm)
    } else if (inherits(fun_data, 'multiFunData')) {
      data_list[[cpt]] <- multiFunData2list(fun_data, norm)
    } else{
      stop('Something went wrong with one of the functional data object!')
    }
    cpt <- cpt + 1
  }
  data_list
}

#' Check the input data and return a list in the right format.
#' 
#' @param data An object from the package \code{funData}. It could be a 
#'  \code{funData::funData} or \code{funData::irregFunData} object.
#' 
#' @return A list, where each element represents a curve. Each curve is defined 
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' 
#' @references C. Happ-Kurz (2020) Object-Oriented Software for Functional Data. 
#'  Journal of Statistical Software, 93(5): 1-38 .
#' @export
checkData <- function(data){
  if (inherits(data, 'funData')) {
    data_ <- funData2list(data)
  } else if (inherits(data, 'irregFunData')) {
    data_ <- irregFunData2list(data)
  } else {
    stop('Wrong data type!')
  }
  data_
}

#' Convert \code{ssfcov} objects into the right format.
#' 
#' @param time Vector, observation times.
#' @param x Vector, observation values.
#' @param subject Vector, observation id.
#' 
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  
#' @references ssfcov, R package. See Cai and Yuan, Nonparametric Covariance
#'  Function Estimation for Functional and Longitudinal Data - 2010 (Technical
#'  Report).
#' @export
cai2list <- function(time, x, subject){
  results <- list()
  for (zz in unique(subject)) {
    results[[zz]] <- list(t = time[subject == zz], x = x[subject == zz])
  }
  results
}

#' Convert comprehensive lists into \code{ssfcov} objects.
#' 
#' @param data A list, where each element represents a curve. Each curve is 
#'  defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'
#' @return A dataframe where the columns are:
#'  \itemize{
#'   \item \strong{$time} The observation times.
#'   \item \strong{x} The observation values.
#'   \item \strong{subject} The observation id.
#'  }
#'  
#' @references ssfcov, R package. See Cai and Yuan, Nonparametric Covariance
#'  Function Estimation for Functional and Longitudinal Data - 2010 (Technical
#'  Report).
#' @export
list2cai <- function(data){
  data |> purrr::map_dfr(~ data.frame(time = .x$t, x = .x$x), .id = 'obs')
}
# ----

# Various utility functions ----

#' Logarithmic sequence generation.
#' 
#' @param from Numeric, the starting values of the sequence.
#' @param to Numeric, the end value of the sequence.
#' @param length.out Numeric, desired length of the sequence.
#' 
#' @return Vector, the logarithmic sequence.
#' @export
lseq <- function(from = 1, to = 100, length.out = 51) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#' Compute the values for some kernels.
#' 
#' @param u Vector, points to estimate the kernel.
#' @param type Integer, used kernel. 
#'  \itemize{
#'   \item \strong{1} Uniform kernel
#'   \item \strong{2} Epanechnikov kernel
#'   \item \strong{3} Biweight kernel
#'  }
#'  
#' @return Vector, the evaluated kernel.
#' @export
kernel <- function(u, type = 1){
  indicator <- function(u) 2 * stats::dunif(u, -1, 1)
  switch(type,
         indicator(u) / 2,
         0.75 * (1 - u**2) * indicator(u),
         0.9375 * (1 - u**2)**2 * indicator(u))
}

#' Test whether the observation points are in the neighborhood.
#' 
#' Test if there are more than \eqn{k_0} points, defined by \eqn{t}, in the 
#' neighborhood of \eqn{t_0}, defined by the bandwidth \eqn{h}.
#' 
#' @param t Vector, the sampling points.
#' @param t0 Numeric, the estimation point.
#' @param h Numeric, the bandwidth.
#' @param k0 Numeric, the number of points to be in the neighborhood.
#' 
#' @return Boolean, true if there are enough points in the neighborhood.
#' @export
neighbors <- function(t, t0, h, k0){
  sum(ifelse(abs(t - t0) <= h, 1, 0)) >= k0
}

#' Perform 2D linear approximation
#' 
#' Return a matrix of points which linearly interpolate given data points.
#' 
#' @param x Vector, giving the coordinates of the points to be interpolated.
#' @param y Matrix, giving the points to be interpolated.
#' @param xout Vector, values specifying where interpolation is to take place.
#' 
#' @return Matrix, interpolated data points.
#' @export
approx_2D <- function(x, y, xout){
  # Linear approximation for each rows
  rows <- matrix(0, ncol = length(xout), nrow = length(x))
  for (i in 1:nrow(y)) {
    rows[i, ] <- stats::approx(x, y[i, ],
                               xout = xout,
                               yleft = y[i, 1], 
                               yright = y[i, length(y[1, ])])$y  
  }
  # Linear approximation for each columns
  cols <- matrix(0, ncol = length(xout), nrow = length(xout))
  for (i in 1:ncol(rows)) {
    cols[, i] <- stats::approx(x, rows[, i], 
                               xout = xout,
                               yleft = rows[1, i], 
                               yright = rows[length(rows[, 1]), i])$y  
  }
  return(cols)
}
# ----
