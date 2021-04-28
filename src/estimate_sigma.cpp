// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// estimate_sigma.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double estimateSigma(
    const List & curves // Curves list ($x and $t)
){
  // Get the number of curves
  arma::uword N = curves.length();
  List mycurve = curves[0];
  
  double sigma = 0;
  for(arma::uword n=0; n<N; n++){
    mycurve = curves[n];
    arma::vec x = mycurve["x"];
    arma::uword M_n = x.n_elem;
    
    double squared_diff = 0;
    for(arma::uword l=1; l<M_n; l++){
      squared_diff += std::pow(x(l) - x(l-1), 2);
    }
    sigma += (squared_diff / (2 * (M_n - 1)));
  }
  return std::pow((sigma / N), 0.5);
}
