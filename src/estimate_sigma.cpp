// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// estimate_sigma.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double estimateSigma(
    const List & curves, // Curves list ($x and $t)
    const double & delta // Neighborhood for the estimation
){
  // Get the number of curves
  arma::uword N = curves.length();
  
  double sigma = 0;
  for(arma::uword n=0; n<N; n++){
    List mycurve = curves[n];
    arma::vec x = mycurve["x"];
    arma::vec t = mycurve["t"];
    arma::uword M_n = x.n_elem;
    
    double squared_diff = 0;
    double cpt = 0;
    for(arma::uword l=1; l<M_n; l++){
      if (std::abs(t(l) - t(l-1)) <= delta){
        squared_diff += std::pow(x(l) - x(l-1), 2);
        cpt += 1;
      }
    }
    sigma += (squared_diff / (2 * cpt));
  }
  return std::pow((sigma / N), 0.5);
}
