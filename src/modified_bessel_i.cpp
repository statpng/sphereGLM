#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>

// [[Rcpp::depends(BH)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector modified_bessel_i(double nu, NumericVector x) {
  int n = x.size();
  NumericVector result(n);
  
  for(int i = 0; i < n; i++) {
    result[i] = boost::math::cyl_bessel_i(nu, x[i]);
  }
  
  return result;
}