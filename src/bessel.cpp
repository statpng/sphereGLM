#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double bessel_i_approx(double x, double nu) {
  // Approximation of modified Bessel function of the first kind
  double t = x / 3.75;
  double y = t * t;
  if (y <= 1) {
    return 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                                                        + y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))));
  } else {
    double z = 1.0 / x;
    return exp(x) / sqrt(x) * (0.39894228 + z * (0.01328592
                                                 + z * (0.00225319 + z * (-0.00157565 + z * (0.00916281
                                                                                             + z * (-0.02057706 + z * (0.02635537 + z * (-0.01647633
                                                                                                                                         + z * 0.00392377))))))));
  }
}