#include <RcppArmadillo.h>
#include <boost/math/special_functions/bessel.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double Cq_cpp(arma::vec theta, bool logarithm = true) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  double result;
  if (logarithm) {
    result = (q/2.0 - 1.0) * log(NORM) - ((q/2.0) * log(2*M_PI) + 
      log(boost::math::cyl_bessel_i(q/2.0 - 1, NORM, boost::math::policies::policy<boost::math::policies::promote_double<false>>())) + 
      NORM);
  } else {
    result = exp((q/2.0 - 1.0) * log(NORM) - ((q/2.0) * log(2*M_PI) + 
      log(boost::math::cyl_bessel_i(q/2.0 - 1, NORM, boost::math::policies::policy<boost::math::policies::promote_double<false>>())) + 
      NORM));
  }
  
  return result;
}

// [[Rcpp::export]]
double Bq_cpp(arma::vec theta) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  return boost::math::cyl_bessel_i(q/2.0, NORM, boost::math::policies::policy<boost::math::policies::promote_double<false>>()) / 
    boost::math::cyl_bessel_i(q/2.0 - 1, NORM, boost::math::policies::policy<boost::math::policies::promote_double<false>>());
}

// [[Rcpp::export]]
double Hq_cpp(arma::vec theta) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  double I0 = boost::math::cyl_bessel_i(q/2.0 - 1, NORM, boost::math::policies::policy<boost::math::policies::promote_double<false>>());
  double I1 = boost::math::cyl_bessel_i(q/2.0, NORM, boost::math::policies::policy<boost::math::policies::promote_double<false>>());
  
  return 1 - pow(I1/I0, 2) - (q-1)/NORM * I1/I0;
}

// [[Rcpp::export]]
arma::vec subgrad(arma::vec theta) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  if (NORM == 0) {
    arma::vec v = arma::randu(q) * 2 - 1;  // Uniform random values between -1 and 1
    return v / (norm(v) * 1.1);
  } else {
    return theta / NORM;
  }
}

// [[Rcpp::export]]
arma::vec b1_vMF(arma::vec theta) {
  return Bq_cpp(theta) * subgrad(theta);
}

// [[Rcpp::export]]
arma::mat b2_vMF(arma::vec theta) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  arma::vec s = subgrad(theta);
  
  return s * s.t() * Hq_cpp(theta) + Bq_cpp(theta) / NORM * (arma::eye(q, q) - theta * theta.t() / pow(NORM, 2));
}

// [[Rcpp::export]]
arma::mat diag_matrix(const arma::mat& X, double lambda) {
  return lambda * arma::eye(X.n_rows, X.n_cols);
}






// [[Rcpp::export]]
arma::mat calculate_b1_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list, const arma::vec& beta, int n) {
  int q = Offset.n_cols;
  arma::mat eta(n, q);
  
  for (int i = 0; i < n; ++i) {
    arma::mat Xt_i = Rcpp::as<arma::mat>(Xt_list[i]);
    arma::vec theta = Offset.row(i).t() + Xt_i.t() * beta;
    
    // b1_vMF 계산
    eta.row(i) = b1_vMF(theta).t();
  }
  
  return eta;
}



// [[Rcpp::export]]
arma::mat calculate_b2_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list, const arma::vec& beta, int n) {
  int q = Offset.n_cols;
  arma::mat W(n * q, n * q, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    arma::mat Xt_i = Rcpp::as<arma::mat>(Xt_list[i]);
    arma::vec theta = Offset.row(i).t() + Xt_i.t() * beta;
    
    // b2_vMF 계산 및 W 행렬에 할당
    W.submat(i*q, i*q, (i+1)*q-1, (i+1)*q-1) = b2_vMF(theta);
  }
  
  return W;
}