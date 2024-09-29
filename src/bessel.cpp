#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Optimized BesselI function
Rcpp::NumericVector besselI_scaled_opt(Rcpp::NumericVector x, double nu, int k_max = 5) {
  int n = x.size();
  Rcpp::NumericVector result(n);
  
  double pi = 3.141592653589793238462643383280;
  
  for(int i = 0; i < n; i++) {
    double z = x[i] / nu;
    double sz = std::sqrt(1 + z*z);
    double t = 1 / sz;
    double eta = 1 / (sz + std::abs(z)) + std::log(z / (1 + sz));
    
    double d = 0;
    if (k_max >= 1) {
      double t2 = t * t;
      double u1_t = (t * (3 - 5 * t2)) / 24;
      
      if (k_max >= 2) {
        double u2_t = t2 * (81 + t2 * (-462 + t2 * 385)) / 1152;
        
        if (k_max >= 3) {
          double u3_t = t * t2 * (30375 + t2 * (-369603 + t2 * (765765 - t2 * 425425))) / 414720;
          
          if (k_max >= 4) {
            double t4 = t2 * t2;
            double u4_t = t4 * (4465125 + t2 * (-94121676 + t2 * (349922430 + t2 * (-446185740 + t2 * 185910725)))) / 39813120;
            
            if (k_max == 5) {
              double u5_t = t * t4 * (1519035525 + t2 * (-49286948607 + t2 * (284499769554 + t2 * (-614135872350 + t2 * (566098157625 - t2 * 188699385875))))) / 6688604160;
              d = (u1_t + (u2_t + (u3_t + (u4_t + u5_t/nu)/nu)/nu)/nu)/nu;
            } else {
              d = (u1_t + (u2_t + (u3_t + u4_t/nu)/nu)/nu)/nu;
            }
          } else {
            d = (u1_t + (u2_t + u3_t/nu)/nu)/nu;
          }
        } else {
          d = (u1_t + u2_t/nu)/nu;
        }
      } else {
        d = u1_t/nu;
      }
    }
    
    result[i] = (1 + d) * std::exp(nu * eta) / std::sqrt(2 * pi * nu * sz);
  }
  
  return result;
}

// [[Rcpp::export]]
double Cq_cpp(arma::vec theta, bool logarithm = true) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  double result;
  if (logarithm) {
    Rcpp::NumericVector bessel_result = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0 - 1);
    result = (q/2.0 - 1.0) * log(NORM) - ((q/2.0) * log(2*M_PI) +
      log(bessel_result[0]) +
      NORM);
  } else {
    Rcpp::NumericVector bessel_result = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0 - 1);
    result = exp((q/2.0 - 1.0) * log(NORM) - ((q/2.0) * log(2*M_PI) +
      log(bessel_result[0]) +
      NORM));
  }
  
  return result;
}

// [[Rcpp::export]]
double Bq_cpp(arma::vec theta) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  Rcpp::NumericVector bessel_result1 = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0);
  Rcpp::NumericVector bessel_result2 = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0 - 1);
  
  return bessel_result1[0] / bessel_result2[0];
}

// [[Rcpp::export]]
double Hq_cpp(arma::vec theta) {
  int q = theta.n_elem;
  double NORM = norm(theta);
  
  Rcpp::NumericVector bessel_result1 = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0);
  Rcpp::NumericVector bessel_result2 = besselI_scaled_opt(Rcpp::NumericVector::create(NORM), q/2.0 - 1);
  
  double I0 = bessel_result2[0];
  double I1 = bessel_result1[0];
  
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