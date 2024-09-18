#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat diag_matrix(const arma::mat& X, double lambda);
arma::mat b2_vMF(arma::vec theta);
arma::vec b1_vMF(arma::vec theta);
arma::vec subgrad(arma::vec theta);
double Hq_cpp(arma::vec theta);
double Bq_cpp(arma::vec theta);
double Cq_cpp(arma::vec theta, bool logarithm = true);
arma::mat calculate_b1_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list, const arma::vec& beta, int n);
arma::mat calculate_b2_vMF(const arma::mat& Offset, const Rcpp::List& Xt_list, const arma::vec& beta, int n);

  
// [[Rcpp::export]]
List sphereGLM_iteration(arma::mat X, arma::mat Y, arma::mat Offset, arma::vec beta, 
                         arma::mat Xt, List Xt_list, double eps, int maxit, 
                         double lambda, bool orthogonal) {
  int n = X.n_rows;
  int p = X.n_cols;
  int q = Y.n_cols;
  
  // Rcpp::Rcout << "1" << std::endl;
  
  int l1 = 1;
  arma::vec beta_old = beta;
  arma::vec beta_new = beta + 1;
  List beta_list;
  List Fn_list;
  std::vector<double> loglik_list;
  std::vector<double> crit_list;
  double crit = 1;
  
  // Rcpp::Rcout << "2" << std::endl;
  
  while (crit > eps && l1 <= maxit) {
    beta_list.push_back(beta);
    beta_old = beta;
    
    // Rcpp::Rcout << "3" << std::endl;
    
    arma::mat eta(n, q);
    arma::mat W(n*q, n*q);
    arma::vec theta(q);
    
    eta = calculate_b1_vMF(Offset, Xt_list, beta, n);
    W = calculate_b2_vMF(Offset, Xt_list, beta, n);
    
    // Rcpp::Rcout << "4" << std::endl;
    
    arma::mat W_inv = inv(W + diag_matrix(W, lambda));
    
    // Rcpp::Rcout << "6 :" << W_inv << std::endl;
    
    arma::vec yt = vectorise(Y.t());
    arma::vec Z = Xt.t() * beta + W_inv * (yt - vectorise(eta.t()));
    
    // beta 업데이트
    arma::vec mu = beta.subvec(0, q-1);
    arma::mat XWX = Xt * W * Xt.t();
    
    if (orthogonal) {
      double gamma = 9999;
      arma::mat penalty = arma::zeros(XWX.n_rows, XWX.n_cols);
      for (int i = 1; i <= p; ++i) {
        penalty.submat(i*q, i*q, (i+1)*q-1, (i+1)*q-1) = gamma * mu * mu.t();
      }
      beta = solve(XWX + penalty, Xt * W * Z);
    } else {
      arma::mat penalty = arma::zeros(XWX.n_rows, XWX.n_cols);
      beta = solve(XWX + penalty, Xt * W * Z);
    }
    
    beta_new = beta;
    
    // loglik 계산
    double loglik = 0;
    for (int i = 0; i < n; ++i) {
      arma::mat Xt_i = Rcpp::as<arma::mat>(Xt_list[i]);
      
      theta = Offset.row(i).t() + Xt_i.t() * beta;
      loglik += arma::dot(theta, Y.row(i).t()) + Cq_cpp(theta, true);
    }
    
    crit = arma::norm(beta_old - beta_new, 2);
    
    loglik_list.push_back(loglik);
    crit_list.push_back(crit);
    Fn_list.push_back(XWX);
    
    l1++;
  }
  
  return List::create(
    Named("beta") = beta,
    Named("beta_list") = beta_list,
    Named("loglik_list") = loglik_list,
    Named("crit_list") = crit_list,
    Named("Fn_list") = Fn_list,
    Named("iterations") = l1 - 1
  );
}