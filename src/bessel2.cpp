
// // [[Rcpp::export]]
// double Cq_cpp(arma::vec theta, bool logarithm = true) {
//   int q = theta.n_elem;
//   double NORM = norm(theta);
// 
//   double result;
//   if (logarithm) {
//     result = (q/2.0 - 1.0) * log(NORM) - ((q/2.0) * log(2*M_PI) +
//                                             log(R::bessel_i(NORM, q/2.0 - 1, 2)) +
//                                             NORM);
//   } else {
//     result = exp((q/2.0 - 1.0) * log(NORM) - ((q/2.0) * log(2*M_PI) +
//                                                 log(R::bessel_i(NORM, q/2.0 - 1, 2)) +
//                                                 NORM));
//   }
// 
//   return result;
// }
// 
// 
// 
// 
// // [[Rcpp::export]]
// double Bq_cpp(arma::vec theta) {
//   int q = theta.n_elem;
//   double NORM = norm(theta);
// 
//   return R::bessel_i(NORM, q/2.0, 1) /
//     R::bessel_i(NORM, q/2.0 - 1, 1);
// }
// 
// 
// 
// 
// // [[Rcpp::export]]
// double Hq_cpp(arma::vec theta) {
//   int q = theta.n_elem;
//   double NORM = norm(theta);
// 
//   double I0 = R::bessel_i(NORM, q/2.0 - 1, 1);
//   double I1 = R::bessel_i(NORM, q/2.0, 1);
// 
//   return 1 - pow(I1/I0, 2) - (q-1)/NORM * I1/I0;
// }
// 
// 
// 