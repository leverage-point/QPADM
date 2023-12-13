#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]] 
vec update_beta (vec old_beta, vec old_eta, double gamma, double lambda, int m) {
  vec ind = old_beta + old_eta / gamma;
  vec pos = old_beta + old_eta / gamma - lambda / (m * gamma);
  vec neg = abs(old_beta + old_eta / gamma) - lambda / (m * gamma);
  mat zero(ind.n_elem,1);
  zero.zeros();
  return((ind>0)%max(pos,zero)-(ind<0)%max(neg,zero));
}

// [[Rcpp::export]]
vec update_r(mat x, vec y, vec old_u, vec old_beta_b, double tau, double gamma) {
  vec pos1 = old_u / gamma + y - x * old_beta_b - tau / gamma;
  vec pos2 = -1 * old_u / gamma - y + x * old_beta_b + (tau - 1) / gamma;
  mat zero(pos1.n_elem,1);
  zero.zeros();
  return(max(pos1,zero)-max(pos2,zero));
}

// [[Rcpp::export]]
vec update_beta_b(mat x, vec y, vec new_r, vec old_u, vec old_eta, vec new_beta, double gamma) {
  return(
    (eye(x.n_cols,x.n_cols) - trans(x) * (inv(eye(x.n_rows,x.n_rows) + x * trans(x))) * x) *
      (trans(x) * (y - new_r + old_u / gamma) - old_eta / gamma + new_beta)
  );
}

// [[Rcpp::export]]
vec update_u(mat x, vec y, vec old_u, vec new_beta_b, vec new_r, double gamma) {
  return(old_u + gamma * (y - x * new_beta_b - new_r));
}

// [[Rcpp::export]]
vec update_eta(vec old_eta, vec new_beta, vec new_beta_b, double gamma) {
  return(old_eta + gamma * (new_beta_b - new_beta));
}





















