#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
rowvec gaussianKernel(rowvec xprime, mat X) {
  
  int n = X.n_rows;
  double distance;
  rowvec output(n);
  
  for (int i = 0; i < n; i++){
    mat diff = xprime - X.row(i);
    distance = pow(norm(diff), 2);
    output(i) = exp(-0.5 * distance);
  }
  return output;
}

// [[Rcpp::export]]
void compute_kernel_z(rowvec & z, mat & K, mat X, int dim, int N){
  
  Rcpp::Environment package_env("package:mvtnorm");
  Rcpp::Function pmvnorm = package_env["pmvnorm"];
  mat identityMatrix(dim, dim, fill::eye);
  NumericVector lowerEnd(dim);
  NumericVector upperEnd(dim);
  upperEnd.fill(1);
  
  for (int i = 0; i < N; i++){
    K.row(i) = gaussianKernel(X.row(i), X);
    Rcpp::List probabilityList = pmvnorm(lowerEnd, upperEnd, X.row(i), identityMatrix);
    NumericVector probability = probabilityList[0];
    z(i) = probability[0] * pow(2.0 * PI, dim/2.0);
  }
}

// [[Rcpp::export]]
void calculate_candidate(rowvec & candidate_prob, rowvec & candidate_Var, rowvec z, NumericMatrix candidateSet, mat K_prime, mat X, double dim, double N, int p){
  
  Rcpp::Environment package_env_1("package:mvtnorm");
  Rcpp::Function pmvnorm = package_env_1["pmvnorm"];
  
  mat identityMatrix(dim, dim, fill::eye);
  NumericVector lowerEnd(dim);
  NumericVector upperEnd(dim);
  upperEnd.fill(1);
  
  for (int i = 0; i < 10; i++){
    Rcpp::List probabilityList = pmvnorm(lowerEnd, upperEnd, candidateSet.row(i), identityMatrix);
    NumericVector probability = probabilityList[0];
    candidate_prob(i) = probability[0] * pow(2.0 * PI, dim/2.0 );
    rowvec kernel = gaussianKernel(candidateSet(i, _), X);
    K_prime(span(0, N+p-2), N+p-1) = trans(kernel);
    K_prime(N+p-1, span(0, N+p-2)) = kernel;
    z(N+p-1) = candidate_prob(i);
    candidate_Var(i) = arma::as_scalar(z * pinv(K_prime) * trans(z));
  }
}