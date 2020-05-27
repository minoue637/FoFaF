//// Cpp functions 2020-2-25 Masaaki Inoue
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
//#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins("cpp11")]]
//// [[Rcpp::plugins(openmp)]]

// create log-likelihood
// [[Rcpp::export("create_LL")]]
double create_LL(const List& field,
                 const int& T_unique,
                 const int& G_k,
                 const int& G_b,
                 const arma::vec& A,
                 const arma::vec& B) {
  arma::cube FoF_cube;
  arma::cube Z;



  double LL = 0;
  double denominator;

  for(int t = 1; t < T_unique; t++){
    FoF_cube = as<arma::cube>(field(t - 1, 0));
    Z        = as<arma::cube>(field(t, 1));
    denominator = 0;
    for(int i = 0; i < G_k; i++){
      for(int j = i; j < G_k; j++){
        for(int l = 0; l < G_b; l++){
          if(Z(i, j, l) != 0){
            LL += Z(i, j, l)*log(A(i)*A(j)*B(l));
          }
          denominator += FoF_cube(i, j, l)*A(i)*A(j)*B(l);
        }
      }
    }
    LL += (-1)*accu(Z)*log(denominator);
  }
  return LL;
}

// update function of A_k
// [[Rcpp::export("create_update_A")]]
double create_update_A(const List& field,
                const int& T_unique,
                const int& G_k,
                const int& G_b,
                const arma::vec& A,
                const arma::vec& B,
                const int& k_0) {
  int k = k_0 - 1;
  arma::cube FoF_cube;
  arma::cube Z;

  double A_k;
  double numerator   = 0;
  double denominator = 0;
  double deno_n; // numerator of denominator
  double deno_d; // denominator of denominator

  for(int t = 1; t < T_unique; t++){
    FoF_cube = as<arma::cube>(field(t - 1, 0));
    Z        = as<arma::cube>(field(t, 1));

    deno_n = 0;
    deno_d = 0;

    for(int i = 0; i < G_k; i++){
      for(int l = 0; l < G_b; l++){
        if(i < k) {
          numerator += Z(i, k, l);
          deno_n += FoF_cube(i, k, l)*A(i)*B(l)/pow(A(k), 3);
        } else if (k < i) {
          numerator += Z(k, i, l);
          deno_n += FoF_cube(k, i, l)*A(i)*B(l)/pow(A(k), 3);
        } else {   // k = i
          numerator += 2*Z(k, k, l);
          deno_n += 2*FoF_cube(k, k, l)*A(k)*B(l)/pow(A(k), 3);
        }
        for(int j = i; j < G_k; j++) {
          deno_d += FoF_cube(i, j, l)*A(i)*A(j)*B(l);
        }
      }
    }
    denominator += accu(Z)*deno_n/deno_d;
  }
  // 0/0 = 1
  if(numerator == 0 || denominator == 0){
    A_k = 0;
  } else {
    A_k = sqrt(sqrt(numerator/denominator));
  }
  return A_k;
}

// update function of B_b
// [[Rcpp::export("create_update_B")]]
double create_update_B(const List& field,
                const int& T_unique,
                const int& G_k,
                const int& G_b,
                const arma::vec& A,
                const arma::vec& B,
                const int& b_0) {
  int b = b_0 - 1;
  arma::cube FoF_cube;
  arma::cube Z;

  double B_b;
  double numerator   = 0;
  double denominator = 0;
  double deno_n; // numerator of denominator
  double deno_d; // denominator of denominator

  for(int t = 1; t < T_unique; t++){
    FoF_cube = as<arma::cube>(field(t - 1, 0));
    Z        = as<arma::cube>(field(t, 1));

    deno_n = 0;
    deno_d = 0;

    for(int i = 0; i < G_k; i++){
      for(int j = i; j < G_k; j++){
        numerator += Z(i, j, b);
        deno_n    += FoF_cube(i, j, b)*A(i)*A(j)/B(b);
        for(int l = 0; l < G_b; l++){
          deno_d += FoF_cube(i, j, l)*A(i)*A(j)*B(l);
        }
      }
    }
    denominator += accu(Z)*deno_n/deno_d;
  }
  if(numerator == 0 || denominator == 0){
    B_b = 0;
  } else {
    B_b = sqrt(numerator/denominator);
  }
  return B_b;
}
