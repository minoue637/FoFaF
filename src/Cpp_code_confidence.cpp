//// Cpp functions 2020-2-25 Masaaki Inoue
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
//#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins("cpp11")]]
//// [[Rcpp::plugins(openmp)]]


// create hessian
// [[Rcpp::export("create_hessian")]]
arma::mat create_hessian(const List& field,
                         const int& T_unique,
                         const arma::vec&A,
                         const arma::vec&B,
                         const arma::vec&A_index,
                         const arma::vec&B_index,
                         const int& G_k,
                         const int& G_b) {

  arma::cube FoF_cube;
  arma::cube Z;

  int n_A = sum(A_index);
  int n_B = sum(B_index);
  int n_hes = n_A + n_B;
  arma::mat hessian(n_hes, n_hes);
  hessian.zeros();

  arma::vec A_est_index(n_A);
  arma::vec B_est_index(n_B);
  int count;

  count = 0;
  for(int i = 0; i < G_k; i++){
    if(A_index(i) == 1){
      A_est_index(count) = i;
      ++count;
    }
  }
  count = 0;
  for(int i = 0; i < G_b; i++){
    if(B_index(i) == 1){
      B_est_index(count) = i;
      ++count;
    }
  }

  arma::mat Z_kk;
  arma::vec Z_b;
  arma::vec C(G_k);
  arma::vec D(G_b);
  double E;
  double F;
  // double G;
  double m;

  for(int t = 1; t < T_unique; t++){
    C.zeros();
    D.zeros();
    E = 0;
    F = 0;
    // G = 0;

    FoF_cube = as<arma::cube>(field(t - 1, 0));
    Z        = as<arma::cube>(field(t, 1));

    Z_kk = sum(Z, 2);
    Z_b  = sum(sum(Z, 0), 1);
    m = accu(Z);

    // set C
    for(int k = 0; k < G_k; k++){
      if(k == 0){
        for(int i = 0; i < G_k; i++){
          for(int j = 0; j < G_b; j++){
            if(i != k){
              C(k) += FoF_cube(k, i, j)*A(i)*B(j);
            } else {
              C(k) += 2*FoF_cube(k, k, j)*A(k)*B(j);
            }
          }
        }
      } else if (k == (G_k-1)) {
        for(int i = 0; i < G_k; i++){
          for(int j = 0; j < G_b; j++){
            if(i != k){
              C(k) += FoF_cube(i, k, j)*A(i)*B(j);
            } else {
              C(k) += 2*FoF_cube(k, k, j)*A(k)*B(j);
            }
          }
        }
      } else {
        for(int i = 0; i < k; i++){
          for(int j = 0; j < G_b; j++){
            C(k) += FoF_cube(i, k, j)*A(i)*B(j);
          }
        }
        for(int i = (k+1); i < G_k; i++){
          for(int j = 0; j < G_b; j++){
            C(k) += FoF_cube(k, i, j)*A(i)*B(j);
          }
        }
        for(int j = 0; j < G_b; j++){
          C(k) += 2*FoF_cube(k, k, j)*A(k)*B(j);
        }
      }
    }

    // set D
    for(int b = 0; b < G_b; b++){
      for(int i = 0; i < G_k; i++){
        for(int j = i; j < G_k; j++){
          D(b) += FoF_cube(i, j, b)*A(i)*A(j);
        }
      }
    }

    // set E
    for(int i = 0; i < G_k; i++){
      for(int j = i; j < G_k; j++){
        for(int l = 0; l < G_b; l++){
          E += FoF_cube(i, j, l)*A(i)*A(j)*B(l);
        }
      }
    }

    int k_0;
    //int k_1;
    int b_0;
    //int b_1;

    for(int i = 0; i < n_hes; i++){
      if(i < n_A){
        k_0 = A_est_index(i);
        hessian(i, i) += (-1)*(sum(Z_kk.col(k_0))  + sum(Z_kk.row(k_0)))/pow(A(k_0), 2); //- m*(2*F*E - pow(C(k_0), 2))/pow(E,2);
      }
      if(i >= n_A){
        b_0 = B_est_index(i-n_A);
        hessian(i, i) += (-1)*Z_b(b_0)/pow(B(b_0), 2); //+ m*D(b_0)*D(b_0)/pow(E,2);
      }
    }

    /*
    for(int i = 0; i < n_hes; i++){
      for(int j = i; j < n_hes; j++){
        if(i < n_A && j < n_A) {
          k_0 = A_est_index(i);
          k_1 = A_est_index(j);
          for(int l = 0; l < G_b; l++){
            F += FoF_cube(k_0, k_1, l)*B(l);
          }

          // hessian
          if(k_0 == k_1){ // i = j
            hessian(i, j) += (-1)*(sum(Z_kk.col(k_0)) + sum(Z_kk.row(k_0)))/pow(A(k_0), 2) - m*(2*F*E - pow(C(k_0), 2))/pow(E,2);
          } else {

            hessian(i, j) += (-1)*m*(F*E - C(k_0)*C(k_1))/pow(E,2);
            hessian(j, i) += (-1)*m*(F*E - C(k_0)*C(k_1))/pow(E,2);

          }

        } else if (i >= n_A && j >= n_A) {
          b_0 = B_est_index(i-n_A);
          b_1 = B_est_index(j-n_A);
          // hessian
          if(b_0 == b_1){ // i = j
            hessian(i, j) += (-1)*Z_b(b_0)/pow(B(b_0), 2) + m*D(b_0)*D(b_0)/pow(E,2);
          } else {

            hessian(i, j) += m*D(b_0)*D(b_1)/pow(E,2);
            hessian(j, i) += m*D(b_0)*D(b_1)/pow(E,2);

          }

        } else {
          k_0 = A_est_index(i);
          b_0 = B_est_index(j-n_A);
          if(k_0 == 0){
            for(int l = (k_0+1); l < G_k; l++){
              G += FoF_cube(k_0, l, b_0)*A(l);
            }
            G += 2*FoF_cube(k_0, k_0, b_0)*A(k_0);
          } else if (k_0 == (G_k-1)) {
            for(int l = 0; l < k_0; l++){
              G += FoF_cube(l, k_0, b_0)*A(l);
            }
            G += 2*FoF_cube(k_0, k_0, b_0)*A(k_0);
          } else {
            for(int l = (k_0+1); l < G_k; l++){
              G += FoF_cube(k_0, l, b_0)*A(l);
            }
            for(int l = 0; l < k_0; l++){
              G += FoF_cube(l, k_0, b_0)*A(l);
            }
            G += 2*FoF_cube(k_0, k_0, b_0)*A(k_0);
          }


          // hessian
          hessian(i, j) += (-1)*m*(G*E - C(k_0)*D(b_0))/pow(E,2);
          hessian(j, i) += (-1)*m*(G*E - C(k_0)*D(b_0))/pow(E,2);
        }
      }
    }
    */
  }

  return hessian;

}
