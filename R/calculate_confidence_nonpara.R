# function to estimate A_k and B_b
calculate_confidence_nonpara <- function(result, nonpara_estim, A_k1)  {

  T_unique <- result$T_unique

  G_k      <- result$bin_result_k$bin_num
  K_center_bin    <- result$bin_result_k$center_bin

  G_b      <- result$bin_result_b$bin_num
  B_center_bin    <- result$bin_result_b$center_bin

  field    <- result$field

  A <- nonpara_estim$A
  B <- nonpara_estim$B


  # initial value of A and B
  A_index <- rep(0, G_k)
  B_index <- rep(0, G_b)
  for(i in 2:G_k){
    if(A[i] != 0){
      A_index[i] = 1
    }
  }
  for(i in 2:G_b){
    if(B[i] != 0){
      B_index[i] = 1
    }
  }

  if(A_k1 == F){
    A_index[2] <- 0
    A[2] <- 0
  }

  A[1] <- 0 # Fixed value
  B[1] <- 0 # Fixed value
  hessian <- create_hessian(field, T_unique, A, B, A_index, B_index, G_k, G_b)
  hessian_inv <- solve(hessian)
  sigma <- sqrt(diag(-hessian_inv))


  count <- 1
  sigma_A <- rep(0, G_k)
  sigma_B <- rep(0, G_b)

  for(i in 1:G_k){
    if(A_index[i] == 1){
      sigma_A[i] <- sigma[count]
      count <- count + 1
    } else{
      sigma_A[i] <- 0
    }
  }
  for(i in 1:G_b){
    if(B_index[i] == 1){
      sigma_B[i] <- sigma[count]
      count <- count + 1
    } else{
      sigma_B[i] <- 0
    }
  }
  A[1] <- nonpara_estim$A[1]
  B[1] <- nonpara_estim$B[1]

  if(A_k1 == F){
    A[2] <- nonpara_estim$A[2]
  }

  res <- list(bin_result_k = result$bin_result_k, bin_result_b = result$bin_result_b, A = A, B = B, sigma_A = sigma_A, sigma_B = sigma_B )
  return(res)
}
