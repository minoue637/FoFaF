# function to estimate A_k and B_b
estimation_k1 <- function(result,
                                     epsilon,
                                     q = 1,
                                     step_size = 0.5,
                                     alpha,
                                     beta)  {
  # a,b : initial value of alpha, beta of optim()
  T_unique <- result$T_unique

  G_k      <- result$bin_result_k$bin_num
  K_center_bin    <- result$bin_result_k$center_bin

  G_b      <- result$bin_result_b$bin_num
  B_center_bin    <- result$bin_result_b$center_bin

  field    <- result$field

  LL <- function (A, B) {
    return(create_LL(field, T_unique, G_k, G_b, A, B))
  }

  # initial value of A and B
  A_initial <- rep(1, G_k)
  B_initial <- rep(1, G_b)
  for(i in 2:G_k){
    A_initial[i] <- (K_center_bin[i]+1)^(alpha)
  }
  for(i in 2:G_b){
    B_initial[i] <- (B_center_bin[i]+1)^beta
  }


  update_A <- function (A, B, i) {
    return(create_update_A(field, T_unique, G_k, G_b, A, B, i))
  }
  update_B <- function (A, B, i) {
    return(create_update_B(field, T_unique, G_k, G_b, A, B, i))
  }

  count        <- 0
  log_likelihood <- LL(A_initial, B_initial)

  A   <- A_initial
  B   <- B_initial

  if (q > 1) {
    parameter_save <- matrix(0,nrow = length(A) + length(B), ncol = q)
    U              <- matrix(0,nrow = length(A) + length(B), ncol = q - 1)
    V              <- matrix(0,nrow = length(A) + length(B), ncol = q - 1)
    candidate      <- rep(0,length(A) + length(B))
  }

  while(1){
    log_likelihood_old <- log_likelihood
    # acceleration
    if(q > 1){
      flag <- 0
      if (count > q) {
        current_pos <- c(A,B)
        U <- parameter_save[,1:(q-1)] - parameter_save[,q]
        V <- parameter_save[,2:q]     - current_pos
        candidate <- tryCatch(current_pos - step_size * ifelse(rep(q > 2,length(current_pos)), V%*%solve(crossprod(U,U) - crossprod(U,V),crossprod(U, V[,q-1])), 1/(sum(U* U) - sum(U*V)) * sum(V*U) * V),
                              error = function(e) {
                                flag <- 1; return(current_pos)})
        positive <- prod(candidate > 0)

        if ((0 == flag) && (1 == positive)) {
          #candidate_ok <- candidate_ok + 1
          A_temp        <- candidate[1:length(A)]
          B_temp        <- candidate[(length(A) + 1):(length(A) + length(B))]
          #offset_temp                <- candidate[length(candidate)]

          LL_candidate <- LL(A_temp, B_temp)
          if (LL_candidate > log_likelihood) {
            print(paste("acceleration:", LL_candidate))
            A      <- candidate[1 : length(A)]
            B      <- candidate[(length(A) + 1) : (length(A) + length(B))]
            #offset <- candidate[length(candidate)]
          }
        }
        parameter_save[, 1:(q-1)]         <- parameter_save[, 2:q]
        parameter_save[1:length(A), q]    <- A
        parameter_save[(length(A) + 1):(length(A) + length(B)), q] <- B
        #parameter_save[dim(parameter_save)[1],q] <- offset
      } else {
        parameter_save[1:length(A),count]    <- A
        parameter_save[(length(A) + 1):(length(A) + length(B)),count] <- B
        #parameter_save[dim(parameter_save)[1],count] <- offset
      }
    }

    # update A_k, B_b
    A_old <- A
    B_old <- B
    for(i in 2:G_k){
      A[i] <- update_A(A_old, B_old, i)
    }
    for(i in 2:G_b){
      B[i] <- update_B(A_old, B_old, i)
    }

    log_likelihood <- LL(A, B)

    count <- count + 1
    if(abs(log_likelihood - log_likelihood_old)/(abs(log_likelihood_old) + 1) <= epsilon){
      break
    }
  }

  res <- list(A = A, B = B,
              bin_result_k = result$bin_result_k, bin_result_b = result$bin_result_b,
              epsilon = epsilon, iteration = count)
  return(res)
}
