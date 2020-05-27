contribution_pa_trans <- function (statistics,
                                        estimation_nonpara,
                                        estimation_para){
  alpha <- estimation_para$par[1]
  beta <- estimation_para$par[2]

  k_center <- statistics$bin_result_k$center_bin
  b_center <- statistics$bin_result_b$center_bin
  k_bin_num <- statistics$bin_result_k$bin_num
  b_bin_num <- statistics$bin_result_b$bin_num

  A <- estimation_nonpara$A
  if(A[length(A)] == 0){
    for(i in (length(A)-1):1){
      if(A[i] != 0){
        A[(i+1):length(A)] <- A[i]
        break
      }
    }
  }
  estimation_nonpara$A <- A

  B <- estimation_nonpara$B
  if(B[length(B)] == 0){
    for(i in (length(B)-1):1){
      if(B[i] != 0){
        B[(i+1):length(B)] <- B[i]
        break
      }
    }
  }
  estimation_nonpara$B <- B

  A_k <- function(n){
    if(estimation_nonpara$A[n] != 0){
      res <- estimation_nonpara$A[n]
    } else {
      res <- (k_center[n] + 1)^alpha
    }
  }
  B_b <- function(n){
    if(estimation_nonpara$B[n] != 0){
      res <- estimation_nonpara$B[n]
    } else {
      res <- (b_center[n] + 1)^beta
    }
  }

  T_max <- dim(statistics$field)[1]

  result <- list()

  for(t in 1:T_max){
    cube_snapshot <- statistics$field[t, 1][[1]]

    A_A <- c()
    B <- c()
    count_A_A <- 0
    count_B <- 0
    tmp <- 0

    ### denominator \Sigma A_{k_{i}(t)}A_{k_{j}(t)}B_{b_{ij}(t)}
    denominator <- 0
    for(i in 1:k_bin_num){
      for(j in i:k_bin_num){
        for(l in 1:b_bin_num){
          tmp <- cube_snapshot[i, j, l]
          if(tmp != 0){
            denominator <- denominator + tmp*A_k(i)*A_k(j)*B_b(l)
          }
        }
      }
    }


    E_pa <- 0
    E_trans <- 0

    for(i in 1:k_bin_num){
      for(j in i:k_bin_num){
        for(l in 1:b_bin_num){
          tmp <- cube_snapshot[i, j, l]
          if(tmp != 0){
            E_pa    <- E_pa + tmp*(A_k(i)*A_k(j)*B_b(l)/denominator)*log2(A_k(i)*A_k(j))
            E_trans <- E_trans + tmp*(A_k(i)*A_k(j)*B_b(l)/denominator)*log2(B_b(l))
          }
        }
      }
    }

    contribution_pa <-  0
    contribution_trans <-  0

    for(i in 1:k_bin_num){
      for(j in i:k_bin_num){
        for(l in 1:b_bin_num){
          tmp <- cube_snapshot[i, j, l]
          if(tmp != 0){
            contribution_pa <- contribution_pa + tmp*(A_k(i)*A_k(j)*B_b(l)/denominator)*(log2(A_k(i)*A_k(j)) - E_pa)^2
            contribution_trans <- contribution_trans + tmp*(A_k(i)*A_k(j)*B_b(l)/denominator)*(log2(B_b(l)) - E_trans)^2
          }
        }
      }
    }

    # print(var(A_A))
    # print(var(B))
    #print(t)

    contribution_pa <- sqrt(contribution_pa)
    contribution_trans <- sqrt(contribution_trans)

    result[[t]] <- list(contribution_pa = contribution_pa, contribution_trans = contribution_trans,
                        E_pa = E_pa, E_trans = E_trans, denominator = denominator)
  }


  contribution <- shape_data_name(result)
  return(contribution)
}


shape_data_name <- function (contribution){
  time_step_num <- length(contribution)

  contribution_A <- rep(0, time_step_num)
  contribution_B <- rep(0, time_step_num)
  ## target contribution
  for(time in 1:time_step_num){
    contribution_A[time] <- contribution[[time]]$contribution_pa
    contribution_B[time] <- contribution[[time]]$contribution_trans
  }


  contribution_result <- list(contribution_A = contribution_A, contribution_B = contribution_B)
}
