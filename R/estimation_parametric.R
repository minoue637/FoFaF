# function to estimate alpha and beta
estimation_parametric <- function(result, a = 1, b = 1, n_thread = 1)  {
  # a,b : initial value of alpha, beta of optim()
  T_unique <- result$T_unique

  G_k      <- result$bin_result_k$bin_num
  K_center_bin    <- result$bin_result_k$center_bin
  #print(paste0("G_k outside:",G_k))

  G_b      <- result$bin_result_b$bin_num
  B_center_bin    <- result$bin_result_b$center_bin
  #print(paste0("G_b outside: ",G_b))

  field    <- result$field

  #print(paste("Size of field outside log_likelihood: " , paste(dim(field[1,1][[1]]))))

  # x : c(alpha, beta)
  LL <- function (x) {
    return(create_log_likelihood(field, T_unique, K_center_bin, B_center_bin, G_k, G_b, n_thread, x[1], x[2]))
  }

  grad <- function (x) {
    return(create_gradient(field, T_unique, K_center_bin, B_center_bin, G_k, G_b, n_thread, x[1], x[2]))
  }


  #res <- optim(par = c(a,b), fn = LL, control=list(fnscale=-1))
  res <- optim(par = c(a,b), fn = LL, gr = grad, control=list(fnscale=-1), method = "BFGS",
               hessian = TRUE)

  return(res)
}
