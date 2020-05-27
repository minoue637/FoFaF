estimation_nonparametric <- function(stats,
                                        A_k1 = T,
                                        epsilon = 10^(-8),
                                        q = 1,
                                        step_size = 0.5,
                                        alpha = 1,
                                        beta = 1,
                                        set_initial_parametric = T){
  
  if(set_initial_parametric == T){
    para <- estimation_parametric(stats, a = 1, b = 1)
    alpha_initial <- para$par[1]
    beta_initial <- para$par[2]
  } else {
    alpha_initial <- alpha
    beta_initial <- beta
  }
  
  if(A_k1 == T){
    nonpara <- estimation_k1(stats,
                                       epsilon,
                                       q,
                                       step_size,
                                       alpha_initial,
                                       beta_initial)
  } else {
    nonpara <- estimation_k2(stats,
                                       epsilon,
                                       q,
                                       step_size,
                                       alpha_initial,
                                       beta_initial)
  }
  result <- calculate_confidence_nonpara(stats, nonpara, A_k1)
  return(result)
}