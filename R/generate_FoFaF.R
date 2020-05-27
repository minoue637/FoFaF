# function to generate FoFaF network 
generate_FoFaF <- function(n_0            = 300,   # number of initial nodes
                           m              = 100,   # number of edges added at each time step
                           n              = 1,     # number of nodes added at each time step
                           T_step         = 300,   # number of time steps
                           alpha          = 1.2,   # parameter of PA function: (k+1)^alpha
                           beta           = 0.6    # parameter of Transitivity function: (b+1)^beta
){
  # adjacency matrix of all nodes
  # n_0 + (T_step - 1)*n: number of all nodes
  u    <- matrix(0,  n_0 + (T_step - 1)*n, n_0 + (T_step - 1)*n)
  
  
  
  ##### loop T_step times
  # T: Number of time-step
  for (t in 1:T_step){           
    if (t == 1){          # add only new_edges to graph edges
      ##### initial
      new_edges  <- matrix(0,n_0,3)
      
      ### create new edges
      for (i in 1:n_0){   # circle graph
        new_edges[i,1]  <- i
        
        if(i < n_0){
          new_edges[i,2]  <- i+1 
        } else {
          new_edges[i,2]  <- 1
        }
        new_edges[i,3]  <- t
      }
      
      ### graph initial
      graph <- new_edges
      
      ### set u_initial 
      u     <- u_from_graph(graph)
      
      
    } else {                       # add new_edges and new_isolated to graph edges
      ###### not initial
      new_edges      <- 0
      new_isolated   <- 0
      u_old          <- 0
      u_old_binary   <- 0
      u_old_square   <- 0
      degree_old     <- 0
      existing_nodes <- 0
      probability    <- 0
      
      u_old          <- u
      
      ### new edges
      n_old               <- dim(u_old)[1]
      u_old_binary        <- matrix(ifelse(u_old >= 1,1,0),nrow = n_old, ncol = n_old)
      u_old_square        <- u_old_binary %*% t(u_old_binary)
      diag(u_old_square)  <- 0
      degree_old          <- rowSums(u_old_binary)
      
      # create new
      # {1 ,..., (n_0 + (t - 2)*n)} are existing node ID at t (omitting new node)
      # probability: probability matrix of which node pair attain a new edge (i<j only)
      # denominator: sum of k_i^alpha*k_j^alpha*b_ij^beta of existing node
      existing_nodes      <- c(1:(n_0 + (t - 2)*n))
      probability         <- matrix(0,length(existing_nodes),length(existing_nodes))
      denominator         <- 0
      for(i in 1:(length(existing_nodes)-1)) {
        for(j in (i+1):length(existing_nodes)) {
          probability[i,j] <- ((degree_old[i]+1)^alpha)*((degree_old[j]+1)^alpha)*((u_old_square[i,j]+1)^beta)
          denominator      <- denominator + ((degree_old[i]+1)^alpha)*((degree_old[j]+1)^alpha)*((u_old_square[i,j]+1)^beta)
        }
      }
      probability  <- probability/denominator
      
      # choose m edges between existing nodes 
      run <- runif(m)
      for(i in 1:m){
        x <- 0
        y <- 0
        a <- run[i]
        b <- 0
        for(j in 1:(length(existing_nodes)-1)){
          for(k in (j+1):length(existing_nodes)){
            b  <- b + probability[j,k]
            if(a <= b){
              x <- j
              y <- k
              break
            }
          }
          if((x !=0) && (y !=0)){
            break
          }
        }
        
        if(i == 1){
          new_edges  <- c(x,y,t)
        } else {
          new_edges  <- rbind(new_edges, c(x,y,t))
        }
      }
      graph          <- rbind(graph, new_edges)
      
      
      ### new isolated node
      new_isolated  <- matrix(0,n,3)
      new_isolated[,1]  <- c((n_0 + (t - 2)*n + 1):(n_0 + (t - 1)*n))
      new_isolated[,2]  <- -1
      new_isolated[,3]  <- t
      graph             <- rbind(graph, new_isolated)
      
      
      u           <- 0
      u_replace   <- 0
      ok_indices  <- 0
      ok_indices  <- graph[,2] != -1
      u           <- matrix(0,  n_0 + (T_step - 1)*n, n_0 + (T_step - 1)*n)
      u_replace   <- u_from_graph(graph)
      u[1:dim(u_replace)[1], 1:dim(u_replace)[1]]   <- u_replace
    }
  }
  return(graph)
}

u_from_graph <- function(graph) {
  edge_indices  <- graph[,2] != -1
  iso_indices   <- graph[,2] == -1
  g_edge        <- as.matrix(graph[edge_indices,1:2, drop = FALSE])
  g_iso         <- as.matrix(graph[iso_indices,1, drop = FALSE])
  
  v  <-  c(c(g_edge),g_iso)
  w  <-  sort(unique(v))
  
  u  <- matrix(0,length(w),length(w))
  
  for(i in 1:dim(g_edge)[1]){
    x       <-  g_edge[i,1]
    y       <-  g_edge[i,2]
    x_new   <-  which(w == x)
    y_new   <-  which(w == y)
    u[x_new,y_new]  <- u[x_new,y_new] + 1
    u[y_new,x_new]  <- u[y_new,x_new] + 1
  }
  return(u)
}
