# R with Cpp function to gather stats from input object

gather_stats <- function(input_object, # object from the graph_from_file function
                         G_k = 10,     # number of bins for the degree
                         G_b = 10      # number of bins for the number of mutual friends
                         ) {
    # graph matrix
    unique_graph <- as.matrix(input_object)

    # unique IDs of all the nodes
    all_nodes        <- unique(c(input_object[,1],input_object[,2][input_object[,2] != -1]))
    names(all_nodes) <- as.character(as.integer(all_nodes))

    # rename all_nodes
    unique_nodes <- as.vector(sort(all_nodes))
    nodes        <- c(1:length(unique_nodes))
    n_node       <- length(unique_nodes)

    # rename nodes of graph
    graph        <- rename_graph(unique_graph, unique_nodes, nodes)

    # number of distinct time-steps
    unique_timestep  <- sort(unique(input_object[,3]))
    T_unique         <- length(unique_timestep)


    #Adjecency matrix
    u           <- u_from_graph_cpp(graph, n_node)
    n_u         <- dim(u)[1]
    u_binary    <- create_u_binary(u)
    u_square    <- tcrossprod(u_binary)
    diag(u_square)  <- 0

    K_max       <- max(colSums(u_binary))[1] # maximum number of friends of a node
    B_max       <- max(u_square)      # maximum number of mutual friends between pairs of nodes

    # create bin vectors from 0 to K_max
    bin_result_k <- binning(K_max, G_k)
    G_k          <- bin_result_k$bin_num
    K_bin        <- bin_result_k$bin

    # create bin vectors from 0 to B_max
    bin_result_b <- binning(B_max, G_b)
    G_b          <- bin_result_b$bin_num
    B_bin        <- bin_result_b$bin

    # eliminate bin label
    B_bin_no_label <- as.vector(B_bin)
    K_bin_no_label <- as.vector(K_bin)

    #print(paste0("G_k in gather_stats before main_loop:",G_k))
    #print(paste0("G_b in gather_stats before main_loop:",G_b))
    ### Main loop ###

    field <- main_loop(graph, T_unique, unique_timestep, nodes, K_bin_no_label, G_k, B_bin_no_label, G_b)

    ###
    edge_indices <- (graph[,2] != -1)
    g <- graph_from_edgelist(graph[edge_indices,1:2], directed = FALSE)

    power_law <- fit_power_law(degree(g))
    clustering_coefficient <- transitivity(g)

    result <- list(T_unique = T_unique, field = field,
                   bin_result_b = bin_result_b, bin_result_k = bin_result_k,
                   K_max = K_max, B_max = B_max,
                   power_law = power_law, clustering_coefficient = clustering_coefficient)
    return(result)
}



