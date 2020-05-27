//// Cpp functions 2020-2-25 Masaaki Inoue
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
//#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins("cpp11")]]
// // [[Rcpp::plugins(openmp)]]

// graph: graph with new node ID (not unique node ID)
// u: Adjacency matrix, u_square: mutual friend matrix, d: degree vector
// u_label: label of u, u_square, degree
// B_bin: bin_vector, G_b: bin_num, K_bin: bin_vector, G_k: bin_num
// new_edges: new edges at t (matrix)
// new_isolated: isolated nodes data at t (vector), isolated_nodes: isolated nodes in graph at t (vector)
// is_new_nodes: whether a node is new or not (vector)
// FoF_cube: friends of friends (arma::cube), Z: new edge (arma::cube)

arma::mat create_new_edges(const int& time, const arma::mat& graph);
arma::vec create_new_isolated(const int& time, const arma::mat&  graph);
arma::vec create_is_appear_initial(const arma::mat& new_edges, const arma::vec& new_isolated, const arma::vec& all_nodes);
arma::vec create_u_label_initial(const arma::mat& new_edges, const arma::vec& all_nodes);
arma::sp_mat create_u_initial(const arma::mat& new_edges, const arma::vec& u_label);
arma::vec create_isolated_nodes(const arma::vec& u_label_initial, const arma::vec& is_appear_initial);
arma::vec update_is_appear(const arma::vec& is_appear_old, const arma::vec& new_isolated, const arma::mat& new_edges);
arma::vec create_is_new_nodes(const arma::vec& is_appear_old, const arma::vec& is_appear_new);
arma::vec update_isolated_nodes(const arma::vec& is_new_nodes, const arma::mat& new_edges, const arma::vec& isolated_nodes_old);
arma::sp_mat create_u_binary(const arma::sp_mat& u);
arma::sp_mat create_u_square(const arma::sp_mat& u_binary);
arma::vec create_degree(const arma::sp_mat& u_binary);
arma::cube create_FoF_cube(
    const arma::sp_mat& u_square, //isn't contain isolated nodes
    const arma::vec& u_label,  //isn't contain isolated nodes
    const arma::vec& degree,   //isn't contain isolated nodes
    const arma::vec& isolated_nodes,
    const arma::vec& K_bin,
    const int& G_k,
    const arma::vec& B_bin,
    const int& G_b);
arma::vec update_u_label(
    const arma::vec& u_label_old,
    const arma::vec& is_appear_old,
    const arma::vec& is_appear_new,
    const arma::vec& isolated_nodes_old,
    const arma::vec& isolated_nodes_new);
arma::sp_mat update_u(const arma::sp_mat& u_old,
                const arma::vec& u_label_new,
                const arma::mat& new_edges);
arma::cube create_Z(const arma::sp_mat& u_square_old,
              const arma::vec& u_label_old,
              const arma::vec& degree_old,
              const arma::mat& new_edges,
              const arma::vec& is_appear_old,
              const arma::vec& K_bin,
              const int& G_k,
              const arma::vec& B_bin,
              const int& G_b);

arma::cube return_cube(const List& field, const int& x, const int& y);



arma::vec create_gradient(const List&     field,
                          const int&      T_unique,
                          const arma::vec&      K_center_bin,
                          const arma::vec&      B_center_bin,
                          const int&      G_k,
                          const int&      G_b,
                          const int&  n_thread,
                          const double&   alpha,
                          const double&   beta);

double create_log_likelihood(SEXP  my_field,
                             const int&   T_unique,
                             const arma::vec&   K_center_bin,
                             const arma::vec&   B_center_bin,
                             const int&   G_k,
                             const int&   G_b,
                             const int&  n_thread,
                             const double&  alpha,
                             const double&  beta);

arma::mat rename_graph(const arma::mat& graph_unique,    // graph with unique node ID
                       const arma::vec& unique_nodes,    // unique node ID
                       const arma::vec& nodes);

arma::sp_mat u_from_graph_cpp(const arma::mat& graph, const int& n_node) ;


///////////////////////////////////////////////////////////////////////////

///////////////////////////////////main////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

// main loop to gather FoF_cube and Z at each timestep

// [[Rcpp::export("main_loop")]]
SEXP main_loop(const arma::mat& graph,
                      const int& T_unique,
                      const arma::vec& unique_timestep,
                      const arma::vec& all_nodes,
                      const arma::vec& K_bin,
                      const int& G_k,
                      const arma::vec& B_bin,
                      const int& G_b) {
  arma::mat new_edges;
  arma::vec new_isolated;

  arma::vec is_appear;
  arma::sp_mat u;
  arma::vec u_label;
  arma::vec isolated_nodes;
  arma::sp_mat u_binary;
  arma::sp_mat u_square;
  arma::vec degree;
  arma::cube FoF_cube;
  arma::cube Z;
  arma::field<arma::cube> F(T_unique,2);

  arma::vec is_appear_old;
  arma::vec isolated_nodes_old;
  arma::sp_mat u_old;
  arma::vec u_label_old;
  arma::sp_mat u_square_old;
  arma::vec degree_old;

  arma::vec is_appear_new;
  arma::vec is_new_nodes;
  arma::vec isolated_nodes_new;
  arma::sp_mat u_new;
  arma::vec u_label_new;
  arma::sp_mat u_square_new;
  arma::vec degree_new;

  //std::cout << "G_k and G_b inside main_loop: " << G_k << " " << G_b << "\n";


  for(int t = 1; t < T_unique + 1; t++){
    new_edges = create_new_edges(unique_timestep(t - 1), graph);
    new_isolated = create_new_isolated(unique_timestep(t - 1), graph);

    if(t == 1){            // initial

      is_appear         = create_is_appear_initial(new_edges, new_isolated, all_nodes);
      u_label           = create_u_label_initial(new_edges, all_nodes);
      u                 = create_u_initial(new_edges, u_label);
      isolated_nodes    = create_isolated_nodes(u_label, is_appear);
      u_binary          = create_u_binary(u);
      u_square          = create_u_square(u_binary);
      degree            = create_degree(u_binary);

      FoF_cube          = create_FoF_cube(u_square, u_label, degree, isolated_nodes, K_bin, G_k, B_bin, G_b);

      F(t - 1,0) = FoF_cube;
      //std::cout << "Size of FOF_cube inside main_loop at t = 1: " << FoF_cube.n_rows << " " << FoF_cube.n_cols << " " << FoF_cube.n_slices << "\n";

      //std::cout << "t = " << t << "\n";
    } else {                //not initial

      // set old
      is_appear_old       = is_appear;
      isolated_nodes_old  = isolated_nodes;
      u_old               = u;
      u_label_old         = u_label;
      u_square_old        = u_square;
      degree_old          = degree;

      // update
      is_appear_new       = update_is_appear(is_appear_old, new_isolated, new_edges);
      is_new_nodes        = create_is_new_nodes(is_appear_old, is_appear_new);
      isolated_nodes_new  = update_isolated_nodes(is_new_nodes, new_edges, isolated_nodes_old);
      u_label_new         = update_u_label(u_label_old, is_appear_old, is_appear_new, isolated_nodes_old,isolated_nodes_new);
      u_new               = update_u(u_old, u_label_new, new_edges);
      u_binary            = create_u_binary(u_new);
      u_square_new        = create_u_square(u_binary);
      degree_new          = create_degree(u_binary);


      // desired statistics
      Z                  = create_Z(u_square_old, u_label_old, degree_old, new_edges, is_appear_old, K_bin, G_k, B_bin, G_b);
      FoF_cube           = create_FoF_cube(u_square_new, u_label_new, degree_new, isolated_nodes_new, K_bin, G_k, B_bin, G_b);

      // set new
      is_appear       = is_appear_new;
      isolated_nodes  = isolated_nodes_new;
      u               = u_new;
      u_label         = u_label_new;
      u_square        = u_square_new;
      degree          = degree_new;

      F(t - 1,0) = FoF_cube;
      F(t - 1,1) = Z;
      //std::cout << "Size of FOF_cube inside main_loop: " << FoF_cube.n_rows << " " << FoF_cube.n_cols << " " << FoF_cube.n_slices << "\n";
      //std::cout << "Size of Z inside main_loop: " << Z.n_rows << " " << Z.n_cols << " " << Z.n_slices << "\n";
      //std::cout << "t = " << t << "\n";
    }

  }
  return Rcpp::wrap(F);
}




///////////////////////////////////////////////////////////////////////////

/////////////////////////////////Functions/////////////////////////////////

///////////////////////////////////////////////////////////////////////////


// all time
// [[Rcpp::export("create_new_edges")]]
arma::mat create_new_edges(const int& time,
                           const arma::mat& graph) {
  arma::mat new_edges;
  arma::vec check1 = graph.cols(1,1);
  arma::vec check2 = graph.cols(2,2);
  int count = 0;
  for(int i = 0; i < graph.n_rows; i++){

    if(check1(i) != -1 && check2(i) == time){
      arma::vec edge = trans(graph.rows(i,i));
      new_edges.insert_rows(count, trans(edge.rows(0, 1)));
      count = count + 1;
    }
  }
  if(count == 0){
    new_edges = zeros(2,2);
  }
  if(count == 1){
    arma::vec x = zeros(2);
    new_edges.insert_rows(1, trans(x));
  }
  return new_edges;
}

// all time
// [[Rcpp::export("create_new_isolated")]]
arma::vec create_new_isolated(const int& time,
                              const arma::mat& graph) {
  arma::vec new_isolated = zeros(2);
  arma::vec check1 = graph.cols(1,1);
  arma::vec check2 = graph.cols(2,2);
  int count = 0;
  for(int i = 0; i < graph.n_rows; i++){
    if(check1(i) == -1 && check2(i) == time){
      int node = graph(i,0);
      if(count >= 2){
        new_isolated.resize(new_isolated.n_elem + 1);
      }
      new_isolated(count) = node;
      count = count + 1;
    }
  }
  return new_isolated;
}


// initial only
// [[Rcpp::export("create_is_appear_initial")]]
arma::vec create_is_appear_initial(const arma::mat& new_edges,
                                   const arma::vec& new_isolated,
                                   const arma::vec& all_nodes) {
  arma::vec is_appear(size(all_nodes));
  is_appear.zeros();
  arma::vec check1 = new_edges.cols(0,0);
  arma::vec check2 = new_edges.cols(1,1);
  arma::vec check3 = new_isolated;
  for(int i = 0; i < new_edges.n_rows; i++){
    int node1 = check1(i);
    int node2 = check2(i);
    if(node1 == 0){
      break;
    }
    if(is_appear(node1 - 1) == 0){
      is_appear(node1 - 1) = 1;
    }
    if(is_appear(node2 - 1) == 0){
      is_appear(node2 - 1) = 1;
    }
  }
  for(int i = 0; i < new_isolated.n_elem; i++){
    int node3 = check3(i);
    if(node3 == 0){
      break;
    }
    if(is_appear(node3 - 1) == 0){
      is_appear(node3 - 1) = 1;
    }
  }
  return is_appear;
}


// initial only
// [[Rcpp::export("create_u_label_initial")]]
arma::vec create_u_label_initial(const arma::mat& new_edges,
                                 const arma::vec& all_nodes) {
  arma::vec node_initial(size(all_nodes));
  node_initial.zeros();
  arma::vec check1 = new_edges.cols(0,0);
  arma::vec check2 = new_edges.cols(1,1);

  // create node_initial
  for(int i = 0; i < new_edges.n_rows; i++){
    int node1 = check1(i);
    int node2 = check2(i);
    if(node1 == 0){
      break;
    }
    if(node_initial(node1 - 1) == 0){
      node_initial(node1 - 1) = 1;
    }
    if(node_initial(node2 - 1) == 0){
      node_initial(node2 - 1) = 1;
    }
  }

  // create u_label
  int num = sum(node_initial);
  arma::vec u_label(num);
  u_label.zeros();

  int count = 0;
  for(int i = 0; i < node_initial.n_elem; i++){       // check all initial nodes
    if(node_initial(i) == 1){
      u_label(count) = i + 1;
      count = count + 1;
    }
  }

  return u_label;
}

// initial only
// [[Rcpp::export("create_u_initial")]]
arma::sp_mat create_u_initial(const arma::mat& new_edges,
                              const arma::vec& u_label) {
  int num = u_label.n_elem;
  arma::sp_mat u_initial(num,num);
  u_initial.zeros();

  // create u_initial
  for(int i = 0; i < new_edges.n_rows; i++){       // count all edges in new_edges

    int node1 = new_edges(i,0);
    int node2 = new_edges(i,1);
    int a;
    int b;
    if(node1 == 0){
      break;
    }
    for(int j = 0; j < num; j++){
      if(u_label(j) == node1){
        a = j;
        for(int k = 0; k < num; k++){
          if(u_label(k) == node2){
            b = k;
            u_initial(a,b) = u_initial(a,b) + 1;
            u_initial(b,a) = u_initial(b,a) + 1;
            break;
          }
        }
        break;
      }
    }
  }

  return u_initial;
}

// initial only
// remove not isolated nodes from is_appear_initial
// [[Rcpp::export("create_isolated_nodes")]]
arma::vec create_isolated_nodes(const arma::vec& u_label_initial,
                                const arma::vec& is_appear_initial) {
  arma::vec isolated_nodes = is_appear_initial;
  for(int i = 0; i < u_label_initial.n_elem; i++){
    int check_node = u_label_initial(i);
    for(int j = 0; j < isolated_nodes.n_elem; j++){
      if(isolated_nodes(check_node - 1) == 1){
        isolated_nodes(check_node - 1) = 0;
      }
    }
  }
  return isolated_nodes;
}



// not initial
// [[Rcpp::export("update_is_appear")]]
arma::vec update_is_appear(const arma::vec& is_appear_old,
                           const arma::vec& new_isolated,
                           const arma::mat& new_edges) {
  arma::vec is_appear_new = is_appear_old;
  arma::vec check1 = new_edges.cols(0,0);
  arma::vec check2 = new_edges.cols(1,1);
  arma::vec check3 = new_isolated;
  for(int i = 0; i < new_edges.n_rows; i++){
    int node1 = check1(i);
    if(node1 == 0){
      break;
    }
    if(is_appear_new(node1 - 1) == 0){
      is_appear_new(node1 - 1) = 1;
    }
  }
  for(int i = 0; i < new_edges.n_rows; i++){
    int node2 = check2(i);
    if(node2 == 0){
      break;
    }
    if(is_appear_new(node2 - 1) == 0){
      is_appear_new(node2 - 1) = 1;
    }
  }
  for(int i = 0; i < new_isolated.n_elem; i++){
    int node3 = check3(i);
    if(node3 == 0){
      break;
    }
    if(is_appear_new(node3 - 1) == 0){
      is_appear_new(node3 - 1) = 1;
    }
  }
  return is_appear_new;
}


// not initial
// [[Rcpp::export("create_is_new_nodes")]]
arma::vec create_is_new_nodes(const arma::vec& is_appear_old,
                             const arma::vec& is_appear_new) {
  arma::vec is_new_nodes = is_appear_new - is_appear_old;
  return is_new_nodes;
}

// not initial
// [[Rcpp::export("update_isolated_nodes")]]
arma::vec update_isolated_nodes(const arma::vec& is_new_nodes,
                          const arma::mat& new_edges,
                          const arma::vec& isolated_nodes_old) {
  arma::vec isolated_nodes_new = isolated_nodes_old;
  for(int i = 0; i < is_new_nodes.n_elem; i++){
    if(is_new_nodes(i) == 1){
      isolated_nodes_new(i) = 1;
    }
  }
  arma::vec check1 = new_edges.cols(0,0);
  arma::vec check2 = new_edges.cols(1,1);
  for(int i = 0; i < new_edges.n_rows; i++){
    double node1 = check1(i);
    double node2 = check2(i);
    if(node1 == 0){
      break;
    }
    if(isolated_nodes_new(node1 - 1) == 1){
      isolated_nodes_new(node1 - 1) = 0;
    }
    if(isolated_nodes_new(node2 - 1) == 1){
      isolated_nodes_new(node2 - 1) = 0;
    }
  }
  return isolated_nodes_new;
}

// all time
// [[Rcpp::export("create_u_binary")]]
arma::sp_mat create_u_binary(const arma::sp_mat& u) {
  arma::sp_mat u_binary = arma::sp_mat(size(u));
  for(int i = 0; i < u.n_rows; i++){
    for(int j = 0; j < u.n_cols; j++){
      if(u(i, j) != 0){
        u_binary(i, j) = 1;
      }
    }
  }
  return u_binary;
}

// all time
// [[Rcpp::export("create_u_square")]]
arma::sp_mat create_u_square(const arma::sp_mat& u_binary) {
  arma::sp_mat u_square = u_binary * trans(u_binary);
  return u_square;
}

// all time
// [[Rcpp::export("create_degree")]]
arma::vec create_degree(const arma::sp_mat& u_binary) {
  arma::vec degree = zeros(u_binary.n_rows);
  for(int i = 0; i < degree.n_rows; i++){
    for(int j = 0; j < u_binary.n_rows; j++){
      degree(i) += u_binary(i, j);
    }
  }
  return degree;
}
// initial only
// [[Rcpp::export("create_FoF_cube")]]
arma::cube create_FoF_cube(
    const arma::sp_mat& u_square, //isn't contain isolated nodes
    const arma::vec& u_label,  //isn't contain isolated nodes
    const arma::vec& degree,   //isn't contain isolated nodes
    const arma::vec& isolated_nodes,
    const arma::vec& K_bin,
    const int& G_k,
    const arma::vec& B_bin,
    const int& G_b) {

  arma::cube FoF_cube = zeros(G_k, G_k, G_b);
  double x;
  double y;
  double z;
  arma::vec w = zeros(2);
  double x_new;
  double y_new;

  for(int i = 0; i < u_square.n_rows; i++){
    if(i+1 < u_square.n_rows){
      for(int j = i+1; j < u_square.n_rows; j++){
        x = K_bin(degree(i));
        y = K_bin(degree(j));
        z = B_bin(u_square(i, j));
        w(0) = x;
        w(1) = y;
        x_new = min(w);
        y_new = max(w);
        FoF_cube(x_new, y_new, z) += 1;
      }
    }
    FoF_cube(K_bin(0), K_bin(degree(i)), B_bin(0)) = FoF_cube(K_bin(0), K_bin(degree(i)), B_bin(0)) + sum(isolated_nodes);
  }
  if(0 < sum(isolated_nodes)){
    FoF_cube(K_bin(0), K_bin(0), B_bin(0)) = FoF_cube(K_bin(0), K_bin(0), B_bin(0)) + (sum(isolated_nodes))*(sum(isolated_nodes) - 1)/2;
  }

  return FoF_cube;
}

// not initial
// [[Rcpp::export("update_u_label")]]
arma::vec update_u_label(
    const arma::vec& u_label_old,
    const arma::vec& is_appear_old,
    const arma::vec& is_appear_new,
    const arma::vec& isolated_nodes_old,
    const arma::vec& isolated_nodes_new) {
  arma::vec u_label_new;
  arma::vec collaborated_nodes_old = is_appear_old - isolated_nodes_old;
  arma::vec collaborated_nodes_new = is_appear_new - isolated_nodes_new;
  arma::vec collaborated_nodes_add = collaborated_nodes_new - collaborated_nodes_old;

  int num = sum(collaborated_nodes_add);

  if(num == 0){
    u_label_new = u_label_old;

  }else if( num > 0 ) {
    u_label_new = u_label_old;
    int prev_size = u_label_new.n_elem;
    u_label_new.resize(sum(collaborated_nodes_new));
    int count = 0;

    for(int i = 0; i < collaborated_nodes_add.n_elem; i++){
      if(collaborated_nodes_add(i) == 1){
        u_label_new(prev_size + count) = i + 1;
        count = count + 1;
      }

      if(count == num){
        break;                                             // need to be checked
      }
    }
  }
  return u_label_new;
}




// not initial
// [[Rcpp::export("update_u")]]
arma::sp_mat update_u(const arma::sp_mat& u_old,
                const arma::vec& u_label_new,
                const arma::mat& new_edges) {
  arma::sp_mat u_new = u_old;

  // update u
  u_new.resize(u_label_new.n_elem, u_label_new.n_elem);
  if(new_edges(0,0) == 0){
    u_new = u_old;
  }else if(new_edges(1,0) == 0){
    int node1 = new_edges(0,0);
    int node2 = new_edges(0,1);
    int a;
    int b;
    for(int i = 0; i < u_label_new.n_elem; i++){
      if(node1 == u_label_new(i) ){
        a = i;
        for(int j = 0; j < u_label_new.n_elem; j++){
          if(node2 == u_label_new(j) ){
            b = j;
            u_new(a,b) = u_new(a,b) + 1;
            u_new(b,a) = u_new(b,a) + 1;
            break;
          }
        }
        break;
      }
    }
  }else{
    for(int i = 0; i < new_edges.n_rows; i++){
      int node1 = new_edges(i,0);
      int node2 = new_edges(i,1);
      int a;
      int b;
      for(int j = 0; j < u_label_new.n_elem; j++){
        if(u_label_new(j) == node1){
          a = j;
          for(int k = 0; k < u_label_new.n_elem; k++){
            if(u_label_new(k) == node2){
              b = k;
              u_new(a,b) = u_new(a,b) + 1;
              u_new(b,a) = u_new(b,a) + 1;
              break;
            }
          }
          break;
        }
      }
    }
  }
  return u_new;
}


// not initial
// doesn't count egdes which contains new node
// [[Rcpp::export("create_Z")]]
arma::cube create_Z(const arma::sp_mat& u_square_old,
              const arma::vec& u_label_old,
              const arma::vec& degree_old,
              const arma::mat& new_edges,
              const arma::vec& is_appear_old,
              const arma::vec& K_bin,
              const int& G_k,
              const arma::vec& B_bin,
              const int& G_b) {
  arma::cube Z = zeros(G_k,G_k,G_b);
  for(int i = 0; i < new_edges.n_rows; i++){
    int node1 = new_edges(i,0);
    int node2 = new_edges(i,1);
    if(node1 == 0){
      break;
    }
    if(is_appear_old(node1 - 1) == 1 && is_appear_old(node2 - 1) == 1){
      int a;
      int b;
      int flag = 0;
      for(int j = 0; j < u_label_old.n_elem; j++){
        if(u_label_old(j) == node1){
          a = j;
          for(int k = 0; k < u_label_old.n_elem; k++){
            if(u_label_old(k) == node2){                  // node1 and node2 are not-isolated
              b = k;
              double a_deg = degree_old(a);
              double b_deg = degree_old(b);
              double a_deg_bin = K_bin(a_deg);
              double b_deg_bin = K_bin(b_deg);
              arma::vec w = zeros(2);
              w(0) = a_deg_bin;
              w(1) = b_deg_bin;
              double x = min(w);
              double y = max(w);
              double z = B_bin(u_square_old(a,b));
              Z(x,y,z) = Z(x,y,z) + 1;
              flag = 1;
              break;
            }
          }
          if(flag == 0){                                      // node1 is not-isolated, but node2 is isolated
            double a_deg = degree_old(a);
            double a_deg_bin = K_bin(a_deg);
            Z(K_bin(0),a_deg_bin,B_bin(0)) = Z(K_bin(0),a_deg_bin,B_bin(0)) + 1;
            flag = 1;
          }
          break;
        }
      }
      if(flag == 0){                                      // node1 is isolated
        for(int k = 0; k < u_label_old.n_elem; k++){
          if(u_label_old(k) == node2){                        // node1 is isolated, but node2 is not-isolated
            b = k;
            double b_deg = degree_old(b);
            double b_deg_bin = K_bin(b_deg);
            Z(K_bin(0),b_deg_bin,B_bin(0)) = Z(K_bin(0),b_deg_bin,B_bin(0)) + 1;
            flag = 1;
            break;
          }
        }
      }
      if(flag == 0){                                     // node1 and node2 are isolated
        Z(K_bin(0),K_bin(0),B_bin(0)) = Z(K_bin(0),K_bin(0),B_bin(0)) + 1;
      }
    }
  }
  return Z;

}

///////////////////////////////////////////////////////////////////////////

//////////////////////////////Rename graph/////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// replace unique node ID with new node ID (rename graph)
// [[Rcpp::export("rename_graph")]]
arma::mat rename_graph(const arma::mat& graph_unique,    // graph with unique node ID
                 const arma::vec& unique_nodes,          // unique node ID
                 const arma::vec& nodes) {               // new node ID
  arma::mat graph;
  graph= graph_unique;
  int i;
  int j;
  int node1;
  int node2;
  int a;
  int b;
//#pragma omp parallel
{
//#pragma omp for private(j, node1, node2, a, b)
  for(i = 0; i < graph_unique.n_rows; i++){
    node1 = graph_unique(i,0);
    node2 = graph_unique(i,1);
    a = 0;
    b = 0;
    for(j = 0; j < unique_nodes.n_elem; j++){
      if(node1 == unique_nodes(j)){
        a = j;
        break;
      }
    }
    if(node2 != -1){
      for(j = 0; j < unique_nodes.n_elem; j++){
        if(node2 == unique_nodes(j)){
          b = j;
          break;
        }
      }
    }
    graph(i,0) = nodes(a);
    if(node2 != -1){
      graph(i,1) = nodes(b);
    }
  }
}
return graph;
}


///////////////////////////////////////////////////////////////////////////

/////////////////////////////u_from_graph_cpp//////////////////////////////

///////////////////////////////////////////////////////////////////////////
// [[Rcpp::export("u_from_graph_cpp")]]
arma::sp_mat u_from_graph_cpp(const arma::mat& graph, const int& n_node) {
  arma::vec indices(graph.n_rows);
  int s;
  int i;
  for(i = 0; i < graph.n_rows; i++){
    graph(i, 1) == -1 ? s = 0 : s = 1;
    indices(i) = s;
  }
  arma::mat graph_edge(sum(indices), 2);
  arma::vec graph_iso(graph.n_rows - sum(indices));
  int a = 0;
  int b = 0;
  for(i = 0; i < graph.n_rows; i++){
    if(indices(i) == 1){
      graph_edge(a,0) = graph(i, 0);
      graph_edge(a,1) = graph(i, 1);
      a += 1;
    } else {
      graph_iso(b) = graph(i,0);
      b += 1;
    }
  }

  arma::sp_mat u(n_node, n_node);
  int x;
  int y;
  for(i = 0; i < graph_edge.n_rows; i++){
    x = graph_edge(i,0) - 1;
    y = graph_edge(i,1) - 1;
    u(x, y)  = u(x, y) + 1;
    u(y, x)  = u(y, x) + 1;
  }
  return u;
}


///////////////////////////////////////////////////////////////////////////

//////////////////////////////log-likelihood///////////////////////////////

///////////////////////////////////////////////////////////////////////////



// [[Rcpp::export("create_log_likelihood")]]
double create_log_likelihood(SEXP my_field,
                             const int&   T_unique,
                             const arma::vec&   K_center_bin,
                             const arma::vec&   B_center_bin,
                             const int&   G_k,
                             const int&   G_b,
                             const int&  n_thread,
                             const double&  alpha,
                             const double&  beta) {
  double l = 0;// log-likelihood
  int t, i, j, k, x, y, z;
  double denominator;
  //arma::field<arma::cube> field(as<arma::field<arma::cube>>(my_field));
   List field(my_field);
  //std::cout << "Size of field at start of log_likelihood: " << arma::size(field) << " " << "\n";
//#pragma omp parallel num_threads(n_thread)
{
//#pragma omp parallel for reduction(+:l) private(i, j, k, x, y, z, denominator)
  for(t = 1; t < T_unique; t++){
    arma::cube * FoF_cube = new arma::cube(as<arma::cube>(field(t - 1, 0)));
    arma::cube * Z        = new arma::cube(as<arma::cube>(field(t, 1)));
    //std::cout << "Log likelihood" << "\n";
    //std::cout << "Size of FOF_cube: " << FoF_cube->n_rows << " " << FoF_cube->n_cols << " " << FoF_cube->n_slices << "\n";
    //std::cout << "Size of Z: " << Z->n_rows << " " << Z->n_cols << " " << Z->n_slices << "\n";
    //std::cout << "G_k and G_b: " << G_k << " " << G_b << "\n";

    //FoF_cube = return_cube(field, t - 1, 0);
    //Z        = return_cube(field, t, 1);

    // calculate denominator of P_{ij}(t) from FoF_cube
    denominator = 0;
    for(i = 0; i < G_k; i++){
      for(j = 0; j < G_k; j++){
        for(k = 0; k < G_b; k++){
          denominator += pow((K_center_bin(i) + 1), alpha)*pow((K_center_bin(j) + 1), alpha)*pow((B_center_bin(k) + 1), beta)*FoF_cube->at(i,j,k);
        }
      }
    }

    // calculate log-likelihood at t
    for(x = 0; x < G_k; x++){
      for(y = 0; y < G_k; y++){
        for(z = 0; z < G_b; z++){
          l += Z->at(x,y,z)*( log(pow((K_center_bin(x) + 1), alpha)*pow((K_center_bin(y) + 1), alpha)*pow((B_center_bin(z) + 1), beta)) - log(denominator) );
        }
      }
    }
    delete Z;
    delete FoF_cube;
  }
}
return l;
}

// [[Rcpp::export("return_cube")]]
arma::cube return_cube(const List& field, const int& x, const int& y){
  arma::cube C = field(x, y);
  return C;
}


// gradient of log-likelihood
// [[Rcpp::export("create_gradient")]]
arma::vec create_gradient(const List&     field,
                    const int&      T_unique,
                    const arma::vec&      K_center_bin,
                    const arma::vec&      B_center_bin,
                    const int&      G_k,
                    const int&      G_b,
                    const int&  n_thread,
                    const double&   alpha,
                    const double&   beta) {
  arma::vec gradient(2);
  double grad_alpha = 0;
  double grad_beta  = 0;
  int t, i, j, k, x, y, z;
  double denominator;
  double numerator_1;
  double numerator_2;
//#pragma omp parallel num_threads(n_thread)
{
//#pragma omp parallel for reduction(+ : grad_alpha, grad_beta) private(i, j, k, x, y, z, denominator, numerator_1, numerator_2)
  for(t = 1; t < T_unique; t++){
    arma::cube * FoF_cube = new arma::cube(as<arma::cube>(field(t - 1, 0)));
    arma::cube * Z        = new arma::cube(as<arma::cube>(field(t, 1)));
   //std::cout << "Gradient" << "\n";
   //std::cout << "Size of FOF_cube: " << FoF_cube->n_rows << " " << FoF_cube->n_cols << " " << FoF_cube->n_slices << "\n";
   //std::cout << "Size of Z: " << Z->n_rows << " " << Z->n_cols << " " << Z->n_slices << "\n";
   //std::cout << "G_k and G_b: " << G_k << " " << G_b << "\n";
    //FoF_cube = return_cube(field, t - 1, 0);
    //Z        = return_cube(field, t, 1);
   // Rcpp::as< arma::Cube<T> > - GitHub
    // culculate denominator and numerator
    denominator = 0;
    numerator_1 = 0;
    numerator_2 = 0;
    for(i = 0; i < G_k; i++){
      for(j = 0; j < G_k; j++){
        for(k = 0; k < G_b; k++){
          denominator += pow((K_center_bin(i) + 1), alpha)*pow((K_center_bin(j) + 1), alpha)*pow((B_center_bin(k) + 1), beta)*FoF_cube->at(i,j,k);
          numerator_1 += pow((K_center_bin(i) + 1), alpha)*pow((K_center_bin(j) + 1), alpha)*pow((B_center_bin(k) + 1), beta)*log((K_center_bin(i) + 1)*(K_center_bin(j) + 1))*FoF_cube->at(i,j,k);
          numerator_2 += pow((K_center_bin(i) + 1), alpha)*pow((K_center_bin(j) + 1), alpha)*pow((B_center_bin(k) + 1), beta)*log(B_center_bin(k) + 1)*FoF_cube->at(i,j,k);
        }
      }
    }

    // culculate gradient at t
    for(x = 0; x < G_k; x++){
      for(y = 0; y < G_k; y++){
        for(z = 0; z < G_b; z++){
          grad_alpha += Z->at(x,y,z)*( log((K_center_bin(x) + 1)*(K_center_bin(y) + 1)) - numerator_1/denominator );
          grad_beta  += Z->at(x,y,z)*( log(B_center_bin(z) + 1) - numerator_2/denominator );
        }
      }
    }
    delete Z;
    delete FoF_cube;
  }
}
gradient(0) = grad_alpha;
gradient(1) = grad_beta;
return gradient;
}

