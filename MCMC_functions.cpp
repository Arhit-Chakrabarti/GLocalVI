#include <RcppDist.h> 
#include <RcppArmadilloExtensions/sample.h> 
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace RcppArmadillo;

// Function to calculate 
// for(l in 1:L){
//  m[l]= sum(x == l)
//  }
// [[Rcpp::export]]
Rcpp::IntegerVector count_my(const arma::vec &x, const int &L){
  int n = x.size(); // Define the size of vector
  
  Rcpp::IntegerVector z(n); // To calculate 0 or 1 for the given vector of labels
  Rcpp::IntegerVector m(L); // To calculate 0 or 1 for the given vector of labels
  for(int l = 0; l < L; l++){
    for(int i = 0; i < n; i++){
      if(x(i) == (l + 1)){  // If x_i = l, then return z_i = 1
        z(i) = 1;
      }else{          // If x_i != l, then return z_i = 0
        z(i) = 0;
      }
    }
    int count = sum(z);
    m(l) = count;
  }
  return m;         // Return the vector of labels
}

// Function to subset the data. Here X is an n x p matrix.   Rcpp equivalent of X[z == h, ].
// Supplied mat x should be n x p matrix. 
// [[Rcpp::export]]
arma::mat data_subset(const arma::mat &x,  const arma::vec &z, const int &h){
  int p = x.n_cols;   // Number of columns correspond to number of variables
  arma::uvec index = find(z == h); // Find which elements of z == h
  
  int n = index.size(); // Calculate the size of the index vector, which stores how many z == h
  arma::mat sub_mat(n, p, arma::fill::zeros); // Define a sub-matrix to store the subsetted matrix data[which(z == h), ]
  
  sub_mat = x.rows(index);  // Calculates data[which(z == h), ]
  return sub_mat; // Return the sub matrix
}

// Rcpp equivalent of if k is a vector of indices and t is a vector of indices then to filter k[t].
// Example: k1.samples is a vector of size L with values 1:L. k1.samples is a vector of size n taking values 1:L.
// Then k_t(k1.samples, t1.samples) is an n component vector k1.samples[t1.samples]
// [[Rcpp::export]]
NumericVector k_t(const IntegerVector k, const IntegerVector t) {
  // NumericVector k_t(const arma::vec k, const IntegerVector t) { 
  // Length of the index vector
  int n_t = t.size();
  // Initialize output vector
  NumericVector out(n_t);
  
  // Loop through index vector and extract values of x at the given positions
  for (int i = 0; i < n_t; i++) {
    out[i] = k[t[i] - 1];
  }
  
  // Return output
  return out;
}
// This function is needed to convert Armadillo vector to Rcpp Numeric vector
// Maybe this function is not needed
NumericVector convertArmadilloToNumeric(arma::vec& x) {
  int n = x.n_elem;
  NumericVector result(n);
  
  // Copy elements from Armadillo vector to NumericVector
  for (int i = 0; i < n; i++) {
    result[i] = x(i);
  }
  
  return result;
}
// This function finds the index Rcpp equivalent of which(z == k)
// Maybe this function is not needed
NumericVector find_index(const arma::vec z, const int k){
  arma::vec index = arma::conv_to<arma::vec>::from(arma::find(z == k));
  NumericVector index_numeric = convertArmadilloToNumeric(index);
  return(index_numeric);
} 

// Rcpp function to calculate the sum of two matrices
arma::mat mat_add(const arma::mat A, const arma::mat B) {
  int m = A.n_rows;
  int n = A.n_cols;
  if (m != B.n_rows | n != B.n_cols) {
    stop("error: matrices dimensions do not match");
  }
  
  arma::mat R(m,n);
  int i = 0, j = 0;
  for (i = 0; i < m; i++) {     // loop over rows
    for (j = 0; j < n; j++) {   // loop over columns
      R(i,j) = A(i,j) + B(i,j); // elementwise addition 
    }
  } 
  return(R);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of 
// for(k in 1:L){
//   if(n.k[k] == 0){
//     x.k.bar[[k]] = c(0, 0)
//   }else{
//     x.k.bar[[k]] = (rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE]) + 
//       rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE]) + 
//       rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE]))
//   }
// }
// This function takes n_k: number of non-empty clusters i.e. n.k in the R program
// k1: k1.samples[[iter - 1]]; t1: t1.samples[[iter - 1]]
// k2: k2.samples[[iter - 1]]; t2: t2.samples[[iter - 1]]
// k3: k3.samples[[iter - 1]]; t3: t3.samples[[iter - 1]]
// x1: X.global[[1]]; x2: X.global[[2]]; x3: X.global[[3]]
// Unfortunately cannot take only the list X.global (MAYBE LOOK AT THIS LATER)
// [[Rcpp::export]]
arma::mat x_k_bar_c(const arma::vec n_k, const int p, const IntegerVector k1, const IntegerVector t1, const arma::mat& x1,
                    const IntegerVector k2, const IntegerVector t2, const arma::mat& x2,
                    const IntegerVector k3, const IntegerVector t3, const arma::mat& x3
) {
  int L = n_k.size(); // Find out the dimension L 
  arma::mat y1 = x1;  // Copy x1 into y1 
  arma::mat y2 = x2;  // Copy x2 into y2
  arma::mat y3 = x3;  // Copy x3 into y3
  // Initiate a matrix x1_k_sum which is pxL dimension to store
  // for(k in 1:L){ x2.k.sum[[k]] = rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE])}
  arma::mat x1_k_sum(p, L, arma::fill::zeros); 
  // Initiate a matrix x2_k_sum which is pxL dimension to store
  // for(k in 1:L){ x2.k.sum[[k]] = rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE])}
  arma::mat x2_k_sum(p, L, arma::fill::zeros);
  // Initiate a matrix x3_k_sum which is pxL dimension to store
  // for(k in 1:L){ x3.k.sum[[k]] = rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE])}
  arma::mat x3_k_sum(p, L, arma::fill::zeros);
  // Initiate a matrix x_k_bar which is pxL dimension to return the final means x.k.bar[[k]] in the form of k th column of x_k_bar
  arma::mat x_k_bar(p, L, arma::fill::zeros);
  for(int k = 1; k <= L; k ++){ //Here k runs from 1 to L
    if(n_k(k-1) > 1){
      arma::vec z1 = k_t(k1, t1); // find out k1.samples[[iter - 1]][t1.samples[[iter - 1]]]
      arma::vec index1 = arma::conv_to<arma::vec>::from(arma::find(z1 == k)); // which(k1.samples[[iter - 1]][t1.samples[[iter - 1]]]==k)
      int num1 = index1.size(); // find out the length(which(k1.samples[[iter - 1]][t1.samples[[iter - 1]]]==k))
      // This num1 is needed to define the dimension of out1 as below
      arma::mat out1(p, num1);
      // This out1 is needed to subset the data X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE]
      out1 = data_subset(y1.t(), z1, k).t();
      // This x1_sum is needed to calculate  rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE])
      arma::mat x1_sum = sum(out1, 1);
      // Store the  rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE]) into the kth column of x1_k_sum
      x1_k_sum.col(k-1) = x1_sum ;
      
      arma::vec z2 = k_t(k2, t2); // find out k2.samples[[iter - 1]][t2.samples[[iter - 1]]]
      arma::vec index2 = arma::conv_to<arma::vec>::from(arma::find(z2 == k)); // which(k2.samples[[iter - 1]][t2.samples[[iter - 1]]]==k)
      
      int num2 = index2.size();
      // This num2 is needed to define the dimension of out2 as below
      arma::mat out2(p, num2);
      // This out2 is needed to subset the data X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE]
      out2 = data_subset(y2.t(), z2, k).t();
      // This x2_sum is needed to calculate  rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE])
      arma::mat x2_sum = sum(out2, 1);
      // Store the  rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE]) into the kth column of x2_k_sum
      x2_k_sum.col(k-1) = x2_sum ;
      
      arma::vec z3 = k_t(k3, t3); // find out k3.samples[[iter - 1]][t3.samples[[iter - 1]]]
      arma::vec index3 = arma::conv_to<arma::vec>::from(arma::find(z3 == k)); // which(k3.samples[[iter - 1]][t3.samples[[iter - 1]]]==k)
      
      int num3 = index3.size();
      // This num3 is needed to define the dimension of out3 as below
      arma::mat out3(p, num3);
      // This out3 is needed to subset the data X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE]
      out3 = data_subset(y3.t(), z3, k).t();
      // This x3_sum is needed to calculate  rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE])
      arma::mat x3_sum = sum(out3, 1);
      // Store the  rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE]) into the kth column of x3_k_sum
      x3_k_sum.col(k-1) = x3_sum;
    }
    
  }
  // Once x1_k_sum, x2_k_sum, x3_k_sum is found out we add them i.e. we need this:
  // (rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE]) + 
  //  rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE]) + 
  //  rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE]))
  arma::mat x1_x2_sum = mat_add(x1_k_sum, x2_k_sum);
  arma::mat x_k_sum = mat_add(x1_x2_sum, x3_k_sum);
  
  // Once we have the sum to calulcate the cluster specific means we divide by the cluster specific sample sizes
  for(int k = 1; k <= L; k ++){
    if(n_k(k-1) > 1){
      x_k_bar.col(k-1) = x_k_sum.col(k-1)/n_k(k-1);
    }
  }
  
  return(x_k_bar);
}

// [[Rcpp::export]]
arma::mat x_k_bar_2_c(const arma::vec n_k, const int p, const IntegerVector k1, const IntegerVector t1, const arma::mat& x1,
                    const IntegerVector k2, const IntegerVector t2, const arma::mat& x2,
                    const IntegerVector k3, const IntegerVector t3, const arma::mat& x3,
                    const IntegerVector k4, const IntegerVector t4, const arma::mat& x4
) {
  int L = n_k.size(); // Find out the dimension L 
  arma::mat y1 = x1;  // Copy x1 into y1 
  arma::mat y2 = x2;  // Copy x2 into y2
  arma::mat y3 = x3;  // Copy x3 into y3
  arma::mat y4 = x4;  // Copy x3 into y3
  // Initiate a matrix x1_k_sum which is pxL dimension to store
  // for(k in 1:L){ x2.k.sum[[k]] = rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE])}
  arma::mat x1_k_sum(p, L, arma::fill::zeros); 
  // Initiate a matrix x2_k_sum which is pxL dimension to store
  // for(k in 1:L){ x2.k.sum[[k]] = rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE])}
  arma::mat x2_k_sum(p, L, arma::fill::zeros);
  // Initiate a matrix x3_k_sum which is pxL dimension to store
  // for(k in 1:L){ x3.k.sum[[k]] = rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE])}
  arma::mat x3_k_sum(p, L, arma::fill::zeros);
  // Initiate a matrix x4_k_sum which is pxL dimension to store
  // for(k in 1:L){ x4.k.sum[[k]] = rowSums(X.global[[4]][, k4.samples[[iter - 1]][t4.samples[[iter - 1]]] == k, drop = FALSE])}
  arma::mat x4_k_sum(p, L, arma::fill::zeros);
  // Initiate a matrix x_k_bar which is pxL dimension to return the final means x.k.bar[[k]] in the form of k th column of x_k_bar
  arma::mat x_k_bar(p, L, arma::fill::zeros);
  for(int k = 1; k <= L; k ++){ //Here k runs from 1 to L
    if(n_k(k-1) > 1){
      arma::vec z1 = k_t(k1, t1); // find out k1.samples[[iter - 1]][t1.samples[[iter - 1]]]
      arma::vec index1 = arma::conv_to<arma::vec>::from(arma::find(z1 == k)); // which(k1.samples[[iter - 1]][t1.samples[[iter - 1]]]==k)
      int num1 = index1.size(); // find out the length(which(k1.samples[[iter - 1]][t1.samples[[iter - 1]]]==k))
      // This num1 is needed to define the dimension of out1 as below
      arma::mat out1(p, num1);
      // This out1 is needed to subset the data X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE]
      out1 = data_subset(y1.t(), z1, k).t();
      // This x1_sum is needed to calculate  rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE])
      arma::mat x1_sum = sum(out1, 1);
      // Store the  rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE]) into the kth column of x1_k_sum
      x1_k_sum.col(k-1) = x1_sum ;
      
      arma::vec z2 = k_t(k2, t2); // find out k2.samples[[iter - 1]][t2.samples[[iter - 1]]]
      arma::vec index2 = arma::conv_to<arma::vec>::from(arma::find(z2 == k)); // which(k2.samples[[iter - 1]][t2.samples[[iter - 1]]]==k)
      
      int num2 = index2.size();
      // This num2 is needed to define the dimension of out2 as below
      arma::mat out2(p, num2);
      // This out2 is needed to subset the data X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE]
      out2 = data_subset(y2.t(), z2, k).t();
      // This x2_sum is needed to calculate  rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE])
      arma::mat x2_sum = sum(out2, 1);
      // Store the  rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE]) into the kth column of x2_k_sum
      x2_k_sum.col(k-1) = x2_sum ;
      
      arma::vec z3 = k_t(k3, t3); // find out k3.samples[[iter - 1]][t3.samples[[iter - 1]]]
      arma::vec index3 = arma::conv_to<arma::vec>::from(arma::find(z3 == k)); // which(k3.samples[[iter - 1]][t3.samples[[iter - 1]]]==k)
      
      int num3 = index3.size();
      // This num3 is needed to define the dimension of out3 as below
      arma::mat out3(p, num3);
      // This out3 is needed to subset the data X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE]
      out3 = data_subset(y3.t(), z3, k).t();
      // This x3_sum is needed to calculate  rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE])
      arma::mat x3_sum = sum(out3, 1);
      // Store the  rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE]) into the kth column of x3_k_sum
      x3_k_sum.col(k-1) = x3_sum;
      
      arma::vec z4 = k_t(k4, t4); // find out k4.samples[[iter - 1]][t4.samples[[iter - 1]]]
      arma::vec index4 = arma::conv_to<arma::vec>::from(arma::find(z4 == k)); // which(k4.samples[[iter - 1]][t4.samples[[iter - 1]]]==k)
      
      int num4 = index4.size();
      // This num3 is needed to define the dimension of out3 as below
      arma::mat out4(p, num4);
      // This out4 is needed to subset the data X.global[[4]][, k4.samples[[iter - 1]][t4.samples[[iter - 1]]] == k, drop = FALSE]
      out4 = data_subset(y4.t(), z4, k).t();
      // This x4_sum is needed to calculate  rowSums(X.global[[4]][, k4.samples[[iter - 1]][t4.samples[[iter - 1]]] == k, drop = FALSE])
      arma::mat x4_sum = sum(out4, 1);
      // Store the  rowSums(X.global[[4]][, k4.samples[[iter - 1]][t4.samples[[iter - 1]]] == k, drop = FALSE]) into the kth column of x4_k_sum
      x4_k_sum.col(k-1) = x4_sum;
    }
    
  }
  // Once x1_k_sum, x2_k_sum, x3_k_sum is found out we add them i.e. we need this:
  // (rowSums(X.global[[1]][, k1.samples[[iter - 1]][t1.samples[[iter - 1]]] == k, drop = FALSE]) + 
  //  rowSums(X.global[[2]][, k2.samples[[iter - 1]][t2.samples[[iter - 1]]] == k, drop = FALSE]) + 
  //  rowSums(X.global[[3]][, k3.samples[[iter - 1]][t3.samples[[iter - 1]]] == k, drop = FALSE]))
  arma::mat x1_x2_sum = mat_add(x1_k_sum, x2_k_sum);
  arma::mat x1_x2_x3_sum = mat_add(x1_x2_sum, x3_k_sum);
  arma::mat x_k_sum = mat_add(x1_x2_x3_sum, x3_k_sum);
  
  // Once we have the sum to calulcate the cluster specific means we divide by the cluster specific sample sizes
  for(int k = 1; k <= L; k ++){
    if(n_k(k-1) > 1){
      x_k_bar.col(k-1) = x_k_sum.col(k-1)/n_k(k-1);
    }
  }
  
  return(x_k_bar);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of 
// for(k in 1:L){
//   phi.samples[[k]] = mvrnorm(n = 1, mu = as.vector(solve(n.k[k] * SigmaG.inv + (lambda) * SigmaG0.inv) %*% (n.k[k] * SigmaG.inv %*% x.k.bar[[k]] + (lambda) * SigmaG0.inv %*% phi0)),
//                              Sigma = solve(n.k[k] * SigmaG.inv + (lambda) * SigmaG0.inv))
// }

// [[Rcpp::export]]
arma::mat draw_phi_samples(const arma::vec n_k, const arma::vec phi_0, const arma::mat x_k_bar, const arma::mat SigmaG_inv,  const arma::mat SigmaG0_inv, const double lambda){
  int L = n_k.size();
  int p = SigmaG_inv.n_cols;
  
  arma::mat phi_samples(p, L);      // Define a matrix to store draws from posterior MVN distribution
  arma::cube Normal_Sigma(p, p, L);
  arma::mat Normal_mean(p, L);
  for(int k = 0; k < L; k++){
    Normal_Sigma.slice(k) = arma::pinv((n_k(k) * SigmaG_inv + (lambda) * SigmaG0_inv), 0.00000000000000001);
    Normal_mean.col(k) = Normal_Sigma.slice(k) * (n_k(k) * SigmaG_inv * x_k_bar.col(k) + (lambda) * SigmaG0_inv * phi_0);
    phi_samples.col(k) = mvnrnd(Normal_mean.col(k), arma::symmatu(Normal_Sigma.slice(k)));
  }
  return(phi_samples);
}


////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of 
// for(t in 1:L){
//   if(n.j.t[[j]][t] == 0){
//     x.j.bar.local[[j]][[t]] = 0
//   }else{
//     x.j.bar.local[[j]][[t]] = sum(X.local[[j]][t.samples[[j]][[iter - 1]] == t])/n.j.t[[j]][t]
//   }
// }
// For j = 1, 2, 3 etc.
// [[Rcpp::export]]
arma::mat x_j_bar_local_c(const arma::vec n_j_t, const int p, const arma::vec t, const arma::mat& x) {
  int L = n_j_t.size(); // Find out the dimension L 
  arma::mat y = x;  // Copy x into y
  arma::vec z = t;
  // Initiate a matrix x_j_t_sum which is pxL dimension to store
  // for(t in 1:L){ x.j.bar.local[[t]] = sum(X.local[[j]][t.samples[[j]][[iter - 1]] == t])
  arma::mat x_j_t_sum(p, L, arma::fill::zeros); 
  
  // Initiate a matrix x_j_t_bar which is pxL dimension to return the final means x.j.bar.local[[j]]in the form of k th column of x_k_bar
  arma::mat x_j_t_bar(p, L, arma::fill::zeros);
  for(int k = 1; k <= L; k ++){ //Here k runs from 1 to L
    if(n_j_t(k-1) >= 1){
      // IntegerVector z1 = t; // copy t.samples[[j]][[iter - 1]] 
      // z = arma::conv_to<arma::vec>::from(Rcpp::z1);
      arma::vec index = arma::conv_to<arma::vec>::from(arma::find(z == k)); // which(t.samples[[j]][[iter - 1]] == k)
      int num = index.size(); // find out the length(which(which(t.samples[[j]][[iter - 1]] == k))
      // This num is needed to define the dimension of out1 as below
      arma::mat out(p, num);
      // This out is needed to subset the data X.local[[j]][, t.samples[[j]][[iter - 1]] == t, drop = FALSE]
      out = data_subset(y.t(), z, k).t();
      // This x_sum is needed to calculate  rowSums(X.local[[j]][, t.samples[[j]][[iter - 1]] == t, drop = FALSE])
      arma::mat x_sum = sum(out, 1);
      // Store the  rowSums(X.local[[j]][, t.samples[[j]][[iter - 1]] == t, drop = FALSE]) into the kth column of x_j_t_sum
      x_j_t_sum.col(k-1) = x_sum ;
    }
  }
  // Once we have the sum to calulcate the cluster specific means we divide by the cluster specific sample sizes
  for(int k = 1; k <= L; k ++){
    if(n_j_t(k-1) >= 1){
      x_j_t_bar.col(k-1) = x_j_t_sum.col(k-1)/n_j_t(k-1);
    }
  }
  
  return(x_j_t_bar);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of (When the data to be drawn in univariate)
// for(t in 1:L){
//   psi.samples[[1]][[t]] = rnorm(n = 1, 
//                                 mean = as.vector(solve(n.j.t[[1]][t] * Sigma1L.inv + (lambda) * Sigma1L0.inv) %*% (n.j.t[[1]][t] * Sigma1L.inv %*% x.j.bar.local[[1]][[t]] + (lambda) * Sigma1L0.inv %*% mu[[1]])),
//                                 sd = sqrt(as.numeric(solve(n.j.t[[1]][t] * Sigma1L.inv + (lambda) * Sigma1L0.inv))))
// }
// [[Rcpp::export]]
arma::mat draw_psi_samples_uni(const arma::vec n_j_t, const double mu, const arma::vec x_j_bar_local, const double SigmajL_inv,  const double SigmajL0_inv, const double lambda){
  int L = n_j_t.size();
  
  arma::mat psi_samples(1, L);      // Define a matrix to store draws from posterior univariate Normal distribution
  arma::vec Normal_Sigma(L);
  arma::vec Normal_mean(L);
  for(int k = 0; k < L; k++){
    Normal_Sigma(k) = sqrt(1/((n_j_t(k) * SigmajL_inv + (lambda) * SigmajL0_inv)));
    Normal_mean(k) = Normal_Sigma(k) * Normal_Sigma(k) * (n_j_t(k) * SigmajL_inv * x_j_bar_local(k) + (lambda) * SigmajL0_inv * mu);
    psi_samples.col(k) = R::rnorm(Normal_mean(k), Normal_Sigma(k));
  }
  return(psi_samples);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of (When the data to be drawn in multivariate)
// for(t in 1:L){
//   psi.samples[[j]][[t]] = mvrnorm(n = 1, 
//                                   mu = as.vector(solve(n.j.t[[j]][t] * Sigma2L.inv + (lambda) * Sigma2L0.inv) %*% (n.j.t[[j]][t] * Sigma2L.inv %*% x.j.bar.local[[j]][[t]] + (lambda) * Sigma2L0.inv %*% mu[[j]])),
//                                   Sigma = solve(n.j.t[[j]][t] * Sigma2L.inv + (lambda) * Sigma2L0.inv))
// }
// For j = 2, 3 etc.

// [[Rcpp::export]]
arma::mat draw_psi_samples_mv(const arma::vec n_j_t, const arma::vec mu, const arma::mat x_j_bar_local, const arma::mat SigmajL_inv,  const arma::mat SigmajL0_inv, const double lambda){
  int L = n_j_t.size();
  int p = SigmajL_inv.n_cols;
  
  arma::mat psi_samples(p, L);      // Define a matrix to store draws from posterior MVN distribution
  arma::cube Normal_Sigma(p, p, L);
  arma::mat Normal_mean(p, L);
  for(int k = 0; k < L; k++){
    Normal_Sigma.slice(k) = arma::pinv((n_j_t(k) * SigmajL_inv + (lambda) * SigmajL0_inv), 0.00000000000000001);
    Normal_mean.col(k) = Normal_Sigma.slice(k) * (n_j_t(k) * SigmajL_inv * x_j_bar_local.col(k) + (lambda) * SigmajL0_inv * mu);
    psi_samples.col(k) = mvnrnd(Normal_mean.col(k), arma::symmatu(Normal_Sigma.slice(k)));
  }
  return(psi_samples);
}
////////////////////////////////////////////////////////////////////////////////
// AUXILLARY FUNCTION NEEDED TO CALCLUATE DENSITY OF MULTIVARIATE NORMAL
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

////////////////////////////////////////////////////////////////////////////////
// Function to calculate the log density of a multivariate Normal
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec log_dmvnorm_my( arma::mat x,  arma::mat mean,  arma::mat sigma){ 
  
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double log2pi = 1.8378770664093454835606594728112352797227949472755668;
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  
  
  return(logretval);
  
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat prob_exponent_no_local(const arma::vec pi_samples, const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k, const arma::mat Sigma_G){
  int L = pi_samples.size();
  int n = X_global.n_rows;
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  arma::mat out3(n, L, arma::fill::zeros);
  arma::mat out(n, L, arma::fill::zeros);
  
  for(int t = 0; t < L; t++){
    int index = k(t) - 1;
    arma::rowvec mean_normal = phi_samples.col(index).t();
    out2.col(t) = arma::as_scalar(log(pi_samples(t))) + out1.col(t);
    out3.col(t) = log_dmvnorm_my(X_global, mean_normal, Sigma_G);
  }
  
  out = mat_add(out2, out3);
  return(out);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat prob_exponent_no_local2(const arma::vec pi_samples, const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k, const arma::cube Sigma_G){
  int L = pi_samples.size();
  int n = X_global.n_rows;
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  arma::mat out3(n, L, arma::fill::zeros);
  arma::mat out(n, L, arma::fill::zeros);
  
  for(int t = 0; t < L; t++){
    int index = k(t) - 1;
    arma::rowvec mean_normal = phi_samples.col(index).t();
    out2.col(t) = arma::as_scalar(log(pi_samples(t))) + out1.col(t);
    out3.col(t) = log_dmvnorm_my(X_global, mean_normal, Sigma_G.slice(index));
  }
  
  out = mat_add(out2, out3);
  return(out);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat prob_exponent_HDP_univariate(const arma::vec pi_samples, const arma::mat phi_samples, const arma::mat X, const IntegerVector k){
  int L = pi_samples.size();
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat Sigma (p, p, arma::fill::eye);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  arma::mat out3(n, L, arma::fill::zeros);
  arma::mat out(n, L, arma::fill::zeros);
  
  for(int t = 0; t < L; t++){
    int index = k(t) - 1;
    arma::rowvec mean_normal = phi_samples.col(index).t();
    out2.col(t) = arma::as_scalar(log(pi_samples(t))) + out1.col(t);
    out3.col(t) = log_dmvnorm_my(X, mean_normal, Sigma);
  }
  
  out = mat_add(out2, out3);
  return(out);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of 
//  for(t in 1:L){
//   # Population 1
//   exponent = log(pi1.samples[[iter]][t]) +  dnorm(as.numeric(X1[c(1), ]), mean = as.numeric(unlist(psi.samples[[1]][t])), sd = Sigma1L, log = TRUE) +
//.             mvtnorm::dmvnorm(t(X1[c(2, 3), ]), mean = phi.samples[[k1.samples[[iter - 1]][t]]], sigma = SigmaG, log = TRUE)
// } 
// THIS FUNCTION CONSIDERS THE FACT THAT THE LOCAL VARIABLE IS UNIVARIATE AND GLOBAL IS ALWAYS MULTIVARIATE
// IN OUR EXAMPLE THIS IS FOR POPULATION 1
// THIS FUNCTION IS USED TO CALCULATE THE PROBABILITIES BEFORE SAMPLING LOCAL LATENT VARIABLES i.e. the t1's
// [[Rcpp::export]]
arma::mat prob_exponent_univ_local(const arma::vec pi_samples, const arma::vec X_local, const arma::vec psi_samples,  const arma::vec Sigma,
                                   const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k, const arma::cube Sigma_G){
  int L = pi_samples.size();
  int n = X_local.size();
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  
  for(int t = 0; t < L; t++){
    int index = k(t) - 1;
    arma::rowvec mean_normal = phi_samples.col(index).t();
    out1.col(t) = log(pi_samples(t)) + arma::log_normpdf(X_local, psi_samples(t), Sigma(t)) ;
    out2.col(t) = log_dmvnorm_my(X_global, mean_normal, Sigma_G.slice(index));
  }
  out = mat_add(out1, out2);
  return(out);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of 
//  for(t in 1:L){
// # Population 2
//   exponent = log(pi2.samples[[iter]][t]) +  mvtnorm::dmvnorm(t(X2[c(1,2), ]), mean = as.numeric(unlist(psi.samples[[2]][t])), sigma = Sigma2L, log = TRUE) + 
//              mvtnorm::dmvnorm(t(X2[c(3, 4), ]), mean = phi.samples[[k2.samples[[iter - 1]][t]]],  sigma = SigmaG, log = TRUE)
// }
// THIS FUNCTION CONSIDERS THE FACT THAT THE LOCAL VARIABLE IS MULTIVARIATE AND GLOBAL IS ALWAYS MULTIVARIATE
// IN OUR EXAMPLE THIS IS FOR POPULATION 2 AND 3
// THIS FUNCTION IS USED TO CALCULATE THE PROBABILITIES BEFORE SAMPLING LOCAL LATENT VARIABLES i.e. the tj's, j = 2, 3
// [[Rcpp::export]]
arma::mat prob_exponent_mv_local(const arma::vec pi_samples, const arma::mat X_local, const arma::mat psi_samples,  const arma::cube Sigma_L,
                                 const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k, const arma::cube Sigma_G){
  int L = pi_samples.size();
  int n = X_local.n_rows;
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  
  for(int t = 0; t < L; t++){
    int index = k(t) - 1;
    arma::rowvec mean_normal_global = phi_samples.col(index).t();
    arma::rowvec mean_normal_local = psi_samples.col(t).t();
    out1.col(t) = log(pi_samples(t)) + log_dmvnorm_my(X_local, mean_normal_local, Sigma_L.slice(t)) ;
    out2.col(t) = log_dmvnorm_my(X_global, mean_normal_global, Sigma_G.slice(index));
  }
  out = mat_add(out1, out2);
  return(out);
}
////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// This function first takes as input the matrix obtained from either the function
// prob_exponent_univ_local or prob_exponent_mv_local then applies the Log-Sum-Exp trick//
// and then normalizes the matrix such the rows add to 1 to return a probability matrix
// Rcpp equivalent of 
// t.prob = t(apply(t.prob, MARGIN = 1, log_sum)), where t.prob is obtained from
// prob_exponent_univ_local or prob_exponent_mv_local and then
// t.prob = t.prob/rowSums(t.prob); 
// [[Rcpp::export]]
arma::mat calc_probability_log_sum_exp_normalized(const arma::mat x){
  int n = x.n_rows;
  int p = x.n_cols;
  arma::mat prob_mat(n, p);
  for(int i = 0; i < n; i++){
    prob_mat.row(i) = exp(x.row(i) - max(x.row(i))); //Log-Sum_Exp trick
    prob_mat.row(i) = prob_mat.row(i) + 0.0000000001; // Add a very small number to stabilize computation
    prob_mat.row(i) = prob_mat.row(i)/sum(prob_mat.row(i)); // Normalize to make it a probability
  }
  return(prob_mat);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to sample observations with a given probability matrix
// The probability matrix is nxL dimension where n is the number of samples to be drawn
// L is the number of categories 
// The R equivalent is
// for(i in 1:nrow(prob.matrix)){
// samples_to_be_drawn[i] = sample(1:ncol(prob.matrix), size = 1, prob = prob.matrix[i, ], replace = TRUE)
//  }
// [[Rcpp::export]]
Rcpp::IntegerVector sample_my(const arma::mat& prob) {
  int n = prob.n_rows; // Defining the number of observations
  Rcpp::IntegerVector z(n);  // Defining the vector z
  int L = prob.n_cols; // Define the number of columns which corresponds to the number of clusters
  arma::vec cluster = arma::regspace(1,  L); //Define a vector to store labels 1, 2, ... , L
  
  for(int i = 0; i < n; i++){
    z(i) = Rcpp::RcppArmadillo::sample(cluster, 1, true, prob.row(i).t())[0];
  }
  return(z);
}
// SOME AUXILLARY FUNCTION
// Rcpp Function to calculate sum(t.samples == t)
// [[Rcpp::export]]
double sum_which(const arma::vec &x, const int &h){
  int n = x.size(); // Define the size of vector
  
  Rcpp::IntegerVector z(n); // To return the vector of labels
  for(int i = 0; i < n; i++){
    if(x(i) == h){  // If x_i = h, then return z_i = 1
      z(i) = 1;
    }else{          // If x_i != h, then return z_i = 0
      z(i) = 0;
    }
  }
  return sum(z);    // Return the vector of labels
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of 
//    for(t in 1:L){
// for(k in 1:L){
// # For population j
//   if(sum(tj.samples[[iter]] == t) > 0){
//     prob[t, k, j] = log(beta.samples[[iter]][k]) + sum(mvtnorm::dmvnorm(t(X.global[[j]][ ,tj.samples[[iter]] == t]), mean = phi.samples[[k]], sigma = SigmaG, log = TRUE))
//   }else{
//     prob[t, k, j] = log(beta.samples[[iter]][k])
//   }
// }
// } For j = 1, 2, 3
// THIS FUNCTION CONSIDERS THE FACT THAT THE GLOBAL IS ALWAYS MULTIVARIATE
// THIS FUNCTION IS USED TO CALCULATE THE PROBABILITIES BEFORE SAMPLING GLOBAL LATENT VARIABLES i.e. the kj's
// [[Rcpp::export]]
arma::mat prob_exponent_mv_global(const arma::vec beta_samples, const arma::mat X_global,
                                  const arma::mat phi_samples, const arma::cube Sigma_G,
                                  const arma::vec t_samples){  
  
  int L = beta_samples.size();
  arma::mat out(L,L, arma::fill::zeros);
  for(int t = 0; t < L; t++){
    for(int k = 0; k < L; k++){
      double cond = sum_which(t_samples, (t + 1));
      if(cond > 0){
        arma::mat log_denisty_MVN = log_dmvnorm_my(data_subset(X_global, t_samples, (t + 1)), phi_samples.col(k).t(), Sigma_G.slice(k));
        double summed_value = log(beta_samples(k)) +  arma::as_scalar(sum(log_denisty_MVN, 0));
        out(t, k) = summed_value;
      }else{
        double summed_value = log(beta_samples(k));
        out(t, k) = summed_value;
      }
      
    }
  }
  return(out);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Function to calculate the Rcpp equivalent of 
//    for(t in 1:L){
// for(k in 1:L){
// # For population j
//   if(sum(tj.samples[[iter]] == t) > 0){
//     prob[t, k, j] = log(beta.samples[[iter]][k]) + sum(mvtnorm::dmvnorm(t(X.global[[j]][ ,tj.samples[[iter]] == t]), mean = phi.samples[[k]], sigma = SigmaG, log = TRUE))
//   }else{
//     prob[t, k, j] = log(beta.samples[[iter]][k])
//   }
// }
// } For j = 1, 2, 3
// THIS FUNCTION CONSIDERS THE FACT THAT THE GLOBAL IS ALWAYS MULTIVARIATE
// THIS FUNCTION IS USED TO CALCULATE THE PROBABILITIES BEFORE SAMPLING GLOBAL LATENT VARIABLES i.e. the kj's
// [[Rcpp::export]]
arma::mat prob_exponent_global_HDP_univariate(const arma::vec beta_samples, const arma::mat X,
                                  const arma::mat phi_samples,
                                  const arma::vec t_samples){  
  
  int L = beta_samples.size();
  arma::mat out(L,L, arma::fill::zeros);
  int p = X.n_cols;
  arma::mat Sigma (p, p, arma::fill::eye);
  
  for(int t = 0; t < L; t++){
    for(int k = 0; k < L; k++){
      double cond = sum_which(t_samples, (t + 1));
      if(cond > 0){
        arma::mat log_denisty_MVN = log_dmvnorm_my(data_subset(X, t_samples, (t + 1)), phi_samples.col(k).t(), Sigma);
        double summed_value = log(beta_samples(k)) +  arma::as_scalar(sum(log_denisty_MVN, 0));
        out(t, k) = summed_value;
      }else{
        double summed_value = log(beta_samples(k));
        out(t, k) = summed_value;
      }
      
    }
  }
  return(out);
}

////////////////////////////////////////////////////////////////////////////////
// AUXILLARY FUNCTION WHICH IS AN INDICATOR FUNCTION. THIS IS NEEDED FOR CALCULATING THE LOG-LIKELIHOOD
////////////////////////////////////////////////////////////////////////////////
double indicator(const int t_sample, const int t){
  double out; 
  if(t_sample == t){
    out = 1;
  }else{
    out = 0;
  }
  return(out);
}


// [[Rcpp::export]]
double logll_no_local(const arma::vec pi_samples,
                      const IntegerVector t_samples,
                      const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k_samples, const arma::mat Sigma_G){
  int L = pi_samples.size();
  int n = X_global.n_rows;
  arma::vec z = k_t(k_samples, t_samples); 
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  for(int i = 0; i < n; i ++){
    for(int t = 0; t < L; t++){
      
      out1(i, t) = arma::as_scalar(indicator(t_samples(i), (t + 1)));
      arma::rowvec mean_normal = phi_samples.col(t).t();
      double z_double = arma::as_scalar(z.row(i));
      int z_int = int(z_double);
      out2(i, t) = arma::as_scalar(indicator(z_int, (t + 1)) * log_dmvnorm_my(X_global.col(i).t(), mean_normal, Sigma_G)) ;
      
    }
  }
  
  out = mat_add(out1, out2);
  double result = arma::accu(out);
  return(result);
}

// [[Rcpp::export]]
double logll_HDP(const arma::vec pi_samples,
                      const IntegerVector t_samples,
                      const arma::mat phi_samples, const arma::mat X, const IntegerVector k_samples){
  int L = pi_samples.size();
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat Sigma (p, p, arma::fill::eye);
  
  arma::vec z = k_t(k_samples, t_samples); 
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  for(int i = 0; i < n; i ++){
    for(int t = 0; t < L; t++){
      
      out1(i, t) = arma::as_scalar(indicator(t_samples(i), (t + 1)));
      arma::rowvec mean_normal = phi_samples.col(t).t();
      double z_double = arma::as_scalar(z.row(i));
      int z_int = int(z_double);
      out2(i, t) = arma::as_scalar(indicator(z_int, (t + 1)) * log_dmvnorm_my(X.row(i).t(), mean_normal, Sigma)) ;
      
    }
  }
  
  out = mat_add(out1, out2);
  double result = arma::accu(out);
  return(result);
}

// [[Rcpp::export]]
double logll_no_local2(const arma::vec pi_samples,
                       const IntegerVector t_samples,
                       const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k_samples, const arma::cube Sigma_G){
  int L = pi_samples.size();
  int n = X_global.n_rows;
  arma::vec z = k_t(k_samples, t_samples); 
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  for(int i = 0; i < n; i ++){
    for(int t = 0; t < L; t++){
      
      out1(i, t) = arma::as_scalar(indicator(t_samples(i), (t + 1)));
      arma::rowvec mean_normal = phi_samples.col(t).t();
      double z_double = arma::as_scalar(z.row(i));
      int z_int = int(z_double);
      out2(i, t) = arma::as_scalar(indicator(z_int, (t + 1)) * log_dmvnorm_my(X_global.col(i).t(), mean_normal, Sigma_G.slice(t))) ;
      
    }
  }
  
  out = mat_add(out1, out2);
  double result = arma::accu(out);
  return(result);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Rcpp equivalent of 
//   for(i in 1:n1){
// for(t in 1:L){
//   logsum = logsum + indic(t1.samples[[iter]][i], t) * dnorm(as.vector(X1[c(1),i]), mean = unlist(Psi1.samples[[iter]][t]), sd = Sigma1L, log = TRUE) + 
//     indic(k1.samples[[iter]][t1.samples[[iter]]][i], t) * mvtnorm::dmvnorm(X1[c(2,3),i], mean = unlist(Phi.samples[[iter]][t]), sigma = SigmaG, log = TRUE) 
//   
// }
// }
// THIS FUNCTION EVALUATES THE SUM OF LOG-LIKELIHOOD WHEN BOTH LOCAL VARIABLE IS UNIVARIATE AND GLOBAL VARIABLE IS MULTIVARIATE
// THIS FUNCTION CONSIDERS THE FACT THAT THE LOCAL VARIABLE IS UNIVARIATE AND GLOBAL IS ALWAYS MULTIVARIATE
// IN OUR EXAMPLE THIS IS FOR POPULATION 1
// [[Rcpp::export]]
double logll_univ_local(const arma::vec pi_samples, const arma::vec X_local, const arma::vec psi_samples,  const double Sigma,
                        const IntegerVector t_samples,
                        const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k_samples, const arma::mat Sigma_G){
  int L = pi_samples.size();
  int n = X_local.size();
  arma::vec z = k_t(k_samples, t_samples); 
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  for(int i = 0; i < n; i ++){
    for(int t = 0; t < L; t++){
      
      out1(i, t) = arma::as_scalar(indicator(t_samples(i), (t + 1)) * arma::log_normpdf(X_local(i), psi_samples(t), Sigma));
      arma::rowvec mean_normal = phi_samples.col(t).t();
      double z_double = arma::as_scalar(z.row(i));
      int z_int = int(z_double);
      out2(i, t) = arma::as_scalar(indicator(z_int, (t + 1)) * log_dmvnorm_my(X_global.col(i).t(), mean_normal, Sigma_G)) ;
      
    }
  }
  
  out = mat_add(out1, out2);
  double result = arma::accu(out);
  return(result);
}

// [[Rcpp::export]]
double logll_univ_local2(const arma::vec pi_samples, const arma::vec X_local, const arma::vec psi_samples,  const arma::vec Sigma,
                        const IntegerVector t_samples,
                        const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k_samples, const arma::cube Sigma_G){
  int L = pi_samples.size();
  int n = X_local.size();
  arma::vec z = k_t(k_samples, t_samples); 
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  for(int i = 0; i < n; i ++){
    for(int t = 0; t < L; t++){
      
      out1(i, t) = arma::as_scalar(indicator(t_samples(i), (t + 1)) * arma::log_normpdf(X_local(i), psi_samples(t), Sigma(t)));
      arma::rowvec mean_normal = phi_samples.col(t).t();
      double z_double = arma::as_scalar(z.row(i));
      int z_int = int(z_double);
      out2(i, t) = arma::as_scalar(indicator(z_int, (t + 1)) * log_dmvnorm_my(X_global.col(i).t(), mean_normal, Sigma_G.slice(t))) ;
      
    }
  }
  
  out = mat_add(out1, out2);
  double result = arma::accu(out);
  return(result);
}

////////////////////////////////////////////////////////////////////////////////
// IMPORTANT FUNCTION
////////////////////////////////////////////////////////////////////////////////
// Rcpp equivalent of 
// for(i in 1:nj){
//   for(t in 1:L){
//     logsum = logsum + indic(tj.samples[[iter]][i], t) * mvtnorm::dmvnorm(X.local[[j]][,i], mean = unlist(Psij.samples[[iter]][t]), sigma = SigmajL, log = TRUE) + 
//       indic(kj.samples[[iter]][tj.samples[[iter]]][i], t) * mvtnorm::dmvnorm(X.glocal[[j]][,i], mean = unlist(Phi.samples[[iter]][t]), sigma = SigmaG, log = TRUE) 
//     
//   }
// } For j = 2, 3
// THIS FUNCTION EVALUATES THE SUM OF LOG-LIKELIHOOD WHEN BOTH LOCAL AND GLOBAL VARIABLE IS MULTIVARIATE
// THIS FUNCTION CONSIDERS THE FACT THAT THE LOCAL VARIABLE IS MULTIVARIATE AND GLOBAL IS ALWAYS MULTIVARIATE
// IN OUR EXAMPLE THIS IS FOR POPULATION 2 AND 3
// [[Rcpp::export]]
double logll_mv_local(const arma::vec pi_samples, const arma::mat X_local, const arma::mat psi_samples,  const arma::mat Sigma,
                      const IntegerVector t_samples,
                      const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k_samples, const arma::mat Sigma_G){
  int L = pi_samples.size();
  int n = X_local.n_cols;
  arma::vec z = k_t(k_samples, t_samples); 
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  for(int i = 0; i < n; i ++){
    for(int t = 0; t < L; t++){
      arma::rowvec mean_normal_local = psi_samples.col(t).t();
      out1(i, t) = arma::as_scalar(indicator(t_samples(i), (t + 1)) * log_dmvnorm_my(X_local.col(i).t(), mean_normal_local, Sigma));
      arma::rowvec mean_normal_global = phi_samples.col(t).t();
      double z_double = arma::as_scalar(z.row(i));
      int z_int = int(z_double);
      out2(i, t) = arma::as_scalar(indicator(z_int, (t + 1)) * log_dmvnorm_my(X_global.col(i).t(), mean_normal_global, Sigma_G)) ;
      
    }
  }
  
  out = mat_add(out1, out2);
  double result = arma::accu(out);
  return(result);
}

// [[Rcpp::export]]
double logll_mv_local2(const arma::vec pi_samples, const arma::mat X_local, const arma::mat psi_samples,  const arma::cube Sigma,
                      const IntegerVector t_samples,
                      const arma::mat phi_samples, const arma::mat X_global, const IntegerVector k_samples, const arma::cube Sigma_G){
  int L = pi_samples.size();
  int n = X_local.n_cols;
  arma::vec z = k_t(k_samples, t_samples); 
  arma::mat out(n, L, arma::fill::zeros);
  arma::mat out1(n, L, arma::fill::zeros);
  arma::mat out2(n, L, arma::fill::zeros);
  for(int i = 0; i < n; i ++){
    for(int t = 0; t < L; t++){
      arma::rowvec mean_normal_local = psi_samples.col(t).t();
      out1(i, t) = arma::as_scalar(indicator(t_samples(i), (t + 1)) * log_dmvnorm_my(X_local.col(i).t(), mean_normal_local, Sigma.slice(t)));
      arma::rowvec mean_normal_global = phi_samples.col(t).t();
      double z_double = arma::as_scalar(z.row(i));
      int z_int = int(z_double);
      out2(i, t) = arma::as_scalar(indicator(z_int, (t + 1)) * log_dmvnorm_my(X_global.col(i).t(), mean_normal_global, Sigma_G.slice(t))) ;
      
    }
  }
  
  out = mat_add(out1, out2);
  double result = arma::accu(out);
  return(result);
}