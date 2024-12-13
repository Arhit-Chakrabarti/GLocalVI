#include <RcppDist.h> 
#include <RcppArmadilloExtensions/sample.h> 
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace RcppArmadillo;
//------------------------------------------------------------------------------
// HELPER FUNCTIONS
//------------------------------------------------------------------------------
// [[Rcpp::export]]
double LogSumExp_cpp(arma::rowvec logX){
  double a = max(logX);
  return(  a + log(accu( exp( logX-a ) )));
}

//------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec reverse_cumsum_cpp(arma::colvec X){
  return( accu(X) - arma::cumsum(X));
}

//------------------------------------------------------------------------------
double log_Const_prod_gamma(int D, double nu){
  double Q = 0.0;
  for( int d=1; d<(D+1); d++){
    Q += lgamma( ( nu + 1.0 - (d) ) * .5);
  }
  return(Q);
}

//----------------------------------------------------------------------------
// [[Rcpp::export]]
double Const_sum_digamma(int D, double nu){
  double Q = 0.0;
  for( int d=1; d<(D+1); d++){
    Q += R::digamma( (nu + 1.0 - (d)) * .5);
  }
  return(Q);
}

//------------------------------------------------------------------------- B.79
double logB_wish(arma::mat W, double nu, double log_Const_prod_gamma){
  
  int D = W.n_cols;
  double p1 = - 0.5 * nu * log(arma::det(W));
  double p2 = (nu * D * .5) * log(2.0) + ( D * (D-1.0) / 4.0 ) * log(arma::datum::pi);
  return( p1 - p2 - log_Const_prod_gamma );
}

//------------------------------------------------------------------------- B.81
// [[Rcpp::export]]
double ElogDetLambda(arma::mat W, double Const_sum_digamma){
  
  int D = W.n_cols;
  return( Const_sum_digamma + log(arma::det(W)) + D*log(2.0) );
  
}

//------------------------------------------------------------------------- B.82
double  H_Lambda(arma::mat W, double nu, 
                 double Const_sum_digamma, 
                 double log_Const_prod_gamma){
  int D = W.n_cols;
  double ElDL = ElogDetLambda(W, Const_sum_digamma);
  double lnB = logB_wish(W, nu, log_Const_prod_gamma);
  return( - lnB - ( nu - D - 1.0) * 0.5 * ElDL + nu * D * 0.5);
  
} 
// [[Rcpp::export]]
arma::colvec E_log_beta(arma::colvec a,
                        arma::colvec b){
  
  int n = a.n_elem;
  arma::colvec res(n);
  for(int i=0; i<n; i++){
    res[i] = R::digamma(a[i]) - R::digamma(a[i] + b[i]);
  }
  return(res);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec E_log_p_Y_Mtheta_cpp_mvt(arma::mat Y,
                                      arma::colvec ml,
                                      double tl,
                                      double cl,
                                      arma::mat Dl){
  
  int p = Y.n_cols;
  int N = Y.n_rows;
  
  double CSDg = Const_sum_digamma(p, cl);
  double ell1 = ElogDetLambda(Dl, CSDg);                   // \ell1 in algorithm
  
  arma::colvec fq(N);
  arma::mat DIFF = Y - arma::repelem(ml.t(), N, 1); // Nj x D
  arma::mat DIFFt = DIFF.t();
  for(int i=0; i<N; i++){
    arma::vec temp = (DIFF.row(i) * Dl * (DIFFt.col(i)));
    fq(i) = temp[0];
  }
  
  arma::colvec ell2 = ( - p * 1.0/(tl) - cl * fq ); // Nj x 1
  
  return(.5 * (ell1 + ell2));
}
//------------------------------------------------------------------------------
// MAIN FUNCTIONS
//------------------------------------------------------------------------------
// --------------------------------------------------------------------------  1
// [[Rcpp::export]]
arma::mat Update_Vk_cpp(double b_bar,
                                arma::field<arma::mat> RHO_jtk,
                                const int J){
  //b_bar is r1/r2, which needs to be updated every time r1 and r2 are updated
  int T = RHO_jtk(0).n_rows;
  int L = RHO_jtk(0).n_cols;
  
  arma::colvec ml(L, arma::fill::zeros);
  
  for(int j = 0; j < J; j++){
    ml += (arma::sum(RHO_jtk(j),0)).t(); // This calculates \sum_{j=1}^{J}\sum_{t=1}^{T}RHO_jtk, k = 1,..,L. ml is Lx1 vector
  }
  arma::colvec a_tilde_vk      = ml        + 1.0; //This is \bar{a}_k = 1 + \sum_{j=1}^{J}\sum_{t=1}^{T}RHO_jtk
  arma::colvec rev_cs_ml = reverse_cumsum_cpp(ml); // This calculates \sum_{l=k+1}^{L}\sum_{j=1}^{J}\sum_{t=1}^{T}RHO_jtk,

  a_tilde_vk[L-1] = 1.0; //Since \bar{a}_k update goes from 1,...,L, We fix \bar{a}_L = 1
  arma::colvec b_tilde_vk      = rev_cs_ml + b_bar; //This is \bar{b}_k = r1/r2 + \sum_{l=k+1}^{L}\sum_{j=1}^{J}\sum_{t=1}^{T}RHO_jtk
  b_tilde_vk[L-1] = 1e-10; //Since \bar{b}_k update goes from 1,...,L, We fix \bar{b}_L to a ver small number
  
  arma::colvec E_ln_Vk    = E_log_beta(a_tilde_vk, b_tilde_vk); //This gives g(\bar{a}_k, \bar{b}_k) for k = 1,...,L as an Lx1 vector
  arma::colvec E_ln_1mVk  = E_log_beta(b_tilde_vk, a_tilde_vk); //This gives g(\bar{b}_k, \bar{a}_k) for k = 1,...,L as an Lx1 vector
  arma::colvec sE_ln_1mVk = shift(E_ln_1mVk, +1); // This shifts g(\bar{b}_k, \bar{a}_k) for k = 1,...,L by 1 cyclically, which yields the vector (g(\bar{b}_L, \bar{a}_L), g(\bar{b}_1, \bar{a}_1), ... , g(\bar{b}_L-1, \bar{a}_L-1))

  sE_ln_1mVk[0] = 0; // Set g(\bar{b}_L, \bar{a}_L) = 0 in the vector (g(\bar{b}_L, \bar{a}_L), g(\bar{b}_1, \bar{a}_1), ... , g(\bar{b}_L-1, \bar{a}_L-1))
  // This is because, we never really consider  g(\bar{b}_H, \bar{a}_H) when calculating \sum{r=1}^{k-1} g(\bar{b}_r, \bar{a}_r) for k = 1,...,L. when r = L, we consider g(\bar{b}_L-1, \bar{a}_L-1)
  arma::colvec CS_E_ln_1mVk = arma::cumsum(sE_ln_1mVk); // This gives\sum{r=1}^{k-1} g(\bar{b}_r, \bar{a}_r) for k = 1,...,L as an Lx1 vector
  arma::mat results(L,3);

  results.col(0) = a_tilde_vk;
  results.col(1) = b_tilde_vk;
  results.col(2) = E_ln_Vk + CS_E_ln_1mVk ; // This gives g(\bar{a}_k, \bar{b}_k) + \sum{r=1}^{k-1} g(\bar{b}_r, \bar{a}_r) for k = 1,...,L as an Lx1 vector, which is needed for ELBO calculation
  return(results);
}
// --------------------------------------------------------------------------- 2
// [[Rcpp::export]]
arma::cube Update_Ujt_cpp(arma::field<arma::mat> XI_jil,
                          double const b_bar,
                          int const T,
                          int const J){
  
  //b_bar is s1/s2, which needs to be updated every time s1 and s2 are updated
  arma::mat N_jt(J,T);
  
  for(int j=0; j<J; j++){
    N_jt.row(j) = arma::sum(XI_jil(j),0); //This calculates \sum_{i=1}^{n_j} XI_jil. So N_jt is a JxT matrix
  }
  
  arma::mat a_bar_Ujt(J,T);
  arma::mat b_bar_Ujt(J,T);
  arma::mat E_lnOmega_jt(J,T);
  
  for(int j = 0; j<J; j++){
  // int j = 2;
    arma::colvec G1 = 1.0 + N_jt.row(j).t();
    G1[T-1] = 1;
    a_bar_Ujt.row(j) = G1.t();
    arma::colvec G2 = reverse_cumsum_cpp(N_jt.row(j).t());
    G2[T-1] = 1e-10;
    b_bar_Ujt.row(j) = G2.t();
    
    arma::colvec E_ln_Ul_k    = E_log_beta(G1, G2); //This gives g(\bar{a}_lk, \bar{b}_lk) for l = 1,...,L as an Lx1 vector for each k = 1,...,H
    arma::colvec E_ln_1mUl_k  = E_log_beta(G2, G1); //This gives g(\bar{b}_lk, \bar{a}_lk) for l = 1,...,L as an Lx1 vector for each k = 1,...,H
    arma::colvec sE_ln_1mUl_k = shift(E_ln_1mUl_k, +1);  //This shifts the vector(g(\bar{b}_lk, \bar{a}_lk)), l=1,...,L by 1 position cyclically
    sE_ln_1mUl_k[0] = 0; //This sets (g(\bar{b}_Lk, \bar{a}_Lk)) = 0
    
    arma::colvec CS_E_ln_1mUL_k = arma::cumsum(sE_ln_1mUl_k); // This gives \sum{r=1}^{l-1} g(\bar{b}_rk, \bar{a}_rk)
    E_lnOmega_jt.row(j) = (CS_E_ln_1mUL_k + E_ln_Ul_k).t(); 
  }
  
  arma::cube results(J,T,3);
  results.slice(0) = a_bar_Ujt;
  results.slice(1) = b_bar_Ujt;
  results.slice(2) = E_lnOmega_jt;
  
  return(results);
}

// --------------------------------------------------------------------------- 3

// [[Rcpp::export]]
arma::colvec Update_gamma_concentration_par(arma::colvec a_tilde_Vk,
                                            arma::colvec b_tilde_Vk,
                                            arma::colvec conc_hyper){
  
  int  L = b_tilde_Vk.n_rows;
  arma::colvec upd_par(2);
  a_tilde_Vk.shed_row(L-1); // Removing the Hth row of a_tilde_Vk
  b_tilde_Vk.shed_row(L-1); // Removing the Hth row of b_tilde_Vk
  //This is needed because we need \sum_{k=1}^{L-1}g(\bar{b}_k, \bar(a)_k)
  
  upd_par[0] = conc_hyper[0] + L - 1.0 ; //This is r1 = a_{\gamma} + (L-1)
  upd_par[1] = conc_hyper[1] - arma::accu(E_log_beta(b_tilde_Vk,a_tilde_Vk)); //This is r2 = b_{\gamma} - \sum_{k=1}^{L-1}g(\bar{b}_k, \bar(a)_k)
  
  return(upd_par);
}

// --------------------------------------------------------------------------- 4
// [[Rcpp::export]]
arma::colvec Update_alpha_concentration_par(arma::mat a_bar_Ulk,
                                           arma::mat b_bar_Ulk,
                                           arma::colvec conc_hyper,
                                           int T,
                                           int J){
  
  arma::colvec upd_par(2);
  a_bar_Ulk.shed_col(T-1); // Removing the Lth row of a_bar_Ulk. a_bar_Ulk is an TxJ matrix originally
  b_bar_Ulk.shed_col(T-1); // Removing the Lth row of b_bar_Ulk. b_bar_Ulk is an TxJ matrix originally
  
  arma::colvec R(J);
  
  for(int j = 0; j < J; j ++){
    
    R[j] =  arma::accu( E_log_beta(b_bar_Ulk.row(j).t(),a_bar_Ulk.row(j).t()) ); // This calculates \sum_{j=1}^{J}g(\bar{b}_{jt}, \bar(a)_{jt})
    
  }
  
  
  upd_par[0] = conc_hyper[0] + J * (T - 1.0); //This is s1 = a_{\alpha} + J(T-1)
  upd_par[1] = conc_hyper[1] - arma::accu(R); //This is s2 = b_{\alpha} - \sum_{t=1}^{T-1}\sum_{j=1}^{J}g(\bar{b}_{jt}, \bar(a)_{jt})
  //This is why we removed the row T-1 a_bar_Ulk and b_bar_Ulk
  
  
  return(upd_par);
}
// Rcpp::List
// --------------------------------------------------------------------------- 5
// [[Rcpp::export]]
Rcpp::List  Update_Phi_k_cpp_mvt(arma::field<arma::mat> Y_grouped,
                                  arma::field<arma::mat> RHO_jtk,
                                  arma::field<arma::mat> XI_jil,
                                  arma::colvec m0,
                                  double tau0,
                                  double nu0,
                                  arma::mat W0,
                                  arma::mat iW0
                                ){
  
  int D = W0.n_rows;
  int J = XI_jil.n_elem;
  int T = XI_jil(0).n_cols;
  int L = RHO_jtk(0).n_cols;
  
  arma::mat SumX_k(D, L, arma::fill::zeros);
  arma::mat Xbar_k(D, L, arma::fill::zeros);
  arma::cube Sk(D,D, L, arma::fill::zeros);
  arma::cube Sk_over_delta_k(D,D, L, arma::fill::zeros);
  
  arma::cube Dk(D,D,L);
  arma::rowvec delta_k(L, arma::fill::zeros);
  
  arma::mat mk(D,L);
  arma::rowvec lambdak(L);
  arma::rowvec ck(L);

  for(int j = 0; j < J; j++){
  // int j = 0;
    arma::rowvec SUM_XI_jil = (arma::sum(XI_jil(j),0)); 
    arma::rowvec PROD_SUM_XI_jil_RHO_jtk = SUM_XI_jil * RHO_jtk(j);
    delta_k += PROD_SUM_XI_jil_RHO_jtk;

    arma::mat subY = Y_grouped[j];
    SumX_k += subY.t() * XI_jil(j) * RHO_jtk(j);
  }
  arma::mat tau0m0  = arma::repelem(tau0 * m0,1,L); // p x L. B0M0 is repeating \lambda_0\mu_0, which is a px1 vector L times along a column to get a pxL matrix. This is to facilitate matrix addition
  lambdak         = delta_k + tau0; // This is t_l = N_l + \lambda0
  ck             = delta_k + nu0 ;    // This is c_l = N_l + \nu_0
  arma::mat ilambdak   = arma::repelem(1.0 / lambdak, D, 1); // p x L. iBL is repeating 1/t_l, which is 1xL vector p times along the rows to get a pxL matrix
  mk              = ( tau0m0 + SumX_k ) % ( ilambdak ); // Element-wise multiplying ( B0M0 + SumY_l ) with iBL gives m_l

  
  // return(subY);
  for(int l=0; l<L; l++){
    if(delta_k(l)>0){
      Xbar_k.col(l) = SumX_k.col(l)/delta_k(l); // This calculates \bar{x}_k
    }
  }

  for(int l=0; l <L; l++){
    arma::mat S_k_temp(D, D);
    for(int j=0; j <J; j++){
      arma::mat subY = Y_grouped[j].t(); // D x Nj
      arma::mat S_j_k_temp(subY.n_rows,subY.n_rows, arma::fill::zeros);
      for(int ii=0; ii<subY.n_cols; ii++){
        double xi = 0.0;
        arma::mat temp = ( subY.col(ii) - Xbar_k.col(l)) *
          (subY.col(ii) - Xbar_k.col(l)).t(); //This calculating (x_ji - \bar{x}_k)(x_ji - \bar{x}_l)^T for a fixed j, i, and l
        for(int t=0; t<T; t++){
          xi += (XI_jil(j))(ii,t) * (RHO_jtk(j))(t,l); // This is XI_{jil} for a fixed j, i, and l
        }
        S_j_k_temp += xi * temp;
      }
      S_k_temp += S_j_k_temp;
    }
    Sk.slice(l) = S_k_temp;
    if(delta_k(l)>0){
        Sk_over_delta_k.slice(l) = Sk.slice(l)/delta_k(l);
        }
    // Dk.slice(l) = (
    //   (iW0)+
    //     arma::symmatu(Sk.slice(l)) +
    //     arma::symmatu(tau0*delta_k(l)/(tau0 + delta_k(l)) *
    //     (Xbar_k.col(l)-m0) * (Xbar_k.col(l)-m0).t())
    // ); // This gives D_k^{-1}
          Dk.slice(l) = arma::inv_sympd(
            (iW0)+
              arma::symmatu(Sk.slice(l)) +
              arma::symmatu(tau0*delta_k(l)/(tau0 + delta_k(l)) *
              (Xbar_k.col(l)-m0) * (Xbar_k.col(l)-m0).t())
          ); // This gives D_k^{-1}
  }

  Rcpp::List results = Rcpp::List::create(
    Rcpp::_["Sk"] = Sk,
    Rcpp::_["mk"]  = mk,
    Rcpp::_["SumX_k"] = SumX_k,
    Rcpp::_["lambdak"] = lambdak.t(),
    Rcpp::_["ck"] = ck.t(),
    Rcpp::_["Dk"] = Dk,
    Rcpp::_["Xbar_k"] = Xbar_k,
    Rcpp::_["deltak"] = delta_k,
    Rcpp::_["Sk_over_delta_k"] = Sk_over_delta_k

  );
  return(results);
  
}
// --------------------------------------------------------------------------- 6
// [[Rcpp::export]]
arma::field<arma::mat> Update_RHOjtk_cpp(arma::field<arma::mat> XI_jil,
                             arma::field<arma::mat> Y_grouped,
                            arma::mat mk,
                            arma::colvec lambdak,
                            arma::colvec ck,
                            arma::cube Dk,
                            arma::colvec ElnPI_k,
                            int const L,
                            int const J,
                            int const T){
  
  
  int D = Y_grouped(0).n_cols; //Dimension of global variables X^G

  arma::cube unn_log_RHO_jtk(T,L,J, arma::fill::zeros);
  arma::cube log_RHO_jtk(T,L,J);
  arma::cube RHO_jtk(T,L,J);

  // ElnPI_k is a Lx1 matrix and ElnPI_k(k) contains g(\bar{a}_k, \bar{b}_k) + \sum_{r=1}^{k-1}g(\bar{b}_r, \bar{a}_r)
  for(int j=0; j < J; j++){
    for(int k = 0; k < L; k++){
      unn_log_RHO_jtk.slice(j).col(k) =  unn_log_RHO_jtk.slice(j).col(k) + ElnPI_k(k);
    }// This gives the matrix for each j = 1,..,J
  } // g(\bar{a}_k, \bar{b}_k) + \sum_{r=1}^{k-1}g(\bar{b}_r, \bar{a}_r)
  
  arma::cube TT(T, L, J);
  for(int j=0; j<J; j++){
    arma::mat subY = Y_grouped[j].t(); // D x Nj
      for(int k=0; k<L; k++){
        arma::rowvec sum_XI_jil_E_log_p_X(T, arma::fill::zeros);
          for(int ii=0; ii<subY.n_cols; ii++){
            sum_XI_jil_E_log_p_X += XI_jil(j).row(ii) * arma::as_scalar(E_log_p_Y_Mtheta_cpp_mvt(subY.col(ii).t(),
              mk.col(k), lambdak(k), ck(k), Dk.slice(k)));
          }
          TT.slice(j).col(k) = sum_XI_jil_E_log_p_X.t();
        }
  }

  for(int j = 0; j < J; j++){
    // Adding the t th row of TT.slice(j) and t th row of unn_log_RHO_jk.slice(j), which gives a TxL matrix
    for(int t = 0; t<T; t++){
      arma::rowvec logunn(L, arma::fill::zeros);
      logunn = unn_log_RHO_jtk.slice(j).row(t) + TT.slice(j).row(t);
      log_RHO_jtk.slice(j).row(t) = logunn - LogSumExp_cpp(logunn); // Normalizing by LogSumExp trick
    }
  }

  arma::field<arma::mat> exp_log_RHO_jtk(J);
  for(int j = 0; j<J; j++){
    exp_log_RHO_jtk(j) = exp(log_RHO_jtk.slice(j));
  }
  // Returning the exp of normalized log_RHO_jtk to gives the probabilities RHO_jk
  return(exp_log_RHO_jtk);
}

// --------------------------------------------------------------------------- 7
// [[Rcpp::export]]
arma::field<arma::mat> Update_XIjil_cpp(arma::field<arma::mat> Y_grouped,
                                        arma::mat mk,
                                        arma::colvec lambdak,
                                        arma::colvec ck,
                                        arma::cube Dk,
                                        
                                        arma::field<arma::mat> Y_j_grouped,
                                        arma::field<arma::mat> m_jt,
                                        arma::field<arma::colvec> lambda_jt,
                                        arma::field<arma::colvec> c_jt,
                                        arma::field<arma::cube> D_jt,
                                        
                                        arma::field<arma::mat> RHO_jtk,
                                        arma::mat ElnOM_lk,
                                      
                                        int const L,
                                        int const J,
                                        int const T){


  arma::field<arma::mat> XI_jil(J); // J different nj x T matrices
  // ElnOM_lk contains the matrix [g(\bar{a}_jt, \bar{b}_jt) + \sum_{l=1}^{t-1}g(\bar{b}_jl, \bar{a}_jl)]
  // For t = 1,...,T; j = 1,...,J
  // So ElnOM_lk is an JxT matrix
  // RHO_jtk is a TxLxJ array

  for(int j = 0; j<J; j++){
    arma::mat subY_j = Y_j_grouped(j).t();
    arma::mat subY = Y_grouped(j).t();
    int N_j = Y_grouped[j].n_rows;
      arma::mat tempres(N_j,T, arma::fill::zeros);
      arma::mat tempres_final(N_j,T, arma::fill::zeros);
        for(int ii = 0; ii<N_j; ii++){
            arma::rowvec temp(T, arma::fill::zeros);
            for(int t = 0; t<T; t++){
                double ell1 = 0.0;
                 ell1 = ElnOM_lk(j,t) + arma::as_scalar(E_log_p_Y_Mtheta_cpp_mvt(subY_j.col(ii).t(),
                                            m_jt(j).col(t), lambda_jt(j)(t), c_jt(j)(t), D_jt(j).slice(t)));
                double sum_RHO_jtk_E_log_p_X = 0.0;
                  for(int k=0; k<L; k++){
                    sum_RHO_jtk_E_log_p_X += RHO_jtk(j)(t,k) * arma::as_scalar(E_log_p_Y_Mtheta_cpp_mvt(subY.col(ii).t(),
                                                    mk.col(k), lambdak(k), ck(k), Dk.slice(k)));
                    }
                  temp(t) = ell1 +  sum_RHO_jtk_E_log_p_X;
          }
            tempres.row(ii) =  temp - LogSumExp_cpp(temp);
        }
        // For a fixed j, returning the exp of normalized XI_jil to gives the probabilities XI_jil
        XI_jil(j) = exp(tempres);
  }
  
  return(XI_jil); // This is essentially a list of J probabilities with jth one being of dimension n_jxT
}

// --------------------------------------------------------------------------- 7
// [[Rcpp::export]]
arma::field<arma::mat> Jitter_XIjil_cpp(const arma::field<arma::mat> XI_jil){
  
  int J = XI_jil.n_elem;
  int T = XI_jil(0).n_cols;
  arma::field<arma::mat> XI_jil_Jitter(J); // J different nj x T matrices

  for(int j = 0; j<J; j++){
    arma::mat subXI = XI_jil(j);
    int N_j = subXI.n_rows;
    arma::mat tempres(N_j,T, arma::fill::zeros);
    for(int ii = 0; ii<N_j; ii++){
      tempres.row(ii) = subXI.row(ii) + 1e-10;
      tempres.row(ii) = tempres.row(ii)/arma::as_scalar(arma::accu(tempres.row(ii)));
    }
    XI_jil_Jitter(j) = tempres;
  }
  
  return(XI_jil_Jitter); // This is essentially a list of J probabilities with jth one being of dimension n_jxT
}
// --------------------------------------------------------------------------- 8
// [[Rcpp::export]]
Rcpp::List Update_Psi_jt_cpp_mvt(arma::field<arma::mat> Y_j_grouped,
                                arma::field<arma::mat> XI_jil,
                                arma::field<arma::colvec> m_j0,
                                arma::colvec tau_j0,
                                arma::colvec nu_j0,
                                arma::field<arma::mat> W_j0,
                                arma::field<arma::mat> iW_j0
){
  int J = XI_jil.n_elem;
  int T = XI_jil(0).n_cols;
  arma::mat delta_jt(J, T);
  arma::field<arma::colvec> lambda_jt(J); 
  arma::field<arma::colvec> c_jt(J);
  arma::field<arma::mat> m_jt(J);
  arma::field<arma::mat> SumX_jt(J);
  arma::field<arma::mat> Xbar_jt(J);
  arma::field<arma::cube> S_jt(J);
  arma::field<arma::cube> S_jt_over_delta_jt(J);
  arma::field<arma::cube> D_jt(J);
  
  for(int j = 0; j<J; j++){
    delta_jt.row(j) = (arma::sum(XI_jil(j),0)); 
    arma::rowvec lambda_jt_temp = delta_jt.row(j) + tau_j0(j);
    arma::rowvec c_jt_temp = delta_jt.row(j) + nu_j0(j);
    c_jt(j) = c_jt_temp.t();
    lambda_jt(j) = lambda_jt_temp.t();
    arma::mat subY_j = Y_j_grouped[j];
    SumX_jt(j) = subY_j.t() * XI_jil(j);
  }
  
  for(int j = 0; j<J; j++){
    int D_j = Y_j_grouped[j].n_cols;
      arma::mat temp_mean(D_j, T);
        for(int t=0; t<T; t++){
          if(delta_jt(j, t)>0){
            temp_mean.col(t) = SumX_jt(j).col(t)/delta_jt(j,t);
          }
      }
      Xbar_jt(j) = temp_mean;
  }
  
  for(int j = 0; j<J; j++){
    arma::mat subY_j = Y_j_grouped[j].t(); // D x Nj
    int D_j = subY_j.n_rows;
    arma::cube temp_St(D_j, D_j, T);
    arma::cube temp_St_over_delta_jt(D_j, D_j, T);
    for(int t = 0; t<T; t++){
      arma::mat temp_S(D_j, D_j);
      for(int ii = 0; ii<subY_j.n_cols; ii++){
        temp_S +=XI_jil(j)(ii, t) * ( subY_j.col(ii) - Xbar_jt(j).col(t)) *
                              (subY_j.col(ii) - Xbar_jt(j).col(t)).t();
      }
      temp_St.slice(t) = temp_S;
      if(delta_jt(j,t) >0){
        temp_St_over_delta_jt.slice(t) = temp_St.slice(t)/delta_jt(j,t);
      }
    }
    S_jt(j) = temp_St;
    S_jt_over_delta_jt(j) = temp_St_over_delta_jt;
  }
  
  for(int j = 0; j<J; j++){
    int D_j = Y_j_grouped[j].n_cols;
    arma::mat B0M0  = arma::repelem(tau_j0(j) * m_j0(j),1,T);
    arma::mat iBL   = arma::repelem(1.0 / lambda_jt(j).t(), D_j, 1); 
    arma::mat temp_mjt  = ( B0M0 + SumX_jt(j) ) % ( iBL );
    m_jt(j) = temp_mjt;
  }
  for(int j = 0; j<J; j++){
    int D_j = Y_j_grouped[j].n_cols;
    arma::cube Wj(D_j, D_j, T);
    for(int t = 0; t<T; t++){
      Wj.slice(t) = arma::inv_sympd(
        (iW_j0(j))+
          arma::symmatu(S_jt(j).slice(t)) +
          arma::symmatu(tau_j0(j)*delta_jt(j,t)/(tau_j0(j) + delta_jt(j,t)) *
          (Xbar_jt(j).col(t) - m_j0(j)) * (Xbar_jt(j).col(t) - m_j0(j)).t())
      ); 
    }
    D_jt(j) = Wj;
  }
  
  
  Rcpp::List results = Rcpp::List::create(
    Rcpp::_["lambda_jt"] = lambda_jt,
    Rcpp::_["c_jt"] = c_jt,
    Rcpp::_["m_jt"] = m_jt,
    Rcpp::_["D_jt"] = D_jt,
    Rcpp::_["SumX_jt"] = SumX_jt,
    Rcpp::_["Xbar_jt"] = Xbar_jt,
    Rcpp::_["delta_jt"] = delta_jt,
    Rcpp::_["S_jt"] = S_jt,
    Rcpp::_["S_jt_over_delta_jt"] = S_jt_over_delta_jt
  );
  return(results);
  
}

