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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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

double lbeta_normconst_cpp(double a, double b){
  double C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - (C_ab));
}

arma::colvec lbeta_normconst_vec_cpp(arma::colvec a, arma::colvec b){
  arma::colvec C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - (C_ab));
}

arma::mat lbeta_normconst_mat_cpp(arma::mat a, arma::mat b){
  arma::mat C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - C_ab);
}

//------------------------------------------------------------------------------
// MAIN FUNCTIONS
//------------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_p_X_Global(const arma::field<arma::mat> XI_ijl,
                       const arma::field<arma::mat> RHO_jtk,
                       const arma::field<arma::mat> Y_grouped,
                       const arma::mat mk,
                       const arma::colvec lambdak,
                       const arma::colvec ck,
                       const arma::cube Dk){

  int J = Y_grouped.n_elem;
  int p = Y_grouped(0).n_cols;
  int L = ck.n_rows;
  int H = RHO_jtk(0).n_rows;
  
  arma::vec CSD(L);
  arma::vec ELDL(L);
  
  for(int k=0; k < L; k++){
    CSD(k)  = Const_sum_digamma(p, ck(k));
    ELDL(k) = ElogDetLambda(Dk.slice(k), CSD(k)); // Corresponds to \ell_k
  }
  arma::vec Z_j(J);
  double Z = 0.0;
  for(int j = 0; j < J; j++){
    int N_j = Y_grouped(j).n_rows;
    arma::mat sub_Y = Y_grouped(j);
      for(int ii = 0; ii < N_j; ii++){
        for(int t = 0; t < H; t++){
          for(int k = 0; k < L; k++){
         
          Z += (XI_ijl(j)(ii, t) * RHO_jtk(j)(t, k))* ( ELDL(k) - 
            (p/lambdak(k)) - 
            (p * log( 2.0 * arma::datum::pi)) - 
            (ck(k) * arma::dot(sub_Y.row(ii).t() - mk.col(k), (Dk.slice(k) *  (sub_Y.row(ii).t() - mk.col(k)))))
          );
        }
      }
    }
  }
  
  return(0.5 * Z);

}

// [[Rcpp::export]]
double elbo_p_X_Global_fast(const arma::field<arma::mat> XI_ijl,
                       const arma::field<arma::mat> RHO_jtk,
                       const arma::mat Xbar_k,
                       const arma::cube Sk,
                       const arma::mat mk,
                       const arma::colvec lambdak,
                       const arma::colvec ck,
                       const arma::cube Dk){
  
  int J = XI_ijl.n_elem;
  int p = Xbar_k.n_rows;
  int L = ck.n_rows;
  int H = RHO_jtk(0).n_rows;
  
  arma::vec CSD(L);
  arma::vec ELDL(L);
  
  for(int k=0; k < L; k++){
    CSD(k)  = Const_sum_digamma(p, ck(k));
    ELDL(k) = ElogDetLambda(Dk.slice(k), CSD(k)); // Corresponds to \ell_k
  }
  
  double p1 = 0.0;
  double p2 = 0.0;
  double p3 = 0.0;
  arma::mat tilde_N_j_tilde_RHO_j_mat(J, L, arma::fill::zeros);
  
  for(int j = 0; j < J; j++){
    arma::rowvec N_jt = (arma::sum(XI_ijl(j),0));
    arma::rowvec tilde_N_j_tilde_RHO_j = N_jt * RHO_jtk(j);
    tilde_N_j_tilde_RHO_j_mat.row(j) = tilde_N_j_tilde_RHO_j;
    arma::rowvec tilde_N_j_tilde_RHO_j_ell = tilde_N_j_tilde_RHO_j % ELDL.t();
    arma::rowvec tilde_N_j_tilde_RHO_j_lambda_k = tilde_N_j_tilde_RHO_j % (1/lambdak).t();
    p1 += arma::accu(tilde_N_j_tilde_RHO_j_ell);
    p2 += arma::accu(tilde_N_j_tilde_RHO_j);
    p3 += arma::accu(tilde_N_j_tilde_RHO_j_lambda_k);
  }
  
  arma::rowvec SUM_j_tilde_N_j_tilde_RHO_j_mat = (arma::sum(tilde_N_j_tilde_RHO_j_mat,0));
  double p4 = 0.0;
  for(int k = 0; k<L; k++){
    p4 += (ck(k) * arma::trace(Dk.slice(k) * Sk.slice(k))) + 
      (ck(k) * SUM_j_tilde_N_j_tilde_RHO_j_mat(k)) * arma::dot( (Xbar_k.col(k) - mk.col(k)), (Dk.slice(k) * (Xbar_k.col(k) - mk.col(k))) );
  }

  double Z = p1 - ((p * log( 2.0 * arma::datum::pi)) * p2 ) - (p * p3) - p4;
  return(0.5 * Z);
}

// [[Rcpp::export]]
double elbo_p_X_Local(const arma::field<arma::mat> XI_ijl,
                      const arma::field<arma::mat> Y_j_grouped,
                      const arma::field<arma::mat> m_jt,
                      const arma::field<arma::colvec> lambda_jt,
                      const arma::field<arma::colvec> c_jt,
                      const arma::field<arma::cube> D_jt){
  
  int J = Y_j_grouped.n_elem;
  int H = c_jt(0).n_rows;
  
  arma::vec Z_j(J);
  for(int j = 0; j < J; j++ ){
    int pj = Y_j_grouped(j).n_cols;
    int N_j = Y_j_grouped(j).n_rows;
    
    arma::vec CSD(H);
    arma::vec ELDL(H);
    
    for(int t=0; t < H; t++){
      CSD(t)  = Const_sum_digamma(pj, c_jt(j)(t));
      ELDL(t) = ElogDetLambda(D_jt(j).slice(t), CSD(t)); // Corresponds to \ell_jt
    }
    
    double Z = 0.0;
    for(int t = 0; t < H; t++){
      for(int ii = 0; ii < N_j; ii++){
        Z += 0.5 * ( (XI_ijl(j)(ii, t) * ELDL(t)) -
          ((pj/lambda_jt(j)(t)) * (XI_ijl(j)(ii, t))) -
          (pj * log( 2.0 * arma::datum::pi) * (XI_ijl(j)(ii, t))) -
          (c_jt(j)(t) * (XI_ijl(j)(ii, t)) * arma::dot(Y_j_grouped(j).row(ii).t() - m_jt(j).col(t), (D_jt(j).slice(t) *  (Y_j_grouped(j).row(ii).t() - m_jt(j).col(t)))))
        );
      }
    }
    Z_j(j) = Z;
  }
  
  return(arma::accu(Z_j));
  
}

// [[Rcpp::export]]
double elbo_p_X_Local_fast(const arma::field<arma::mat> XI_ijl,
                           const arma::field<arma::cube> S_jt,
                           const arma::field<arma::mat> Xbar_jt,
                           const arma::field<arma::mat> m_jt,
                           const arma::field<arma::colvec> lambda_jt,
                           const arma::field<arma::colvec> c_jt,
                           const arma::field<arma::cube> D_jt){
  
  int J = XI_ijl.n_elem;
  int H = c_jt(0).n_rows;
  
  arma::vec p1(J);
  arma::vec p2(J);
  arma::vec p3(J);
  arma::vec p4(J);
  
  for(int j = 0; j < J; j++ ){
    int pj = D_jt(j).n_cols;
    
    arma::vec CSD(H);
    arma::vec ELDL(H);
    
    for(int t=0; t < H; t++){
      CSD(t)  = Const_sum_digamma(pj, c_jt(j)(t));
      ELDL(t) = ElogDetLambda(D_jt(j).slice(t), CSD(t)); // Corresponds to \ell_jt
    }
    
    arma::rowvec N_jt = (arma::sum(XI_ijl(j),0));
    
    p1(j) = arma::accu(N_jt % ELDL.t());
    p2(j) = log( 2.0 * arma::datum::pi) * pj * arma::accu(N_jt);
    p3(j) = pj * arma::accu(N_jt / lambda_jt(j).t());
    
    double temp = 0.0;
    for(int t = 0; t < H; t++){
      temp += c_jt(j)(t) * (arma::trace(D_jt(j).slice(t) * S_jt(j).slice(t)) +
        N_jt(t) * arma::dot( (Xbar_jt(j).col(t) - m_jt(j).col(t)),  D_jt(j).slice(t) * (Xbar_jt(j).col(t) - m_jt(j).col(t))  )
      );
    }
    p4(j) = temp;
  }
  
  double Z = arma::accu(p1) - arma::accu(p2) - arma::accu(p3) - arma::accu(p4); 
  return(0.5 * Z);
  
}


// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_diff_v(arma::colvec a_tilde_k,
                   arma::colvec b_tilde_k,
                   const arma::colvec R_concDP){
  // HERE S_concDP contains r1 and r2
  int L = a_tilde_k.n_rows;
  
  a_tilde_k.shed_row(L-1);
  b_tilde_k.shed_row(L-1);
  
  arma::colvec Y =
    E_log_beta(b_tilde_k,a_tilde_k) * ( (R_concDP[0]/R_concDP[1]) - 1.0 );
  
  double p1 = (L-1.0)  * ( R::digamma(R_concDP[0]) - log(R_concDP[1]) ) + arma::accu(Y);  
  double p2 = arma::accu(lbeta_normconst_vec_cpp(a_tilde_k,b_tilde_k) +
                         E_log_beta(a_tilde_k,b_tilde_k) % (a_tilde_k - 1.0 )+
                         E_log_beta(b_tilde_k,a_tilde_k) % (b_tilde_k - 1.0 ));
  
  return(p1-p2);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_diff_u(arma::mat a_bar_Ult,
                   arma::mat b_bar_Ult,
                   const arma::colvec S_concDP){
  
  int J = a_bar_Ult.n_rows;
  int H = a_bar_Ult.n_cols;
  
  a_bar_Ult.shed_col(H-1);
  b_bar_Ult.shed_col(H-1);
  
  arma::colvec R(J);
  
  for(int j = 0; j < J; j++){
    R[j] =  arma::accu( E_log_beta(b_bar_Ult.row(j).t(),a_bar_Ult.row(j).t()) );
    
  }
  
  double p1 = (J * (H-1.0))  * ( R::digamma(S_concDP[0]) - log(S_concDP[1]) ) + ( (S_concDP[0]/S_concDP[1]) - 1.0 ) * arma::accu(R);

  arma::colvec S(J);

  for(int j = 0; j < J; j++){

    S[j] =  arma::accu(lbeta_normconst_vec_cpp(a_bar_Ult.row(j).t(),b_bar_Ult.row(j).t()) +
      E_log_beta(a_bar_Ult.row(j).t(),b_bar_Ult.row(j).t()) % (a_bar_Ult.row(j).t() - 1.0 )+
      E_log_beta(b_bar_Ult.row(j).t(),a_bar_Ult.row(j).t()) % (b_bar_Ult.row(j).t() - 1.0 ));

  }
  double p2 = arma::accu(S);
  return(p1-p2);
}

// [[Rcpp::export]]
double elbo_diff_alpha(const arma::colvec conc_hyper,
                       const arma::colvec S_concDP){
  double pa =
    conc_hyper[0] * log(conc_hyper[1]) - lgamma(conc_hyper[0]) +
    (conc_hyper[0] - 1.0) * ( R::digamma(S_concDP[0]) - log(S_concDP[1]) ) -
    conc_hyper[1] * S_concDP[0]/S_concDP[1];
  
  double qa = S_concDP[0] * log(S_concDP[1]) - lgamma(S_concDP[0]) +
    (S_concDP[0] - 1.0) * (R::digamma(S_concDP[0])-log(S_concDP[1])) -
    S_concDP[0];
  
  return(pa - qa);
}

// [[Rcpp::export]]
double elbo_diff_gamma(const arma::colvec conc_hyper,
                       const arma::colvec R_concDP){
  double pa =
    conc_hyper[0] * log(conc_hyper[1]) - lgamma(conc_hyper[0]) +
    (conc_hyper[0] - 1.0) * ( R::digamma(R_concDP[0]) - log(R_concDP[1]) ) -
    conc_hyper[1] * R_concDP[0]/R_concDP[1];
  
  double qa = R_concDP[0] * log(R_concDP[1]) - lgamma(R_concDP[0]) +
    (R_concDP[0] - 1.0) * (R::digamma(R_concDP[0])-log(R_concDP[1])) -
    R_concDP[0];
  
  return(pa - qa);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_diff_k(const arma::field<arma::mat> RHO_jtk,
                   const arma::colvec ElnPI){
  
  int J = RHO_jtk.n_elem;
  arma::colvec Z1_Z2(J);
  
  for(int j = 0; j < J; j++){
    arma::colvec mdot_k = arma::sum(RHO_jtk(j), 0).t();
    double Z1 = arma::accu(mdot_k % ElnPI);
    
    double Z2 = arma::accu(RHO_jtk(j) % log(RHO_jtk(j) + 1e-12));
    Z1_Z2[j] = Z1 - Z2;
  }

  double Z = arma::accu(Z1_Z2);
  return(Z);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_diff_t(const arma::field<arma::mat> XI_jil,
                   const arma::mat ElnOM_lk){
  
  int J = ElnOM_lk.n_rows;
  
  arma::colvec Z1_Z2(J);
  
  for(int j = 0; j < J; j++){
    arma::colvec mdot_k = arma::sum(XI_jil(j), 0).t();
    double Z1 = arma::accu(mdot_k % ElnOM_lk.row(j).t());
    
    double Z2 = arma::accu(XI_jil(j) % log(XI_jil(j) + 1e-12));
    Z1_Z2[j] = Z1 - Z2;
  }
  
  double Z = arma::accu(Z1_Z2);
  return(Z);
}

// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_p_Psi(const arma::field<arma::colvec> m_j0,
                  const arma::colvec tau_j0,
                  const arma::colvec nu_j0,
                  const arma::field<arma::mat> W_j0,
                  const arma::field<arma::mat> iW_j0,
                  const arma::vec lCpl0, arma::vec HlogB0,
                  const arma::field<arma::mat> m_jt, arma::field<arma::colvec> lambda_jt,
                  const arma::field<arma::colvec> c_jt, arma::field<arma::cube> D_jt){
  
  int J = m_jt.n_elem;
  arma::vec Z_j(J);
  int H = lambda_jt(0).n_rows;
  
  for(int j = 0; j < J; j ++){
    int pj = W_j0(j).n_rows;
    
    
    arma::vec CSD(H);
    arma::vec LCPL(H);
    
    double p1   = 0.0;
    double ell1 = 0.0;
    double p4   = 0.0;
    
    
    for(int t = 0; t<H; t++){
      
      CSD(t)  = Const_sum_digamma(pj, c_jt(j)(t));
      LCPL(t) = log_Const_prod_gamma(pj, c_jt(j)(t));
      
      ell1   += ElogDetLambda(D_jt(j).slice(t), CSD(t));
      p1     += c_jt(j)(t) * arma::dot((m_jt(j).col(t)-m_j0(j)) ,
                     D_jt(j).slice(t) * ((m_jt(j).col(t)-m_j0(j))));
      p4     += c_jt(j)(t) * arma::trace( iW_j0(j) * D_jt(j).slice(t) );
      
    }
    
    double p0 = .5* (H * pj * log( tau_j0(j)/(2.0*arma::datum::pi) ) +
                     ell1 -
                     pj * tau_j0(j) * arma::accu(1.0/(lambda_jt(j))) -
                     tau_j0(j) * p1 );
    
    
    double  Z =  HlogB0(j) + p0 + (nu_j0(j) - pj - 1.0) * 0.5 * ell1 - 0.5 * p4;
    Z_j(j) = Z;
  }
  
  
  return(arma::accu(Z_j));
}

// [[Rcpp::export]]
double elbo_p_Phi(const arma::colvec m0,
                  const double tau0,
                  const double nu0,
                  const arma::mat W0,
                  const arma::mat iW0,
                  const double lCpl0, const double LlogB0,
                  const arma::mat mk, const arma::colvec lambdak,
                  const arma::colvec ck, const arma::cube Dk){
  
  int p = W0.n_rows;
  int L = lambdak.n_rows;
  
  arma::vec CSD(L);
  arma::vec LCPL(L);
  
  double p1   = 0.0;
  double ell1 = 0.0;
  double p4   = 0.0;
  
  
  for(int k = 0; k<L; k++){
    
    CSD(k)  = Const_sum_digamma(p, ck(k));
    LCPL(k) = log_Const_prod_gamma(p, ck(k));
    
    ell1   += ElogDetLambda(Dk.slice(k), CSD(k));
    p1     += ck(k) * arma::dot((mk.col(k)-m0) ,
                 Dk.slice(k) * ((mk.col(k)-m0)));
    p4     += ck(k) * arma::trace( iW0 * Dk.slice(k) );
    
  }
  
  double p0 = .5* (L * p * log( tau0/(2.0*arma::datum::pi) ) +
                   ell1 -
                   p * tau0 * arma::accu(1.0/(lambdak)) -
                   tau0 * p1 );
  
  
  double  Z =  LlogB0 + p0 + (nu0 - p - 1.0) * 0.5 * ell1 - 0.5 * p4;
  
  return(Z);
}

// ----------------------------------------------------------------------- 10.77
// [[Rcpp::export]]
double elbo_q_Phi_k(const arma::mat mk, const arma::colvec lambdak,
                    const arma::colvec ck, const arma::cube Dk){
  
  int p = Dk.slice(0).n_rows;
  int L = ck.n_rows;
  arma::vec CSD(L);
  arma::vec ELDL(L);
  arma::vec H(L);
  
  for(int k = 0; k<L; k++){
    CSD(k)  = Const_sum_digamma(p, ck(k));
    ELDL(k) = ElogDetLambda(Dk.slice(k), CSD(k));
    H(k)    = H_Lambda(Dk.slice(k),
      ck(k), CSD(k),
      log_Const_prod_gamma(p,ck(k)));
  }
  
  arma::vec Z = .5 * ELDL +
    p * 0.5 * ( log( lambdak / (2 * arma::datum::pi)) - 1.0) -
    H; 
  
  return(arma::accu(Z));
}

// ----------------------------------------------------------------------- 10.77
// [[Rcpp::export]]
double elbo_q_Psi_k(const arma::field<arma::mat> m_jt, const arma::field<arma::colvec> lambda_jt,
                    const arma::field<arma::colvec> c_jt, const arma::field<arma::cube> D_jt){
  
  int J = m_jt.n_elem;
  arma::vec E_q(J);
  
  for(int j = 0; j < J; j++){
    int pj = D_jt(j).slice(0).n_rows;
    int H = c_jt(j).n_rows;
    
    arma::vec CSD(H);
    arma::vec ELDL(H);
    arma::vec H_Wishart(H);
    
    for(int t = 0; t < H; t++){
      CSD(t)  = Const_sum_digamma(pj, c_jt(j)(t));
      ELDL(t) = ElogDetLambda(D_jt(j).slice(t), CSD(t));
      H_Wishart(t)    = H_Lambda(D_jt(j).slice(t),
                c_jt(j)(t), CSD(t),
                log_Const_prod_gamma(pj,c_jt(j)(t)));
    }
    
    arma::vec Z = .5 * ELDL +
      pj * 0.5 * ( log( lambda_jt(j) / (2 * arma::datum::pi)) - 1.0) -
      H_Wishart; 
    
    E_q[j] = arma::accu(Z);
  }
  
  
  return(arma::accu(E_q));
}
