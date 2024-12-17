library(Rcpp)
sourceCpp("VI_functions.cpp")
sourceCpp("ELBO_functions.cpp")
#########################################################################################
H = 30; L = 30; # Set TRUNCATION LEVELS
n = unlist(lapply(X, ncol)) # Set SAMPLE SIZES
J = length(X) # Set Size of NUMBER OF GROUPS
p.local = unlist(lapply(X.local, ncol)) # COLLATE THE DIMENSION OF LOCAL VARIABLES


m0_global = matrix(0, nrow = p.global); tau0_global = 0.01; nu0_global = p.global + 5; W0_global = diag(p.global); iW0_global = solve(W0_global)

tau_j0 = matrix(0.01, nrow = J, ncol = 1)
nu_j0 = matrix(c(p.local + 5), nrow = J, ncol = 1)
m_j0 = list()
W_j0 = list()
iW_j0 = list()
for(j in 1:J){
  m_j0[[j]] = matrix(0, nrow = p.local[j], ncol = 1)
  W_j0[[j]] = diag(1, nrow = p.local[j])
  iW_j0[[j]] = solve(W_j0[[j]])
}

s1 = 1; s2 = 2; b_jt_bar = s1/s2
r1 = 1; r2 = 2; b_k_bar = r1/r2  

a_alpha = 0.1; b_alpha = 0.1
a_gamma = 0.1; b_gamma = 0.1

out = list()

nReplicate = 20

for(r in 1:nReplicate){
  
  cat(paste0("\nRunning VI Replicate: ", r, "\n\n"))
  ###############################################################################################
  XI_jil = list()
  for(j in 1:J){
    log.XI_il  <- array(stats::rbeta( n[j] * H, 1, 1),dim = c( n[j], H))
    Z           <- apply(log.XI_il, c(1), function(x) matrixStats::logSumExp(x))
    XI_jil[[j]] <- exp(sapply(1:H, function(qq) log.XI_il[,qq]-Z,simplify = "array"))
  }
  ###############################################################################################
  RHO_jtk = list()
  for(j in 1:J){
    log.RHO_jtk  <- array(stats::rbeta( H * L, 1, 1),dim = c( H, L))
    Z           <- apply(log.RHO_jtk, c(1), function(x) matrixStats::logSumExp(x))
    RHO_jtk[[j]] <- exp(sapply(1:L, function(qq) log.RHO_jtk[,qq]-Z,simplify = "array"))
  }
  
  maxIter = 500
  elbo = 0
  epsilon = 1e-5
  diff_ELBO = 1
  # time at the beginning
  Time.start = Sys.time()
  for(iter in 1:maxIter){
    # iter = 1
    if(iter == 1){
      cat(paste0("Iteration: ", iter, "\n"))
    }
    if(iter %% floor((10/100)*(maxIter)) == 0) {
      cat(paste0("Iteration: ", iter, "\n"))
    }
    
    ###############################################################################################
    # First Update v_k ~ Beta(\bar{a}_k, \bar{b}_k), for k = 1,.., L-1
    ###############################################################################################
    var_par_v = Update_Vk_cpp(b_bar = b_k_bar,
                              RHO_jtk = RHO_jtk,
                              J = J)
    
    a_vk     = var_par_v[, 1, drop = FALSE];
    b_vk     = var_par_v[, 2, drop = FALSE];
    E_ln_PIk = var_par_v[, 3, drop = FALSE];
    
    ###############################################################################################
    # Second Update u_jt ~ Beta(\bar{a}_jt, \bar{b}_jt), for j = 1,..,J; t = 1,.., T-1
    ###############################################################################################
    XI_jil = Jitter_XIjil_cpp(XI_jil = XI_jil)
    
    UU = Update_Ujt_cpp(XI_jil = XI_jil,
                        b_bar = b_jt_bar,
                        T = H,
                        J = J)
    
    a_bar_Ujt = UU[ , ,1] 
    b_bar_Ujt = UU[ , ,2]
    ElnOM_jt  = UU[ , ,3]
    
    # ElnOM_jt
    ###############################################################################################
    # Third Update alpha ~ Gamma(s1, s2)
    ###############################################################################################
    Update_alpha = Update_alpha_concentration_par(a_bar_Ulk = a_bar_Ujt,
                                                  b_bar_Ulk = b_bar_Ujt,
                                                  conc_hyper = matrix(c(a_alpha, b_alpha), ncol = 1),
                                                  T = H,
                                                  J = J)
    
    b_jt_bar = Update_alpha[1,]/Update_alpha[2,]
    
    ###############################################################################################
    # Fourth Update gamma ~ Gamma(r1, r2)
    ###############################################################################################
    Update_gamma = Update_gamma_concentration_par(a_tilde_Vk = a_vk,
                                                  b_tilde_Vk = b_vk,
                                                  conc_hyper = matrix(c(a_gamma, b_gamma), ncol = 1))
    
    b_k_bar = Update_gamma[1,]/Update_gamma[2,]
    
    ###############################################################################################
    # Fifth Update phi_k ~ NW(m_k, lambda_k, c_k, D_k)
    ###############################################################################################
    out.global = Update_Phi_k_cpp_mvt(Y_grouped = X.global,
                                      RHO_jtk = RHO_jtk,
                                      XI_jil = XI_jil,
                                      m0 = m0_global,
                                      tau0 = tau0_global,
                                      nu0 = nu0_global,
                                      W0 = W0_global,
                                      iW0 = iW0_global)
    
    mk = out.global$mk; lambdak =  out.global$lambdak; ck = out.global$ck; Dk = out.global$Dk
    Xbar_k = out.global$Xbar_k; Sk = out.global$Sk
    ###############################################################################################
    # Sixth Update k_jt ~ Multinomial(RHO_jtk), for j=1,...,J; t = 1,..,T; k = 1,...,L
    ###############################################################################################
    RHO_jtk = Update_RHOjtk_cpp(XI_jil = XI_jil,
                                Y_grouped = X.global, 
                                mk = mk,
                                lambdak =  lambdak,
                                ck = ck,
                                Dk = Dk,
                                ElnPI_k = E_ln_PIk, 
                                L = L,
                                J = J,
                                T = H)
    
    ###############################################################################################
    # Seventh Update Psi_jt ~ NW(m_jt, lambda_jt, c_jt, D_jt)
    ###############################################################################################
    out_Psi = Update_Psi_jt_cpp_mvt(Y_j_grouped = X.local,
                                    XI_jil = XI_jil,
                                    m_j0 = m_j0,
                                    tau_j0 = tau_j0 ,
                                    nu_j0 = nu_j0,
                                    W_j0 = W_j0,
                                    iW_j0 = iW_j0
    )
    
    lambda_jt = out_Psi$lambda_jt
    c_jt = out_Psi$c_jt
    m_jt = out_Psi$m_jt
    D_jt = out_Psi$D_jt
    S_jt = out_Psi$S_jt
    Xbar_jt = out_Psi$Xbar_jt
    ###############################################################################################
    # Eigth Update t_ji ~ Multinomial(XI_jit), for j=1,...,J; i = 1,..,n_j; t = 1,...,T
    ###############################################################################################
    XI_jil = Update_XIjil_cpp(Y_grouped = X.global,
                              mk = mk,
                              lambdak = lambdak,
                              ck = ck,
                              Dk = Dk,
                              
                              Y_j_grouped = X.local,
                              m_jt = m_jt,
                              lambda_jt = lambda_jt,
                              c_jt = c_jt,
                              D_jt = D_jt,
                              
                              RHO_jtk = RHO_jtk,
                              ElnOM_lk = ElnOM_jt,
                              
                              L = L,
                              J = J,
                              T = H)
    
    ###############################################################################################
    # CALCULATION OF ELBO
    ###############################################################################################
    elbo1 = elbo_diff_v(a_tilde_k = a_vk,
                        b_tilde_k = b_vk,
                        R_concDP = Update_gamma)
    
    elbo2 = elbo_diff_u(a_bar_Ult = a_bar_Ujt,
                        b_bar_Ult = b_bar_Ujt,
                        S_concDP = Update_alpha)
    
    elbo3 = elbo_diff_gamma(conc_hyper = matrix(c(a_gamma, b_gamma), ncol = 1),
                            R_concDP = Update_gamma)
    
    elbo4 = elbo_diff_alpha(conc_hyper = matrix(c(a_alpha, b_alpha), ncol = 1),
                            S_concDP = Update_alpha)
    
    elbo5 = elbo_diff_k(RHO_jtk = RHO_jtk,
                        ElnPI = E_ln_PIk)
    
    elbo6 = elbo_diff_t(XI_jil = XI_jil,
                        ElnOM_lk = ElnOM_jt)
    
    elbo7 = elbo_q_Phi_k(mk = mk, lambdak = lambdak,
                         ck = ck, Dk = Dk)
    
    elbo8 = elbo_q_Psi_k(m_jt = m_jt, lambda_jt = lambda_jt,
                         c_jt = c_jt, D_jt = D_jt)
    
    
    lCpl0  = log_Const_prod_gamma(D = p.global, nu = nu0_global)
    LlogB0 = L * logB_wish(W = W0_global , nu = nu0_global, log_Const_prod_gamma = lCpl0);
    
    elbo9 = elbo_p_Phi(m0 = m0_global,
                       tau0 = tau0_global,
                       nu0 = nu0_global,
                       W0 = W0_global,
                       iW0 = iW0_global,
                       lCpl0 = lCpl0, LlogB0 = LlogB0,
                       mk = mk, lambdak = lambdak,
                       ck = ck, Dk = Dk)
    
    lCpl0 = 0
    HlogB0 = 0
    for(j in 1:J){
      lCpl0[j]  = log_Const_prod_gamma(D = p.local[j], nu = nu_j0[j])
      HlogB0[j] = H * logB_wish(W = W_j0[[j]] , nu = nu_j0[j], log_Const_prod_gamma = lCpl0[j]);
    }
    
    elbo10 = elbo_p_Psi(m_j0 = m_j0,
                        tau_j0 = tau_j0,
                        nu_j0 = nu_j0,
                        W_j0 = W_j0,
                        iW_j0 = iW_j0,
                        lCpl0 = lCpl0, HlogB0 = HlogB0,
                        m_jt = m_jt, lambda_jt = lambda_jt,
                        c_jt = c_jt, D_jt = D_jt )
    
    elbo11 = elbo_p_X_Global_fast(XI_ijl = XI_jil,
                                  RHO_jtk = RHO_jtk,
                                  Xbar_k = Xbar_k,
                                  Sk = Sk,
                                  mk = mk,
                                  lambdak = lambdak,
                                  ck = ck,
                                  Dk = Dk)
    
    elbo12 = elbo_p_X_Local_fast(XI_ijl = XI_jil,
                                 S_jt = S_jt,
                                 Xbar_jt = Xbar_jt,
                                 m_jt = m_jt,
                                 lambda_jt = lambda_jt,
                                 c_jt = c_jt,
                                 D_jt = D_jt)
    
    elbo_value = elbo1 +
      elbo2 + 
      elbo3 + 
      elbo4 + 
      elbo5 + 
      elbo6 + 
      elbo7 +
      elbo8 + 
      elbo9 + 
      elbo10 + 
      elbo11 + 
      elbo12 
    ## STORE THE ELBO VALUE
    elbo[iter] = elbo_value
    if(iter == 1){diff_ELBO = 1}
    if(iter > 1){
      diff_ELBO = (elbo[iter]-elbo[iter - 1])}
    
    if( (abs(diff_ELBO) > epsilon) & (iter == maxIter) ) {
      Time.not.converged = Sys.time()
      Time.taken = round(Time.not.converged - Time.start, 3)
      Status = "Not Converged"
      cat(paste0("Warning! Maximum Number of Iterations have reached and ELBO has not yet converged\nThe difference in ELBO at stopping time is ", abs(diff_ELBO), "\nEBLO value at stop = ", elbo[iter], "\nTime taken by algorithm: ", Time.taken," ", attr(Time.taken, "units")[1], "\n"))
    }
    if((abs(diff_ELBO) < epsilon) & (iter < maxIter)){
      # time at the convergence
      Time.converged = Sys.time()
      Time.taken = round(Time.converged - Time.start, 3)
      Status = "Converged"
      cat(paste0("Success! Increase in ELBO < ", epsilon, ". Algorithm has converged in ",iter," iterations\n", "EBLO value at convergence = ", elbo[iter], "\nTime taken to converge: ", Time.taken," ", attr(Time.taken, "units")[1], "\n"))
      break
    }
    
  }
  out[[r]] = list(ELBO_val_conv = elbo[iter],
                  ELBO = elbo,
                  
                  a_vk     = a_vk,
                  b_vk     = b_vk,
                  E_ln_PIk = E_ln_PIk,
                  
                  a_bar_Ujt = a_bar_Ujt,
                  b_bar_Ujt = b_bar_Ujt,
                  ElnOM_jt  = ElnOM_jt,
                  
                  s1 = Update_alpha[1,],
                  s2 = Update_alpha[2,],
                  
                  r1 = Update_gamma[1,],
                  r2 = Update_gamma[2,],
                  
                  mk = mk, lambdak = lambdak, ck = ck, Dk = Dk,
                  m_jt = m_jt, lambda_jt = lambda_jt, c_jt = c_jt, D_jt = D_jt,
                  RHO_jtk = RHO_jtk,
                  XI_jil = XI_jil,
                  Time.taken = Time.taken,
                  Status = Status)
}

ELBO = 0
Time = 0

for(i in 1:length(out)){
  ELBO[i] = out[[i]]$ELBO_val_conv
  Time[i] = out[[i]]$Time.taken
}

best_outputs = out[[which.max(ELBO)]]

est.k_jt = list()
for(j in 1:J){
  est.k_jt[[j]] = apply(best_outputs$RHO_jtk[[j]], 1, which.max)
}

est.t_ji = list()
for(j in 1:J){
  est.t_ji[[j]] = apply(best_outputs$XI_jil[[j]], 1, which.max)
}

DATA_Est = list()
for(j in 1:J){
  DATA_Est[[j]] = data.frame(X.global[[j]],
                             Cluster.true = factor(k.true[[j]][t.true[[j]]]),
                             Cluster.est = factor(est.k_jt[[j]][est.t_ji[[j]]]),
                             Population.org = paste0("Population ", j),
                             Population = paste0("Population ", j, "\n", "ARI = ", round(aricode::ARI(k.true[[j]][t.true[[j]]], est.k_jt[[j]][est.t_ji[[j]]]), 4)))
}


DATA.global_Est = NULL
for(j in 1:J){
  DATA.global_Est = rbind(DATA.global_Est, DATA_Est[[j]])
}


library(pals)
myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (L - 22))))

names(myvalues) = 1:L

library(tidyverse)
True_plot <- DATA.global_Est %>% ggplot(aes(x = X1, y = X2, col = Cluster.true)) + geom_point(size = 3) + labs(title = paste0("True global-level clusters")) + scale_color_manual(values = myvalues) + facet_grid(~Population.org) + 
  theme_minimal() +  
  theme(
    # LABLES APPEARANCE
    # panel.grid.major = element_blank(), 
    # panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3), 
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) + guides(color = guide_legend(title = "Clusters"))

Est_plot <- DATA.global_Est %>% ggplot(aes(x = X1, y = X2, col = Cluster.est)) + geom_point(size = 3) + labs(title = paste0("Estimated global-level clusters")) +  
                                                                                                             # subtitle = paste0("ELBO Value at convergence = ", round(best_outputs$ELBO_val_conv, 4))) 
   scale_color_manual(values = myvalues) + facet_grid(~Population) + 
  theme_minimal() +  
  theme(
    # LABLES APPEARANCE
    # panel.grid.major = element_blank(), 
    # panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3), 
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) + guides(color = guide_legend(title = "Clusters"))

gridExtra::grid.arrange(True_plot, Est_plot)

ARI_Local = 0
ARI_Global = 0
for(j in 1:J){
  ARI_Local[j] = aricode::ARI(est.t_ji[[j]], t.true[[j]])
  ARI_Global[j] = aricode::ARI(k.true[[j]][t.true[[j]]], est.k_jt[[j]][est.t_ji[[j]]])
}

ARI_Global
ARI_Local

aricode::ARI(unlist(lapply(1:J, function(j){ c(k.true[[j]][t.true[[j]]]) })),
             unlist(lapply(1:J, function(j){ c(est.k_jt[[j]][est.t_ji[[j]]]) }))
)

DATA_Est_Local = list()
for(j in 1:J){
  DATA_Est_Local[[j]] = data.frame(X.local[[j]][ , c(1:2)],
                                   Cluster.true = factor(t.true[[j]]),
                                   Cluster.est = factor(est.t_ji[[j]]),
                                   Population.org = paste0("Population ", j),
                                   Population = paste0("Population ", j, "\n", "ARI = ", round(aricode::ARI(t.true[[j]], est.t_ji[[j]]), 4)))
}


DATA.local_Est = NULL
for(j in 1:J){
  DATA.local_Est = rbind(DATA.local_Est, DATA_Est_Local[[j]])
}


library(pals)
myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (L - 22))))

names(myvalues) = 1:L

library(tidyverse)
True_plot_Local <- DATA.local_Est %>% ggplot(aes(x = X1, y = X2, col = Cluster.true)) + geom_point(size = 3) + labs(title = paste0("True local-level clusters")) + scale_color_manual(values = myvalues) + facet_grid(~Population.org) + 
  theme_minimal() +  
  theme(
    # LABLES APPEARANCE
    # panel.grid.major = element_blank(), 
    # panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3), 
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) + guides(color = guide_legend(title = "Clusters"))

Est_plot_Local <- DATA.local_Est %>% ggplot(aes(x = X1, y = X2, col = Cluster.est)) + geom_point(size = 3) + labs(title = paste0("Estimated local-level clusters")) + scale_color_manual(values = myvalues) + facet_grid(~Population) + 
  theme_minimal() +  
  theme(
    # LABLES APPEARANCE
    # panel.grid.major = element_blank(), 
    # panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3), 
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) + guides(color = guide_legend(title = "Clusters"))

plot_Global <- gridExtra::grid.arrange(True_plot, Est_plot)
plot_Local <- gridExtra::grid.arrange(True_plot_Local, Est_plot_Local)
library(grid)
gridExtra::grid.arrange(plot_Local, plot_Global, ncol = 2,
                        top = textGrob(paste0("Local variable dimensions = ", p.local[1],", ", p.local[2], ", ", p.local[3], ". Global variable dimension = ", p.global), ,gp=gpar(fontsize=20,font=3)))

