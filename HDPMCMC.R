################################################################################
# LOAD THE FUNCTIONS NEEDED FOR MCMC
################################################################################
library(Rcpp)
sourceCpp("MCMC_functions.cpp")

# Define Some Auxillary Functions
log_sum <- function(x){
  exp(x - max(x))
}

indic <- function(x, t){
  ifelse(x == t, 1, 0)
}

# @pars argument should be a named list with a0, b0, J, pi1, pi2, pi3
# @ycurrent corresponds to the previous value of alpha
sample_alpha <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  pi1 = as.numeric(pars$pi1); pi2 = as.numeric(pars$pi2); pi3 = as.numeric(pars$pi3)
  a0 = pars$a0; b0 = pars$b0
  J = pars$J
  L = length(pi1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = - b0 * ycand + (a0 - 1)*log(ycand) + sum(((ycand/L) - 1) * log(pi1)) + sum((ycand/L - 1) * log(pi2)) + sum((ycand/L - 1) * log(pi3)) + J * lgamma(ycand) - (J * L * lgamma(ycand/L)) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = - b0 * ycurrent + (a0 - 1)*log(ycurrent) + sum(((ycurrent/L) - 1) * log(pi1)) + sum((ycurrent/L - 1) * log(pi2)) + sum((ycurrent/L - 1) * log(pi3)) + J * lgamma(ycurrent) - (J * L * lgamma(ycurrent/L))  + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(list(out = out, accept = accept))
}

# @pars argument should be a named list with a1, b1 and beta
# @ycurrent corresponds to the previous value of m
sample_m <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  a1 = pars$a1; b1 = pars$b1
  beta = pars$beta
  L = length(beta) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a1, rate = b1) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = - b1 * ycand + (a1 - 1)*log(ycand) + sum(((ycand/L) - 1) * log(beta)) + lgamma(ycand) - (L * lgamma(ycand/L)) + dgamma(ycurrent, shape = a1, rate = b1, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = - b1 * ycurrent + (a1 - 1)*log(ycurrent) + sum(((ycurrent/L) - 1) * log(beta)) + lgamma(ycurrent) - (L * lgamma(ycurrent/L)) + dgamma(ycand, shape = a1, rate = b1, log = TRUE)
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(list(out = out, accept = accept))
}

L = 10
for(j in 1:J){
  X.global[[j]] = t(X.global[[j]])
}

# Prior hyper-parameters corresponding to global parameters
lambda0 = 0.01; nu0 = p.global + 5; m0 = rep(0, p.global); Psi0 = diag(1, p.global); iPsi0 = solve(Psi0)  

alpha.start = 2
m.start = 2

library(extraDistr)
# True weights
beta.start = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting beta
pi.start = replicate(n = J, list())
t.start = replicate(n = J, list())
k.start = replicate(n = J, list())

for(j in 1:J){
  pi.start[[j]] = as.numeric(rdirichlet(n = 1, alpha = rep(1/L, L))) # Starting pi_j
  t.start[[j]] = sample(1:L, size = n[j], prob = pi.start[[j]], replace = TRUE)
  k.start[[j]] = sample(1:L, size = L, prob = beta.start, replace = TRUE)
}

num_iter = 50000 # Number of MCMC iterations

m = replicate(n = J, list(0))

t.samples <- replicate(n = J, list(replicate(n = num_iter, list(0))))
k.samples <- replicate(n = J, list(replicate(n = num_iter, list(0))))
pi.samples <- replicate(n = J, list(replicate(n = num_iter, list())))

for(j in 1:J){
  t.samples[[j]][[1]] = t.start[[j]]
  k.samples[[j]][[1]] = k.start[[j]]
  pi.samples[[j]][[1]] = pi.start[[j]]
}

beta.samples <- list(); beta.samples[[1]] <- beta.start

M.k.samples <- list(); Tau2.k.samples <- list() # List to store global parameter samples (means and variances)

alpha.samples <- list(); alpha.samples[[1]] = alpha.start
m.samples <- list(); m.samples[[1]] = m.start

################################################################################
# SAMPLING 
# ################################################################################
time.start = Sys.time()
for(iter in 2:num_iter){
  # iter = 2
  # Printing the iterations
  if(iter == 2){
    cat(paste0("Iteration: ", iter-1, "\n"))
  }
  if(iter %% floor((10/100)*(num_iter + 1)) == 0) {
    cat(paste0("Iteration: ", iter, "\n"))
  }
  
  
  for(j in 1:J){
    for(l in 1:L){
      m[[j]][l] = sum(t.samples[[j]][[iter - 1]] == l)
    }
  }
  
  
  for(j in 1:J){
    pi.samples[[j]][[iter]] = as.numeric(rdirichlet(n = 1, alpha = m[[j]] + alpha.samples[[iter - 1]]/L))
    pi.samples[[j]][[iter]] = (pi.samples[[j]][[iter]] + 1e-8)/sum(pi.samples[[j]][[iter]] + 1e-8)
  }
  
  d = replicate(n = J, list(0))
  for(j in 1:J){
    for(k in 1:L){
      d[[j]][k] = sum(k.samples[[j]][[iter - 1]] == k)
      
    } 
  }
  
  d.all = Reduce(`+`, d)
  
  beta.samples[[iter]] = as.numeric(rdirichlet(n = 1, alpha = d.all + m.samples[[iter - 1]]/L))
  beta.samples[[iter]] = (beta.samples[[iter]] + 1e-8)/sum(beta.samples[[iter]] + 1e-8)
  
  n.k = replicate(n = J, list(0))
  for(j in 1:J){
    for(k in 1:L){
      n.k[[j]][k] = sum(k.samples[[j]][[iter - 1]][t.samples[[j]][[iter - 1]]] == k)
    }
  }
  
  n.k.all = Reduce(`+`, n.k)
  
  tau2.k.samples <- array(0, dim = c(p.global, p.global, L))
  m.k.samples <- matrix(0, p.global, L)
  
  x.k.bar = list()
  
  for(k in 1:L){
    Sum_x_k_temp = matrix(0, p.global, 1)
    if(n.k.all[k] == 0){
      x.k.bar[[k]] = rep(0, nrow(X.global[[1]]))
    }else{
      for(j in 1:J){
        Sum_x_k_temp = Sum_x_k_temp + rowSums(X.global[[j]][, k.samples[[j]][[iter - 1]][t.samples[[j]][[iter - 1]]] == k, drop = FALSE])
      }
      x.k.bar[[k]] = as.numeric(Sum_x_k_temp)/n.k.all[k]
    }
  }
  
  S.k = array(0, dim = c(p.global, p.global, L))
  Z.k = array(0, dim = c(p.global, p.global, L))
  
  lambda.hat = rep(0, L); nu.hat = rep(0, L); m.k.hat = matrix(0, nrow = p.global, ncol = L)
  
  nu.hat     = nu0 + n.k.all
  lambda.hat = lambda0 + n.k.all
  
  for(k in 1:L){
    if(n.k.all[k] != 0){
      S.k.j = array(0, dim = c(p.global, p.global, J))
      for(j in 1:J){
        X.global.subset = X.global[[j]][, k.samples[[j]][[iter - 1]][t.samples[[j]][[iter - 1]]] == k, drop = FALSE]
        n.global.subset = ncol(X.global.subset)
        S.k.j_temp = matrix(0, nrow = p.global, ncol = p.global)
        if(n.global.subset > 0){
          for(ii in 1:n.global.subset){
            S.k.j_temp = S.k.j_temp +  (X.global.subset[, ii, drop = FALSE] - matrix(x.k.bar[[k]], ncol = 1)) %*% t( (X.global.subset[, ii, drop = FALSE] - matrix(x.k.bar[[k]], ncol = 1))  )
          }
        }
        S.k.j[ , ,j] = S.k.j_temp
      }
      S.k_temp = matrix(0, nrow = p.global, ncol = p.global)
      for(j in 1:J){
        S.k_temp = S.k_temp + S.k.j[ , ,j]
      }
      S.k[, , k] = S.k_temp
      Z.k[, , k] = ((lambda0 * n.k.all[k])/(lambda0 + n.k.all[k])) * (matrix(x.k.bar[[k]], ncol = 1) - matrix(m0, ncol = 1)) %*% t( (matrix(x.k.bar[[k]], ncol = 1) - matrix(m0, ncol = 1)) )
    }
    
    m.k.hat[, k] = ((lambda0 * m0) + (n.k.all[k] * x.k.bar[[k]]))/lambda.hat[k]
  }
  
  
  for(k in 1:L){
    tau2.k.samples[, , k] = stats::rWishart(n = 1, df = nu.hat[k], Sigma = solve(iPsi0 + S.k[ , , k] + Z.k[, , k]))
    
    m.k.samples[ , k] =  MASS::mvrnorm(n = 1, mu = m.k.hat[ , k], Sigma = solve(lambda.hat[k] * tau2.k.samples[, , k]) )
  }
  
  Tau2.k.samples[[iter]] = tau2.k.samples
  M.k.samples[[iter]] = m.k.samples

  
  t.prob = replicate(n = J, list())
  for(j in 1:J){
    t.prob[[j]] = matrix(0, nrow = n[j], ncol = L)
  }
  
  iTau2.j.samples = replicate(J, list())
  iTau2.k.samples = array(0, dim = dim(Tau2.k.samples[[iter]]))
  
  for(l in 1:dim(Tau2.k.samples[[iter]])[3]){
    iTau2.k.samples[, , l] = solve(Tau2.k.samples[[iter]][, , l])
  }


  
  for(j in 1:J){
    t.prob[[j]] =   prob_exponent_no_local2(pi_samples = pi.samples[[j]][[iter]], 
                                            phi_samples = M.k.samples[[iter]],
                                            X_global = t(X.global[[j]]), 
                                            k =  k.samples[[j]][[iter - 1]], 
                                            Sigma_G = iTau2.k.samples)
  }
  
  
  
  # Log-sum-exp trick and normalization
  for(j in 1:J){
    t.prob[[j]] = calc_probability_log_sum_exp_normalized(t.prob[[j]])
  }
  
  for(j in 1:J){
    t.samples[[j]][[iter]] = sample_my(t.prob[[j]])
  }
  
  
  # Calculate the probabilities of dish indices and sample
  prob = array(0, dim = c(L, L, J))
  
  for(j in 1:J){
    prob[, , j] = prob_exponent_mv_global(beta_samples = beta.samples[[iter]], 
                                          X_global = t(X.global[[j]]), 
                                          phi_samples = matrix(m.k.samples, nrow = p.global, L), 
                                          Sigma_G = iTau2.k.samples,
                                          t_samples = t.samples[[j]][[iter]])}
  
  # Log-sum-exp trick and normalization
  for(j in 1:J){prob[, ,j] = calc_probability_log_sum_exp_normalized(prob[, ,j])}
  
  for(j in 1:J){
    k.samples[[j]][[iter]] = sample_my(prob[, , j])}
  
  sampled.alpha = sample_alpha(ycurrent = alpha.samples[[iter - 1]], 
                               pars = list(a0 = 0.1, b0 = 0.1,
                                           pi1 = pi.samples[[1]][[iter]],
                                           pi2 = pi.samples[[2]][[iter]],
                                           pi3 = pi.samples[[3]][[iter]],
                                           J = 3))
  alpha.samples[[iter]] = sampled.alpha$out
  
  sampled.m = sample_m(ycurrent = m.samples[[iter - 1]], 
                       pars = list(a1 = 0.1, b1 = 0.1,
                                   beta = beta.samples[[iter]]))
  m.samples[[iter]] = sampled.m$out
  
}
time.end = Sys.time()
time.end - time.start

################################################################################
# CALCULATE LOG-LIKELIHOOD TO PLOT TRACEPLOT OF LOG-LIKELIHOOD
################################################################################
ll = 0

for(iter in 2:num_iter){
  # Printing the iterations
  if(iter == 2){
    cat(paste0("Iteration: ", iter-1, "\n"))
  }
  if(iter %% floor((10/100)*(num_iter + 1)) == 0) {
    cat(paste0("Iteration: ", iter, "\n"))
  }
  
  logsum = 0
  
  iTau2.k.samples = array(0, dim = dim(Tau2.k.samples[[iter]]))
  
  for(l in 1:dim(Tau2.k.samples[[iter]])[3]){
    iTau2.k.samples[, , l] = solve(Tau2.k.samples[[iter]][, , l])
  }
  
  
  
  for(j in 1:J){
    logsum = logsum + logll_no_local2(pi_samples = pi.samples[[j]][[iter]],
                                      t_samples = t.samples[[j]][[iter]],
                                      phi_samples = M.k.samples[[iter]], 
                                      X_global = X.global[[j]],
                                      k_samples = k.samples[[j]][[iter]], 
                                      Sigma_G = iTau2.k.samples) 
  }
  
  m.j.tilde = replicate(J, list)
  for(j in 1:J){
    m.j.tilde[[j]] = count_my(t.samples[[j]][[iter]], L)
  }
  
  sum_m.j_log_pi.j = 0
  for(j in 1:J){
    sum_m.j_log_pi.j = sum_m.j_log_pi.j + sum(m.j.tilde[[j]] * log(pi.samples[[j]][[iter]]))
  }
  
  logsum = logsum + sum_m.j_log_pi.j
  m.tilde = rowSums(sapply(1:J, function(j) { count_my(k.samples[[j]][[iter]], L) }))
  logsum = logsum + sum(m.tilde * log(beta.samples[[iter]]))
  
  ll[iter] = logsum
}

################################################################################
# PLOT TRACEPLOT OF LOG-LIKELIHOOD AND ACF
################################################################################
burn = 30000 # Burn-in 
samples.thin = seq((burn + 1), (num_iter), by = 20)

log_like <- ll[samples.thin]

library(tidyverse)
# RUN TrueLL_MVNIG.R BEFORE THIS PLOT
ll_plot <- data.frame(x = 1:length(log_like), y = log_like) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Traceplot of log-likelihood", x = "Iteration post burn-in", y = "") +
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
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
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) 


library(forecast)
ACF_plot <-  ggAcf(x = log_like, lag.max = 40) + ggtitle("ACF of log-likelihood") + labs(y = "") +
  # ylim(c(-0.35, 0.65)) + 
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
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
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )


ll_plot
ACF_plot

gridExtra::grid.arrange(ll_plot, ACF_plot, ncol = 2)

library(latex2exp)
alpha_plot <- data.frame(alpha = unlist(alpha.samples[samples.thin]),
                         Iteration = 1:length(samples.thin)) %>% ggplot(aes(x = Iteration, y = alpha)) + geom_line() +
  labs(
    x = "Iteration post burn-in", y = TeX("$\\alpha$")) +
  theme_classic() +
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
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
  )


gamma_plot <- data.frame(m = unlist(m.samples[samples.thin]),
                         Iteration = 1:length(samples.thin)) %>% ggplot(aes(x = Iteration, y = m)) + geom_line() +
  labs(
    x = "Iteration post burn-in", y = TeX("$\\gamma$")) +
  theme_classic() +
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
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
  )

gridExtra::grid.arrange(alpha_plot, gamma_plot, ncol = 2)
################################################################################
# CLUSTERING PERFORMANCE
################################################################################
x.limit.lower <- min(unlist(lapply(1:J, function(j){ min(X.global[[j]][1, ])})))
x.limit.upper <- max(unlist(lapply(1:J, function(j){ max(X.global[[j]][1, ])})))

y.limit.lower <- min(unlist(lapply(1:J, function(j){ min(X.global[[j]][2, ])})))
y.limit.upper <- max(unlist(lapply(1:J, function(j){ max(X.global[[j]][2, ])})))

index <- list()
for(iter in 2:num_iter){
  index[[iter]] = as.numeric(sapply(1:J, function(j){ c(k.samples[[j]][[iter]][t.samples[[j]][[iter]]])}))
}

samples <- samples.thin
posterior_samples <- matrix(0, nrow = length(samples), ncol = length(index[samples][[1]]))

for(i in 1:length(samples)){
  posterior_samples[i, ] = index[samples][[i]]
}

library(mcclust.ext)
sim_mat = comp.psm(posterior_samples)
par(mfrow = c(1, 1))
plotpsm(sim_mat)
################################################################################
# GLOBAL LEVEL CLUSTERING
################################################################################
################################################################################
# CLUSTERING BY MINIMIZING VI
################################################################################
best_mcclust2 = minVI(sim_mat, posterior_samples, method = "all")

summary(best_mcclust2)
singleton_mcclust2 = as.numeric(names(which(table(best_mcclust2$cl[1, ]) <= 0.0 * length(best_mcclust2$cl[1, ]))))

singleton.index_mmclust2 = which(best_mcclust2$cl[1, ] %in% singleton_mcclust2)

z.estimated_mcclust2 = best_mcclust2$cl[1, ]

cluster.global_mcclust <- data.frame(x = unlist(lapply(seq_len(J), function(j){c(X.global[[j]][1, ])})),
                                     y = unlist(lapply(seq_len(J), function(j){c(X.global[[j]][2, ])})),
                                     cluster = factor(z.estimated_mcclust2),
                                     cancer = factor(unlist(lapply(seq_len(J), function(j){ c(rep(paste0("Population ", j), n[j]))}))))

if(length(singleton.index_mmclust2) > 0){
  z.estimated_mcclust2 = z.estimated_mcclust2[-singleton.index_mmclust2]
  cluster.global_mcclust = cluster.global_mcclust[-singleton.index_mmclust2, ]
}else{
  z.estimated_mcclust2 = z.estimated_mcclust2
  cluster.global_mcclust = cluster.global_mcclust
}


library(pals)
L.max = 30
myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (L.max - 22))))


library(patchwork)
library(tidyverse)
library(latex2exp)
ARI_pop = 0
for(j in 1:J){
  ARI_pop[j] = aricode::ARI(k.true[[j]][t.true[[j]]],
                            cluster.global_mcclust %>% filter(cancer == paste0("Population ",j)) %>% select(cluster) %>% pull())
}



cluster.global_mcclust2 <- data.frame(x = unlist(lapply(seq_len(J), function(j){c(X.global[[j]][1, ])})),
                                      y = unlist(lapply(seq_len(J), function(j){c(X.global[[j]][2, ])})),
                                      cluster = factor(z.estimated_mcclust2),
                                      cancer = factor(unlist(lapply(seq_len(J), function(j){ c(rep(paste0("Population ", j, " \nARI = ", round(ARI_pop[j], 3)), n[j]))}))))


plot.global.VI <- cluster.global_mcclust2 %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 3) + labs(x = "X1", y = "X2") + facet_wrap(~cancer, ncol = J) + 
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
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
    legend.title=element_text(size=18),
    legend.text=element_text(size=16)
  ) + labs(title = "Global level Clustering")


plot.global.VI

