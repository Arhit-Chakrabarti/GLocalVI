rm(list = ls())
# Generate the data
if (!require(extraDistr)) install.packages("extraDistr", dependencies = TRUE); suppressPackageStartupMessages(library(extraDistr))
if (!require(MASS)) install.packages("MASS", dependencies = TRUE); suppressPackageStartupMessages(library(MASS))

J = 3
L.local.true = c(3, 5, 4)
L.global.true = 6 # true number of global groups in population
alpha0 = 25
m0 = 25

set.seed(2024)
alpha.true = rgamma(n = 1, shape = alpha0, rate = 1) # True alpha1 for generating data
m.true = rgamma(n = 1, shape = m0, rate = 1) # True alpha2 for generating data

# True weights
beta.true = as.numeric(rdirichlet(n = 1, alpha = rep(m.true/L.global.true, L.global.true))) # True beta
pi.true = lapply(1:J, function(j) {as.numeric(rdirichlet(n = 1, alpha = rep(alpha.true/L.local.true[j], L.local.true[j])))})

# Sample sizes
n = rep(100, J)

t.true = lapply(1:J, function(j){ sample(1:L.local.true[j], size = n[j], prob = pi.true[[j]], replace = TRUE)})

k.true = lapply(1:J, function(j){ sample(1:L.global.true, size = L.global.true, prob = beta.true, replace = TRUE)})
# Draw data from Normal populations
p.global = 3 # Dimension of global variables
p.local = c(2, 3, 4) # Dimension of local variables

# True parameters of the hyperprior for the Global variables
nu0.global.true = 5 * p.global;  Psi0.global.true = diag(1, p.global)

# True parameters of the hyperprior for the Local variables
nu0.local.true = sapply(1:J, function(j) 5 * p.local[j])
Psi0.local.true = sapply(1:J, function(j) diag(1, p.local[j]))

SigmaG0 = simplify2array(lapply(1:L.global.true, function(j){LaplacesDemon::rinvwishart(nu = nu0.global.true, S = Psi0.global.true)}))

SigmaL0 = replicate(J, list())
for(i in 1:J){
  SigmaL0[[i]] = simplify2array(lapply(1:L.local.true[i], function(j){LaplacesDemon::rinvwishart(nu = nu0.local.true[i], S = Psi0.local.true[[i]])}), except = FALSE)
}

phi0 = rep(0, p.global)
lambda.local = 0.05
lambda.local.inv = 1/lambda.local

lambda.global = 0.1
lambda.global.inv = 1/lambda.global

mu = lapply(1:J, function(j) rep(0, p.local[j]))

mean.local.true = replicate(J, list()) 
for(i in 1:J){
  mean.local.true[[i]] = t(matrix(sapply(1:L.local.true[i], function(j){mvrnorm(n = 1, mu = mu[[i]], Sigma = lambda.local.inv * SigmaL0[[i]][, , j])}), nrow = p.local[i]))
}

mean.global.true = t(sapply(1:L.global.true, function(j){mvrnorm(n = 1, mu = phi0, Sigma = lambda.global.inv * SigmaG0[,,j])}))# True mean global

X = replicate(J, list())
for(i in 1:J){
  X[[i]] = sapply(1:n[i], function(j){
    c(MASS::mvrnorm(n = 1, mu = mean.local.true[[i]][t.true[[i]][j], ], Sigma = SigmaL0[[i]][ , ,t.true[[i]][j]]),
      MASS::mvrnorm(n = 1, mu = mean.global.true[k.true[[i]][t.true[[i]]][j], ], Sigma = SigmaG0[, ,k.true[[i]][t.true[[i]]][j]]))
  })
}

X.local = list()
X.global = list()

for(j in 1:J){
  X.local[[j]] = t(X[[j]][c(seq_len(p.local[j])), , drop = FALSE])
  X.global[[j]] = t(X[[j]][-c(seq_len(p.local[j])), , drop = FALSE])
}

DATA = list()
for(j in 1:J){
  DATA[[j]] = data.frame(X.global[[j]], 
                         Cluster.true = factor(k.true[[j]][t.true[[j]]]),
                         Population.org = paste0("Population ", j))
}


DATA.global = NULL
for(j in 1:J){
  DATA.global = rbind(DATA.global, DATA[[j]])
}

if (!require(pals)) install.packages("pals", dependencies = TRUE); suppressPackageStartupMessages(library(pals))

L.max = 30
myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (L.max - 22))))

names(myvalues) = 1:L.max

if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE); suppressPackageStartupMessages(library(tidyverse))

plot1 <- DATA.global %>% ggplot(aes(x = X1, y = X2, col = Cluster.true)) + geom_point(size = 3) + labs(title = paste0("True global-level clusters")) + scale_color_manual(values = myvalues) + facet_grid(~Population.org) + 
  theme_minimal() +  
  theme(
    # LABLES APPEARANCE
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

plot2 <- DATA.global %>% ggplot(aes(x = X2, y = X3, col = Cluster.true)) + geom_point(size = 3) + labs(title = paste0("True global-level clusters")) + scale_color_manual(values = myvalues) + facet_grid(~Population.org) + 
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

gridExtra::grid.arrange(plot1, plot2)



