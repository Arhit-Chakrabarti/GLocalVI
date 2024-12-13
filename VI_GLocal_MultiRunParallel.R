source("GLocalVI.R") # Load the GLocal DP VI function
# Load the necessary library for parallel execution
if (!require(parallel)) install.packages("parallel", dependencies = TRUE); suppressPackageStartupMessages(library(parallel))
# Define the number of cores to use
numberOfCores <- detectCores() - 1  # Use one less than the total cores
# Run the function in parallel
numReplicate = 20

out <- mclapply(1:numReplicate, function(i) {
  GLocalVI(X.local = X.local, X.global = X.global,
           H = 30,
           L = 30,
           maxIter = 500,
           epsilon = 1e-5)}
  , mc.cores = numberOfCores)


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


if (!require(pals)) install.packages("pals", dependencies = TRUE); suppressPackageStartupMessages(library(pals))

myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (8))))

names(myvalues) = 1:30

if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE); suppressPackageStartupMessages(library(tidyverse))

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

if (!require(aricode)) install.packages("aricode", dependencies = TRUE); suppressPackageStartupMessages(library(aricode))

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


if (!require(pals)) install.packages("pals", dependencies = TRUE); suppressPackageStartupMessages(library(pals))

myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (8))))

names(myvalues) = 1:30

if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE); suppressPackageStartupMessages(library(tidyverse))

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

if (!require(gridExtra)) install.packages("gridExtra", dependencies = TRUE); suppressPackageStartupMessages(library(gridExtra))

plot_Global <- gridExtra::grid.arrange(True_plot, Est_plot)
plot_Local <- gridExtra::grid.arrange(True_plot_Local, Est_plot_Local)

if (!require(grid)) install.packages("grid", dependencies = TRUE); suppressPackageStartupMessages(library(grid))
if (!require(gridExtra)) install.packages("gridExtra", dependencies = TRUE); suppressPackageStartupMessages(library(gridExtra))

gridExtra::grid.arrange(plot_Local, plot_Global, ncol = 2,
                        top = textGrob(paste0("Local variable dimensions = ", p.local[1],", ", p.local[2], ", ", p.local[3], ". Global variable dimension = ", p.global), ,gp=gpar(fontsize=20,font=3)))

