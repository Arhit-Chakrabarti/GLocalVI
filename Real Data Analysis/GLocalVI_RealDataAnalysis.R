if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE); suppressPackageStartupMessages(library(tidyverse))
if (!require(data.table)) install.packages("data.table", dependencies = TRUE); suppressPackageStartupMessages(library(data.table))
################################################################################
# READING THE PHENOTYPE DATA
################################################################################
colon_cov = fread("./Real Data Analysis/Data/Colon.gz")
rectal_cov = fread("./Real Data Analysis/Data/Rectal.gz")
stomach_cov = fread("./Real Data Analysis/Data/Stomach.gz")
eoso_cov = fread("./Real Data Analysis/Data/Eoso.gz")

################################################################################
# READING THE GENE-EXPRESSION DATA
################################################################################
colon = fread("./Real Data Analysis/Data/Colon_counts.tsv.gz")
rectal = fread("./Real Data Analysis/Data/Rectal_counts.tsv.gz")
stomach = fread("./Real Data Analysis/Data/Stomach_counts.tsv.gz")
eoso = fread("./Real Data Analysis/Data/Eoso_counts.tsv.gz")

common.genes <- as.vector(intersect(intersect(intersect(colon[c(1:60483), 1], 
                                                        rectal[c(1:60483), 1]),
                                              stomach[c(1:60483), 1]),
                                    eoso[c(1:60483), 1]))

attributes(common.genes) <- NULL
common.genes = unlist(common.genes)

################################################################################
# KEEP ONLY THE GENES THAT ARE COMMON ACROSS THE FOUR CANCERS
################################################################################
colon_data <- t(colon %>% filter(Ensembl_ID %in% common.genes) %>%
                  select(-Ensembl_ID))

rectal_data <- t(rectal %>% filter(Ensembl_ID %in% common.genes) %>%
                   select(-Ensembl_ID))

stomach_data <- t(stomach %>% filter(Ensembl_ID %in% common.genes) %>%
                    select(-Ensembl_ID))
eoso_data <- t(eoso %>% filter(Ensembl_ID %in% common.genes) %>%
                 select(-Ensembl_ID))

################################################################################
# MERGE THE GENE EXPRESSION DATA FOR THE FOUR CANCERS
################################################################################
all_data_merged <- rbind(colon_data,
                         rectal_data,
                         stomach_data,
                         eoso_data)


################################################################################
# PERFORM PCA ON THE COMBINED GENE EXPRESSION DATA
################################################################################
if (!require(irlba)) install.packages("irlba", dependencies = TRUE); suppressPackageStartupMessages(library(irlba))

pca.all = prcomp_irlba(all_data_merged, n = 10) # Top 10 PCs
summary(pca.all) # Look at the summary of PCA

# The Rotated Data form the Global variables
X = pca.all$x
# Rename the rows of X. This is needed for data sub-setting
rownames(X) = rownames(all_data_merged) 

################################################################################
# SEPERATE THE PCA DATA BY THE FOUR CANCERS
################################################################################
pca.colon = X[which(rownames(X) %in% colon_cov$submitter_id.samples), ]
pca.rectal = X[which(rownames(X) %in% rectal_cov$submitter_id.samples), ]
pca.stomach = X[which(rownames(X) %in% stomach_cov$submitter_id.samples), ]
pca.eoso = X[which(rownames(X) %in% eoso_cov$submitter_id.samples), ]


################################################################################
# PLOT THE PC DATA
################################################################################

x.limit.lower <- min(X[, 1])
x.limit.upper <- max(X[, 1])
y.limit.lower <- min(X[, 2])
y.limit.upper <- max(X[, 2])

plot1 = data.frame(pca.colon) %>% ggplot(aes(x = PC1, y = PC2)) + geom_point() + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + ggtitle("Colon cancer") + 
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
  )

plot2 = data.frame(pca.rectal) %>% ggplot(aes(x = PC1, y = PC2)) + geom_point() + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper)  + ggtitle("Rectal cancer") + 
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
  )

plot3 = data.frame(pca.stomach) %>% ggplot(aes(x = PC1, y = PC2)) + geom_point() + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper)  + ggtitle("Stomach cancer") + 
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
  )

plot4 = data.frame(pca.eoso) %>% ggplot(aes(x = PC1, y = PC2)) + geom_point() + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper)  + ggtitle("Esophageal cancer") + 
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
  )

if (!require(gridExtra)) install.packages("gridExtra", dependencies = TRUE); suppressPackageStartupMessages(library(gridExtra))
gridExtra::grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)

################################################################################
# FINAL DATA GENERATION BEFORE RUNNING VI 
################################################################################
# REMOVE THE OBSERVATIONS WITH MISSING CLINICAL VARIABLES
################################################################################
stomach_PCA <- data.frame(submitter_id.samples = rownames(pca.stomach), pca.stomach) %>% na.omit()
rownames(stomach_PCA) <- NULL # REMOVE THE ROWNAMES (NOT NECESSARY)

colon_PCA <- data.frame(submitter_id.samples = rownames(pca.colon), pca.colon) %>% na.omit() 
rownames(colon_PCA) <- NULL # REMOVE THE ROWNAMES (NOT NECESSARY)

rectal_PCA <- data.frame(submitter_id.samples = rownames(pca.rectal), pca.rectal) %>% na.omit() 
rownames(rectal_PCA) <- NULL # REMOVE THE ROWNAMES (NOT NECESSARY)

eoso_PCA <- data.frame(submitter_id.samples = rownames(pca.eoso), pca.eoso) %>% na.omit() 
rownames(eoso_PCA) <- NULL # REMOVE THE ROWNAMES (NOT NECESSARY)

################################################################################
# MERGE THE CLINICAL VARIABLES (LOCAL VARIABLES) WITH THE PCs (GLOBAL VARIABLES) 
################################################################################
colon_merged <- left_join(colon_PCA, colon_cov, by = "submitter_id.samples")
rectal_merged <- left_join(rectal_PCA, rectal_cov, by = "submitter_id.samples")
stomach_merged <- left_join(stomach_PCA, stomach_cov, by = "submitter_id.samples")
eoso_merged <- left_join(eoso_PCA, eoso_cov, by = "submitter_id.samples")

################################################################################
# LOCAL VARIABLE FOR STOMACH CANCER:  NUMBER OF POSITIVE LYMPH NODES
# GLOBAL VARIABLES SHARED ACROSS CANCERS: TOP 10 PCs
################################################################################
stomach_full <- stomach_merged %>% dplyr::select(number_of_lymphnodes_positive_by_he,
                                                 PC1:PC10) %>% na.omit() 

dim(stomach_full)

################################################################################
# LOCAL VARIABLE FOR RECTAL CANCER:  CEA LEVEL
# GLOBAL VARIABLES SHARED ACROSS CANCERS: TOP 10 PCs
################################################################################
rectal_full <- rectal_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                               PC1:PC10) %>% na.omit() 

dim(rectal_full)

################################################################################
# LOCAL VARIABLE FOR COLON CANCER:  CEA LEVEL AND BMI
# GLOBAL VARIABLES SHARED ACROSS CANCERS: TOP 10 PCs
################################################################################
colon_full <- colon_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                             bmi.exposures,
                                             PC1:PC10) %>% na.omit() 
dim(colon_full)

################################################################################
# LOCAL VARIABLE FOR ESOPHAGEAL CANCER:  NUMBER OF CIGARETTES PER DAY
# GLOBAL VARIABLES SHARED ACROSS CANCERS: TOP 10 PCs
################################################################################
eoso_full <- eoso_merged %>% dplyr::select(cigarettes_per_day.exposures,
                                           PC1:PC10) %>% na.omit() 

dim(eoso_full)

X1 <- t(as.matrix(stomach_full))
X2 <- t(as.matrix(rectal_full))
X3 <- t(as.matrix(colon_full))
X4 <- t(as.matrix(eoso_full))

X <- list(X1, X2, X3, X4)

################################################################################
# STORE THE GLOBAL VARIABLES IN A LIST
################################################################################
X.global = list(t(X1[which(colnames(stomach_full) %in% paste0("PC",1:10)), ,drop = F]),
                t(X2[which(colnames(rectal_full) %in% paste0("PC",1:10)), ,drop = F]), 
                t(X3[which(colnames(colon_full) %in% paste0("PC",1:10)), ,drop = F]), 
                t(X4[which(colnames(eoso_full) %in% paste0("PC",1:10)), ,drop = F])) 

################################################################################
# STORE THE LOCAL VARIABLES IN A LIST
################################################################################
X.local = list(t(X1[-which(colnames(stomach_full) %in% paste0("PC",1:10)), ,drop = F]),
               t(X2[-which(colnames(rectal_full) %in% paste0("PC",1:10)), ,drop = F]), 
               t(X3[-which(colnames(colon_full) %in% paste0("PC",1:10)), ,drop = F]), 
               t(X4[-which(colnames(eoso_full) %in% paste0("PC",1:10)), ,drop = F])) 



################################################################################
## CENTER THE LOCAL COVARIATES
################################################################################
X.local = lapply(X.local, FUN = scale, scale = FALSE)


################################################################################
# RUN VI IN PARALLEL FOR 20 INDEPENDENT RUNS
################################################################################
source("GLocalVI.R") # Load the GLocal DP VI function
# Load the necessary library for parallel execution
if (!require(parallel)) install.packages("parallel", dependencies = TRUE); suppressPackageStartupMessages(library(parallel))
# Define the number of cores to use
numberOfCores <- detectCores() - 1  # Use one less than the total cores
# Run the function in parallel
numReplicate = 20 # NUMBER OF PARALLEL RUNS OF VI

out <- mclapply(1:numReplicate, function(i) {
  GLocalVI(X.local = X.local, X.global = X.global,
           H = 30,
           L = 30,
           maxIter = 2000,
           epsilon = 1e-5)}
  , mc.cores = numberOfCores)

################################################################################
# TO EXTRACT THE RUN OF VI WITH HIGHEST ELBO AT CONVERGENCE
################################################################################
ELBO = 0
Time = 0

for(i in 1:length(out)){
  ELBO[i] = out[[i]]$ELBO_val_conv
  Time[i] = out[[i]]$Time.taken
}
 # EXTRACT THE OUTPUT WITH HIGHEST ELBO TO DRAW INFERENCE
best_outputs = out[[which.max(ELBO)]]

J = length(X.global)
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
                             Cluster.est = factor(est.k_jt[[j]][est.t_ji[[j]]]),
                             Population = paste0("Population ", j))
}


DATA.global_Est = NULL
for(j in 1:J){
  DATA.global_Est = rbind(DATA.global_Est, DATA_Est[[j]])
}

DATA.global_Est <- DATA.global_Est %>%
  mutate(Population = case_when(Population == "Population 1" ~ 'Stomach Cancer',
                                Population == "Population 2" ~ 'Rectal Cancer',
                                Population == "Population 3" ~ 'Colon Cancer',
                                Population == "Population 4" ~ 'Esophageal Cancer'))

if (!require(pals)) install.packages("pals", dependencies = TRUE); suppressPackageStartupMessages(library(pals))
myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (8))))

names(myvalues) = 1:30

if (!require(latex2exp)) install.packages("latex2exp", dependencies = TRUE); suppressPackageStartupMessages(library(latex2exp))

Est_plot.global <- DATA.global_Est %>% ggplot(aes(x = PC1, y = PC2, col = Cluster.est)) + geom_point(size = 3) +  
  scale_color_manual(values = myvalues) + facet_grid(~Population) + 
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
  ) + guides(color = guide_legend(title = "Clusters")) + labs(x = TeX("$PC_1$", bold = TRUE), y = TeX("$PC_2$", bold = TRUE))

Est_plot.global


DATA_Est_Local = list()
for(j in 1:J){
  DATA_Est_Local[[j]] = data.frame(X.global[[j]],
                                   Cluster.est = factor(est.t_ji[[j]]),
                                   Population = paste0("Population ", j))
}


DATA.local_Est = NULL
for(j in 1:J){
  DATA.local_Est = rbind(DATA.local_Est, DATA_Est_Local[[j]])
}

DATA.local_Est <- DATA.local_Est %>%
  mutate(Population = case_when(Population == "Population 1" ~ 'Stomach Cancer',
                                Population == "Population 2" ~ 'Rectal Cancer',
                                Population == "Population 3" ~ 'Colon Cancer',
                                Population == "Population 4" ~ 'Esophageal Cancer'))

DATA.local_Est$Cluster.est.global <- DATA.global_Est$Cluster.est


characters <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
######################################################################## 
### RELABELLING LOCAL LEVEL CLUSTERS
######################################################################## 
### STOMACH CANCER
#####################################################
Local.cluster <- DATA.local_Est %>% filter(Population == "Stomach Cancer") %>% select(Cluster.est) %>% pull() %>% as.character()

Global.cluster <- DATA.local_Est %>% filter(Population == "Stomach Cancer") %>% select(Cluster.est.global) %>% pull() %>% as.character()

Local.cluster_unique_sorted <- sort(as.numeric(unique(Local.cluster)))

New.Local.cluster = 0
for(l in 1:length(Local.cluster)){
  for(m in 1:length(Local.cluster_unique_sorted)){
    if(Local.cluster[l] == as.character(Local.cluster_unique_sorted[m]))
      New.Local.cluster[l] = paste0(Global.cluster[l], characters[Local.cluster_unique_sorted[m]])
  }
  
}

Local.cluster.stomach_relabelled <- factor(New.Local.cluster)
#####################################################
### RECTAL CANCER
#####################################################
Local.cluster <- DATA.local_Est %>% filter(Population == "Rectal Cancer") %>% select(Cluster.est) %>% pull() %>% as.character()

Global.cluster <- DATA.local_Est %>% filter(Population == "Rectal Cancer") %>% select(Cluster.est.global) %>% pull() %>% as.character()

Local.cluster_unique_sorted <- sort(as.numeric(unique(Local.cluster)))

New.Local.cluster = 0
for(l in 1:length(Local.cluster)){
  for(m in 1:length(Local.cluster_unique_sorted)){
    if(Local.cluster[l] == as.character(Local.cluster_unique_sorted[m]))
      New.Local.cluster[l] = paste0(Global.cluster[l], characters[Local.cluster_unique_sorted[m]])
  }
  
}

Local.cluster.rectal_relabelled <- factor(New.Local.cluster)

#####################################################
### COLON CANCER
#####################################################
Local.cluster <- DATA.local_Est %>% filter(Population == "Colon Cancer") %>% select(Cluster.est) %>% pull() %>% as.character()

Global.cluster <- DATA.local_Est %>% filter(Population == "Colon Cancer") %>% select(Cluster.est.global) %>% pull() %>% as.character()

Local.cluster_unique_sorted <- sort(as.numeric(unique(Local.cluster)))

New.Local.cluster = 0
for(l in 1:length(Local.cluster)){
  for(m in 1:length(Local.cluster_unique_sorted)){
    if(Local.cluster[l] == as.character(Local.cluster_unique_sorted[m]))
      New.Local.cluster[l] = paste0(Global.cluster[l], characters[Local.cluster_unique_sorted[m]])
  }
  
}

Local.cluster.colon_relabelled <- factor(New.Local.cluster)

#####################################################
### ESOPHAGEAL CANCER
#####################################################
Local.cluster <- DATA.local_Est %>% filter(Population == "Esophageal Cancer") %>% select(Cluster.est) %>% pull() %>% as.character()

Global.cluster <- DATA.local_Est %>% filter(Population == "Esophageal Cancer") %>% select(Cluster.est.global) %>% pull() %>% as.character()

Local.cluster_unique_sorted <- sort(as.numeric(unique(Local.cluster)))

New.Local.cluster = 0
for(l in 1:length(Local.cluster)){
  for(m in 1:length(Local.cluster_unique_sorted)){
    if(Local.cluster[l] == as.character(Local.cluster_unique_sorted[m]))
      New.Local.cluster[l] = paste0(Global.cluster[l], characters[Local.cluster_unique_sorted[m]])
  }
  
}

Local.cluster.eoso_relabelled <- factor(New.Local.cluster)

DATA.local_Est$Cluster.est.relabelled <- c(Local.cluster.stomach_relabelled,
                                           Local.cluster.rectal_relabelled,
                                           Local.cluster.colon_relabelled,
                                           Local.cluster.eoso_relabelled)

myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (26)),
                    alphabet(n = (10))))

names(myvalues) = unique(DATA.local_Est$Cluster.est.relabelled)


Est_plot.local <- DATA.local_Est %>% ggplot(aes(x = PC1, y = PC2, col = Cluster.est.relabelled)) + geom_point(size = 3) +  
  scale_color_manual(values = myvalues) + facet_grid(~Population) + 
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
  ) + guides(color = guide_legend(title = "Clusters")) + guides(color = guide_legend(title = "Clusters")) + labs(x = TeX("$PC_1$", bold = TRUE), y = TeX("$PC_2$", bold = TRUE))

Est_plot.local

if (!require(gridExtra)) install.packages("gridExtra", dependencies = TRUE); suppressPackageStartupMessages(library(gridExtra))
if (!require(grid)) install.packages("grid", dependencies = TRUE); suppressPackageStartupMessages(library(grid))

p.local = unlist(lapply(X.local, ncol))
p.global = ncol(X.global[[1]])

gridExtra::grid.arrange(Est_plot.global, Est_plot.local, nrow = 2,
                        top = textGrob(paste0("Local variables. Stomach: # Postive Lymph nodes, Rectal: Pre-operative CEA, Colon: Pre-operative CEA and BMI, Esophageal: # Cigarrets/day\nLocal variable dimensions = ", p.local[1],", ", p.local[2], ", ", p.local[3], ", ", p.local[4], ". Global variable dimension (number of PCs) = ", p.global), ,gp=gpar(fontsize=15,font=3)))


################################################################################
## PERFORM UMAP ON THE COMBINED GENE EXPRESSION DATA TO VISUALIZE THE CLUSTERS
################################################################################
if (!require(uwot)) install.packages("uwot", dependencies = TRUE); suppressPackageStartupMessages(library(uwot))

umap.all = uwot::umap(all_data_merged, n_neighbors = 30,
                      n_components = 2,
                      metric = "euclidean",
                      min_dist = 1.5,
                      seed = 2024)

################################################################################
# THIS ARE TO MATCH THE UMAP EMBEDDINGS WITH THE PCS
################################################################################
stomach_full2 <- stomach_merged %>% dplyr::select(number_of_lymphnodes_positive_by_he,
                                                  PC1:PC10,
                                                  submitter_id.samples) %>% na.omit() 

rectal_full2 <- rectal_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                                PC1:PC10,
                                                submitter_id.samples) %>% na.omit() 

colon_full2 <- colon_merged %>% dplyr::select(preoperative_pretreatment_cea_level,
                                              bmi.exposures,
                                              PC1:PC10,
                                              submitter_id.samples) %>% na.omit() 

eoso_full2 <- eoso_merged %>% dplyr::select(cigarettes_per_day.exposures,
                                            PC1:PC10,
                                            submitter_id.samples) %>% na.omit() 

################################################################################
# GET THE UMAP DATA BY CANCERS
################################################################################
umap.colon = umap.all[which(rownames(umap.all) %in% colon_full2$submitter_id.samples),]
umap.rectal = umap.all[which(rownames(umap.all) %in% rectal_full2$submitter_id.samples),]
umap.stomach = umap.all[which(rownames(umap.all) %in% stomach_full2$submitter_id.samples),]
umap.eoso = umap.all[which(rownames(umap.all) %in% eoso_full2$submitter_id.samples),]


UMAP_DATA_Est <- rbind(data.frame(umap.stomach, Cluster.est = factor(est.k_jt[[1]][est.t_ji[[1]]]), Cancer = "Stomach Cancer"),
                       data.frame(umap.rectal, Cluster.est = factor(est.k_jt[[2]][est.t_ji[[2]]]), Cancer = "Rectal Cancer"),
                       data.frame(umap.colon, Cluster.est = factor(est.k_jt[[3]][est.t_ji[[3]]]), Cancer = "Colon Cancer"),
                       data.frame(umap.eoso, Cluster.est = factor(est.k_jt[[4]][est.t_ji[[4]]]), Cancer = "Esophageal Cancer")
)

################################################################################
# GLOBAL VARIABLES BY GLOBAL-LEVEL CLUSTERS
################################################################################
myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (8))))

names(myvalues) = 1:30

Global_UMAP_Plot <- UMAP_DATA_Est %>% ggplot(aes(x = X1, y = X2, col = Cluster.est)) + geom_point(size = 3) +   
  scale_color_manual(values = myvalues) + facet_grid(~Cancer) + 
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
  ) + guides(color = guide_legend(title = "Clusters")) + labs(x = TeX("$UMAP_1$", bold = TRUE), y = TeX("$UMAP_2$", bold = TRUE))

################################################################################
# GLOBAL VARIABLES BY LOCAL-LEVEL CLUSTERS
################################################################################
myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (26)),
                    alphabet(n = (10))))


names(myvalues) = unique(DATA.local_Est$Cluster.est.relabelled)

UMAP_DATA_Est_Local <- rbind(data.frame(umap.stomach, Cluster.est = Local.cluster.stomach_relabelled, Cancer = "Stomach Cancer"),
                       data.frame(umap.rectal, Cluster.est = Local.cluster.rectal_relabelled, Cancer = "Rectal Cancer"),
                       data.frame(umap.colon, Cluster.est = Local.cluster.colon_relabelled, Cancer = "Colon Cancer"),
                       data.frame(umap.eoso, Cluster.est = Local.cluster.eoso_relabelled, Cancer = "Esophageal Cancer")
)


Local_UMAP_Plot <- UMAP_DATA_Est_Local %>% ggplot(aes(x = X1, y = X2, col = Cluster.est)) + geom_point(size = 3) +  
  scale_color_manual(values = myvalues) + facet_grid(~Cancer) + 
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
  ) + guides(color = guide_legend(title = "Clusters")) + labs(x = TeX("$UMAP_1$", bold = TRUE), y = TeX("$UMAP_2$", bold = TRUE))

gridExtra::grid.arrange(Global_UMAP_Plot, Local_UMAP_Plot, ncol = 1)


################################################################################
# LOCAL VARIABLES KERNEL DENSITY/SCATTER PLOTS BY LOCAL-LEVEL CLUSTERS
################################################################################
myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (26)),
                    alphabet(n = (10))))

names(myvalues) = unique(DATA.local_Est$Cluster.est.relabelled)

Stomach_Cancer_LocalVarPlot <- data.frame(y = X.local[[1]][, 1] + as.numeric(attr(X.local[[1]], "scaled:center")), Cluster = Local.cluster.stomach_relabelled) %>% group_by(Cluster) %>% ggplot(aes(x = y, col = Cluster)) + geom_density(size = 1) + labs(x = "Number of positive lymph nodes", y = "Density", title = "Stomach Cancer")  + scale_color_manual(values = myvalues)  + theme_classic() +  
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

Rectal_Cancer_LocalVarPlot <- data.frame(y = X.local[[2]][, 1] + as.numeric(attr(X.local[[2]], "scaled:center")), Cluster = Local.cluster.rectal_relabelled) %>% group_by(Cluster) %>% ggplot(aes(x = y, col = Cluster)) + geom_density(size = 1) + labs(x = "CEA level", y = "Density", title = "Rectal Cancer")  + scale_color_manual(values = myvalues)  + theme_classic() +  
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
  ) + xlim(c(0, 200))

Colon_Cancer_LocalVarPlot <- data.frame(x = X.local[[3]][, 1] + as.numeric(attr(X.local[[3]], "scaled:center")[1]),
                                        y = X.local[[3]][, 2] + as.numeric(attr(X.local[[3]], "scaled:center")[2]), Cluster = Local.cluster.colon_relabelled) %>% ggplot(aes(x = x, y = y, col = Cluster)) + geom_point(size = 3) + labs(x = "CEA level", y = "BMI", title = "Colon Cancer")  + scale_color_manual(values = myvalues)  + theme_classic() +  
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
  ) + xlim(c(0, 300)) + ylim(c(15, 55))

Esophageal_Cancer_LocalVarPlot <- data.frame(y = X.local[[4]][, 1] + as.numeric(attr(X.local[[4]], "scaled:center")), Cluster = Local.cluster.eoso_relabelled) %>% group_by(Cluster) %>% ggplot(aes(x = y, col = Cluster)) + geom_density(size = 1) + labs(x = "Number of cigarettes per day", y = "Density", title = "Esophageal Cancer")  + scale_color_manual(values = myvalues)  + theme_classic() +  
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

gridExtra::grid.arrange(Colon_Cancer_LocalVarPlot,
                        Esophageal_Cancer_LocalVarPlot,
                        Rectal_Cancer_LocalVarPlot,
                        Stomach_Cancer_LocalVarPlot, ncol = 2)

################################################################################
## LOCAL LEVEL KM ANALYSIS
################################################################################
################################################################################
# READING THE PHENOTYPE DATA GAIN AS THE ORIGINAL DATA WAS MODIFED PREVIOSULY
################################################################################
colon_cov = fread("./Real Data Analysis/Data/Colon.gz")
rectal_cov = fread("./Real Data Analysis/Data/Rectal.gz")
stomach_cov = fread("./Real Data Analysis/Data/Stomach.gz")
eoso_cov = fread("./Real Data Analysis/Data/Eoso.gz")
################################################################################
# READING THE SURVIVAL DATA
################################################################################
colon_surv = fread("./Real Data Analysis/Data/Colon_Survival.txt")
rectal_surv = fread("./Real Data Analysis/Data/Rectal_Survival.txt")
stomach_surv = fread("./Real Data Analysis/Data/Stomach_Survival.txt")
eoso_surv = fread("./Real Data Analysis/Data/Eoso_Survival.txt")

mytheme <- theme_minimal() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=24, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=20, face="bold", colour = "black"),    
    axis.title.y = element_text(size=20, face="bold", colour = "black"),    
    axis.text.x = element_text(size=20, face="bold", colour = "black"), 
    axis.text.y = element_text(size=20, face="bold", colour = "black"),
    strip.text.x = element_text(size = 16, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 16, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=18),
    legend.text=element_text(size=17)
  )
######################################################################## 
### ESOPHAGEAL CANCER
######################################################################## 
eoso.cluster <- DATA.local_Est %>% filter(Population == "Esophageal Cancer")
eoso_merged_new <- left_join(eoso.cluster, eoso_merged, by = paste0("PC",1:10))

if (!require(stringr)) install.packages("stringr", dependencies = TRUE); suppressPackageStartupMessages(library(stringr))

eoso_merged_new$submitter_id.samples <- str_sub(eoso_merged_new$submitter_id.samples, end = -2)
eoso_surv <- rename(eoso_surv, submitter_id.samples = sample)
eoso_surv$submitter_id.samples <- str_sub(eoso_surv$submitter_id.samples, end = -2)

eoso_merged_surv_new <- left_join(eoso_merged_new, eoso_surv, by = "submitter_id.samples")
eoso_merged_surv_new <- eoso_merged_surv_new %>% dplyr::select(OS, OS.time, Cluster.est.relabelled)

if (!require(survival)) install.packages("survival", dependencies = TRUE); suppressPackageStartupMessages(library(survival))

fit.eoso.local <- survfit(Surv(OS.time, OS) ~ Cluster.est.relabelled, data = eoso_merged_surv_new)

if (!require(survminer)) install.packages("survminer", dependencies = TRUE); suppressPackageStartupMessages(library(survminer))

myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (10))))

names(myvalues) = as.character(unique(DATA.local_Est$Cluster.est.relabelled))

os.eoso.local = ggsurvplot(fit.eoso.local, data = eoso_merged_surv_new, conf.int = F,
                           legend.title = "Cluster",
                           legend.labs = as.character(levels(Local.cluster.eoso_relabelled)),
                           palette = myvalues,
                           ggtheme = mytheme) + ggtitle("Overall Survival: Esophageal Cancer")

######################################################################## 
### STOMACH CANCER
######################################################################## 
stomach.cluster <- DATA.local_Est %>% filter(Population == "Stomach Cancer")
stomach_merged_new <- left_join(stomach.cluster, stomach_merged, by = paste0("PC",1:10))

stomach_merged_new$submitter_id.samples <- str_sub(stomach_merged_new$submitter_id.samples, end = -2)
stomach_surv <- rename(stomach_surv, submitter_id.samples = sample)
stomach_surv$submitter_id.samples <- str_sub(stomach_surv$submitter_id.samples, end = -2)

stomach_merged_surv_new <- left_join(stomach_merged_new, stomach_surv, by = "submitter_id.samples")
stomach_merged_surv_new <- stomach_merged_surv_new %>% dplyr::select(OS, OS.time, Cluster.est.relabelled)

fit.stomach.local <- survfit(Surv(OS.time, OS) ~ Cluster.est.relabelled, data = stomach_merged_surv_new)

myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (10))))

names(myvalues) = as.character(unique(DATA.local_Est$Cluster.est.relabelled))

os.stomach.local = ggsurvplot(fit.stomach.local, data = stomach_merged_surv_new, conf.int = F,
                              legend.title = "Cluster",
                              legend.labs = as.character(levels(Local.cluster.stomach_relabelled)),
                              palette = myvalues,
                              ggtheme = mytheme) + ggtitle("Overall Survival: Stomach Cancer")

######################################################################## 
### RECTAL CANCER
######################################################################## 
rectal.cluster <- DATA.local_Est %>% filter(Population == "Rectal Cancer")
rectal_merged_new <- left_join(rectal.cluster, rectal_merged, by = paste0("PC",1:10))

rectal_merged_new$submitter_id.samples <- str_sub(rectal_merged_new$submitter_id.samples, end = -2)
rectal_surv <- rename(rectal_surv, submitter_id.samples = sample)

rectal_merged_surv_new <- left_join(rectal_merged_new, rectal_surv, by = "submitter_id.samples")
rectal_merged_surv_new <- rectal_merged_surv_new %>% dplyr::select(OS, OS.time, Cluster.est.relabelled)

fit.rectal.local <- survfit(Surv(OS.time, OS) ~ Cluster.est.relabelled, data = rectal_merged_surv_new)

myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (10))))

names(myvalues) = as.character(unique(DATA.local_Est$Cluster.est.relabelled))

os.rectal.local = ggsurvplot(fit.rectal.local, data = rectal_merged_surv_new, conf.int = F,
                             legend.title = "Cluster",
                             legend.labs = as.character(levels(Local.cluster.rectal_relabelled)),
                             palette = myvalues,
                             ggtheme = mytheme) + ggtitle("Overall Survival: Rectal Cancer")

######################################################################## 
### COLON CANCER
######################################################################## 
colon.cluster <- DATA.local_Est %>% filter(Population == "Colon Cancer")
colon_merged_new <- left_join(colon.cluster, colon_merged, by = paste0("PC",1:10))

colon_merged_new$submitter_id.samples <- str_sub(colon_merged_new$submitter_id.samples, end = -2)
colon_surv <- rename(colon_surv, submitter_id.samples = sample)

colon_merged_surv_new <- left_join(colon_merged_new, colon_surv, by = "submitter_id.samples")
colon_merged_surv_new <- colon_merged_surv_new %>% dplyr::select(OS, OS.time, Cluster.est.relabelled)

fit.colon.local <- survfit(Surv(OS.time, OS) ~ Cluster.est.relabelled, data = colon_merged_surv_new)

myvalues = unname(c(kelly(n = 22),
                    alphabet2(n = (10))))

names(myvalues) = as.character(unique(DATA.local_Est$Cluster.est.relabelled))

os.colon.local = ggsurvplot(fit.colon.local, data = colon_merged_surv_new, conf.int = F,
                            legend.title = "Cluster",
                            legend.labs = as.character(levels(Local.cluster.colon_relabelled)),
                            palette = myvalues,
                            ggtheme = mytheme) + ggtitle("Overall Survival: Colon Cancer")


if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)
splots <- list()
splots[[1]] <- os.colon.local
splots[[2]] <- os.rectal.local
splots[[3]] <- os.eoso.local
splots[[4]] <- os.stomach.local


res <- arrange_ggsurvplots(splots, print = TRUE,
                           ncol = 2, nrow = 2)
