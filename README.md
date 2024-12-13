
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Global-Local (GLocal) Dirichlet Process: Variational Bayes Posterior Inference

<!-- badges: start -->
<!-- badges: end -->

<span style="font-size:24px;">Code repository for “**Global-Local
Dirichlet Processes for Clustering Grouped Data in the Presence of
Group-Specific Idiosyncratic Variables**”</span>.

**Description of R scripts**:

- The R script ***GenerateData.R*** is used to generate data for
  simulations.

- The R script ***GLocalVI.R*** contains the main function that
  implements the Variational Inference (VI) algorithm for GLocal
  Dirichlet Process. This R script is not to be run by the user.

- The R script ***VI_GLocal_MultiRunParallel.R*** performs multiple
  parallel runs of the proposed VI algorithm. This script is run after
  generating the data using the R script ***GenerateData.R***.

- The R script ***GLocalMCMC.R*** performs MCMC-based posterior
  inference for the GLocal DP with the proposed Metropolis Within Gibbs
  sampler. This script is run after generating the data using the R
  script ***GenerateData.R***.

- The R script ***HDPMCMC.R*** performs MCMC-based posterior inference
  for the hierarchical Dirichlet Process (HDP). HDP uses only the global
  variables whereas GLocal DP used both global and local variables. This
  script is run after generating the data using the R script
  ***GenerateData.R*** and takes into account only the global variables.

Description of **Folders**:

- The folder **Simulation Plots** contains all simulation plots used in
  the paper.
- The folder **Real Data Analysis** contains the data, code, and plots
  from the Real Data Analysis Section of the paper.
  - The sub-folder **Data** contains the gene expression, clinical, and
    survival data from the TCGA database corresponding to the four GI
    cancers.
  - The sub-folder **Plots** contains the plots used in the Real Data
    Analysis Section of the paper.
  - The R script ***GLocalVI_RealDataAnalysis.R*** is used in the
    analysis of the real data. This script contains all pre-processing
    steps before running the proposed VI algorithm, followed by steps to
    generate the plots as reported in the paper and stored in the
    **Plots** sub-folder.
