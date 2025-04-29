# 3BNC117_mPBPK_Model
R code and data used for simulation and parameter estimation of a minimal physicologically based pharmacokinetic (mPBPK) model for 3BNC117 PK in  
"Factors Influencing Monoclonal Antibody Pharmacokinetics Across Varying Immune Perturbations" by DeBonis & Davis et al., in preparation.

## Contents
Below are descriptions of the code and files found in this document.  
Note, nlmixr2 (https://nlmixr2.org/), rxode2 (https://nlmixr2.github.io/rxode2/), and xpose (https://uupharmacometrics.github.io/xpose/)  
packages are required for running these files.

### mPBPK_model.R
R file containing mPBPK model in nlmixr2 format for parameter estimation

### mPBPK_model_ODE.R
R file containing mPBPK model in rxode2 format for visualization

### estimate_mPBPK_model_parameters.R
R file for estimating parameters using nonlinear mixed effects modeling with the nlmixr2 package.

### plot_estimation_results.R
R file for visualization of parameter estiamtion results including concentration v. time profiles, EBEs, and model assessment plots using the xpose package

### population_PK_data.xlsx
Excel file containing population PK dataset to be used in parameter estimation.
