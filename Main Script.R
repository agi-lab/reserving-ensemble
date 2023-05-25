##### Notes ####################################################################################################
##Please put all the R scripts within the same folder in order to let the main script run successfully#######
################################################################################################################

### Loading the required packages: -----------------------------------------------------

source("Loading Required Packages.R")

### Import the relevant functions: -----------------------------------------------------
#### Import the functions to build the ensembles
source("Defining Relevant Functions.R")

### Generate simulating triangles and fit individual models -----------------------------------------------------

##### Specify the number of simulations that the user wants to run:
ntri <- 10
### Recommend choosing a number greater than or equal to 10; the number of simulations used in this paper is 100

##### Run the script:
source("Fitting Component Models.R")

### Build ensembles --------------------------------------------------------------------------------------------

#### Fit the standard linear pools:

source("Building Standard Linear Pool.R")

#### Fit the ADLP ensembles:

source("Building ADLP Ensembles.R")


### Reserve calculation -----------------------------------------------------

### Calculate the central estimation of reserve
source("Calculate Central Reserve Bias.R")

### Simulate the reserve quantiles for Equally weighted ensemble, BMV and ADLP:
source("Simulating Reserve Quantiles.R")

### Simulate the reserve quantiles for SLP:
source("Simulating Quantiles Par0.R")

### CRPS calculations -----------------------------------------------------

#### Import the CRPS calculation functions
source("CRPS_Calculation_Functions.R")

#### Perform the calculations:
source("CRPS_calculation_script.R")

### Statistical tests -----------------------------------------------------
source("Statistical Tests.R")

### Plotting of key results -----------------------------------------------------
### (This script is optional to run)
source("Plottings.R")

### Getting additional results in the appendix (Optional to run) -----------------------------------------------------

#### Testing ADLP ensembles with three subsets
source("Appendix_Building Three Subsets ADLP ensembles.R")

#### Testing GAMLSS and smoothing spline
source("Appendix_GAMLSS with Smoothing Spline.R")
