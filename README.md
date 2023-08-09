# reserving-ensemble
## Introduction

The code is written for the paper "Stochastic Ensemble Loss Reserving", with the following authors: Benjamin Avanzi, Yanfeng Li, Bernard Wong, Alan Xian. [paper on arxiv](https://arxiv.org/abs/2206.08541)

The scripts for reserving ensembles can produce the standard linear pool ensemble and the Accident Development periods dependent on Linear Pools (ADLP) ensembles tailored to the general insurance loss data.  

Please run through each of the chunks in the file 'reserving-ensemble.Rmd', which also include the plotting of the key results in the paper.

## Code Overview

The code contains three main directories used for analysis and a `plotting` directory to store the plotted outputs from the .Rmd file. This section provides an overview of each file and also 
a general framework for how the analysis is performed.

### Code Structure

* reserving-ensemble  
    * reserving-ensemble.Rmd  
    * utils  
        * `defining_scoring_functions.R`
        * `defining_helper_functions.R`
        * `functions_for_trisize.R`
        * `load_required_packages.R`
    * ensemble  
        * `defining_component_functions.R`
        * `fitting_component_models.R`
        * `fitting_ensemble_models.R`
        * `defining_partitions_triangle_40.R`
        * `defining_partitions_triangle_20.R`
    * simulation
        * triangle_40-data
            * *Contains 100 simulations of size 40 triangles and empirical simulations of total reserves from fitted ADLP ensembles*
            * `simX-full-data.csv`
            * `simX-past-data.csv`
            * `empirical_simulations_ensemble-40.csv`
        * triangle_20-data
            * *Contains 100 simulations of size 20 triangles and empirical simulations of total reserves from fitted ADLP ensembles*
            * `simX-full-data.csv`
            * `simX-past-data.csv`
            * `empirical_simulations_ensemble-20.csv`
        * `simulate_claims.R`
        * `simulate_ensembles.R`
    * plotting


The following table describes the purpose and type of functions each `.R` file 
contains.

| R file | Purpose |
|-|-|
| `load_required_packages.R`            | Defines packages to be loaded in, and installs such packages if not available on user device. |
| `defining_helper_functions.R`         | Defines miscellaneous helper functions that are commonly used across all files. |
| `functions_for_trisize.R`             | Allocate specific variables/functions to general variables/function names, whose behaviour would vary depending on triangle size |
| `defining_partitions_triangle_20.R`   | Defines the ADLP partitions and train/validation/test set partitions to be used for size 20 triangles  <br /><br /> This creates the `simX-full-data.csv` and `simX-past-data.csv` files |
| `defining_partitions_triangle_40.R`   | Defines the ADLP partitions and train/validation/test set partitions to be used for size 40 triangles  <br /><br /> This creates the `simX-full-data.csv` and `simX-past-data.csv` files |
| `simulate_claims.R`                   | Simulates claims triangles from the SynthETIC package |
| `defining_component_functions.R`      | Defines the density, distribution and prediction functions for varying component models to be used in ensembles |
| `fitting_component_models.R`          | Fits all the component models against partitioned train/validation/test datasets |
| `fitting_ensemble_models.R`           | Uses the density and train/validation/test datasets to fit the component model weights for relevant ADLP ensemble models with different partitions |
| `defining_scoring_functions.R`        | Defines functions to calculate the Log Score and Continuously Ranked Probability Score (CRPS) for ADLP ensembles models |
| `simulate_ensembles.R`                | Uses fitted ensemble models and underlying component models to simulate aggregate reserves for each ensemble model <br /><br /> This creates the `empirical_simulations_ensemble-X.csv` files |

### Analysis Framework

As this code is provided to complement the aforementioned "Stochastic Ensemble Loss Reserving" paper, the code follows a similar framework in analysis shown in the paper. It is encouraged that the reader uses the paper to complement their understanding of the code and vice versa.

#### 1) Function definitions and Data Preparation 

To prepare the data and functions for analysis, the following R files with helper functions are sourced.   

* `load_required_packages.R`  
* `defining_helper_functions.R`  
* `functions_for_trisize.R`  
* `defining_partitions_triangle_20.R`  
* `defining_partitions_triangle_40.R`  

If non-existent, past and full triangles are simulated from the `SynthETIC` package by using `simulate_claims.R`. 
Past triangles denote the dataset containing upper left claims triangle only, whilst full triangles denote the full 
dataset where all the claims in each accident year have been fully developed. Note that these triangles will mostly contain the same values with the exception of information along the diagonal (where 'future' settlement data would not have been received).

#### 2) Fitting of Component Models

Whilst the paper uses simulated claims, the following step can also be performed using real claims data.
The following R files are used:  

* `defining_component_functions.R`  
* `fitting_component_models.R`   

The claims data would be partitioned into a train/validation/test datasets and then both
in-sample and out-of-sample models are fitted for each of the component model. The code 
calculates the relevant densities of each datapoint for each model, which is to be used in
the next section.

#### 3) Fitting of Ensemble Models

Similar to the previous step, this step can be performed using real claims data instead of the
simulated ones. Note that the datasets (and train/validation/test split) should be the same
in Step 2) and Step 3).

The following R file is used:

* `fitting_ensemble_models.R`

This section uses the calculated densities of the component models to calculate the
required model weights for each component model for each ADLP partition. The estimation 
of the 'model weights for each ADLP parition' is the fitting process of the ensemble model.

#### 4) Analysis of Ensemble Models

After fitting ensemble model, we compare the predictive performance of the ensembles using 
two methods. The first is through statistical tests and the second is comparing the predictive
performance on aggregate reserve level. 

The following R files are used:

* `defining_scoring_functions.R`
* `simulate_ensembles.R` 

Further statistical tests and graphical results are then performed, which is shown in the `.Rmd` file.
