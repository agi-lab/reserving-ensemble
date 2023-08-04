################################################################################
### Function for fitting ensemble models (automatically calculates weights)
###
### Function accepts three 5 parameters:
### - component_dens_valid: densities of component models from validation set
### - component_dens_out: densities of component models from out of sample set
### - valid_partitions: list containing indices for the paritions of the validation set
### - out_partitions: list containing indicies for the partitions of the out of sample set
################################################################################

fit_ensemble_model <- function(
    component_dens_valid, 
    component_dens_out, 
    partition_func
){
  
    # number of simulations
    n.sims <- length(component_dens_valid)
    # number of components (take away origin and dev columns)
    n.components <- length(component_dens_valid[[1]]) - 2
    
    w_init_1_3<-rep(1/n.components,n.components)
    
    n.partitions <- length(partition_func(component_dens_valid[[1]]))
    
    OW_out_dens_sims <- list()
    model_weights_simul_sims <- list()
    optim_MM_sims <- list()
    
    for (sim in 1:n.sims) {
    
    OW_out_dens <-list()
    OW_out_CDF <- list()
    model_weights_simul<-list()
    optim_MM_par <- list()
    
    valid_partitions <- partition_func(component_dens_valid[[sim]])
    out_partitions <- partition_func(component_dens_out[[sim]])
    
    
    for (i in 1:n.partitions) {
        partition_in <- do.call('rbind', valid_partitions[1:i])
        parition_out <- valid_partitions[[i]]
        
        dens_valid <- partition_in[, 3:(n.components+2)]
        dens_out <- parition_out[, 3:(n.components+2)]

        #Train the model weights using the MM Algorithm
        mat_valid<-as.matrix(dens_valid)
        mat_test<-as.matrix(dens_out)
        optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
        finalw_MM_aug_new1<-optim_MM$finalparams
        
        #Calculate the predictive density by the ensemble
        model_weights_simul[[i]]<-finalw_MM_aug_new1
        optim_MM_par[[i]] <- optim_MM
    }
    model_weights_simul_sims[[sim]] <- model_weights_simul
    optim_MM_sims[[sim]] <- optim_MM_par
    }
    
    return (list(model_weights_simul_sims=model_weights_simul_sims, partition_func = partition_func, optim_MM=optim_MM_sims))
}

################################################################################
### Function for fitting ensemble models for different partitions for size 40 triangle
################################################################################


fit_all_partition_ensembles_40 <- function(meta_dens_valid, meta_dens_outsample) {
    source('ensemble/defining_partitions_triangle_40.R')
  
    ADLP_par0 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par0)
    
    ADLP_par1 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par1_2_40)
    
    ADLP_par2 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par2_2_40)
    
    ADLP_par3 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par3_2_40)
    
    ADLP_par4 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par4_2_40)
    
    ADLP_par5 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par5_2_40)
    
    ADLP_par6 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par6_2_40)
    
    ADLP_par7 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par7_2_40)
    
    ADLP_par8 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par8_2_40)
    
    ADLP_par9 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par9_2_40)
    
    ADLP_par10 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par10_2_40)
    
    ADLP_par11 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par11_2_40)
    
    ADLP_par12 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par12_2_40)
    
    ADLP_par13 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par13_2_40)
    
    ADLP_par14 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par14_2_40)
    
    ADLP_par15 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par15_2_40)
    
    ADLP_par16 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par16_2_40)
    
    ADLP_par17 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par17_2_40)
    
    ADLP_par18 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par18_2_40)
    
    all_ensembles <- list(
      ADLP_par0=ADLP_par0,
      ADLP_par1=ADLP_par1,
      ADLP_par2=ADLP_par2,
      ADLP_par3=ADLP_par3,
      ADLP_par4=ADLP_par4,
      ADLP_par5=ADLP_par5,
      ADLP_par6=ADLP_par6,
      ADLP_par7=ADLP_par7,
      ADLP_par8=ADLP_par8,
      ADLP_par9=ADLP_par9,
      ADLP_par10=ADLP_par10,
      ADLP_par11=ADLP_par11,
      ADLP_par12=ADLP_par12,
      ADLP_par13=ADLP_par13,
      ADLP_par14=ADLP_par14,
      ADLP_par15=ADLP_par15,
      ADLP_par16=ADLP_par16,
      ADLP_par17=ADLP_par17,
      ADLP_par18=ADLP_par18
    )
    return (all_ensembles)
}


fit_all_partition_ensembles_40_3par <- function(meta_dens_valid, meta_dens_outsample) {
  source('ensemble/defining_partitions_triangle_40.R')
  
  ADLP_sub3_par1 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par1_3_40)
  
  ADLP_sub3_par2 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par2_3_40)
  
  ADLP_sub3_par3 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par3_3_40)
  
  ADLP_sub3_par4 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par4_3_40)
  
  all_ensembles <- list(
    ADLP_sub3_par1=ADLP_sub3_par1,
    ADLP_sub3_par2=ADLP_sub3_par2,
    ADLP_sub3_par3=ADLP_sub3_par3,
    ADLP_sub3_par4=ADLP_sub3_par4
  )
  return (all_ensembles)
}

################################################################################
### Function for fitting ensemble models for different partitions for size 20 triangle
################################################################################


fit_all_partition_ensembles_20 <- function(meta_dens_valid, meta_dens_outsample) {
    source('ensemble/defining_partitions_triangle_20.R')
    
    ADLP_par0 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par0)
    
    ADLP_par1 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par1_2_20)
    
    ADLP_par2 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par2_2_20)
    
    ADLP_par3 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par3_2_20)
    
    ADLP_par4 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par4_2_20)
    
    ADLP_par5 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par5_2_20)
    
    ADLP_par6 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par6_2_20)
    
    ADLP_par7 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par7_2_20)
    
    ADLP_par8 <- fit_ensemble_model(meta_dens_valid, meta_dens_outsample, par8_2_20)
    
    all_ensembles <- list(
        ADLP_par0=ADLP_par0,
        ADLP_par1=ADLP_par1,
        ADLP_par2=ADLP_par2,
        ADLP_par3=ADLP_par3,
        ADLP_par4=ADLP_par4,
        ADLP_par5=ADLP_par5,
        ADLP_par6=ADLP_par6,
        ADLP_par7=ADLP_par7,
        ADLP_par8=ADLP_par8
    )
    return (all_ensembles)
}

