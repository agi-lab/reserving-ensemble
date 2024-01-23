################################################################################
### This module defines two main scoring methods for the ensemble models. This 
###     is composed of the log score (LS) and CRPS
### 
### To calculate the log score, the density of the ensemble is firstly calculated
### 
### To calculate the CRPS, the distribution function of each component model
###     is firstly evaluated given the out of sample dataset. Then, for each 
###     parition, the CRPS of each data point is calculated and summed together.
################################################################################

library("coda")
if (!exists('tri.size')) {stop("Size of triangles 'tri.size' should be defined")}

################################################################################
### Module for calculating component densities
################################################################################

calc_dens_ensemble <- function(meta_dens_components, ensemble_models) {
    
    n.ensembles <- length(ensemble_models)
    all_dens_ensemble <- list()
    
    for (i in 1:n.ensembles){
        ensemble <- ensemble_models[[i]]
        dens_ensemble <- list()
        
        for (sim in 1:n.sims) {
            dens_partitions <- ensemble$partition_func(meta_dens_components[[sim]])
            n.partitions <- length(dens_partitions)
            ensemble_w <- ensemble$model_weights_simul_sims[[sim]]
            dens <- c()
            dens_ij <- c()
            
            for (j in 1:n.partitions) {
                
                meta_partition <- dens_partitions[[j]][, -c(1, 2)]
                
                dens_predict <- as.matrix(meta_partition) %*% ensemble_w[[j]]
                dens <- c(dens, dens_predict) 
                dens_ij <- c(dens_ij, paste(dens_partitions[[j]]$origin, "-", dens_partitions[[j]]$dev, sep = ""))
            }
            
            names(dens) <- dens_ij
            dens_ensemble[[sim]] <- dens[order(names(dens))]
        }
        all_dens_ensemble[[i]] <- dens_ensemble
    }
    
    names(all_dens_ensemble) <-  names(ensemble_models)
    
    return (all_dens_ensemble)
}

################################################################################
### Module for calculating CRPS metrics
################################################################################

calc_crps_ensemble_40 <- function(n.sims, sample.interval, component_models, ensemble_models) {
    
    I<-function(y,z) ifelse(y<=z,1,0)
    
    n.ensembles <- length(ensemble_models)
    
    #tau_Ga <- 5
    #tau_LN <- 5
    
    crps_ensembles <- list()
    for (j in 1:n.ensembles){
        crps_ensembles[[j]] <- list()
    }
    
    for (s in 1:n.sims) {
        set.seed(20200130+s)
        
        fit_ODP_GLM = component_models[[s]]$fit_ODP_GLM
        fit_GAGLM = component_models[[s]]$fit_GAGLM
        fit_LNGLM = component_models[[s]]$fit_LNGLM
        fit_ZAGA = component_models[[s]]$fit_ZAGA
        fit_ZALN = component_models[[s]]$fit_ZALN
        fit_ODPHo = component_models[[s]]$fit_ODPHo
        fit_GaHo = component_models[[s]]$fit_GaHo
        fit_LNHo = component_models[[s]]$fit_LNHo
        fit_ODPCal = component_models[[s]]$fit_ODPCal
        fit_GaCal = component_models[[s]]$fit_GaCal
        fit_LNCal = component_models[[s]]$fit_LNCal
        fit_SpNormal = component_models[[s]]$fit_SpNormal
        fit_SpGamma = component_models[[s]]$fit_SpGamma
        fit_SpLN = component_models[[s]]$fit_SpLN
        fit_GaGAMLSS = component_models[[s]]$fit_GaGAMLSS
        fit_LNGAMLSS = component_models[[s]]$fit_LNGAMLSS
        fit_PPCI = component_models[[s]]$fit_PPCI
        fit_PPCF = component_models[[s]]$fit_PPCF
        
        full_data <- read.csv(sprintf('simulation/triangle_%s-data/sim%s-full-data.csv', tri.size, s))
        in_sample <- claims_df_in_out(full_data)$train
        out_sample <- claims_df_in_out(full_data)$test
        
        z_l<-1
        z_u <- 2*round(max(out_sample$aggregate_claims),0)
        z<- seq(z_l, z_u, by = sample.interval)
        ## Define a scale factor:
        scale_factor <- (z_u-z_l)/(length(z)-1)
        
        data <- in_sample
        newdata <- out_sample
        
        # This should be a list of length 18 for each component, each with a 
        # length(z) x nrow(newdata) matrix
        
        sink(nullfile())
        pred_CDF <- list(
            cal_CDF_ODP(z, fit_param_ODP(fit_ODP_GLM, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GAGLM, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNGLM, data, newdata, tau_LN), new_y = T),
            cal_CDF_ZAGA(z, fit_param_ZAGA(fit_ZAGA, data, newdata), new_y = T),
            cal_CDF_ZALN(z, fit_param_ZALN(fit_ZALN, data, newdata), new_y = T),
            cal_CDF_ODP(z, fit_param_ODP(fit_ODPHo, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GaHo, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNHo, data, newdata, tau_LN), new_y = T),
            cal_CDF_ODP(z, fit_param_ODP(fit_ODPCal, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GaCal, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNCal, data, newdata, tau_LN), new_y = T),
            cal_CDF_Normal(z, fit_param_NO(fit_SpNormal, data, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_SpGamma, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_SpLN, data, newdata, tau_LN), new_y = T),
            cal_CDF_GA_Gamlss(z, fit_param_GAGamlss(fit_GaGAMLSS, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN_Gamlss(z, fit_param_LNGamlss(fit_LNGAMLSS, data, newdata, tau_LN), new_y = T),
            cal_CDF_PPCI(z, fit_param_PPCI(fit_PPCI$model, fit_PPCI$N, newdata), new_y = T),
            cal_CDF_PPCF(z, fit_param_PPCF(fit_PPCF$model_subCount, fit_PPCF$model_subPayments, fit_PPCF$N, data, newdata), new_y = T)
        )
        sink()
        
        for (j in 1:n.ensembles){
            ensemble <- ensemble_models[[j]]
            
            out_partitions <- ensemble$partition_func(out_sample)
            n.partitions <- length(out_partitions)
            
            ensemble_w <- ensemble$model_weights_simul_sims[[s]]
            crps <- c()
            crps_ij <- c()
            
            for (k in 1:n.partitions) {
                y = out_partitions[[k]][ , "aggregate_claims"]
                partition_ind <- rownames(newdata) %in% rownames(out_partitions[[k]])
                w <- ensemble_w[[k]]
                
                pred_CDF_ensemble <- matrix(0, nrow = length(z), ncol = length(y))
                for (cdf_i in 1:length(w)) {
                    pred_CDF_ensemble <- pred_CDF_ensemble + pred_CDF[[cdf_i]][, partition_ind] * w[cdf_i]
                }
                
                diff <- colSums((pred_CDF_ensemble-t(outer(y, z, I)))^2)
                
                crps <- c(crps, diff)
                crps_ij <- c(crps_ij, paste(out_partitions[[k]]$origin, "-", out_partitions[[k]]$dev, sep = ""))
            }
            
            names(crps) = crps_ij
            crps_ensembles[[j]][[s]] <- scale_factor*crps[order(names(crps))]
        }
        
        if (floor(s %% (n.sims/4)) == 0) {
            print(sprintf('... Completed %s of %s simulations...', s, n.sims))
        }
    }
    
    names(crps_ensembles) <- names(ensemble_models)
    
    return(crps_ensembles)
}





## Define a function to calculate weighted crps: 
calc_weighted_crps_ensemble_40 <- function(n.sims, sample.interval, component_models, ensemble_models, focus) {
    
    I<-function(y,z) ifelse(y<=z,1,0)
    
    n.ensembles <- length(ensemble_models)
    ## Define a weight function:
    weights_func <- function(z, mu, sigma, focus){
        if(focus == "Center"){
            w <- dnorm(z, mu, sigma)
        }else if(focus == "Tails"){
            w <- 1-dnorm(z, mu, sigma)/dnorm(mu, mu, sigma)
        }else if(focus == "Right tail"){
            w <- pnorm(z, mu, sigma)
        }else if(focus == "Left tail"){
            w <- 1-pnorm(z, mu, sigma)
        }
        return(w)
    }
    
    #tau_Ga <- 5
    #tau_LN <- 5
    
    crps_ensembles <- list()
    for (j in 1:n.ensembles){
        crps_ensembles[[j]] <- list()
    }
    
    for (s in 1:n.sims) {
        set.seed(20200130+s)
        
        fit_ODP_GLM = component_models[[s]]$fit_ODP_GLM
        fit_GAGLM = component_models[[s]]$fit_GAGLM
        fit_LNGLM = component_models[[s]]$fit_LNGLM
        fit_ZAGA = component_models[[s]]$fit_ZAGA
        fit_ZALN = component_models[[s]]$fit_ZALN
        fit_ODPHo = component_models[[s]]$fit_ODPHo
        fit_GaHo = component_models[[s]]$fit_GaHo
        fit_LNHo = component_models[[s]]$fit_LNHo
        fit_ODPCal = component_models[[s]]$fit_ODPCal
        fit_GaCal = component_models[[s]]$fit_GaCal
        fit_LNCal = component_models[[s]]$fit_LNCal
        fit_SpNormal = component_models[[s]]$fit_SpNormal
        fit_SpGamma = component_models[[s]]$fit_SpGamma
        fit_SpLN = component_models[[s]]$fit_SpLN
        fit_GaGAMLSS = component_models[[s]]$fit_GaGAMLSS
        fit_LNGAMLSS = component_models[[s]]$fit_LNGAMLSS
        fit_PPCI = component_models[[s]]$fit_PPCI
        fit_PPCF = component_models[[s]]$fit_PPCF
        
        full_data <- read.csv(sprintf('simulation/triangle_%s-data/sim%s-full-data.csv', tri.size, s))
        in_sample <- claims_df_in_out(full_data)$train
        out_sample <- claims_df_in_out(full_data)$test
        
        
        z_l<-1
        z_u <- 2*round(max(out_sample$aggregate_claims),0)
        z<- seq(z_l, z_u, by = sample.interval)
        ## Define a scale factor:
        scale_factor <- (z_u-z_l)/(length(z)-1)
        
        data <- in_sample
        newdata <- out_sample
        
        # This should be a list of length 18 for each component, each with a 
        # length(z) x nrow(newdata) matrix
        
        sink(nullfile())
        pred_CDF <- list(
            cal_CDF_ODP(z, fit_param_ODP(fit_ODP_GLM, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GAGLM, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNGLM, data, newdata, tau_LN), new_y = T),
            cal_CDF_ZAGA(z, fit_param_ZAGA(fit_ZAGA, data, newdata), new_y = T),
            cal_CDF_ZALN(z, fit_param_ZALN(fit_ZALN, data, newdata), new_y = T),
            cal_CDF_ODP(z, fit_param_ODP(fit_ODPHo, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GaHo, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNHo, data, newdata, tau_LN), new_y = T),
            cal_CDF_ODP(z, fit_param_ODP(fit_ODPCal, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GaCal, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNCal, data, newdata, tau_LN), new_y = T),
            cal_CDF_Normal(z, fit_param_NO(fit_SpNormal, data, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_SpGamma, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_SpLN, data, newdata, tau_LN), new_y = T),
            cal_CDF_GA_Gamlss(z, fit_param_GAGamlss(fit_GaGAMLSS, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN_Gamlss(z, fit_param_LNGamlss(fit_LNGAMLSS, data, newdata, tau_LN), new_y = T),
            cal_CDF_PPCI(z, fit_param_PPCI(fit_PPCI$model, fit_PPCI$N, newdata), new_y = T),
            cal_CDF_PPCF(z, fit_param_PPCF(fit_PPCF$model_subCount, fit_PPCF$model_subPayments, fit_PPCF$N, data, newdata), new_y = T)
        )
        sink()
        
        for (j in 1:n.ensembles){
            ensemble <- ensemble_models[[j]]
            
            out_partitions <- ensemble$partition_func(out_sample)
            n.partitions <- length(out_partitions)
            
            ensemble_w <- ensemble$model_weights_simul_sims[[s]]
            crps <- c()
            crps_ij <- c()
            
            for (k in 1:n.partitions) {
                y = out_partitions[[k]][ , "aggregate_claims"]
                partition_ind <- rownames(newdata) %in% rownames(out_partitions[[k]])
                w <- ensemble_w[[k]]
                
                pred_CDF_ensemble <- matrix(0, nrow = length(z), ncol = length(y))
                for (cdf_i in 1:length(w)) {
                    pred_CDF_ensemble <- pred_CDF_ensemble + pred_CDF[[cdf_i]][, partition_ind] * w[cdf_i]
                }
                
                ## Calculate the mean and sigma used for the weights:
                mu_w <- mean(y)
                sigma_w <- sd(y)
                ### Calculate the weights: 
                crps_w <- weights_func(z, mu_w, sigma_w, focus)
                
                diff <- colSums(sweep((pred_CDF_ensemble-t(outer(y, z, I)))^2, MARGIN = 1, crps_w, FUN = "*"))
                
                crps <- c(crps, diff)
                crps_ij <- c(crps_ij, paste(out_partitions[[k]]$origin, "-", out_partitions[[k]]$dev, sep = ""))
            }
            
            names(crps) = crps_ij
            crps_ensembles[[j]][[s]] <- scale_factor*crps[order(names(crps))]
        }
        
        if (floor(s %% (n.sims/4)) == 0) {
            print(sprintf('... Completed %s of %s simulations...', s, n.sims))
        }
    }
    
    names(crps_ensembles) <- names(ensemble_models)
    
    return(crps_ensembles)
}










### We exclude ZALN and ZAGA as there are no zero incremental claims for 10x10 and 20x20 

calc_crps_ensemble_20 <- function(n.sims, component_models, ensemble_models) {
    
    I<-function(y,z) ifelse(y<=z,1,0)
    
    n.ensembles <- length(ensemble_models)
    
    #tau_Ga <- 5
    #tau_LN <- 5
    
    crps_ensembles <- list()
    for (j in 1:n.ensembles){
        crps_ensembles[[j]] <- list()
    }
    
    for (s in 1:n.sims) {
        set.seed(20200130+s)
        
        fit_ODP_GLM = component_models[[s]]$fit_ODP_GLM
        fit_GAGLM = component_models[[s]]$fit_GAGLM
        fit_LNGLM = component_models[[s]]$fit_LNGLM
        fit_ODPHo = component_models[[s]]$fit_ODPHo
        fit_GaHo = component_models[[s]]$fit_GaHo
        fit_LNHo = component_models[[s]]$fit_LNHo
        fit_ODPCal = component_models[[s]]$fit_ODPCal
        fit_GaCal = component_models[[s]]$fit_GaCal
        fit_LNCal = component_models[[s]]$fit_LNCal
        fit_SpNormal = component_models[[s]]$fit_SpNormal
        fit_SpGamma = component_models[[s]]$fit_SpGamma
        fit_SpLN = component_models[[s]]$fit_SpLN
        fit_GaGAMLSS = component_models[[s]]$fit_GaGAMLSS
        fit_LNGAMLSS = component_models[[s]]$fit_LNGAMLSS
        fit_PPCI = component_models[[s]]$fit_PPCI
        fit_PPCF = component_models[[s]]$fit_PPCF
        
        full_data <- read.csv(sprintf('simulation/triangle_%s-data/sim%s-full-data.csv', tri.size, s))
        in_sample <- claims_df_in_out(full_data)$train
        out_sample <- claims_df_in_out(full_data)$test
        
        z_l<-1
        z_u <- 2*round(max(out_sample$aggregate_claims),0)
        z<-z_l:z_u
        
        
        
        data <- in_sample
        newdata <- out_sample
        
        # This should be a list of length 18 for each component, each with a 
        # length(z) x nrow(newdata) matrix
        
        sink(nullfile())
        pred_CDF <- list(
            cal_CDF_ODP(z, fit_param_ODP(fit_ODP_GLM, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GAGLM, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNGLM, data, newdata, tau_LN), new_y = T),
            cal_CDF_ODP(z, fit_param_ODP(fit_ODPHo, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GaHo, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNHo, data, newdata, tau_LN), new_y = T),
            cal_CDF_ODP(z, fit_param_ODP(fit_ODPCal, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_GaCal, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_LNCal, data, newdata, tau_LN), new_y = T),
            cal_CDF_Normal(z, fit_param_NO(fit_SpNormal, data, newdata), new_y = T),
            cal_CDF_GA(z, fit_param_GA(fit_SpGamma, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN(z, fit_param_LN(fit_SpLN, data, newdata, tau_LN), new_y = T),
            cal_CDF_GA_Gamlss(z, fit_param_GAGamlss(fit_GaGAMLSS, data, newdata, tau_Ga), new_y = T),
            cal_CDF_LN_Gamlss(z, fit_param_LNGamlss(fit_LNGAMLSS, data, newdata, tau_LN), new_y = T),
            cal_CDF_PPCI(z, fit_param_PPCI(fit_PPCI$model, fit_PPCI$N, newdata), new_y = T),
            cal_CDF_PPCF(z, fit_param_PPCF(fit_PPCF$model_subCount, fit_PPCF$model_subPayments, fit_PPCF$N, data, newdata), new_y = T)
        )
        sink()
        
        
        for (j in 1:n.ensembles){
            ensemble <- ensemble_models[[j]]
            
            out_partitions <- ensemble$partition_func(out_sample)
            n.partitions <- length(out_partitions)
            
            ensemble_w <- ensemble$model_weights_simul_sims[[s]]
            crps <- c()
            crps_ij <- c()
            
            for (k in 1:n.partitions) {
                y = out_partitions[[k]][ , "aggregate_claims"]
                partition_ind <- rownames(newdata) %in% rownames(out_partitions[[k]])
                w <- ensemble_w[[k]]
                
                pred_CDF_ensemble <- matrix(0, nrow = length(z), ncol = length(y))
                for (cdf_i in 1:length(w)) {
                    pred_CDF_ensemble <- pred_CDF_ensemble + pred_CDF[[cdf_i]][, partition_ind] * w[cdf_i]
                }
                
                diff <- colSums((pred_CDF_ensemble-t(outer(y, z, I)))^2)
                
                crps <- c(crps, diff)
                crps_ij <- c(crps_ij, paste(out_partitions[[k]]$origin, "-", out_partitions[[k]]$dev, sep = ""))
            }
            
            names(crps) = crps_ij
            crps_ensembles[[j]][[s]] <- crps[order(names(crps))]
        }
        
        if (floor(s %% (n.sims/4)) == 0) {
            print(sprintf('... Completed %s of %s simulations...', s, n.sims))
        }
    }
    
    names(crps_ensembles) <- names(ensemble_models)
    
    return(crps_ensembles)
}

################################################################################
### Module for calculating DM and adj DM test stats
################################################################################

calc_DM_test_stat <- function(model1_scores, model2_scores) {
    # Assume dim of scores is (n.sample, n.sim)
    n <- ncol(model1_scores)
    
    F_bar <- apply(model1_scores, MARGIN = 1, FUN = mean)
    S_bar <- apply(model2_scores, MARGIN = 1, FUN = mean)
    sigma2_1 <- sqrt(1/n * rowSums((model1_scores - model2_scores)^2))
    
    test_stat <- (sqrt(n) * (F_bar - S_bar)) / (sigma2_1)
    return (test_stat)
}


calc_adj_DM_test_stat <- function(model1_scores, model2_scores) {
    n <- ncol(model1_scores)
    d_bar <- apply(model1_scores - model2_scores, MARGIN = 1, FUN = mean)
    ##Calculate the Spectral Density at Zero
    fd_0 <- apply(model1_scores - model2_scores, MARGIN = 1, FUN = function(x) spectrum0.ar(x)$spec)
    
    test_stat <- d_bar / sqrt(fd_0 / n)
}

sim_data_to_matrix <- function(model_scores) {
    # Assume model_scores is in the format (list of lists):
    # - Simulation
    # - Data
    # We want to get the following list
    #   - simulation 
    #     - data
    # and turn it into (data, simulation) matrix
    
    return (t(matrix(unlist(model_scores), ncol = n.sims)))
}

calc_DM_test_models <- function(ensemble_scores, model1_ind, model2_ind) {
    # Assume ensemble_scores is in the format:
    # - Ensemble
    # - Simulation
    # - Data
    # We want to get the following list
    #   - simulation 
    #     - data
    # and turn it into (data, simulation) matrix
    
    model1_scores <- sim_data_to_matrix(ensemble_scores[[model1_ind]])
    model2_scores <- sim_data_to_matrix(ensemble_scores[[model2_ind]])
    
    return (na.omit(calc_DM_test_stat(model1_scores, model2_scores)))
}

calc_adj_DM_test_models <- function(ensemble_scores, model1_ind, model2_ind) {
    # Assume ensemble_scores is in the format:
    # - Ensemble
    # - Simulation
    # - Data
    # We want to get the following list
    #   - simulation 
    #     - data
    # and turn it into (data, simulation) matrix
    
    model1_scores <- sim_data_to_matrix(ensemble_scores[[model1_ind]])
    model2_scores <- sim_data_to_matrix(ensemble_scores[[model2_ind]])
    
    return (na.omit(calc_adj_DM_test_stat(model1_scores, model2_scores)))
}