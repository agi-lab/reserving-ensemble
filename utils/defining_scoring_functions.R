library("coda")
if (!exists('tri.size')) {stop("Size of triangles 'tri.size' should be defined")}

################################################################################
### Module for calculating component densities
################################################################################

calc_dens_ensemble <- function(meta_dens_components, ensemble_models) {
    
    #ensemble_models <- fitted_components_outsample
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

calc_crps_ensemble <- function(n.sims, component_models, ensemble_models) {
    
    I<-function(y,z) ifelse(y<=z,1,0)
    
    n.ensembles <- length(ensemble_models)
    #component_models <- fitted_components_outsample
    #ensemble_models <- all_ensembles
    
    tau_Ga <- 5
    tau_LN <- 5
    
    crps_ensembles <- list()
    for (j in 1:n.ensembles){
        crps_ensembles[[j]] <- list()
    }
    
    for (s in 1:n.sims) {
        # s <- 1
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
        #JL
        z_u <- 2*round(max(out_sample$aggregate_claims),0)
        z<-z_l:z_u
        
        
        
        data <- in_sample
        newdata <- out_sample
        
        # This should be a list of length 18 for each component, each with a 
        # length(z) x nrow(newdata) matrix
        
        pred_CDF <- list(
            cal_CDF_ODP(z, fit_ODP_GLM, newdata, new_y = T),
            cal_CDF_GA(z, fit_GAGLM, tau_Ga, data, newdata, new_y = T),
            cal_CDF_LN(z, fit_LNGLM, tau_LN, data, newdata, new_y = T),
            cal_CDF_ZAGA(z, fit_ZAGA, data, newdata, new_y = T),
            cal_CDF_ZALN(z, fit_ZALN, data, newdata, new_y = T),
            cal_CDF_ODP(z, fit_ODPHo, newdata, new_y = T),
            cal_CDF_GA(z, fit_GaHo, tau_Ga, data, newdata, new_y = T),
            cal_CDF_LN(z, fit_LNHo, tau_LN, data, newdata, new_y = T),
            cal_CDF_ODP(z, fit_ODPCal, newdata, new_y = T),
            cal_CDF_GA(z, fit_GaCal, tau_Ga, data, newdata, new_y = T),
            cal_CDF_LN(z, fit_LNCal, tau_LN, data, newdata, new_y = T),
            cal_CDF_Normal(z, fit_SpNormal, data, newdata, new_y = T),
            cal_CDF_GA(z, fit_SpGamma, tau_Ga, data, newdata, new_y = T),
            cal_CDF_LN(z, fit_SpLN, tau_LN, data, newdata, new_y = T),
            cal_CDF_GA_Gamlss(z, fit_GaGAMLSS, tau_Ga, data, newdata, new_y = T),
            cal_CDF_LN_Gamlss(z, fit_LNGAMLSS, tau_LN, data, newdata, new_y = T),
            cal_CDF_PPCI(z, fit_PPCI$model, fit_PPCI$N, newdata, new_y = T),
            cal_CDF_PPCF(z, fit_PPCF$model_subCount, fit_PPCF$model_subPayments, fit_PPCF$N, data, newdata, new_y = T)
        )
        
        
        
        for (j in 1:n.ensembles){
            # j <- 2
            ensemble <- ensemble_models[[j]]
            
            out_partitions <- ensemble$partition_func(out_sample)
            n.partitions <- length(out_partitions)
            
            ensemble_w <- ensemble$model_weights_simul_sims[[s]]
            crps <- c()
            crps_ij <- c()
            
            for (k in 1:n.partitions) {
                # k <- 1
                y = out_partitions[[k]][ , "aggregate_claims"]
                partition_ind <- rownames(newdata) %in% rownames(out_partitions[[k]])
                w <- ensemble_w[[k]]
                
                #partition_ind <- 1
                pred_CDF_ensemble <- matrix(0, nrow = length(z), ncol = length(y))
                for (cdf_i in 1:length(w)) {
                    pred_CDF_ensemble <- pred_CDF_ensemble + pred_CDF[[cdf_i]][, partition_ind] * w[cdf_i]
                }
                
                # Testing code: 
                # pred_CDF_all_models <- matrix(NA, nrow = 702, ncol = 18)
               # for(i in 1:18){
            #        pred_CDF_all_models[, i] <- pred_CDF[[i]][, 1]
            #    }
                # View(pred_CDF_all_models%*%w)
                
                diff <- colSums((pred_CDF_ensemble-t(outer(y, z, I)))^2)
                
                crps <- c(crps, diff)
                crps_ij <- c(crps_ij, paste(out_partitions[[k]]$origin, "-", out_partitions[[k]]$dev, sep = ""))
            }
            
            names(crps) = crps_ij
            crps_ensembles[[j]][[s]] <- crps[order(names(crps))]
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
    # Example: Checks
    #model1_scores <- sim_data_to_matrix(ensemble_logS$ADLP_par0)
    #model2_scores <- sim_data_to_matrix(ensemble_logS$EW)
    n <- ncol(model1_scores)
    
    F_bar <- apply(model1_scores, MARGIN = 1, FUN = mean)
    S_bar <- apply(model2_scores, MARGIN = 1, FUN = mean)
    sigma2_1 <- sqrt(1/n * rowSums((model1_scores - model2_scores)^2))
    
    test_stat <- (sqrt(n) * (F_bar - S_bar)) / (sigma2_1)
    return (test_stat)
}

# Check
#sqrt(sum((model1_scores[2, ]-model2_scores[2, ])^2)/n)


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
    
    #model1_ind <- SLP_ind 
    #model2_ind <- EW_ind
    #ensemble_scores <- ensemble_logS
    
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