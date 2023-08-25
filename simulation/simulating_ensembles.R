################################################################################
## This module simulations claims triangles for analysis assuming the
##  fitted ensemble weights and fitted component models
################################################################################

if (!exists('tri.size')) {stop("Size of triangles 'tri.size' should be defined")}

################################################################################
### Function for simulating claim values from ensembles
################################################################################
 
simulate_reserves_40 <- function(n.sims, sims, component_models, ensemble_models) {
    sink(nullfile())
    
    U <- runif(sims)
    
    n.ensembles <- length(ensemble_models)
    
    all_reserve_simulations <- list()
    for (j in 1:n.ensembles){
        all_reserve_simulations[[j]] <- matrix(NA,nrow=sims,ncol=n.sims)
    }
    for (s in 1:n.sims) {
        set.seed(20200130+s)
        
        full_data <- read.csv(sprintf('simulation/triangle_%s-data/sim%s-full-data.csv', tri.size, s))
        in_sample <- claims_df_in_out(full_data)$train
        out_sample <- claims_df_in_out(full_data)$test
        
        data <- in_sample
        
        for (j in 1:n.ensembles){
            ensemble <- ensemble_models[[j]]
            out_sample_partition <- ensemble$partition_func(out_sample)
            n.partitions <- length(out_sample_partition)
            ensemble_w <- ensemble$model_weights_simul_sims[[s]]
            ensemble_reserve <- c()
            
            for (k in 1:n.partitions) {
                mix_ind_subset<-findInterval(U, cumsum(unlist(ensemble_w[[k]])))+1
                
                newdata <- out_sample_partition[[k]]
                
                component_partition_param <- predict_component_partitions(component_models, s, data, newdata)
                
                individual_claim_simulations<-matrix(NA,nrow=nrow(newdata),ncol=sims)
                for (i in 1:sims){
                    if(mix_ind_subset[i]==1){individual_claim_simulations[,i]<-as.vector(simulate_ODP(component_partition_param$param_ODP_GLM))}
                    else if(mix_ind_subset[i]==2){individual_claim_simulations[,i]<-as.vector(simulate_GA(component_partition_param$param_GAGLM))}
                    else if(mix_ind_subset[i]==3){individual_claim_simulations[,i]<-as.vector(simulate_LN(component_partition_param$param_LNGLM))}
                    else if(mix_ind_subset[i]==4){individual_claim_simulations[,i]<-as.vector(simulate_ZAGA(component_partition_param$param_ZAGA))}
                    else if(mix_ind_subset[i]==5){individual_claim_simulations[,i]<-as.vector(simulate_ZALN(component_partition_param$param_ZALN))}
                    else if(mix_ind_subset[i]==6){individual_claim_simulations[,i]<-as.vector(simulate_ODP(component_partition_param$param_ODPHo))}
                    else if(mix_ind_subset[i]==7){individual_claim_simulations[,i]<-as.vector(simulate_GA(component_partition_param$param_GaHo))}
                    else if(mix_ind_subset[i]==8){individual_claim_simulations[,i]<-as.vector(simulate_LN(component_partition_param$param_LNHo))}
                    else if(mix_ind_subset[i]==9){individual_claim_simulations[,i]<-as.vector(simulate_ODP(component_partition_param$param_ODPCal))}
                    else if(mix_ind_subset[i]==10){individual_claim_simulations[,i]<-as.vector(simulate_GA(component_partition_param$param_GaCal))}
                    else if(mix_ind_subset[i]==11){individual_claim_simulations[,i]<-as.vector(simulate_LN(component_partition_param$param_LNCal))}
                    else if(mix_ind_subset[i]==12){individual_claim_simulations[,i]<-as.vector(simulate_NO(component_partition_param$param_SpNormal))}
                    else if(mix_ind_subset[i]==13){individual_claim_simulations[,i]<-as.vector(simulate_GA(component_partition_param$param_SpGamma))}
                    else if(mix_ind_subset[i]==14){individual_claim_simulations[,i]<-as.vector(simulate_LN(component_partition_param$param_SpLN))}
                    else if(mix_ind_subset[i]==15){individual_claim_simulations[,i]<-as.vector(simulate_GAGamlss(component_partition_param$param_GaGAMLSS))}
                    else if(mix_ind_subset[i]==16){individual_claim_simulations[,i]<-as.vector(simulate_LNGamlss(component_partition_param$param_LNGAMLSS))}
                    else if(mix_ind_subset[i]==17){individual_claim_simulations[,i]<-as.vector(simulate_PPCI(component_partition_param$param_PPCI))}
                    else{individual_claim_simulations[,i]<-as.vector(simulate_PPCF(component_partition_param$param_PPCF))}
                }       
                
                ensemble_reserve <- rbind(ensemble_reserve, individual_claim_simulations) 
            }
            
            all_reserve_simulations[[j]][, s] <- apply(ensemble_reserve,FUN=function(x) sum(x*10000),MARGIN=2)
        }
    }
    
    sink()
    return (all_reserve_simulations)
}

predict_component_partitions_40 <- function(component_models, s, data, newdata) {
    tau_Ga <- 5
    tau_LN <- 5
    
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
    
    param_ODP_GLM <- fit_param_ODP(fit_ODP_GLM, newdata=newdata)
    param_GAGLM <- fit_param_GA(fit_GAGLM, data=data, newdata=newdata,tau=tau_Ga)
    param_LNGLM <- fit_param_LN(fit_LNGLM, data=data, newdata=newdata,tau=tau_LN)
    param_ZAGA <- fit_param_ZAGA(fit_ZAGA, data=data, newdata=newdata)
    param_ZALN <- fit_param_ZALN(fit_ZALN, data=data, newdata=newdata)
    param_ODPHo <- fit_param_ODP(fit_ODPHo, newdata=newdata)
    param_GaHo <- fit_param_GA(fit_GaHo, data=data, newdata=newdata,tau=tau_Ga)
    param_LNHo <- fit_param_LN(fit_LNHo, data=data, newdata=newdata,tau=tau_LN)
    param_ODPCal <- fit_param_ODP(fit_ODPCal, newdata)
    param_GaCal <- fit_param_GA(fit_GaCal, data=data, newdata,tau=tau_Ga)
    param_LNCal <- fit_param_LN(fit_LNCal, data=data, newdata,tau=tau_LN)
    param_SpNormal <- fit_param_NO(fit_SpNormal, data=data, newdata=newdata)
    param_SpGamma <- fit_param_GA(fit_SpGamma, data=data, newdata,tau=tau_Ga)
    param_SpLN <- fit_param_LN(fit_SpLN, data=data, newdata,tau=tau_LN)
    param_GaGAMLSS <- fit_param_GAGamlss(fit_GaGAMLSS, data=data, newdata=newdata,tau=tau_Ga)
    param_LNGAMLSS <- fit_param_LNGamlss(fit_LNGAMLSS, data=data, newdata=newdata,tau=tau_LN)
    param_PPCI <- fit_param_PPCI(fit_PPCI$model, fit_PPCI$N, newdata = newdata)
    param_PPCF <- fit_param_PPCF(fit_PPCF$model_subCount, fit_PPCF$model_subPayments, N=fit_PPCF$N, data=data, newdata=newdata)
    
    return (
        list(
            param_ODP_GLM = param_ODP_GLM,
            param_GAGLM = param_GAGLM,
            param_LNGLM = param_LNGLM,
            param_ZAGA = param_ZAGA,
            param_ZALN = param_ZALN,
            param_ODPHo = param_ODPHo,
            param_GaHo = param_GaHo,
            param_LNHo = param_LNHo,
            param_ODPCal = param_ODPCal,
            param_GaCal = param_GaCal,
            param_LNCal = param_LNCal,
            param_SpNormal = param_SpNormal,
            param_SpGamma = param_SpGamma,
            param_SpLN = param_SpLN,
            param_GaGAMLSS = param_GaGAMLSS,
            param_LNGAMLSS = param_LNGAMLSS,
            param_PPCI = param_PPCI,
            param_PPCF = param_PPCF
        )
    )
}

################################################################################
### Function for simulating claim values from ensembles for size 20 triangles
################################################################################

simulate_reserves_20 <- function(n.sims, sims, component_models, ensemble_models) {
    sink(nullfile())
    
    U <- runif(sims)
    
    n.ensembles <- length(ensemble_models)
    
    all_reserve_simulations <- list()
    for (j in 1:n.ensembles){
        all_reserve_simulations[[j]] <- matrix(NA,nrow=sims,ncol=n.sims)
    }
    for (s in 1:n.sims) {
        set.seed(20200130+s)
        
        full_data <- read.csv(sprintf('simulation/triangle_%s-data/sim%s-full-data.csv', tri.size, s))
        in_sample <- claims_df_in_out(full_data)$train
        out_sample <- claims_df_in_out(full_data)$test
        
        data <- in_sample
        
        for (j in 1:n.ensembles){
            ensemble <- ensemble_models[[j]]
            out_sample_partition <- ensemble$partition_func(out_sample)
            n.partitions <- length(out_sample_partition)
            ensemble_w <- ensemble$model_weights_simul_sims[[s]]
            ensemble_reserve <- c()
            
            for (k in 1:n.partitions) {
                mix_ind_subset<-findInterval(U, cumsum(unlist(ensemble_w[[k]])))+1
                
                newdata <- out_sample_partition[[k]]
                
                component_partition_param <- predict_component_partitions(component_models, s, data, newdata)
                
                individual_claim_simulations<-matrix(NA,nrow=nrow(newdata),ncol=sims)
                for (i in 1:sims){
                    if(mix_ind_subset[i]==1){individual_claim_simulations[,i]<-as.vector(simulate_ODP(component_partition_param$param_ODP_GLM))}
                    else if(mix_ind_subset[i]==2){individual_claim_simulations[,i]<-as.vector(simulate_GA(component_partition_param$param_GAGLM))}
                    else if(mix_ind_subset[i]==3){individual_claim_simulations[,i]<-as.vector(simulate_LN(component_partition_param$param_LNGLM))}
                    # else if(mix_ind_subset[i]==4){individual_claim_simulations[,i]<-as.vector(simulate_ZAGA(component_partition_param$param_ZAGA))}
                    # else if(mix_ind_subset[i]==5){individual_claim_simulations[,i]<-as.vector(simulate_ZALN(component_partition_param$param_ZALN))}
                    else if(mix_ind_subset[i]==4){individual_claim_simulations[,i]<-as.vector(simulate_ODP(component_partition_param$param_ODPHo))}
                    else if(mix_ind_subset[i]==5){individual_claim_simulations[,i]<-as.vector(simulate_GA(component_partition_param$param_GaHo))}
                    else if(mix_ind_subset[i]==6){individual_claim_simulations[,i]<-as.vector(simulate_LN(component_partition_param$param_LNHo))}
                    else if(mix_ind_subset[i]==7){individual_claim_simulations[,i]<-as.vector(simulate_ODP(component_partition_param$param_ODPCal))}
                    else if(mix_ind_subset[i]==8){individual_claim_simulations[,i]<-as.vector(simulate_GA(component_partition_param$param_GaCal))}
                    else if(mix_ind_subset[i]==9){individual_claim_simulations[,i]<-as.vector(simulate_LN(component_partition_param$param_LNCal))}
                    else if(mix_ind_subset[i]==10){individual_claim_simulations[,i]<-as.vector(simulate_NO(component_partition_param$param_SpNormal))}
                    else if(mix_ind_subset[i]==11){individual_claim_simulations[,i]<-as.vector(simulate_GA(component_partition_param$param_SpGamma))}
                    else if(mix_ind_subset[i]==12){individual_claim_simulations[,i]<-as.vector(simulate_LN(component_partition_param$param_SpLN))}
                    else if(mix_ind_subset[i]==13){individual_claim_simulations[,i]<-as.vector(simulate_GAGamlss(component_partition_param$param_GaGAMLSS))}
                    else if(mix_ind_subset[i]==14){individual_claim_simulations[,i]<-as.vector(simulate_LNGamlss(component_partition_param$param_LNGAMLSS))}
                    else if(mix_ind_subset[i]==15){individual_claim_simulations[,i]<-as.vector(simulate_PPCI(component_partition_param$param_PPCI))}
                    else{individual_claim_simulations[,i]<-as.vector(simulate_PPCF(component_partition_param$param_PPCF))}
                }       
                
                ensemble_reserve <- rbind(ensemble_reserve, individual_claim_simulations) 
            }
            
            all_reserve_simulations[[j]][, s] <- apply(ensemble_reserve,FUN=function(x) sum(x*10000),MARGIN=2)
        }
    }
    
    sink()
    return (all_reserve_simulations)
}

predict_component_partitions_20 <- function(component_models, s, data, newdata) {
    tau_Ga <- 5
    tau_LN <- 5
    
    fit_ODP_GLM = component_models[[s]]$fit_ODP_GLM
    fit_GAGLM = component_models[[s]]$fit_GAGLM
    fit_LNGLM = component_models[[s]]$fit_LNGLM
    # fit_ZAGA = component_models[[s]]$fit_ZAGA
    # fit_ZALN = component_models[[s]]$fit_ZALN
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
    
    param_ODP_GLM <- fit_param_ODP(fit_ODP_GLM, newdata=newdata)
    param_GAGLM <- fit_param_GA(fit_GAGLM, data=data, newdata=newdata,tau=tau_Ga)
    param_LNGLM <- fit_param_LN(fit_LNGLM, data=data, newdata=newdata,tau=tau_LN)
    # param_ZAGA <- fit_param_ZAGA(fit_ZAGA, data=data, newdata=newdata)
    # param_ZALN <- fit_param_ZALN(fit_ZALN, data=data, newdata=newdata)
    param_ODPHo <- fit_param_ODP(fit_ODPHo, newdata=newdata)
    param_GaHo <- fit_param_GA(fit_GaHo, data=data, newdata=newdata,tau=tau_Ga)
    param_LNHo <- fit_param_LN(fit_LNHo, data=data, newdata=newdata,tau=tau_LN)
    param_ODPCal <- fit_param_ODP(fit_ODPCal, newdata)
    param_GaCal <- fit_param_GA(fit_GaCal, data=data, newdata,tau=tau_Ga)
    param_LNCal <- fit_param_LN(fit_LNCal, data=data, newdata,tau=tau_LN)
    param_SpNormal <- fit_param_NO(fit_SpNormal, data=data, newdata=newdata)
    param_SpGamma <- fit_param_GA(fit_SpGamma, data=data, newdata,tau=tau_Ga)
    param_SpLN <- fit_param_LN(fit_SpLN, data=data, newdata,tau=tau_LN)
    param_GaGAMLSS <- fit_param_GAGamlss(fit_GaGAMLSS, data=data, newdata=newdata,tau=tau_Ga)
    param_LNGAMLSS <- fit_param_LNGamlss(fit_LNGAMLSS, data=data, newdata=newdata,tau=tau_LN)
    param_PPCI <- fit_param_PPCI(fit_PPCI$model, fit_PPCI$N, newdata = newdata)
    param_PPCF <- fit_param_PPCF(fit_PPCF$model_subCount, fit_PPCF$model_subPayments, N=fit_PPCF$N, data=data, newdata=newdata)
    
    return (
        list(
            param_ODP_GLM = param_ODP_GLM,
            param_GAGLM = param_GAGLM,
            param_LNGLM = param_LNGLM,
            # param_ZAGA = param_ZAGA,
            # param_ZALN = param_ZALN,
            param_ODPHo = param_ODPHo,
            param_GaHo = param_GaHo,
            param_LNHo = param_LNHo,
            param_ODPCal = param_ODPCal,
            param_GaCal = param_GaCal,
            param_LNCal = param_LNCal,
            param_SpNormal = param_SpNormal,
            param_SpGamma = param_SpGamma,
            param_SpLN = param_SpLN,
            param_GaGAMLSS = param_GaGAMLSS,
            param_LNGAMLSS = param_LNGAMLSS,
            param_PPCI = param_PPCI,
            param_PPCF = param_PPCF
        )
    )
}