
################################################################################
## This file defines functions required for component models and their relevant
## distributions, including:
## - Density
## - Predictive parameters
## - Predictive density
## - Simulations
##    Matrix with as many rows as newdata and k columns, where k is the number of simulations
################################################################################

# The component models used in this code is discussed in the paper

################################################################################
# ODP GLMcc
################################################################################

dODP<-function(y,lambda,phi, new_y = F) {
     
     if (!new_y) {
         ifelse(rep(phi, length(y)) < 1e-6, dpois(floor(y), lambda), dpois(floor(y/phi), floor(lambda/phi))/phi)
     } else {
         ifelse(matrix(rep(phi, length(y)), ncol = length(phi)), 
                mapply(function(lambda) dpois(floor(y), lambda), lambda),
                mapply(function(lambda, phi) dpois(floor(y/phi), floor(lambda/phi))/phi, lambda, phi))
     }
 }

# dODP<-function(y, lambda,phi, new_y = F) {
#     if (!new_y) {
#         dtweedie(y = y, mu = lambda, phi = phi, power = 1)
#     } else {
#         mapply(function(lambda, phi) dtweedie(x = y, mu = lambda, phi = phi, power = 1), lambda, phi)
#     }
# }


pODP<-function(y, lambda,phi, new_y = F) {
    if (!new_y) {
        ptweedie(q = y, mu = lambda, phi = phi, power = 1)
    } else {
        mapply(function(lambda, phi) ptweedie(q = y, mu = lambda, phi = phi, power = 1), lambda, phi)
    }
}

fit_param_ODP<-function(model,newdata){
    pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
    pred_ODP_mu<-pred_ODP$fit
    pred_ODP_phi<-(pred_ODP$residual.scale)^2
    return(list(pred_ODP_mu=pred_ODP_mu, pred_ODP_phi=pred_ODP_phi))
}

##Calculate ODP density from ODP GLM
cal_dens_ODP<-function(y, param){
  pred_ODP_mu<-param$pred_ODP_mu
  pred_ODP_phi<-param$pred_ODP_phi
  return(dODP(y,lambda=pred_ODP_mu,phi=pred_ODP_phi))
}

##Calculate ODP mean for ODP GLM
cal_mu_ODP<-function(param){
  pred_ODP_mu<-param$pred_ODP_mu
  return(pred_ODP_mu)
}

##Calculate ODP mean for ODP GLM
cal_sigma_ODP<-function(param){
  pred_ODP_mu<-param$pred_ODP_mu
  pred_ODP_phi<-param$pred_ODP_phi
  pred_ODP_sigma<-pred_ODP_phi*pred_ODP_mu
  return(pred_ODP_sigma)
}

##Calculate ODP CDF from ODP GLM
cal_CDF_ODP<-function(y, param, new_y = F){
  pred_ODP_mu<-param$pred_ODP_mu
  pred_ODP_phi<-param$pred_ODP_phi
  
  if (!new_y) {return(pODP(y,lambda=pred_ODP_mu,phi=pred_ODP_phi))}
  else {return(mapply(
    function(lambda, phi)
      pODP(y,lambda=lambda,phi=phi, new_y=new_y), pred_ODP_mu, pred_ODP_phi)
    )}
}

simulate_ODP<-function(param){
  k=1
  simy<-replicate(k,rtweedie(length(param$pred_ODP_mu),xi=1,mu=param$pred_ODP_mu,phi=param$pred_ODP_phi))
  return(simy)
  # replicate the simulations process k times; the ODP distribution corresponding to tweedie distribution with power=1
}


################################################################################
# LN GLM
################################################################################

fit_param_LN<-function(model,data,newdata,tau){
    pred_mu<-predict(model,what="mu",data=data,newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",data=data,newdata=newdata,type="response")
    return(list(pred_mu=pred_mu, pred_sigma=pred_sigma, tau=tau))
}

##Calculate Log-Normal density from Log-Normal GLM
cal_dens_LN<-function(y, param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    return(dLNO(x=y+tau,mu=pred_mu,sigma=pred_sigma))
}

##Calculate LN mean for LN GLM
cal_mu_LN<-function(param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    mean<-exp(pred_mu+pred_sigma^2/2)
    return(mean)
}

##Calculate LN mean for LN GLM
cal_sigma_LN<-function(param){
    pred_sigma<-param$pred_sigma
    return(pred_sigma)
}

##Calculate Log-Normal CDF from Log-Normal GLM
cal_CDF_LN<-function(y, param, new_y = F){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    if (!new_y) {return(pLNO(q=y+tau,mu=pred_mu,sigma=pred_sigma))}
    else {return(mapply(
        function(mu, sigma)
            pLNO(q=y+tau,mu=mu,sigma=sigma), pred_mu, pred_sigma)
    )}
}

simulate_LN<-function(param){
    k=1
    simy<-replicate(k, rLNO(n=length(param$pred_mu),mu=param$pred_mu,sigma=param$pred_sigma)-param$tau)
    simy[simy<0]<-0
    return(simy)
}

################################################################################
# Gamma GLM
################################################################################

fit_param_GA<-function(model,data,newdata,tau){
    pred_mu<-predict(model,what="mu",data=data,newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",data=data,newdata=newdata,type="response")
    return(list(pred_mu=pred_mu, pred_sigma=pred_sigma, tau=tau))
}

##Calculate Gamma density from Gamma GLM
cal_dens_GA<-function(y, param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    return(dGA(x=y+tau,mu=pred_mu,sigma=pred_sigma))
}

##Calculate GA mean for GA GLM
cal_mu_GA<-function(param){
    pred_mu<-param$pred_mu
    return(pred_mu)
}

##Calculate GA mean for GA GLM
cal_sigma_GA<-function(param){
    pred_sigma<-param$pred_sigma
    return(pred_sigma)
}


##Calculate Gamma CDF from Gamma GLM
cal_CDF_GA<-function(y, param, new_y = F){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    if (!new_y) {return(pGA(q=y+tau,mu=pred_mu,sigma=pred_sigma))}
    else {return(mapply(
        function(mu, sigma)
            pGA(q=y+tau,mu=mu,sigma=sigma), pred_mu, pred_sigma)
    )}
}

simulate_GA<-function(param){
    k=1
    simy<-replicate(k, rGA(n=length(param$pred_mu),mu=param$pred_mu,sigma=param$pred_sigma)-param$tau)
    ###Assign the negative value to be zero
    simy[simy<0]<-0
    return(simy)
    # replicate the simulations process k times
}

################################################################################
# Normal
################################################################################

fit_param_NO<-function(model, data, newdata){
    pred_mu<-predict(model,what="mu",data=data,newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",data=data,newdata=newdata,type="response")
    return(list(pred_mu=pred_mu, pred_sigma=pred_sigma))
}

##Calculate Normal density from Normal GLM
cal_dens_Normal<-function(y, param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    return(dNO(x=y,mu=pred_mu,sigma=pred_sigma))
}

##Calculate Normal mean for Normal GLM
cal_mu_Normal<-function(param){
    pred_mu<-param$pred_mu
    return(pred_mu)
}

##Calculate Normal mean for Normal GLM
cal_sigma_Normal<-function(param){
    pred_sigma<-param$pred_sigma
    return(pred_sigma)
}

##Calculate Normal CDF from Normal GLM
cal_CDF_Normal<-function(y, param, new_y = F){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    if (!new_y) {return(pNO(q=y,mu=pred_mu,sigma=pred_sigma))}
    else {return(mapply(
        function(mu, sigma)
            pNO(q=y,mu=mu,sigma=sigma), pred_mu, pred_sigma)
    )}
}

simulate_NO<-function(param){
    k=1
    simy<-replicate(k,rNO(n=length(param$pred_mu),mu=param$pred_mu,sigma=param$pred_sigma))
    simy[simy<0]<-0
    return(simy)
}

################################################################################
# ZALN
################################################################################

#Create density function for ZALN
dZALN<-Zadj.d(family="LOGNO")

#Create CDF function for ZALN
pZALN<-Zadj.p(family="LOGNO")

rZALN<-Zadj.r(family="LOGNO")

fit_param_ZALN<-function(model,data,newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,parameter="mu",data=data,newdata=newdata,type="response")
    pred_sigma<-predict(model,parameter="sigma",data=data,newdata=newdata,type="response")
    pred_nu<-predict(model,parameter="xi0",data=data,newdata=newdata_nu,type="response")
    return(list(pred_sigma=pred_sigma, pred_mu=pred_mu, pred_nu=pred_nu))
}


##Calculate ZALN density from ZALN GLM
cal_dens_ZALN<-function(y, param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    pred_nu<-param$pred_nu
    return(dZALN(x=y,mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))
}

##Calculate ZALN mean for ZALN GLM
cal_mu_ZALN<-function(param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    pred_nu<-param$pred_nu
    mean<-(1-pred_nu)*exp(pred_mu+pred_sigma^2/2)
    return(mean)
}

##Calculate ZALN mean for ZALN GLM
cal_sigma_ZALN<-function(param){
    pred_sigma<-param$pred_sigma
    return(pred_sigma)
}

##Calculate ZALN CDF from ZALN GLM
cal_CDF_ZALN<-function(y, param, new_y = F){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    pred_nu<-param$pred_nu
    if (!new_y) {return(pZALN(q=y,mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))}
    else {return(mapply(
        function(mu, sigma, xi0)
            pZALN(q=y,mu=mu,sigma=sigma, xi0 = xi0), pred_mu, pred_sigma, pred_nu)
    )}
}

simulate_ZALN<-function(param){
    k=1
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    pred_nu<-mean(param$pred_nu)
    simy<-replicate(k,rZALN(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))
    return(simy)
}

################################################################################
# ZAGA
################################################################################

fit_param_ZAGA<-function(model,data,newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,what="mu",data=data,newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",data=data,newdata=newdata,type="response")
    pred_nu<-predict(model,what="nu",data=data,newdata=newdata_nu,type="response")
    return(list(pred_mu=pred_mu, pred_sigma=pred_sigma, pred_nu=pred_nu)) 
}

##Calculate ZAGA density from ZAGA GLM
cal_dens_ZAGA<-function(y, param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    pred_nu<-param$pred_nu
    return(dZAGA(x=y,mu=pred_mu,sigma=pred_sigma,nu=pred_nu))
}

##Calculate ZAGA mean for ZAGA GLM
cal_mu_ZAGA<-function(param){
    pred_mu<-param$pred_mu
    pred_nu<-param$pred_nu
    mean<-(1-pred_nu)*pred_mu
    return(mean)
}

##Calculate ZAGA mean for ZAGA GLM
cal_sigma_ZAGA<-function(param){
    pred_sigma<-param$pred_sigma
    return(pred_sigma)
}

##Calculate ZAGA CDF from ZAGA GLM
cal_CDF_ZAGA<-function(y, param, new_y = F){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    pred_nu<-param$pred_nu
    if (!new_y) {return(pZAGA(q=y,mu=pred_mu,sigma=pred_sigma,nu=pred_nu))}
    else {return(mapply(
        function(mu, sigma, nu)
            pZAGA(q=y,mu=mu,sigma=sigma, nu = nu), pred_mu, pred_sigma, pred_nu)
    )}
}

simulate_ZAGA<-function(param){
    k=1
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    pred_nu<-mean(param$pred_nu)
    simy<-replicate(k, rZAGA(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma,nu=pred_nu))
    return(simy)
}

################################################################################
# GLMppci
################################################################################

# Fit PPCI submodel
calc_PPCI <- function(N, train.data) {
    ##Fit a second model for paid loss
    ##Payment per Notified Claim
    ###Repeat each element in N_rep by the number of entries in each accident period
    N_rep<-vlookup_N_i(N,train.data)
    PPCI=train.data[order(as.numeric(train.data$origin)),]$aggregate_claims/N_rep
    
    return (PPCI)
}

fit_param_PPCI<-function(model,N,newdata){
    N_rep<-vlookup_N_i(N,newdata)
    pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
    mu<-pred_ODP$fit*N_rep
    phi<-(pred_ODP$residual.scale)^2*(N_rep)^2
    return(list(mu=mu, phi=phi))
}

##Calculate ODP density from PPCI model
cal_dens_PPCI<-function(y, param){
  mu<-param$mu
  phi<-param$phi
  dens<-c()
  for (i in 1:length(y)){
    dens[i]<-dODP(y=y[i],lambda=mu[i],phi=phi[i])
  }
  return(dens)
}

##Calculate ODP mean for PPCI model
cal_mu_PPCI<-function(param){
  mu<-param$mu
  return(mu)
}

##Calculate ODP mean for PPCI model
cal_sigma_PPCI<-function(param){
  mu<-param$mu
  phi<-param$phi
  return(mu*phi)
}

cal_CDF_PPCI<-function(y, param, new_y = F){
  mu<-param$mu
  phi<-param$phi
  if (!new_y) {return (sapply(1:length(y), function(i) pODP(y=y[i],lambda=mu[i],phi=phi[i])))}
  else {return(mapply(
    function(mu, phi)
      sapply(y, function(x) pODP(y=x,lambda=mu,phi=phi)), mu, phi)
  )}
}

simulate_PPCI<-function(param){
  k=1
  simy<-replicate(k, rtweedie(length(param$mu),xi=1,mu=param$mu,phi=param$phi))
  return(simy)
}

################################################################################
# GLMppcf
################################################################################

fit_param_PPCF<-function(model_subCount,model_subPayments,N, data,newdata){
    pred_F<-predict(model_subCount,newdata=newdata,type="response")
    pred_OT <- calc_pred_OT(model_subCount, N, data, newdata)
    
    pred_payment<-predict(model_subPayments,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
    mu<-pred_payment$fit*pred_F
    phi<-(pred_payment$residual.scale)^2*(pred_F)^2
    return(list(mu=mu, phi=phi))
}

# Fit PPCF submodel
calc_PPCF <- function(N, data) {
    PPCF=data$aggregate_claims/data$settle_count
    #Create a column for operation time
    ##Fit a second model for paid loss
    ##Payment per Notified Claim
    ###Repeat each element in N_rep by the number of entries in each accident period
    N_rep<-vlookup_N_i(N,data)
    OT=data$cum_settle_count/N_rep
    #Remove the cells with zero finalized claim count but with positive payments; 
    PPCF_data <- data.frame(list(PPCF = PPCF, OT=OT))
    PPCF_data<-na.omit(PPCF_data[PPCF_data$PPCF!=Inf,])
    
    return (PPCF_data)
}

calc_pred_OT <- function(model_subCount, N, data, newdata) {
    N_rep<-vlookup_N_i(N,newdata)
    # Predicted Finalised Claims for Each dataset
    # Assumes model only uses dev
    pred_F<-predict(model_subCount,newdata=data.frame(dev=1:tri.size),type="response")
    
    cum_F_count_train <-  matrix(0, nrow = tri.size, ncol = tri.size)
    for (i in 1:tri.size) {
      for (j in 1:tri.size) {
        # If cumulative settled count does not exist in train dataset
        # Add the most recent cumulative settled count from previous dev year and the predicted settle count for the year
        # Assumes that there is always a value for dev year 1 
        cum_F_count_train[i, j] <- ifelse(nrow(data[(data$origin==i)&(data$dev==j), ]) == 0, 
                                          cum_F_count_train[i, j-1] + pred_F[j], 
                                          data[(data$origin==i)&(data$dev==j), 'cum_settle_count'])
      }
    }
    
    pred_F_cum <- unlist(mapply(function(i, j) cum_F_count_train[i, j], 
                                as.numeric(as.character(newdata$origin)),
                                as.numeric(as.character(newdata$dev))
    ))
    
    pred_OT<-pred_F_cum/N_rep
    
    return(pred_OT)
}

##Calculate ODP density from PPCF GLM
cal_dens_PPCF<-function(y, param){
    mu<-param$mu
    phi<-param$phi

    dens_PPCF<-c()
    for (i in 1:length(y)){
      dens_PPCF[i]<-dODP(y=y[i],lambda=mu[i],phi=phi[i])
    }
    dens_PPCF[dens_PPCF==0]=min(dens_PPCF[dens_PPCF!=0])
    #assign a small number to the zero density in order to prevent -Inf Log Score
    return(dens_PPCF)
}

##Calculate ODP mean for PPCF model
cal_mu_PPCF<-function(param){
    mu<-param$mu
    return(mu)
}

##Calculate ODP mean for PPCF model
cal_sigma_PPCF<-function(param){
    mu<-param$mu
    phi<-param$phi
    return(phi*mu)
}

cal_CDF_PPCF<-function(y, param, new_y = F){
    mu<-param$mu
    phi<-mean(param$phi)
    
    if (!new_y) {
      CDF_PPCF <- pODP(y=y,lambda=mu,phi=phi)
    }
    else {
      CDF_PPCF <- mapply(
        function(mu, phi)
        pODP(y=y,lambda=mu,phi=phi, new_y=new_y), mu, phi)
    }
    return (ifelse(CDF_PPCF == 0, min(CDF_PPCF[CDF_PPCF!=0]), CDF_PPCF))
}

simulate_PPCF<-function(param){
    k=1
    simy<-replicate(k, rtweedie(length(param$mu), xi=1,mu=param$mu,phi=param$phi))
    return(simy)
}

################################################################################
# GAMLSS Gamma
################################################################################

fit_param_GAGamlss<-function(model,data,newdata,tau){
    pred_mu<-predict(model,what="mu",data=data,newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",data=data,newdata=newdata,type="response")
    return(list(pred_mu=pred_mu, pred_sigma=pred_sigma, tau=tau))
}

##Calculate Gamma density for Gamma GAMLSS
cal_dens_GA_Gamlss<-function(y, param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    return(dGA(x=y+tau,mu=pred_mu,sigma=pred_sigma))
}

##Calculate GA mean for GA GAMLSS
cal_mu_GA_Gamlss<-function(param){
    pred_mu<-param$pred_mu
    return(pred_mu)
}

##Calculate GA mean for GA GAMLSS
cal_sigma_GA_Gamlss<-function(param){
    pred_sigma<-param$pred_sigma
    return(pred_sigma)
}

##Calculate Gamma CDF for Gamma GAMLSS
cal_CDF_GA_Gamlss<-function(y, param, new_y = F){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    if (!new_y) {return(pGA(q=y+tau,mu=pred_mu,sigma=pred_sigma))}
    else {return(mapply(
        function(mu, sigma)
            pGA(q=y+tau,mu=mu,sigma=sigma), pred_mu, pred_sigma)
    )}
}

simulate_GAGamlss<-function(param){
    k=1
    simy<-replicate(k, rGA(n=length(param$pred_mu),mu=param$pred_mu,sigma=param$pred_sigma)-param$tau)
    simy[simy<0]<-0
    return(simy)
}

################################################################################
# GAMLSS Lognormal
################################################################################

fit_param_LNGamlss<-function(model,data,newdata,tau){
    pred_mu<-predict(model,what="mu",data=data,newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",data=data,newdata=newdata,type="response")
    return(list(pred_mu=pred_mu, pred_sigma=pred_sigma, tau=tau))
}

##Calculate Log-Normal density for Log-Normal GAMLSS
cal_dens_LN_Gamlss<-function(y, param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    return(dLNO(x=y+tau,mu=pred_mu,sigma=pred_sigma))
}

##Calculate LN mean for LN GAMLSS
cal_mu_LN_Gamlss<-function(param){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    mean<-exp(pred_mu+1/2*pred_sigma^2)
    return(mean)
}

##Calculate LN mean for LN GAMLSS
cal_sigma_LN_Gamlss<-function(param){
    pred_sigma<-param$pred_sigma
    return(pred_sigma)
}

##Calculate Log-Normal CDF for Log-Normal GAMLSS
cal_CDF_LN_Gamlss<-function(y, param, new_y = F){
    pred_mu<-param$pred_mu
    pred_sigma<-param$pred_sigma
    tau<-param$tau
    if (!new_y) {return(pLNO(q=y+tau,mu=pred_mu,sigma=pred_sigma))}
    else {return(mapply(
        function(mu, sigma)
            pLNO(q=y+tau,mu=mu,sigma=sigma), pred_mu, pred_sigma)
    )}
}

simulate_LNGamlss<-function(param){
    k=1
    simy<-replicate(k, rLNO(n=length(param$pred_mu),mu=param$pred_mu,sigma=param$pred_sigma)-param$tau)
    simy[simy<0]<-0
    return(simy)
}