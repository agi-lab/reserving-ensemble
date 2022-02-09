#Create functions that can create simulated results:

#return a matrix with as many rows as newdata and k columns, where k is the number of simulations
simulate_ODP<-function(model,newdata,k){
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu<-pred_ODP$fit
  pred_ODP_phi<-(pred_ODP$residual.scale)^2
  simy<-replicate(k, rtweedie(length(pred_ODP_mu), xi=1,mu=pred_ODP_mu,phi=pred_ODP_phi))
  return(simy)
  # replicate the simulations process k times; the ODP distribution corresponding to tweedie distribution with power=1
}


simulate_GA<-function(model,newdata,k,tau){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  simy<-replicate(k, rGA(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma)-tau)
  ###Assign the negative value to be zero
  simy[simy<0]<-0
  return(simy)
  # replicate the simulations process k times
}

simulate_LN<-function(model,newdata,k,tau){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  simy<-replicate(k, rLNO(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma)-tau)
  simy[simy<0]<-0
  return(simy)
}

simulate_ZAGA<-function(model,newdata,k){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  pred_nu<-mean(predict(model,what="nu",newdata=newdata_nu,type="response"))
  simy<-replicate(k, rZAGA(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma,nu=pred_nu))
  return(simy) 
}
rZALN<-Zadj.r(family="LOGNO")
simulate_ZALN<-function(model,newdata,k){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu<-predict(model,parameter="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,parameter="sigma",newdata=newdata,type="response")
  pred_nu<-mean(predict(model,parameter="xi0",newdata=newdata_nu,type="response"))
  simy<-replicate(k,rZALN(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))
  return(simy)
}

simulate_NO<-function(model,newdata,k){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  simy<-replicate(k,rNO(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma))
  simy[simy<0]<-0
  return(simy)
}

simulate_GAGamlss<-function(model,newdata_mu,newdata_sigma,k,tau){
  pred_mu<-predict(model,what="mu",newdata=newdata_mu,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  simy<-replicate(k, rGA(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma)-tau)
  simy[simy<0]<-0
  return(simy)
}

simulate_LNGamlss<-function(model,newdata_mu,newdata_sigma,k,tau){
  pred_mu<-predict(model,what="mu",newdata=newdata_mu,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  simy<-replicate(k, rLNO(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma)-tau)
  simy[simy<0]<-0
  return(simy)
}

simulate_PPCI<-function(model,N,newdata,index_start,index_end,k){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  mu<-pred_ODP$fit*N_rep_valid
  phi<-(pred_ODP$residual.scale)^2*(N_rep_valid)^2
  simy<-replicate(k, rtweedie(length(mu),xi=1,mu=mu,phi=phi))
  return(simy)
}


simulate_PPCF<-function(model_subCount,model_subPayments,train_cumF,newdataFC,newdata_Pay,N,index_start,index_end,k){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
  pred_F<-predict(model_subCount,newdata=newdataFC,type="response")
  pred_F_cum<-list()
  for(i in index_start:index_end){
    pf<-predict(model_subCount,newdata=newdataFC[newdataFC$origin==i,],type="response")
    pf_cum<-round(cumsum(c(train_cumF[train_cumF$origin==i,]$value[nrow(train_cumF[train_cumF$origin==i,])],pf))[-1],0)
    pred_F_cum[[i]]<-pf_cum
  }
  pred_F_cum<-as.vector(unlist(pred_F_cum[index_start:index_end]))
  
  pred_OT<-pred_F_cum/N_rep_valid
  pred_payment<-predict.glm(model_subPayments,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
  mu<-pred_payment$fit*pred_F
  phi<-(pred_payment$residual.scale)^2*(pred_F)^2
  simy<-replicate(k, rtweedie(length(mu), xi=1,mu=mu,phi=phi))
  return(simy)
}

reserve_par1_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)
reserve_par2_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)
reserve_par3_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)
reserve_par4_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)
reserve_par5_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)
reserve_par6_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)
reserve_par7_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)
reserve_par8_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)


##Equally Weighted Ensemble and BMV
reserve_BMV_newSimulScheme<-matrix(NA,nrow=1000,ncol=ntri)
reserve_par1_newEW<-matrix(NA,nrow=1000,ncol=ntri)

for (D in 1:ntri){
  
  set.seed(20200130+D)
  U<-runif(1000)
  # Package-wise global parameters
  set_parameters(ref_claim = 200000, time_unit = 1/4)
  ref_claim = return_parameters()[1]
  time_unit = return_parameters()[2]
  
  # Claim occurrence
  # Example 1.1: Claim occurrence: Constant exposure and frequency
  
  # Input parameters
  years = 10
  I = years / time_unit
  E = c(rep(12000, I)) # effective annual exposure rates
  lambda = c(rep(0.03, I))
  
  # Number of claims occurring for each period i
  n_vector = claim_frequency(I = I, E = E, freq = lambda)
  
  # Occurrence time of each claim r, for each period i
  occurrence_times = claim_occurrence(frequency_vector = n_vector)
  
  ##############Claim size: Default power normal
  # S^0.2 ~ N(9.5, 3), left truncated at 30
  
  S_df <- function(s) {
    # truncate and rescale
    if (s < 1) {
      return(0)
    } else {
      p_trun <- pnorm(s^0.2, 9.5, 3) - pnorm(1^0.2, 9.5, 3)
      p_rescaled <- p_trun/(1 - pnorm(1^0.2, 9.5, 3))
      return(p_rescaled)
    }
  }
  
  claim_sizes <- claim_size(frequency_vector = n_vector, simfun = S_df, type = "p", range = c(0, 1e24))
  
  ###Notification Delay
  notidel_param <- function(claim_size, occurrence_period){
    mean <- 2.47; cv <- 1.54
    shape <- get_Weibull_parameters(mean, cv)[1, ]
    scale <- get_Weibull_parameters(mean, cv)[2, ]
    c(shape = shape, scale = scale)
  }
  
  ## output
  notidel <- claim_notification(n_vector, claim_size(n_vector),
                                paramfun = notidel_param)
  #By changing claim notification to transformed Gamma, the first development periods contain no zero value claims 
  
  
  # Claim closure: Default Weibull
  # Claim closure: Default Weibull
  setldel_param <- function(claim_size, occurrence_period) {
    target_mean <- 11.74
    
    # specify the target Weibull coefficient of variation
    target_cv <- 0.61
    
    c(shape = get_Weibull_parameters(target_mean, target_cv)[1, ],
      scale = get_Weibull_parameters(target_mean, target_cv)[2, ])
  }
  
  ## output
  # simulate the settlement delays from the Weibull with parameters above
  setldel <- claim_closure(n_vector, claim_sizes, rfun = rweibull, paramfun = setldel_param)
  
  # Claim partial payment: Default mixture distribution
  
  no_payments <- claim_payment_no(n_vector, claim_sizes)
  
  # Claim Partial payment (Without Inflation): Default Distribution
  
  payment_sizes <- claim_payment_size(n_vector, claim_sizes,no_payments)
  
  # Claim payment time: Default Distribution
  
  ## input
  r_pmtdel <- function(n, claim_size, setldel, setldel_mean) {
    result <- c(rep(NA, n))
    
    # First simulate the unnormalised values of d, sampled from a Weibull distribution
    if (n >= 4) {
      # 1) Simulate the last payment delay
      unnorm_d_mean <- (1 / 4) / time_unit
      unnorm_d_cv <- 0.20
      parameters <- get_Weibull_parameters(target_mean = unnorm_d_mean, target_cv = unnorm_d_cv)
      result[n] <- stats::rweibull(1, shape = parameters[1], scale = parameters[2])
      
      # 2) Simulate all the other payment delays
      for (i in 1:(n - 1)) {
        unnorm_d_mean <- setldel_mean / n
        unnorm_d_cv <- 0.35
        parameters <- get_Weibull_parameters(target_mean = unnorm_d_mean, target_cv = unnorm_d_cv)
        result[i] <- stats::rweibull(1, shape = parameters[1], scale = parameters[2])
      }
      
    } else {
      for (i in 1:n) {
        unnorm_d_mean <- setldel_mean / n
        unnorm_d_cv <- 0.35
        parameters <- get_Weibull_parameters(target_mean = unnorm_d_mean, target_cv = unnorm_d_cv)
        result[i] <- stats::rweibull(1, shape = parameters[1], scale = parameters[2])
      }
    }
    
    # Normalise d such that sum(inter-partial delays) = settlement delay
    # To make sure that the pmtdels add up exactly to setldel, we treat the last one separately
    result[1:n-1] <- (setldel/sum(result)) * result[1:n-1]
    result[n] <- setldel - sum(result[1:n-1])
    
    return(result)
  }
  
  param_pmtdel <- function(claim_size, setldel, occurrence_period) {
    # mean settlement delay
    if (claim_size < (0.10 * ref_claim) & occurrence_period >= 21) {
      a <- min(0.85, 0.65 + 0.02 * (occurrence_period - 21))
    } else {
      a <- max(0.85, 1 - 0.0075 * occurrence_period)
    }
    mean_quarter <- a * min(25, max(1, 6 + 4*log(claim_size/(0.10 * ref_claim))))
    target_mean <- mean_quarter / 4 / time_unit
    
    c(claim_size = claim_size,
      setldel = setldel,
      setldel_mean = target_mean)
  }
  
  ## output
  payment_delays <- claim_payment_delay(
    n_vector, claim_sizes, no_payments, setldel,
    rfun = r_pmtdel, paramfun = param_pmtdel,
    occurrence_period = rep(1:I, times = n_vector))
  
  # payment times on a continuous time scale
  payment_times <- claim_payment_time(n_vector, occurrence_times, notidel, payment_delays)
  # payment times in periods
  payment_periods <- claim_payment_time(n_vector, occurrence_times, notidel, payment_delays,discrete = TRUE)
  
  
  # Claim inflation
  # Base inflation: a vector of quarterly rates
  # In this data set we set base inflation to be at 2% p.a. constant for both past and future
  
  demo_rate <- (1 + 0.02)^(1/4) - 1
  base_inflation_past <- rep(demo_rate, times = 40)
  base_inflation_future <- rep(demo_rate, times = 40)
  base_inflation_vector <- c(base_inflation_past, base_inflation_future)
  
  # Remove the Superimposed inflation
  
  SI_occurrence <- function(occurrence_time, claim_size){1}
  # 2) With respect to payment "time" (continuous scale)
  # -> compounding by user-defined time unit
  SI_payment <-  function(payment_time, claim_size){1}
  
  payment_inflated <- claim_payment_inflation(
    n_vector,
    payment_sizes,
    payment_times,
    occurrence_times,
    claim_sizes,
    base_inflation_vector,
    SI_occurrence,
    SI_payment
  )
  agg_sqrt_withTail<-claim_output(n_vector, payment_times, payment_inflated, incremental = TRUE,future=TRUE,adjust=FALSE)
  agg_sqrt<-agg_sqrt_withTail[,-ncol(agg_sqrt_withTail)]
  colnames(agg_sqrt)<-c(1:40)
  rownames(agg_sqrt)<-c(1:40)
  #Generate a Claim Triangle
  agg_tri_withInflation <- claim_output(n_vector, payment_times, payment_inflated,incremental = TRUE)
  colnames(agg_tri_withInflation)<-c(1:40)
  rownames(agg_tri_withInflation)<-c(1:40)
  full_dat<-as.data.frame(as.triangle(agg_sqrt))
  full_dat$value=full_dat$value/10000
  full_dat$Calendar=as.numeric(full_dat$origin)+as.numeric(full_dat$dev)
  full_dat$origin=as.factor(full_dat$origin)
  full_dat$dev=as.factor(full_dat$dev)
  in_sample<-full_dat[full_dat$Calendar<=41,]
  out_sample<-full_dat[full_dat$Calendar>41,]
  
  tri_withInflation<-as.triangle(agg_tri_withInflation)
  
  # convert the aggregate data generated by SynthETIC into triangle format: value, development year and accident year
  dt_withInflation<-as.data.frame(tri_withInflation)
  # Convert the triangle into dataframe format
  dt_withInflation_past<-as.data.frame(tri_withInflation,na.rm=TRUE)
  # only contains the past observations (the upper triangle) as the lower triangle are coded as N.A
  dt_withInflation_past$Calendar=as.character(as.numeric(dt_withInflation_past$origin)+as.numeric(dt_withInflation_past$dev))
  #add calendar periods in the data frame
  
  #Convert the square into data frame
  dt_sqrt<-as.data.frame(as.triangle(agg_sqrt))
  
  #Rescale the Claims Value in the unit of 0000's
  dt_withInflation_past$value=dt_withInflation_past$value/10000
  tau_LN<-5
  tau_Ga<-5
  ####PPCI
  all_claims <- claims(
    frequency_vector = n_vector,
    occurrence_list = occurrence_times,
    claim_size_list = claim_sizes,
    notification_list = notidel,
    settlement_list = setldel,
    no_payments_list = no_payments,
    payment_size_list = payment_sizes,
    payment_delay_list = payment_delays,
    payment_time_list = payment_times,
    payment_inflated_list = payment_inflated
  )
  transaction_dataset <- generate_transaction_dataset(all_claims,adjust=TRUE)
  transaction_dataset$ReportingPeriod=ceiling(transaction_dataset$notidel)
  transaction_dataset$Calendar=transaction_dataset$ReportingPeriod+transaction_dataset$occurrence_period
  matri_NotiCount<-matrix(NA,nrow=40,ncol=40)
  
  for (i in 1:40){
    for (j in 1:40){
      if (i+j<=41){matri_NotiCount[i,j]=nrow(transaction_dataset[transaction_dataset$occurrence_period==i&transaction_dataset$ReportingPeriod==j&transaction_dataset$Calendar<=41,])}
      else{matri_NotiCount[i,j]='NA'}
    }
  }
  
  #Using the Chain-Ladder package to generate a triangular object:
  NotiCount_tri<-as.triangle(matri_NotiCount)
  #Convert it to a dataframe format
  NotiCount_dt<-as.data.frame(NotiCount_tri,na.rm=TRUE)
  NotiCount_dt$Calendar=as.numeric(NotiCount_dt$origin)+as.numeric(NotiCount_dt$dev)
  #Construct a full claim count square
  matri_NotiCount_sqrt<-matrix(NA,nrow=40,ncol=40)
  
  for (i in 1:40){
    for (j in 1:40){
      matri_NotiCount_sqrt[i,j]=nrow(transaction_dataset[transaction_dataset$occurrence_period==i&transaction_dataset$ReportingPeriod==j,])
    }
  }
  NotiCount_full_dt<-as.data.frame(as.triangle(matri_NotiCount_sqrt))
  NotiCount_full_dt$Calendar=as.numeric(NotiCount_full_dt$origin)+as.numeric(NotiCount_full_dt$dev)
  insample_NC<-NotiCount_full_dt[NotiCount_full_dt$Calendar<=41,]
  test_NC<-NotiCount_full_dt[NotiCount_full_dt$Calendar>41,]
  
  ###Fitting of all the models in the in-sample 
  #####Models with Basic Structures
  #Fitting of ODP GLM
  ODP_GLM_in<-glm(formula=value~factor(origin)+factor(dev),family=quasipoisson(link="log"),data=in_sample)
  #Fitting of ZAGA
  ZAGA_in<-gamlss(formula=value~factor(origin)+factor(dev),nu.formula=~as.numeric(as.character(dev)),data=in_sample,family=ZAGA(mu.link="log",sigma.link = "log", nu.link = "logit"),trace=FALSE)
  #Fitting of ZALN
  LN_in<-gamlssZadj(y=value,mu.formula =~factor(origin)+factor(dev),xi0.formula=~as.numeric(as.character(dev)),data=in_sample,family=LOGNO(mu.link="identity",sigma.link="log"))
  #Fitting of opt_LN:
  tau_LN<-5
  LN_optimTau<-gamlss(formula=(value+tau_LN)~factor(origin)+factor(dev),data=in_sample,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  tau_Ga<-5
  Ga_optimTau<-gamlss(formula=(value+tau_Ga)~factor(origin)+factor(dev),data=in_sample,family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
  
  ###Hoerl Curve
  #Under ODP Assumption
  in_sample_Ho<-in_sample
  in_sample_Ho$dev<-as.numeric(as.character(in_sample$dev))
  glm_ODP_Ho_In<-glm(formula=value~factor(origin)+log(as.numeric(as.character(dev)))+as.numeric(as.character(dev)),family=quasipoisson(link="log"),data=in_sample_Ho,trace=FALSE)
  
  #Density for Hoerl curve
  out_sample_Ho<-out_sample
  out_sample_Ho$dev<-as.numeric(as.character(out_sample$dev))
  
  
  #Under Gamma Assumption
  
  glm_Ga_Ho_In<-gamlss(formula=(value+tau_Ga)~factor(origin)+log(dev)+dev,data=in_sample_Ho,family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
  
  #Under Log-Normal Assumption
  glm_LN_Ho_In<-gamlss(formula=(value+tau_LN)~factor(origin)+log(dev)+dev,data=in_sample_Ho,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  
  ###GLM with Calendar Periods Effects
  
  #Under ODP Assumption
  
  glm_ODP_Cal_In<-glm(formula=value~factor(dev)+Calendar,family=quasipoisson(link="log"),data=in_sample)
  
  #Under Gamma Assumption
  
  glm_Ga_Cal_In<-gamlss(formula=(value+tau_Ga)~factor(dev)+Calendar,data=in_sample,family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
  
  #Under Log-Normal Assumption
  glm_LN_Cal_In<-gamlss(formula=(value+tau_LN)~factor(dev)+Calendar,data=in_sample,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  in_sample_numeric<-in_sample
  in_sample_numeric$origin=as.numeric(in_sample$origin)
  in_sample_numeric$dev=as.numeric(in_sample$dev)
  in_sample_numeric$Calendar=as.numeric(in_sample$Calendar)
  
  ###Smoothing Spline
  sp_Normal_In<-gamlss(formula=value~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=in_sample_numeric,family=NO(),trace=FALSE)
  sp_Gamma_In<-gamlss(formula=(value+tau_Ga)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=in_sample_numeric,family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
  sp_LN_In<-gamlss(formula=(value+tau_LN)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=in_sample_numeric,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  
  ###GAMLSS
  gamlss2_GA_In<-gamlss(formula=(value+tau_Ga)~factor(origin)+factor(dev),data=in_sample,sigma.formula=~cs(as.numeric(as.character(dev))),family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
  gamlss2_LN_In<-gamlss(formula=(value+tau_LN)~factor(origin)+factor(dev),data=in_sample,sigma.formula=~cs(as.numeric(as.character(dev))),family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  
  
  
  
  
  
  
  
  
  
  
  ###Predictive density for out-sample data
  out_sample<-out_sample[order(as.numeric(as.character(out_sample$origin))),]
  
  
  ###simulated Results for Out-Sample Data
  out_sample<-out_sample[order(as.numeric(as.character(out_sample$origin))),]
  
  out_sample_numeric<-out_sample
  out_sample_numeric$origin=as.numeric(out_sample$origin)
  out_sample_numeric$dev=as.numeric(out_sample$dev)
  out_sample_numeric$Calendar=as.numeric(out_sample$Calendar)
  
  #Simulations for Hoerl curve
  out_sample_Ho<-out_sample
  out_sample_Ho$dev<-as.numeric(as.character(out_sample$dev))
  
  
  
  #quantile for GAMLSS
  out_sample_mu<-out_sample
  out_sample_sigma<-out_sample_numeric
  
  #Fit a PPCI model
  
  ##Fit a ODP/Chain-Ladder model to predict total claim notification count at each accident period 
  fit_nc_In<-glm(value~factor(origin)+factor(dev),data=insample_NC,family=quasipoisson(link="log"))
  ##Obtain Estimated notification claims count for each accident period: Nk
  N_in<-c()
  N_in[1]<-sum(insample_NC[insample_NC$origin==1,]$value)
  for (i in 2:40){
    N_in[i]<-sum(insample_NC[insample_NC$origin==i,]$value)+sum(round(predict(fit_nc_In,newdata=test_NC[test_NC$origin==i,],type="response"),0))
  }
  
  ##Fit a second model for paid loss
  
  ##Payment per Notified Claim
  ###Create a function that counts the number of entries in each accident period
  count_by_origin<-function(data,start,end){
    cbo<-c()
    for (i in start:end){
      cbo[i]<-nrow(data[data$origin==i,]) 
    }
    return(cbo)
  }
  ###Sort the data by accident period
  in_sample<-in_sample[order(as.numeric(as.character(in_sample$origin))),]
  ###Repete each element in N_rep by the number of entries in each accident period
  N_rep_In<-rep(N_in,times=count_by_origin(in_sample,1,40))
  in_sample$PPCI=in_sample$value/N_rep_In
  
  ##Fitting a ODP model on PPCI
  fit_ODP_ppci_In<-glm(PPCI~factor(dev),family=quasipoisson(link="log"),data=in_sample)
  
  ###For PPCF model
  transaction_dataset$SettlementPeriod=ceiling(transaction_dataset$setldel)
  transaction_dataset$ReportingPeriod=ceiling(transaction_dataset$notidel)
  transaction_dataset$Calendar=transaction_dataset$ReportingPeriod+transaction_dataset$occurrence_period
  matri_SettlementCount<-matrix(NA,nrow=40,ncol=40)
  transaction_dataset$Calendar=as.numeric(transaction_dataset$occurrence_period)+as.numeric(transaction_dataset$SettlementPeriod)
  for (i in 1:40){
    for (j in 1:40){
      if (i+j<=41){matri_SettlementCount[i,j]=nrow(transaction_dataset[transaction_dataset$occurrence_period==i&transaction_dataset$SettlementPeriod==j&transaction_dataset$Calendar<=41,])}
      else{matri_SettlementCount[i,j]='NA'}
    }
  }
  
  sqrt_SettlementCount<-matrix(NA,nrow=40,ncol=40)
  
  for (i in 1:40){
    for (j in 1:40){
      sqrt_SettlementCount[i,j]=nrow(transaction_dataset[transaction_dataset$occurrence_period==i&transaction_dataset$SettlementPeriod==j,])
    }
  }
  
  #Using the Chain-Ladder package to generate a triangular object:
  SettlCount_tri<-as.triangle(matri_SettlementCount)
  SettlCount_sqrt<-as.triangle(sqrt_SettlementCount)
  
  #Convert it into a dataframe
  SettlCount_dt<-as.data.frame(SettlCount_tri,na.rm=TRUE)
  SettlCount_dt$Calendar=as.numeric(SettlCount_dt$origin)+as.numeric(SettlCount_dt$dev)
  SettlCount_in<-SettlCount_dt
  
  SettlCount_full<-as.data.frame(SettlCount_sqrt)
  SettlCount_full$Calendar=as.numeric(SettlCount_full$origin)+as.numeric(SettlCount_full$dev)
  SettlCount_out<-SettlCount_full[SettlCount_full$Calendar>41,]
  
  
  #Create a triangle for Full Cummulative Claims Count
  cum_Count_tri_full<-incr2cum(as.triangle(matri_NotiCount_sqrt))
  cum_Count_dt_full<-as.data.frame(cum_Count_tri_full,na.rm=TRUE)
  cum_Count_dt_full<-cum_Count_dt_full[order(as.numeric(cum_Count_dt_full$origin)),]
  cum_Count_dt_full$Calendar=as.numeric(as.character(cum_Count_dt_full$origin))+as.numeric(cum_Count_dt_full$dev)
  
  #Create a triangle for In-Sample Cummulative finalised claims count
  cum_F_tri<-incr2cum(SettlCount_tri)
  cum_F_dt<-as.data.frame(cum_F_tri,na.rm=TRUE)
  cum_F_dt<-cum_F_dt[order(as.numeric(cum_F_dt$origin)),]
  cum_F_dt$Calendar=as.numeric(as.character(cum_F_dt$origin))+as.numeric(cum_F_dt$dev)
  
  #Create a triangle for Out-Sample Cummulative finalised claims count
  cum_F_tri_full<-incr2cum(SettlCount_sqrt)
  cum_F_dt_full<-as.data.frame(cum_F_tri_full)
  cum_F_dt_full$Calendar=as.numeric(as.character(cum_F_dt_full$origin))+as.numeric(as.character(cum_F_dt_full$dev))
  cum_F_dt_full<-cum_F_dt_full[order(as.numeric(as.character(cum_F_dt_full$origin))),]
  cum_F_dt_out<-cum_F_dt_full[cum_F_dt_full$Calendar>41,]
  
  #Create a triangle for Unfinalized claim count triangle (full)
  cum_U_dt_full<-data.frame(origin=cum_F_dt_full$origin,dev=cum_F_dt_full$dev,Calendar=cum_F_dt_full$Calendar,value=cum_Count_dt_full$value-cum_F_dt_full$value)
  cum_U_dt_full<-cum_U_dt_full[order(as.numeric(as.character(cum_U_dt_full$origin))),]
  cum_U_dt_In<-cum_U_dt_full[cum_U_dt_full$Calendar<=41,]
  cum_U_dt_Out<-cum_U_dt_full[cum_U_dt_full$Calendar>41,]
  
  #Generate a full data frame containing finalized claim count, unfinalized claim count and notification claim count
  
  NotiCount_full_dt<-NotiCount_full_dt[order(as.numeric(NotiCount_full_dt$origin)),]
  UplusN_full<-list()
  SettlCount_full<-SettlCount_full[order(as.numeric(as.character(SettlCount_full$origin))),]
  
  
  for (i in 1:40){
    UplusN_full[[i]]<-cum_U_dt_full[cum_U_dt_full$origin==i,]$value[1:(nrow(cum_U_dt_full[cum_U_dt_full$origin==i,])-1)]+NotiCount_full_dt[NotiCount_full_dt$origin==i,]$value[2:nrow(NotiCount_full_dt[NotiCount_full_dt$origin==i,])]
  }
  
  UplusN_full<-as.vector(unlist(UplusN_full))
  FCount_full<-SettlCount_full[SettlCount_full$dev!=1,]$value
  
  FUN_dat_full<-data.frame(origin=SettlCount_full[SettlCount_full$dev!=1,]$origin,dev=SettlCount_full[SettlCount_full$dev!=1,]$dev,UplusN=UplusN_full,FCount=FCount_full)
  
  FUN_dat_full$Calendar=as.numeric(as.character(FUN_dat_full$origin))+as.numeric(FUN_dat_full$dev)
  FUN_dat_full<-FUN_dat_full[order(as.numeric(as.character(FUN_dat_full$origin))),]
  addrow<-data.frame(origin=as.factor(40),dev=1,UplusN=UplusN_full[length(UplusN_full)],FCount=SettlCount_dt[SettlCount_dt$origin==40,]$value,Calendar=41)
  FUN_dat_full<-rbind(FUN_dat_full,addrow)
  
  FUN_dat_In<-FUN_dat_full[FUN_dat_full$Calendar<=41,]
  FUN_dat_Out<-FUN_dat_full[FUN_dat_full$Calendar>41,]
  
  in_sample$PPCF=in_sample$value/SettlCount_in$value
  odp_FC_In<-glm(FCount~factor(dev),data=FUN_dat_In,family=quasipoisson(link="log"))
  #Create a column for operation time
  in_sample$OT=cum_F_dt$value/N_rep_In
  #Remove the cells with zero finalized claim count but with positive payments; 
  new_in_PPCF<-na.omit(in_sample[in_sample$PPCF!=Inf,])
  ODP_PPCF_In<-glm(PPCF~OT,family=quasipoisson(link="log"),data=new_in_PPCF)
  
  ####Simulate Quantile for Ensemble 1
  out_sample_numeric<-out_sample
  out_sample_numeric$origin=as.numeric(out_sample$origin)
  out_sample_numeric$dev=as.numeric(out_sample$dev)
  out_sample_numeric$Calendar=as.numeric(out_sample$Calendar)
  
  out_sample_mu<-out_sample
  out_sample_sigma<-out_sample_numeric
  
  out_sample_Ho<-out_sample
  out_sample_Ho$dev<-as.numeric(as.character(out_sample$dev))
  
  
  mix_ind_subset1<-findInterval(U,cumsum(unlist(model_weights_simul_par1_new[[D]][1])))+1
  OW_par1_new_Var1<-matrix(NA,nrow=nrow(out_sample[ind_subset1_par1_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_subset1[i]==1){OW_par1_new_Var1[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset1_par1_new,],k=1)}
    else if(mix_ind_subset1[i]==2){OW_par1_new_Var1[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset1_par1_new,],k=1)}
    else if(mix_ind_subset1[i]==3){OW_par1_new_Var1[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset1_par1_new,],k=1)}
    else if(mix_ind_subset1[i]==4){OW_par1_new_Var1[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset1_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset1[i]==5){OW_par1_new_Var1[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset1_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset1[i]==6){OW_par1_new_Var1[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset1_par1_new,],k=1)}
    else if(mix_ind_subset1[i]==7){OW_par1_new_Var1[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset1_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset1[i]==8){OW_par1_new_Var1[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset1_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset1[i]==9){OW_par1_new_Var1[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset1_par1_new,],k=1)}
    else if(mix_ind_subset1[i]==10){OW_par1_new_Var1[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset1_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset1[i]==11){OW_par1_new_Var1[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset1_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset1[i]==12){OW_par1_new_Var1[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset1_par1_new,],k=1)}
    else if(mix_ind_subset1[i]==13){OW_par1_new_Var1[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset1_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset1[i]==14){OW_par1_new_Var1[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset1_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset1[i]==15){OW_par1_new_Var1[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset1_par1_new,],newdata_sigma=out_sample_sigma[ind_subset1_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset1[i]==16){OW_par1_new_Var1[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset1_par1_new,],newdata_sigma =out_sample_sigma[ind_subset1_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset1[i]==17){OW_par1_new_Var1[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset1_par1_new,],2,40,k=1)}
    else{OW_par1_new_Var1[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset1_par1_new,],out_sample[ind_subset1_par1_new,],N_in,2,40,k=1)}
  }
  
  
  mix_ind_subset2<-findInterval(U,cumsum(unlist(model_weights_simul_par1_new[[D]][2])))+1
  OW_par1_new_Var2<-matrix(NA,nrow=nrow(out_sample[ind_subset2_par1_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_subset2[i]==1){OW_par1_new_Var2[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset2_par1_new,],k=1)}
    else if(mix_ind_subset2[i]==2){OW_par1_new_Var2[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset2_par1_new,],k=1)}
    else if(mix_ind_subset2[i]==3){OW_par1_new_Var2[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset2_par1_new,],k=1)}
    else if(mix_ind_subset2[i]==4){OW_par1_new_Var2[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset2_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset2[i]==5){OW_par1_new_Var2[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset2_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset2[i]==6){OW_par1_new_Var2[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset2_par1_new,],k=1)}
    else if(mix_ind_subset2[i]==7){OW_par1_new_Var2[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset2_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset2[i]==8){OW_par1_new_Var2[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset2_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset2[i]==9){OW_par1_new_Var2[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset2_par1_new,],k=1)}
    else if(mix_ind_subset2[i]==10){OW_par1_new_Var2[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset2_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset2[i]==11){OW_par1_new_Var2[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset2_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset2[i]==12){OW_par1_new_Var2[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset2_par1_new,],k=1)}
    else if(mix_ind_subset2[i]==13){OW_par1_new_Var2[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset2_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset2[i]==14){OW_par1_new_Var2[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset2_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset2[i]==15){OW_par1_new_Var2[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset2_par1_new,],newdata_sigma=out_sample_sigma[ind_subset2_par1_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_subset2[i]==16){OW_par1_new_Var2[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset2_par1_new,],newdata_sigma =out_sample_sigma[ind_subset2_par1_new,],k=1,tau=tau_LN)}
    else if(mix_ind_subset2[i]==17){OW_par1_new_Var2[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset2_par1_new,],2,40,k=1)}
    else{OW_par1_new_Var2[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset2_par1_new,],out_sample[ind_subset2_par1_new,],N_in,2,40,k=1)}
  }
  OW_par1_new_Var<-rbind(OW_par1_new_Var1,OW_par1_new_Var2)
  
  
  ###Simulate Quantile for Ensemble 2
  
  mix_ind_par2_subset1<-findInterval(U,cumsum(unlist(model_weights_simul_par2_new[[D]][1])))+1
  OW_par2_new_Var1<-matrix(NA,nrow=nrow(out_sample[ind_subset1_par2_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par2_subset1[i]==1){OW_par2_new_Var1[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset1_par2_new,],k=1)}
    else if(mix_ind_par2_subset1[i]==2){OW_par2_new_Var1[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset1_par2_new,],k=1)}
    else if(mix_ind_par2_subset1[i]==3){OW_par2_new_Var1[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset1_par2_new,],k=1)}
    else if(mix_ind_par2_subset1[i]==4){OW_par2_new_Var1[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset1_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset1[i]==5){OW_par2_new_Var1[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset1_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset1[i]==6){OW_par2_new_Var1[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset1_par2_new,],k=1)}
    else if(mix_ind_par2_subset1[i]==7){OW_par2_new_Var1[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset1_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset1[i]==8){OW_par2_new_Var1[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset1_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset1[i]==9){OW_par2_new_Var1[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset1_par2_new,],k=1)}
    else if(mix_ind_par2_subset1[i]==10){OW_par2_new_Var1[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset1_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset1[i]==11){OW_par2_new_Var1[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset1_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset1[i]==12){OW_par2_new_Var1[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset1_par2_new,],k=1)}
    else if(mix_ind_par2_subset1[i]==13){OW_par2_new_Var1[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset1_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset1[i]==14){OW_par2_new_Var1[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset1_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset1[i]==15){OW_par2_new_Var1[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset1_par2_new,],newdata_sigma=out_sample_sigma[ind_subset1_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset1[i]==16){OW_par2_new_Var1[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset1_par2_new,],newdata_sigma =out_sample_sigma[ind_subset1_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset1[i]==17){OW_par2_new_Var1[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset1_par2_new,],2,40,k=1)}
    else{OW_par2_new_Var1[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset1_par2_new,],out_sample[ind_subset1_par2_new,],N_in,2,40,k=1)}
  }
  
  
  mix_ind_par2_subset2<-findInterval(U,cumsum(unlist(model_weights_simul_par2_new[[D]][2])))+1
  OW_par2_new_Var2<-matrix(NA,nrow=nrow(out_sample[ind_subset2_par2_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par2_subset2[i]==1){OW_par2_new_Var2[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset2_par2_new,],k=1)}
    else if(mix_ind_par2_subset2[i]==2){OW_par2_new_Var2[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset2_par2_new,],k=1)}
    else if(mix_ind_par2_subset2[i]==3){OW_par2_new_Var2[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset2_par2_new,],k=1)}
    else if(mix_ind_par2_subset2[i]==4){OW_par2_new_Var2[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset2_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset2[i]==5){OW_par2_new_Var2[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset2_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset2[i]==6){OW_par2_new_Var2[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset2_par2_new,],k=1)}
    else if(mix_ind_par2_subset2[i]==7){OW_par2_new_Var2[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset2_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset2[i]==8){OW_par2_new_Var2[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset2_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset2[i]==9){OW_par2_new_Var2[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset2_par2_new,],k=1)}
    else if(mix_ind_par2_subset2[i]==10){OW_par2_new_Var2[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset2_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset2[i]==11){OW_par2_new_Var2[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset2_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset2[i]==12){OW_par2_new_Var2[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset2_par2_new,],k=1)}
    else if(mix_ind_par2_subset2[i]==13){OW_par2_new_Var2[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset2_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset2[i]==14){OW_par2_new_Var2[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset2_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset2[i]==15){OW_par2_new_Var2[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset2_par2_new,],newdata_sigma=out_sample_sigma[ind_subset2_par2_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par2_subset2[i]==16){OW_par2_new_Var2[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset2_par2_new,],newdata_sigma =out_sample_sigma[ind_subset2_par2_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par2_subset2[i]==17){OW_par2_new_Var2[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset2_par2_new,],2,40,k=1)}
    else{OW_par2_new_Var2[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset2_par2_new,],out_sample[ind_subset2_par2_new,],N_in,2,40,k=1)}
  }
  OW_par2_new_Var<-rbind(OW_par2_new_Var1,OW_par2_new_Var2)
  
  
  ##Ensemble 3
  mix_ind_par3_subset1<-findInterval(U,cumsum(unlist(model_weights_simul_par3_new[[D]][1])))+1
  OW_par3_new_Var1<-matrix(NA,nrow=nrow(out_sample[ind_subset1_par3_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par3_subset1[i]==1){OW_par3_new_Var1[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset1_par3_new,],k=1)}
    else if(mix_ind_par3_subset1[i]==2){OW_par3_new_Var1[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset1_par3_new,],k=1)}
    else if(mix_ind_par3_subset1[i]==3){OW_par3_new_Var1[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset1_par3_new,],k=1)}
    else if(mix_ind_par3_subset1[i]==4){OW_par3_new_Var1[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset1_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset1[i]==5){OW_par3_new_Var1[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset1_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset1[i]==6){OW_par3_new_Var1[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset1_par3_new,],k=1)}
    else if(mix_ind_par3_subset1[i]==7){OW_par3_new_Var1[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset1_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset1[i]==8){OW_par3_new_Var1[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset1_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset1[i]==9){OW_par3_new_Var1[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset1_par3_new,],k=1)}
    else if(mix_ind_par3_subset1[i]==10){OW_par3_new_Var1[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset1_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset1[i]==11){OW_par3_new_Var1[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset1_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset1[i]==12){OW_par3_new_Var1[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset1_par3_new,],k=1)}
    else if(mix_ind_par3_subset1[i]==13){OW_par3_new_Var1[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset1_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset1[i]==14){OW_par3_new_Var1[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset1_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset1[i]==15){OW_par3_new_Var1[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset1_par3_new,],newdata_sigma=out_sample_sigma[ind_subset1_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset1[i]==16){OW_par3_new_Var1[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset1_par3_new,],newdata_sigma =out_sample_sigma[ind_subset1_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset1[i]==17){OW_par3_new_Var1[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset1_par3_new,],2,40,k=1)}
    else{OW_par3_new_Var1[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset1_par3_new,],out_sample[ind_subset1_par3_new,],N_in,2,40,k=1)}
  }
  
  
  mix_ind_par3_subset2<-findInterval(U,cumsum(unlist(model_weights_simul_par3_new[[D]][2])))+1
  OW_par3_new_Var2<-matrix(NA,nrow=nrow(out_sample[ind_subset2_par3_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par3_subset2[i]==1){OW_par3_new_Var2[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset2_par3_new,],k=1)}
    else if(mix_ind_par3_subset2[i]==2){OW_par3_new_Var2[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset2_par3_new,],k=1)}
    else if(mix_ind_par3_subset2[i]==3){OW_par3_new_Var2[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset2_par3_new,],k=1)}
    else if(mix_ind_par3_subset2[i]==4){OW_par3_new_Var2[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset2_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset2[i]==5){OW_par3_new_Var2[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset2_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset2[i]==6){OW_par3_new_Var2[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset2_par3_new,],k=1)}
    else if(mix_ind_par3_subset2[i]==7){OW_par3_new_Var2[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset2_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset2[i]==8){OW_par3_new_Var2[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset2_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset2[i]==9){OW_par3_new_Var2[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset2_par3_new,],k=1)}
    else if(mix_ind_par3_subset2[i]==10){OW_par3_new_Var2[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset2_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset2[i]==11){OW_par3_new_Var2[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset2_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset2[i]==12){OW_par3_new_Var2[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset2_par3_new,],k=1)}
    else if(mix_ind_par3_subset2[i]==13){OW_par3_new_Var2[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset2_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset2[i]==14){OW_par3_new_Var2[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset2_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset2[i]==15){OW_par3_new_Var2[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset2_par3_new,],newdata_sigma=out_sample_sigma[ind_subset2_par3_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par3_subset2[i]==16){OW_par3_new_Var2[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset2_par3_new,],newdata_sigma =out_sample_sigma[ind_subset2_par3_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par3_subset2[i]==17){OW_par3_new_Var2[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset2_par3_new,],2,40,k=1)}
    else{OW_par3_new_Var2[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset2_par3_new,],out_sample[ind_subset2_par3_new,],N_in,2,40,k=1)}
  }
  OW_par3_new_Var<-rbind(OW_par3_new_Var1,OW_par3_new_Var2)
  
  ##Ensemble 4
  
  mix_ind_par4_subset1<-findInterval(U,cumsum(unlist(model_weights_simul_par4_new[[D]][1])))+1
  OW_par4_new_Var1<-matrix(NA,nrow=nrow(out_sample[ind_subset1_par4_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par4_subset1[i]==1){OW_par4_new_Var1[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset1_par4_new,],k=1)}
    else if(mix_ind_par4_subset1[i]==2){OW_par4_new_Var1[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset1_par4_new,],k=1)}
    else if(mix_ind_par4_subset1[i]==3){OW_par4_new_Var1[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset1_par4_new,],k=1)}
    else if(mix_ind_par4_subset1[i]==4){OW_par4_new_Var1[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset1_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset1[i]==5){OW_par4_new_Var1[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset1_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset1[i]==6){OW_par4_new_Var1[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset1_par4_new,],k=1)}
    else if(mix_ind_par4_subset1[i]==7){OW_par4_new_Var1[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset1_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset1[i]==8){OW_par4_new_Var1[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset1_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset1[i]==9){OW_par4_new_Var1[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset1_par4_new,],k=1)}
    else if(mix_ind_par4_subset1[i]==10){OW_par4_new_Var1[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset1_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset1[i]==11){OW_par4_new_Var1[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset1_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset1[i]==12){OW_par4_new_Var1[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset1_par4_new,],k=1)}
    else if(mix_ind_par4_subset1[i]==13){OW_par4_new_Var1[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset1_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset1[i]==14){OW_par4_new_Var1[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset1_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset1[i]==15){OW_par4_new_Var1[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset1_par4_new,],newdata_sigma=out_sample_sigma[ind_subset1_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset1[i]==16){OW_par4_new_Var1[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset1_par4_new,],newdata_sigma =out_sample_sigma[ind_subset1_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset1[i]==17){OW_par4_new_Var1[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset1_par4_new,],2,40,k=1)}
    else{OW_par4_new_Var1[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset1_par4_new,],out_sample[ind_subset1_par4_new,],N_in,2,40,k=1)}
  }
  
  
  mix_ind_par4_subset2<-findInterval(U,cumsum(unlist(model_weights_simul_par4_new[[D]][2])))+1
  OW_par4_new_Var2<-matrix(NA,nrow=nrow(out_sample[ind_subset2_par4_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par4_subset2[i]==1){OW_par4_new_Var2[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset2_par4_new,],k=1)}
    else if(mix_ind_par4_subset2[i]==2){OW_par4_new_Var2[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset2_par4_new,],k=1)}
    else if(mix_ind_par4_subset2[i]==3){OW_par4_new_Var2[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset2_par4_new,],k=1)}
    else if(mix_ind_par4_subset2[i]==4){OW_par4_new_Var2[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset2_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset2[i]==5){OW_par4_new_Var2[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset2_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset2[i]==6){OW_par4_new_Var2[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset2_par4_new,],k=1)}
    else if(mix_ind_par4_subset2[i]==7){OW_par4_new_Var2[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset2_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset2[i]==8){OW_par4_new_Var2[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset2_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset2[i]==9){OW_par4_new_Var2[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset2_par4_new,],k=1)}
    else if(mix_ind_par4_subset2[i]==10){OW_par4_new_Var2[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset2_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset2[i]==11){OW_par4_new_Var2[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset2_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset2[i]==12){OW_par4_new_Var2[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset2_par4_new,],k=1)}
    else if(mix_ind_par4_subset2[i]==13){OW_par4_new_Var2[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset2_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset2[i]==14){OW_par4_new_Var2[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset2_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset2[i]==15){OW_par4_new_Var2[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset2_par4_new,],newdata_sigma=out_sample_sigma[ind_subset2_par4_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par4_subset2[i]==16){OW_par4_new_Var2[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset2_par4_new,],newdata_sigma =out_sample_sigma[ind_subset2_par4_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par4_subset2[i]==17){OW_par4_new_Var2[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset2_par4_new,],2,40,k=1)}
    else{OW_par4_new_Var2[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset2_par4_new,],out_sample[ind_subset2_par4_new,],N_in,2,40,k=1)}
  }
  OW_par4_new_Var<-rbind(OW_par4_new_Var1,OW_par4_new_Var2)
  
  
  ##Ensemble 5
  mix_ind_par5_subset1<-findInterval(U,cumsum(unlist(model_weights_simul_par5_new[[D]][1])))+1
  OW_par5_new_Var1<-matrix(NA,nrow=nrow(out_sample[ind_subset1_par5_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par5_subset1[i]==1){OW_par5_new_Var1[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset1_par5_new,],k=1)}
    else if(mix_ind_par5_subset1[i]==2){OW_par5_new_Var1[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset1_par5_new,],k=1)}
    else if(mix_ind_par5_subset1[i]==3){OW_par5_new_Var1[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset1_par5_new,],k=1)}
    else if(mix_ind_par5_subset1[i]==4){OW_par5_new_Var1[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset1_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset1[i]==5){OW_par5_new_Var1[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset1_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset1[i]==6){OW_par5_new_Var1[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset1_par5_new,],k=1)}
    else if(mix_ind_par5_subset1[i]==7){OW_par5_new_Var1[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset1_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset1[i]==8){OW_par5_new_Var1[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset1_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset1[i]==9){OW_par5_new_Var1[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset1_par5_new,],k=1)}
    else if(mix_ind_par5_subset1[i]==10){OW_par5_new_Var1[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset1_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset1[i]==11){OW_par5_new_Var1[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset1_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset1[i]==12){OW_par5_new_Var1[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset1_par5_new,],k=1)}
    else if(mix_ind_par5_subset1[i]==13){OW_par5_new_Var1[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset1_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset1[i]==14){OW_par5_new_Var1[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset1_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset1[i]==15){OW_par5_new_Var1[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset1_par5_new,],newdata_sigma=out_sample_sigma[ind_subset1_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset1[i]==16){OW_par5_new_Var1[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset1_par5_new,],newdata_sigma =out_sample_sigma[ind_subset1_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset1[i]==17){OW_par5_new_Var1[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset1_par5_new,],2,40,k=1)}
    else{OW_par5_new_Var1[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset1_par5_new,],out_sample[ind_subset1_par5_new,],N_in,2,40,k=1)}
  }
  
  
  mix_ind_par5_subset2<-findInterval(U,cumsum(unlist(model_weights_simul_par5_new[[D]][2])))+1
  OW_par5_new_Var2<-matrix(NA,nrow=nrow(out_sample[ind_subset2_par5_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par5_subset2[i]==1){OW_par5_new_Var2[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset2_par5_new,],k=1)}
    else if(mix_ind_par5_subset2[i]==2){OW_par5_new_Var2[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset2_par5_new,],k=1)}
    else if(mix_ind_par5_subset2[i]==3){OW_par5_new_Var2[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset2_par5_new,],k=1)}
    else if(mix_ind_par5_subset2[i]==4){OW_par5_new_Var2[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset2_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset2[i]==5){OW_par5_new_Var2[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset2_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset2[i]==6){OW_par5_new_Var2[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset2_par5_new,],k=1)}
    else if(mix_ind_par5_subset2[i]==7){OW_par5_new_Var2[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset2_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset2[i]==8){OW_par5_new_Var2[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset2_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset2[i]==9){OW_par5_new_Var2[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset2_par5_new,],k=1)}
    else if(mix_ind_par5_subset2[i]==10){OW_par5_new_Var2[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset2_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset2[i]==11){OW_par5_new_Var2[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset2_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset2[i]==12){OW_par5_new_Var2[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset2_par5_new,],k=1)}
    else if(mix_ind_par5_subset2[i]==13){OW_par5_new_Var2[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset2_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset2[i]==14){OW_par5_new_Var2[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset2_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset2[i]==15){OW_par5_new_Var2[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset2_par5_new,],newdata_sigma=out_sample_sigma[ind_subset2_par5_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par5_subset2[i]==16){OW_par5_new_Var2[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset2_par5_new,],newdata_sigma =out_sample_sigma[ind_subset2_par5_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par5_subset2[i]==17){OW_par5_new_Var2[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset2_par5_new,],2,40,k=1)}
    else{OW_par5_new_Var2[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset2_par5_new,],out_sample[ind_subset2_par5_new,],N_in,2,40,k=1)}
  }
  OW_par5_new_Var<-rbind(OW_par5_new_Var1,OW_par5_new_Var2)
  
  
  ## Ensemble 6
  mix_ind_par6_subset1<-findInterval(U,cumsum(unlist(model_weights_simul_par6_new[[D]][1])))+1
  OW_par6_new_Var1<-matrix(NA,nrow=nrow(out_sample[ind_subset1_par6_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par6_subset1[i]==1){OW_par6_new_Var1[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset1_par6_new,],k=1)}
    else if(mix_ind_par6_subset1[i]==2){OW_par6_new_Var1[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset1_par6_new,],k=1)}
    else if(mix_ind_par6_subset1[i]==3){OW_par6_new_Var1[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset1_par6_new,],k=1)}
    else if(mix_ind_par6_subset1[i]==4){OW_par6_new_Var1[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset1_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset1[i]==5){OW_par6_new_Var1[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset1_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset1[i]==6){OW_par6_new_Var1[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset1_par6_new,],k=1)}
    else if(mix_ind_par6_subset1[i]==7){OW_par6_new_Var1[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset1_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset1[i]==8){OW_par6_new_Var1[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset1_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset1[i]==9){OW_par6_new_Var1[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset1_par6_new,],k=1)}
    else if(mix_ind_par6_subset1[i]==10){OW_par6_new_Var1[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset1_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset1[i]==11){OW_par6_new_Var1[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset1_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset1[i]==12){OW_par6_new_Var1[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset1_par6_new,],k=1)}
    else if(mix_ind_par6_subset1[i]==13){OW_par6_new_Var1[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset1_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset1[i]==14){OW_par6_new_Var1[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset1_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset1[i]==15){OW_par6_new_Var1[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset1_par6_new,],newdata_sigma=out_sample_sigma[ind_subset1_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset1[i]==16){OW_par6_new_Var1[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset1_par6_new,],newdata_sigma =out_sample_sigma[ind_subset1_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset1[i]==17){OW_par6_new_Var1[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset1_par6_new,],2,40,k=1)}
    else{OW_par6_new_Var1[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset1_par6_new,],out_sample[ind_subset1_par6_new,],N_in,2,40,k=1)}
  }
  
  
  mix_ind_par6_subset2<-findInterval(U,cumsum(unlist(model_weights_simul_par6_new[[D]][2])))+1
  OW_par6_new_Var2<-matrix(NA,nrow=nrow(out_sample[ind_subset2_par6_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par6_subset2[i]==1){OW_par6_new_Var2[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset2_par6_new,],k=1)}
    else if(mix_ind_par6_subset2[i]==2){OW_par6_new_Var2[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset2_par6_new,],k=1)}
    else if(mix_ind_par6_subset2[i]==3){OW_par6_new_Var2[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset2_par6_new,],k=1)}
    else if(mix_ind_par6_subset2[i]==4){OW_par6_new_Var2[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset2_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset2[i]==5){OW_par6_new_Var2[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset2_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset2[i]==6){OW_par6_new_Var2[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset2_par6_new,],k=1)}
    else if(mix_ind_par6_subset2[i]==7){OW_par6_new_Var2[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset2_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset2[i]==8){OW_par6_new_Var2[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset2_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset2[i]==9){OW_par6_new_Var2[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset2_par6_new,],k=1)}
    else if(mix_ind_par6_subset2[i]==10){OW_par6_new_Var2[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset2_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset2[i]==11){OW_par6_new_Var2[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset2_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset2[i]==12){OW_par6_new_Var2[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset2_par6_new,],k=1)}
    else if(mix_ind_par6_subset2[i]==13){OW_par6_new_Var2[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset2_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset2[i]==14){OW_par6_new_Var2[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset2_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset2[i]==15){OW_par6_new_Var2[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset2_par6_new,],newdata_sigma=out_sample_sigma[ind_subset2_par6_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par6_subset2[i]==16){OW_par6_new_Var2[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset2_par6_new,],newdata_sigma =out_sample_sigma[ind_subset2_par6_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par6_subset2[i]==17){OW_par6_new_Var2[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset2_par6_new,],2,40,k=1)}
    else{OW_par6_new_Var2[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset2_par6_new,],out_sample[ind_subset2_par6_new,],N_in,2,40,k=1)}
  }
  OW_par6_new_Var<-rbind(OW_par6_new_Var1,OW_par6_new_Var2)
  
  
  ##Ensemble 7
  
  mix_ind_par7_subset1<-findInterval(U,cumsum(unlist(model_weights_simul_par7_new[[D]][1])))+1
  OW_par7_new_Var1<-matrix(NA,nrow=nrow(out_sample[ind_subset1_par7_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par7_subset1[i]==1){OW_par7_new_Var1[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset1_par7_new,],k=1)}
    else if(mix_ind_par7_subset1[i]==2){OW_par7_new_Var1[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset1_par7_new,],k=1)}
    else if(mix_ind_par7_subset1[i]==3){OW_par7_new_Var1[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset1_par7_new,],k=1)}
    else if(mix_ind_par7_subset1[i]==4){OW_par7_new_Var1[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset1_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset1[i]==5){OW_par7_new_Var1[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset1_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset1[i]==6){OW_par7_new_Var1[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset1_par7_new,],k=1)}
    else if(mix_ind_par7_subset1[i]==7){OW_par7_new_Var1[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset1_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset1[i]==8){OW_par7_new_Var1[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset1_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset1[i]==9){OW_par7_new_Var1[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset1_par7_new,],k=1)}
    else if(mix_ind_par7_subset1[i]==10){OW_par7_new_Var1[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset1_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset1[i]==11){OW_par7_new_Var1[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset1_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset1[i]==12){OW_par7_new_Var1[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset1_par7_new,],k=1)}
    else if(mix_ind_par7_subset1[i]==13){OW_par7_new_Var1[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset1_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset1[i]==14){OW_par7_new_Var1[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset1_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset1[i]==15){OW_par7_new_Var1[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset1_par7_new,],newdata_sigma=out_sample_sigma[ind_subset1_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset1[i]==16){OW_par7_new_Var1[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset1_par7_new,],newdata_sigma =out_sample_sigma[ind_subset1_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset1[i]==17){OW_par7_new_Var1[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset1_par7_new,],2,40,k=1)}
    else{OW_par7_new_Var1[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset1_par7_new,],out_sample[ind_subset1_par7_new,],N_in,2,40,k=1)}
  }
  
  
  mix_ind_par7_subset2<-findInterval(U,cumsum(unlist(model_weights_simul_par7_new[[D]][2])))+1
  OW_par7_new_Var2<-matrix(NA,nrow=nrow(out_sample[ind_subset2_par7_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par7_subset2[i]==1){OW_par7_new_Var2[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset2_par7_new,],k=1)}
    else if(mix_ind_par7_subset2[i]==2){OW_par7_new_Var2[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset2_par7_new,],k=1)}
    else if(mix_ind_par7_subset2[i]==3){OW_par7_new_Var2[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset2_par7_new,],k=1)}
    else if(mix_ind_par7_subset2[i]==4){OW_par7_new_Var2[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset2_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset2[i]==5){OW_par7_new_Var2[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset2_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset2[i]==6){OW_par7_new_Var2[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset2_par7_new,],k=1)}
    else if(mix_ind_par7_subset2[i]==7){OW_par7_new_Var2[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset2_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset2[i]==8){OW_par7_new_Var2[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset2_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset2[i]==9){OW_par7_new_Var2[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset2_par7_new,],k=1)}
    else if(mix_ind_par7_subset2[i]==10){OW_par7_new_Var2[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset2_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset2[i]==11){OW_par7_new_Var2[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset2_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset2[i]==12){OW_par7_new_Var2[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset2_par7_new,],k=1)}
    else if(mix_ind_par7_subset2[i]==13){OW_par7_new_Var2[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset2_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset2[i]==14){OW_par7_new_Var2[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset2_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset2[i]==15){OW_par7_new_Var2[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset2_par7_new,],newdata_sigma=out_sample_sigma[ind_subset2_par7_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par7_subset2[i]==16){OW_par7_new_Var2[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset2_par7_new,],newdata_sigma =out_sample_sigma[ind_subset2_par7_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par7_subset2[i]==17){OW_par7_new_Var2[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset2_par7_new,],2,40,k=1)}
    else{OW_par7_new_Var2[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset2_par7_new,],out_sample[ind_subset2_par7_new,],N_in,2,40,k=1)}
  }
  OW_par7_new_Var<-rbind(OW_par7_new_Var1,OW_par7_new_Var2)
  
  
  ##Ensemble 8
  
  mix_ind_par8_subset1<-findInterval(U,cumsum(unlist(model_weights_simul_par8_new[[D]][1])))+1
  OW_par8_new_Var1<-matrix(NA,nrow=nrow(out_sample[ind_subset1_par8_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par8_subset1[i]==1){OW_par8_new_Var1[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset1_par8_new,],k=1)}
    else if(mix_ind_par8_subset1[i]==2){OW_par8_new_Var1[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset1_par8_new,],k=1)}
    else if(mix_ind_par8_subset1[i]==3){OW_par8_new_Var1[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset1_par8_new,],k=1)}
    else if(mix_ind_par8_subset1[i]==4){OW_par8_new_Var1[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset1_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset1[i]==5){OW_par8_new_Var1[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset1_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset1[i]==6){OW_par8_new_Var1[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset1_par8_new,],k=1)}
    else if(mix_ind_par8_subset1[i]==7){OW_par8_new_Var1[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset1_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset1[i]==8){OW_par8_new_Var1[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset1_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset1[i]==9){OW_par8_new_Var1[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset1_par8_new,],k=1)}
    else if(mix_ind_par8_subset1[i]==10){OW_par8_new_Var1[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset1_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset1[i]==11){OW_par8_new_Var1[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset1_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset1[i]==12){OW_par8_new_Var1[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset1_par8_new,],k=1)}
    else if(mix_ind_par8_subset1[i]==13){OW_par8_new_Var1[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset1_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset1[i]==14){OW_par8_new_Var1[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset1_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset1[i]==15){OW_par8_new_Var1[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset1_par8_new,],newdata_sigma=out_sample_sigma[ind_subset1_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset1[i]==16){OW_par8_new_Var1[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset1_par8_new,],newdata_sigma =out_sample_sigma[ind_subset1_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset1[i]==17){OW_par8_new_Var1[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset1_par8_new,],2,40,k=1)}
    else{OW_par8_new_Var1[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset1_par8_new,],out_sample[ind_subset1_par8_new,],N_in,2,40,k=1)}
  }
  
  
  mix_ind_par8_subset2<-findInterval(U,cumsum(unlist(model_weights_simul_par8_new[[D]][2])))+1
  OW_par8_new_Var2<-matrix(NA,nrow=nrow(out_sample[ind_subset2_par8_new,]),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_par8_subset2[i]==1){OW_par8_new_Var2[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample[ind_subset2_par8_new,],k=1)}
    else if(mix_ind_par8_subset2[i]==2){OW_par8_new_Var2[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample[ind_subset2_par8_new,],k=1)}
    else if(mix_ind_par8_subset2[i]==3){OW_par8_new_Var2[,i]<-simulate_ZALN(LN_in,newdata=out_sample[ind_subset2_par8_new,],k=1)}
    else if(mix_ind_par8_subset2[i]==4){OW_par8_new_Var2[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample[ind_subset2_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset2[i]==5){OW_par8_new_Var2[,i]<-simulate_LN(LN_optimTau,newdata=out_sample[ind_subset2_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset2[i]==6){OW_par8_new_Var2[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho[ind_subset2_par8_new,],k=1)}
    else if(mix_ind_par8_subset2[i]==7){OW_par8_new_Var2[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho[ind_subset2_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset2[i]==8){OW_par8_new_Var2[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho[ind_subset2_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset2[i]==9){OW_par8_new_Var2[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample[ind_subset2_par8_new,],k=1)}
    else if(mix_ind_par8_subset2[i]==10){OW_par8_new_Var2[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample[ind_subset2_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset2[i]==11){OW_par8_new_Var2[,i]<-simulate_LN(glm_LN_Cal_In,out_sample[ind_subset2_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset2[i]==12){OW_par8_new_Var2[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric[ind_subset2_par8_new,],k=1)}
    else if(mix_ind_par8_subset2[i]==13){OW_par8_new_Var2[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric[ind_subset2_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset2[i]==14){OW_par8_new_Var2[,i]<-simulate_LN(sp_LN_In,out_sample_numeric[ind_subset2_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset2[i]==15){OW_par8_new_Var2[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu[ind_subset2_par8_new,],newdata_sigma=out_sample_sigma[ind_subset2_par8_new,],k=1,tau=tau_Ga)}
    else if(mix_ind_par8_subset2[i]==16){OW_par8_new_Var2[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu[ind_subset2_par8_new,],newdata_sigma =out_sample_sigma[ind_subset2_par8_new,],k=1,tau=tau_LN)}
    else if(mix_ind_par8_subset2[i]==17){OW_par8_new_Var2[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample[ind_subset2_par8_new,],2,40,k=1)}
    else{OW_par8_new_Var2[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out[ind_subset2_par8_new,],out_sample[ind_subset2_par8_new,],N_in,2,40,k=1)}
  }
  OW_par8_new_Var<-rbind(OW_par8_new_Var1,OW_par8_new_Var2)
  
 
  
  
  ###Equally Weighted Ensemble
  
  mix_ind_EW<-findInterval(U,cumsum(rep(1/18,18)))+1
  OW_newEW<-matrix(NA,nrow=nrow(out_sample),ncol=1000)
  for (i in 1:1000){
    if(mix_ind_EW[i]==1){OW_newEW[,i]<-simulate_ODP(ODP_GLM_in,newdata=out_sample,k=1)}
    else if(mix_ind_EW[i]==2){OW_newEW[,i]<-simulate_ZAGA(ZAGA_in,newdata=out_sample,k=1)}
    else if(mix_ind_EW[i]==3){OW_newEW[,i]<-simulate_ZALN(LN_in,newdata=out_sample,k=1)}
    else if(mix_ind_EW[i]==4){OW_newEW[,i]<-simulate_GA(Ga_optimTau,newdata=out_sample,k=1,tau=tau_Ga)}
    else if(mix_ind_EW[i]==5){OW_newEW[,i]<-simulate_LN(LN_optimTau,newdata=out_sample,k=1,tau=tau_LN)}
    else if(mix_ind_EW[i]==6){OW_newEW[,i]<-simulate_ODP(glm_ODP_Ho_In,newdata=out_sample_Ho,k=1)}
    else if(mix_ind_EW[i]==7){OW_newEW[,i]<-simulate_GA(glm_Ga_Ho_In,newdata=out_sample_Ho,k=1,tau=tau_Ga)}
    else if(mix_ind_EW[i]==8){OW_newEW[,i]<-simulate_LN(glm_LN_Ho_In,newdata=out_sample_Ho,k=1,tau=tau_LN)}
    else if(mix_ind_EW[i]==9){OW_newEW[,i]<-simulate_ODP(glm_ODP_Cal_In,out_sample,k=1)}
    else if(mix_ind_EW[i]==10){OW_newEW[,i]<-simulate_GA(glm_Ga_Cal_In,out_sample,k=1,tau=tau_Ga)}
    else if(mix_ind_EW[i]==11){OW_newEW[,i]<-simulate_LN(glm_LN_Cal_In,out_sample,k=1,tau=tau_LN)}
    else if(mix_ind_EW[i]==12){OW_newEW[,i]<-simulate_NO(sp_Normal_In,newdata=out_sample_numeric,k=1)}
    else if(mix_ind_EW[i]==13){OW_newEW[,i]<-simulate_GA(sp_Gamma_In,out_sample_numeric,k=1,tau=tau_Ga)}
    else if(mix_ind_EW[i]==14){OW_newEW[,i]<-simulate_LN(sp_LN_In,out_sample_numeric,k=1,tau=tau_LN)}
    else if(mix_ind_EW[i]==15){OW_newEW[,i]<-simulate_GAGamlss(gamlss2_GA_In,newdata_mu=out_sample_mu,newdata_sigma=out_sample_sigma,k=1,tau=tau_Ga)}
    else if(mix_ind_EW[i]==16){OW_newEW[,i]<-simulate_LNGamlss(gamlss2_LN_In,newdata_mu=out_sample_mu,newdata_sigma =out_sample_sigma,k=1,tau=tau_LN)}
    else if(mix_ind_EW[i]==17){OW_newEW[,i]<-simulate_PPCI(fit_ODP_ppci_In,N_in,out_sample,2,40,k=1)}
    else{OW_newEW[,i]<-simulate_PPCF(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out,out_sample,N_in,2,40,k=1)}
  }
  ## BMV
  BMV_var1<-simulate_LN(sp_LN_In,out_sample_numeric,k=1000,tau=tau_LN)
  
  ###Calculate the Reserve
  
  ####Log-Score Optimized Weighted Ensemble
  reserve_par1_newSimulScheme[,D]<-apply(OW_par1_new_Var,FUN=function(x) sum(x*10000),MARGIN=2)
  reserve_par2_newSimulScheme[,D]<-apply(OW_par2_new_Var,FUN=function(x) sum(x*10000),MARGIN=2)
  reserve_par3_newSimulScheme[,D]<-apply(OW_par3_new_Var,FUN=function(x) sum(x*10000),MARGIN=2)
  reserve_par4_newSimulScheme[,D]<-apply(OW_par4_new_Var,FUN=function(x) sum(x*10000),MARGIN=2)
  reserve_par5_newSimulScheme[,D]<-apply(OW_par5_new_Var,FUN=function(x) sum(x*10000),MARGIN=2)
  reserve_par6_newSimulScheme[,D]<-apply(OW_par6_new_Var,FUN=function(x) sum(x*10000),MARGIN=2)
  reserve_par7_newSimulScheme[,D]<-apply(OW_par7_new_Var,FUN=function(x) sum(x*10000),MARGIN=2)
  reserve_par8_newSimulScheme[,D]<-apply(OW_par8_new_Var,FUN=function(x) sum(x*10000),MARGIN=2)
  

  ##Equally Weighted Ensemble and BMV
  reserve_par1_newEW[,D]<-apply(OW_newEW,FUN=function(x) sum(x*10000),MARGIN=2)
  reserve_BMV_newSimulScheme[,D]<-apply(BMV_var1,FUN=function(x) sum(x*10000),MARGIN=2)
  
}