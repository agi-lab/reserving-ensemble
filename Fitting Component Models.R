
index_best_validmod2<-c()
valid_dens_list<-list()
out_dens_list<-list()
out_mu_list<-list()
valid_mu_list<-list()
model_weights_simul<-list()
future_claims<-matrix(NA,nrow=780,ncol=ntri)



for (D in 1:ntri){
  
  set.seed(20200130+D)
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
  
  # Construction of Training Set:
  train_1<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)==1,]
  train_2<-dt_withInflation_past[as.numeric(dt_withInflation_past$Calendar)<=34&as.numeric(dt_withInflation_past$origin)<=32&as.numeric(dt_withInflation_past$origin)>=2,]
  train_3<-dt_withInflation_past[dt_withInflation_past$dev==1&as.numeric(dt_withInflation_past$origin)>=33&as.numeric(dt_withInflation_past$origin)<=40,]
  train_4<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)>=33&as.numeric(dt_withInflation_past$origin)<=34&as.numeric(dt_withInflation_past$dev)==2,]
  train_5<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)==40&as.numeric(dt_withInflation_past$dev)==1,]
  train<-rbind(train_1,train_2,train_3,train_4,train_5)
  # Construction of Validation set:
  valid_1<-dt_withInflation_past[as.numeric(dt_withInflation_past$Calendar)>=35&as.numeric(dt_withInflation_past$Calendar)<=38&as.numeric(dt_withInflation_past$origin)>=2&as.numeric(dt_withInflation_past$origin)<=32,]
  valid_2<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)==33&as.numeric(dt_withInflation_past$dev)>=3&as.numeric(dt_withInflation_past$dev)<=5,]
  valid_3<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)==34&as.numeric(dt_withInflation_past$dev)>=3&as.numeric(dt_withInflation_past$dev)<=4,]
  valid_4<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)>=35&as.numeric(dt_withInflation_past$origin)<=37&as.numeric(dt_withInflation_past$dev)>=2&as.numeric(dt_withInflation_past$dev)<=3,]
  valid_5<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)>=38&as.numeric(dt_withInflation_past$origin)<=39&as.numeric(dt_withInflation_past$dev)==2,]
  valid_6<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)==40&as.numeric(dt_withInflation_past$dev)==1,]
  valid<-rbind(valid_1,valid_2,valid_3,valid_4,valid_5,valid_6)
  # Construction of Test Data
  test_1<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)>=2&as.numeric(dt_withInflation_past$origin)<=35&as.numeric(dt_withInflation_past$Calendar)>=39&as.numeric(dt_withInflation_past$Calendar)<=41,]
  test_2<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)==36&as.numeric(dt_withInflation_past$dev)>=4&as.numeric(dt_withInflation_past$dev)<=5,]
  test_3<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)==37&as.numeric(dt_withInflation_past$dev)==4,]
  test_4<-dt_withInflation_past[as.numeric(dt_withInflation_past$origin)==38&as.numeric(dt_withInflation_past$dev)==3,]
  test<-rbind(test_1,test_2,test_3,test_4)
  test$origin=as.factor(test$origin)
  test$dev=as.factor(test$dev)
  test<-test[order(as.numeric(as.character(test$origin))),]
  train$origin=as.factor(train$origin)
  train$dev=as.factor(train$dev)
  train$Calendar=as.numeric(train$Calendar)
  #######Generate Predictive Density for Validation Sets
  test<-test[order(as.numeric(as.character(test$origin))),]
  aug_valid<-rbind(valid,test)
  aug_valid<-aug_valid[order(as.numeric(as.character(aug_valid$origin))),]
  aug_valid$origin<-as.factor(aug_valid$origin)
  aug_valid$dev<-as.factor(aug_valid$dev)
  aug_valid$Calendar<-as.numeric(aug_valid$Calendar)
  aug_valid_numeric<-aug_valid
  aug_valid_numeric$origin=as.numeric(as.character(aug_valid$origin))
  aug_valid_numeric$dev=as.numeric(as.character(aug_valid$dev))
  aug_valid_numeric$Calendar=as.numeric(as.character(aug_valid$Calendar))
  ###Fitting of all the models in the training Sets
  ####Basic Models 
  #ODP
  ODP_GLM_train<-glm(formula=value~factor(origin)+factor(dev),family=quasipoisson(link="log"),data=train)
  #ZAGA
  gamma_1_train<-gamlss(formula=value~factor(origin)+factor(dev),nu.formula=~as.numeric(as.character(dev)),data=train,family=ZAGA(mu.link="log",sigma.link = "log", nu.link = "logit"))
  #ZALN
  LN_1_train<-gamlssZadj(y=value,mu.formula = ~factor(origin)+factor(dev),xi0.formula=~as.numeric(as.character(dev)),data=train,family=LOGNO(mu.link="identity",sigma.link="log"))
  tau_LN<-5
  tau_Ga<-5
  #Gamma 
  Ga_optimTau_train<-gamlss(formula=(value+tau_Ga)~factor(origin)+factor(dev),data=train,family=GA(mu.link="log", sigma.link ="log"))
  #Log-Normal
  LN_optimTau_train<-gamlss(formula=(value+tau_LN)~factor(origin)+factor(dev),data=train,family=LOGNO(mu.link="identity",sigma.link="log"))
  dens_ODPGLM2<-cal_dens_ODP(ODP_GLM_train,aug_valid)
  dens_GAGLM2<-cal_dens_GA(Ga_optimTau_train,tau_Ga,aug_valid)
  dens_LNGLM2<-cal_dens_LN(LN_optimTau_train,tau_LN,aug_valid)
  dens_ZAGA2<-cal_dens_ZAGA(gamma_1_train,aug_valid)
  dens_ZALN2<-cal_dens_ZALN(LN_1_train,aug_valid)
  
  train_new<-train
  train_new$origin=as.numeric(as.character(train$origin))
  train_new$dev=as.numeric(as.character(train$dev))
  
  ####With Hoerl Curve
  #Under ODP Assumption
  train_Ho<-train
  train_Ho$dev=as.numeric(as.character(train$dev))
  glm_ODP_Ho_tr1<-glm(formula=value~factor(origin)+log(dev)+dev,family=quasipoisson(link="log"),data=train_Ho)
  #R retun warning for rank deficient fit for Hoerl Curve ODP
  #Under Gamma Assumption
  glm_Ga_Ho_tr1<-gamlss(formula=(value+tau_Ga)~factor(origin)+log(dev)+dev,data=train_Ho,family=GA(mu.link="log", sigma.link ="log"))
  #Under Log-Normal Assumption
  glm_LN_Ho_tr1<-gamlss(formula=(value+tau_LN)~factor(origin)+log(dev)+dev,data=train_Ho,family=LOGNO(mu.link="identity",sigma.link="log"))
  #Calculate the density for Hoerl Curve
  #Density for Hoerl curve
  aug_valid_Ho<-aug_valid
  aug_valid_Ho$dev<-as.numeric(as.character(aug_valid$dev))
  dens_ODPHo2<-cal_dens_ODP(glm_ODP_Ho_tr1,aug_valid_Ho)
  dens_GaHo2<-cal_dens_GA(glm_Ga_Ho_tr1,tau_Ga,aug_valid_Ho)
  dens_LNHo2<-cal_dens_LN(glm_LN_Ho_tr1,tau_LN,aug_valid_Ho)
  
  ####With Calendar Periods
  #Under ODP Assumption
  glm_ODP_Cal_tr1<-glm(formula=value~factor(dev)+Calendar,family=quasipoisson(link="log"),data=train)
  #Under Gamma Assumption
  glm_Ga_Cal_tr1<-gamlss(formula=(value+tau_Ga)~factor(dev)+Calendar,data=train,family=GA(mu.link="log", sigma.link ="log"))
  #Under Log-Normal Assumption
  glm_LN_Cal_tr1<-gamlss(formula=(value+tau_LN)~factor(dev)+Calendar,data=train,family=LOGNO(mu.link="identity",sigma.link="log"))
  #Density for GLM that includes Calendar period effects
  dens_ODPCal2<-cal_dens_ODP(glm_ODP_Cal_tr1,aug_valid)
  dens_GaCal2<-cal_dens_GA(glm_Ga_Cal_tr1,tau_Ga,aug_valid)
  dens_LNCal2<-cal_dens_LN(glm_LN_Cal_tr1,tau_LN,aug_valid)
  
  ####Smoothing Spline
  sp_Normal2<-gamlss(formula=value~scs(origin)+scs(dev),data=train_new,family=NO(),trace=FALSE)
  sp_Gamma2<-gamlss(formula=(value+tau_Ga)~scs(origin)+scs(dev),data=train_new,family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
  sp_LN2<-gamlss(formula=(value+tau_LN)~scs(origin)+scs(dev),data=train_new,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  
  #Density for Smoothing Spline
  dens_SpNormal2<-cal_dens_Normal(sp_Normal2,aug_valid_numeric)
  dens_SpGamma2<-cal_dens_GA(sp_Gamma2,tau_Ga,aug_valid_numeric)
  dens_SpLN2<-cal_dens_LN(sp_LN2,tau_LN,aug_valid_numeric)
  
  ####GAMLSS
  #GAMLSS 2: Smooth Effects on the predictor for sigma term
  gamlss2_GA<-gamlss(formula=(value+tau_Ga)~factor(origin)+factor(dev),data=train,sigma.formula=~cs(as.numeric(as.character(dev))),family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
  gamlss2_LN<-gamlss(formula=(value+tau_Ga)~factor(origin)+factor(dev),data=train,sigma.formula=~cs(as.numeric(as.character(dev))),family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  #Density for GAMLSS
  aug_valid_mu<-aug_valid
  aug_valid_sigma<-aug_valid_numeric
  dens_GaGAMLSS2<-cal_dens_GA_Gamlss(gamlss2_GA,tau_Ga,aug_valid_mu,aug_valid_sigma)
  dens_LNGAMLSS2<-cal_dens_LN_Gamlss(gamlss2_LN,tau_LN,aug_valid_mu,aug_valid_sigma)
  
  
  
  
  
  
  ##Mean for Basic Models
  mu_ODPGLM2<-cal_mu_ODP(ODP_GLM_train,aug_valid)
  mu_GAGLM2<-ifelse(cal_mu_GA(Ga_optimTau_train,tau_Ga,aug_valid)>tau_Ga,cal_mu_GA(Ga_optimTau_train,tau_Ga,aug_valid)-tau_Ga,0)
  mu_LNGLM2<-ifelse(cal_mu_LN(LN_optimTau_train,tau_LN,aug_valid)>tau_LN,cal_mu_LN(LN_optimTau_train,tau_LN,aug_valid)-tau_LN,0)
  mu_ZAGA2<-cal_mu_ZAGA(gamma_1_train,aug_valid)
  mu_ZALN2<-cal_mu_ZALN(LN_1_train,aug_valid)
  
  
  
  
  #Mean for Hoerl curve
  aug_valid_Ho<-aug_valid
  aug_valid_Ho$dev<-as.numeric(as.character(aug_valid$dev))
  mu_ODPHo2<-cal_mu_ODP(glm_ODP_Ho_tr1,aug_valid_Ho)
  mu_GaHo2<-ifelse(cal_mu_GA(glm_Ga_Ho_tr1,tau_Ga,aug_valid_Ho)>tau_Ga,cal_mu_GA(glm_Ga_Ho_tr1,tau_Ga,aug_valid_Ho)-tau_Ga,0)
  mu_LNHo2<-ifelse(cal_mu_LN(glm_LN_Ho_tr1,tau_LN,aug_valid_Ho)>tau_LN,cal_mu_LN(glm_LN_Ho_tr1,tau_LN,aug_valid_Ho)-tau_LN,0)
  
  
  
  #Mean for Calendar Periods
  mu_ODPCal2<-cal_mu_ODP(glm_ODP_Cal_tr1,aug_valid)
  mu_GaCal2<-ifelse(cal_mu_GA(glm_Ga_Cal_tr1,tau_Ga,aug_valid)>tau_Ga,cal_mu_GA(glm_Ga_Cal_tr1,tau_Ga,aug_valid)-tau_Ga,0)
  mu_LNCal2<-ifelse(cal_mu_LN(glm_LN_Cal_tr1,tau_LN,aug_valid)>tau_LN,cal_mu_LN(glm_LN_Cal_tr1,tau_LN,aug_valid)-tau_LN,0)
  
  
  
  #Mean for Smoothing Spline
  mu_SpNormal2<-cal_mu_Normal(sp_Normal2,aug_valid_numeric)
  mu_SpGamma2<-ifelse(cal_mu_GA(sp_Gamma2,tau_Ga,aug_valid_numeric)>tau_Ga,cal_mu_GA(sp_Gamma2,tau_Ga,aug_valid_numeric)-tau_Ga,0)
  mu_SpLN2<-ifelse(cal_mu_LN(sp_LN2,tau_LN,aug_valid_numeric)>tau_LN,cal_mu_LN(sp_LN2,tau_LN,aug_valid_numeric)-tau_LN,0)
  
  #Mean for GAMLSS
  aug_valid_mu<-aug_valid
  aug_valid_sigma<-aug_valid_numeric
  #Mean for GAMLSS
  mu_GaGAMLSS2<-ifelse(cal_mu_GA_Gamlss(gamlss2_GA,tau_Ga,aug_valid_mu,aug_valid_sigma)>tau_Ga,cal_mu_GA_Gamlss(gamlss2_GA,tau_Ga,aug_valid_mu,aug_valid_sigma)-tau_Ga,0)
  mu_LNGAMLSS2<-ifelse(cal_mu_LN_Gamlss(gamlss2_LN,tau_LN,aug_valid_mu,aug_valid_sigma)>tau_LN,cal_mu_LN_Gamlss(gamlss2_LN,tau_LN,aug_valid_mu,aug_valid_sigma)-tau_LN,0)
  
  
  
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
  # Construction of Training Set for Notification Counts(to avoid data leakage):
  train_1_N<-NotiCount_dt[as.numeric(NotiCount_dt$origin)==1,]
  train_2_N<-NotiCount_dt[as.numeric(NotiCount_dt$Calendar)<=34&as.numeric(NotiCount_dt$origin)<=32&as.numeric(NotiCount_dt$origin)>=2,]
  train_3_N<-NotiCount_dt[NotiCount_dt$dev==1&as.numeric(NotiCount_dt$origin)>=33&as.numeric(NotiCount_dt$origin)<=40,]
  train_4_N<-NotiCount_dt[as.numeric(NotiCount_dt$origin)>=33&as.numeric(NotiCount_dt$origin)<=34&as.numeric(NotiCount_dt$dev)==2,]
  train_5_N<-NotiCount_dt[as.numeric(NotiCount_dt$origin)==40&as.numeric(NotiCount_dt$dev)==1,]
  train_NC<-rbind(train_1_N,train_2_N,train_3_N,train_4_N,train_5_N)
  #Construct of Validation set(identical to Paid Loss to avoid data leakage)
  valid_1_N<-NotiCount_dt[as.numeric(NotiCount_dt$Calendar)>=35&as.numeric(NotiCount_dt$Calendar)<=38&as.numeric(NotiCount_dt$origin)>=2&as.numeric(NotiCount_dt$origin)<=32,]
  valid_2_N<-NotiCount_dt[as.numeric(NotiCount_dt$origin)==33&as.numeric(NotiCount_dt$dev)>=3&as.numeric(NotiCount_dt$dev)<=5,]
  valid_3_N<-NotiCount_dt[as.numeric(NotiCount_dt$origin)==34&as.numeric(NotiCount_dt$dev)>=3&as.numeric(NotiCount_dt$dev)<=4,]
  valid_4_N<-NotiCount_dt[as.numeric(NotiCount_dt$origin)>=35&as.numeric(NotiCount_dt$origin)<=37&as.numeric(NotiCount_dt$dev)>=2&as.numeric(NotiCount_dt$dev)<=3,]
  valid_5_N<-NotiCount_dt[as.numeric(NotiCount_dt$origin)>=38&as.numeric(NotiCount_dt$origin)<=39&as.numeric(NotiCount_dt$dev)==2,]
  valid_6_N<-NotiCount_dt[as.numeric(NotiCount_dt$origin)==40&as.numeric(NotiCount_dt$dev)==1,]
  valid_N<-rbind(valid_1_N,valid_2_N,valid_3_N,valid_4_N,valid_5_N,valid_6_N)
  ##Fit a ODP/Chain-Ladder model to predict total claim notification count at each accident period 
  fit_nc<-glm(value~factor(origin)+factor(dev),data=train_NC,family=quasipoisson(link="log"))
  ##Obtain Estimated notification claims count for each accident period: Nk
  N<-c()
  N[1]<-sum(train_NC[train_NC$origin==1,]$value)
  for (i in 2:40){
    N[i]<-sum(train_NC[train_NC$origin==i,]$value)+sum(round(predict(fit_nc,newdata=test_NC[test_NC$origin==i,],type="response"),0))+sum(round(predict(fit_nc,newdata=valid_N[valid_N$origin==i,],type="response"),0))
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
  new_train<-train[order(as.numeric(train$origin)),]
  ###Repete each element in N_rep by the number of entries in each accident period
  N_rep<-rep(N,times=count_by_origin(train,1,40))
  new_train$PPCI=new_train$value/N_rep
  
  ##Fitting a ODP model on PPCI
  fit_ODP_ppci<-glm(PPCI~factor(dev),family=quasipoisson(link="log"),data=new_train)
  
  #####PPCF
  transaction_dataset$SettlementPeriod=ceiling(transaction_dataset$setldel)
  matri_SettlementCount<-matrix(NA,nrow=40,ncol=40)
  
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
  # Construction of Training Set for Finalized Counts(to avoid data leakage):
  train_1_F<-SettlCount_dt[as.numeric(SettlCount_dt$origin)==1,]
  train_2_F<-SettlCount_dt[as.numeric(SettlCount_dt$Calendar)<=34&as.numeric(SettlCount_dt$origin)<=32&as.numeric(SettlCount_dt$origin)>=2,]
  train_3_F<-SettlCount_dt[SettlCount_dt$dev==1&as.numeric(SettlCount_dt$origin)>=33&as.numeric(SettlCount_dt$origin)<=40,]
  train_4_F<-SettlCount_dt[as.numeric(SettlCount_dt$origin)>=33&as.numeric(SettlCount_dt$origin)<=34&as.numeric(SettlCount_dt$dev)==2,]
  train_5_F<-SettlCount_dt[as.numeric(SettlCount_dt$origin)==40&as.numeric(SettlCount_dt$dev)==1,]
  train_F<-rbind(train_1_F,train_2_F,train_3_F,train_4_F,train_5_F)
  train_F<-train_F[order(as.numeric(train_F$origin)),]
  
  
  #Construct of Validation set(identical to Paid Loss to avoid data leakage)
  valid_1_F<-SettlCount_dt[as.numeric(SettlCount_dt$Calendar)>=35&as.numeric(SettlCount_dt$Calendar)<=38&as.numeric(SettlCount_dt$origin)>=2&as.numeric(SettlCount_dt$origin)<=32,]
  valid_2_F<-SettlCount_dt[as.numeric(SettlCount_dt$origin)==33&as.numeric(SettlCount_dt$dev)>=3&as.numeric(SettlCount_dt$dev)<=5,]
  valid_3_F<-SettlCount_dt[as.numeric(SettlCount_dt$origin)==34&as.numeric(SettlCount_dt$dev)>=3&as.numeric(SettlCount_dt$dev)<=4,]
  valid_4_F<-SettlCount_dt[as.numeric(SettlCount_dt$origin)>=35&as.numeric(SettlCount_dt$origin)<=37&as.numeric(SettlCount_dt$dev)>=2&as.numeric(SettlCount_dt$dev)<=3,]
  valid_5_F<-SettlCount_dt[as.numeric(SettlCount_dt$origin)>=38&as.numeric(SettlCount_dt$origin)<=39&as.numeric(SettlCount_dt$dev)==2,]
  valid_6_F<-SettlCount_dt[as.numeric(SettlCount_dt$origin)==40&as.numeric(SettlCount_dt$dev)==1,]
  valid_F<-rbind(valid_1_F,valid_2_F,valid_3_F,valid_4_F,valid_5_F,valid_6_F)
  valid_F<-valid_F[order(as.numeric(valid_F$origin)),]
  #Create a triangle for Cummulative Claims Count
  cum_Count_tri<-incr2cum(NotiCount_tri)
  cum_Count_dt<-as.data.frame(cum_Count_tri,na.rm=TRUE)
  cum_Count_dt<-cum_Count_dt[order(as.numeric(cum_Count_dt$origin)),]
  cum_Count_dt$Calendar=as.numeric(as.character(cum_Count_dt$origin))+as.numeric(cum_Count_dt$dev)
  #Create a triangle for Cummulative finalised claims count
  cum_F_tri<-incr2cum(SettlCount_tri)
  cum_F_dt<-as.data.frame(cum_F_tri,na.rm=TRUE)
  cum_F_dt<-cum_F_dt[order(as.numeric(cum_F_dt$origin)),]
  cum_F_dt$Calendar=as.numeric(as.character(cum_F_dt$origin))+as.numeric(cum_F_dt$dev)
  # Construction of Training Set for Cummulative Notified Claim Counts(to avoid data leakage):
  train_1_cumN<-cum_Count_dt[as.numeric(as.character(cum_Count_dt$origin))==1,]
  train_2_cumN<-cum_Count_dt[as.numeric(cum_Count_dt$Calendar)<=34&as.numeric(as.character(cum_Count_dt$origin))<=32&as.numeric(as.character(cum_Count_dt$origin))>=2,]
  train_3_cumN<-cum_Count_dt[cum_Count_dt$dev==1&as.numeric(as.character(cum_Count_dt$origin))>=33&as.numeric(as.character(cum_Count_dt$origin))<=40,]
  train_4_cumN<-cum_Count_dt[as.numeric(as.character(cum_Count_dt$origin))>=33&as.numeric(as.character(cum_Count_dt$origin))<=34&as.numeric(cum_Count_dt$dev)==2,]
  train_5_cumN<-cum_Count_dt[as.numeric(as.character(cum_Count_dt$origin))==40&as.numeric(cum_Count_dt$dev)==1,]
  train_cumN<-rbind(train_1_cumN,train_2_cumN,train_3_cumN,train_4_cumN,train_5_cumN)
  train_cumN<-train_cumN[order(as.numeric(train_cumN$origin)),]
  # Construction of Training Set for Cummulative Finalized Counts(to avoid data leakage):
  train_1_cumF<-cum_F_dt[as.numeric(cum_F_dt$origin)==1,]
  train_2_cumF<-cum_F_dt[as.numeric(cum_F_dt$Calendar)<=34&as.numeric(cum_F_dt$origin)<=32&as.numeric(cum_F_dt$origin)>=2,]
  train_3_cumF<-cum_F_dt[cum_F_dt$dev==1&as.numeric(cum_F_dt$origin)>=33&as.numeric(cum_F_dt$origin)<=40,]
  train_4_cumF<-cum_F_dt[as.numeric(cum_F_dt$origin)>=33&as.numeric(cum_F_dt$origin)<=34&as.numeric(cum_F_dt$dev)==2,]
  train_5_cumF<-cum_F_dt[as.numeric(cum_F_dt$origin)==40&as.numeric(cum_F_dt$dev)==1,]
  train_cumF<-rbind(train_1_cumF,train_2_cumF,train_3_cumF,train_4_cumF,train_5_cumF)
  train_cumF<-train_cumF[order(as.numeric(train_cumF$origin)),]
  cum_U_dt<-data.frame(origin=cum_F_dt$origin,dev=cum_F_dt$dev,value=cum_Count_dt$value-cum_F_dt$value)
  #Generate a full data frame containing finalized claim count, unfinalized claim count and notification claim count
  NotiCount_dt<-NotiCount_dt[order(as.numeric(NotiCount_dt$origin)),]
  SettlCount_dt<-SettlCount_dt[order(as.numeric(SettlCount_dt$origin)),]
  UplusN<-list()
  
  
  #UplusN[1]<-NA
  #UplusN[2:nrow(cum_U_dt)]<-cum_U_dt$value[1:(nrow(cum_U_dt)-1)]+NotiCount_dt$value[2:nrow(NotiCount_dt)]
  
  for (i in 1:40){
    UplusN[[i]]<-cum_U_dt[cum_U_dt$origin==i,]$value[1:(nrow(cum_U_dt[cum_U_dt$origin==i,])-1)]+NotiCount_dt[NotiCount_dt$origin==i,]$value[2:nrow(NotiCount_dt[NotiCount_dt$origin==i,])]
  }
  
  UplusN<-as.vector(unlist(UplusN))
  FCount<-SettlCount_dt[SettlCount_dt$dev!=1,]$value
  FUN_dat<-data.frame(origin=SettlCount_dt[SettlCount_dt$dev!=1,]$origin,dev=SettlCount_dt[SettlCount_dt$dev!=1,]$dev,UplusN=UplusN[1:(length(UplusN)-2)],FCount=FCount)
  FUN_dat$Calendar=as.numeric(as.character(FUN_dat$origin))+as.numeric(FUN_dat$dev)
  #full_FUN<-data.frame(origin=NotiCount_dt$origin,dev=NotiCount_dt$dev,U=cum_U_dt$value,N=NotiCount_dt$value,UplusN=UplusN,F=SettlCount_dt$value)
  #View(full_FUN)
  addrow<-data.frame(origin=as.factor(40),dev=1,UplusN=UplusN[length(UplusN)],FCount=SettlCount_dt[SettlCount_dt$origin==40,]$value,Calendar=41)
  FUN_dat<-rbind(FUN_dat,addrow)
  train_1_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))==1,]
  train_2_FC<-FUN_dat[as.numeric(FUN_dat$Calendar)<=34&as.numeric(as.character(FUN_dat$origin))<=32&as.numeric(as.character(FUN_dat$origin))>=2,]
  train_3_FC<-FUN_dat[FUN_dat$dev==1&as.numeric(as.character(FUN_dat$origin))>=33&as.numeric(as.character(FUN_dat$origin))<=40,]
  train_4_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))>=33&as.numeric(as.character(FUN_dat$origin))<=34&as.numeric(FUN_dat$dev)==2,]
  train_5_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))==40&as.numeric(FUN_dat$dev)==1,]
  train_FC<-rbind(train_1_FC,train_2_FC,train_3_FC,train_4_FC,train_5_FC)
  train_FC<-train_FC[order(as.numeric(as.character(train_FC$origin))),]
  
  valid_1_FC<-FUN_dat[as.numeric(FUN_dat$Calendar)>=35&as.numeric(FUN_dat$Calendar)<=38&as.numeric(as.character(FUN_dat$origin))>=2&as.numeric(as.character(FUN_dat$origin))<=32,]
  valid_2_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))==33&as.numeric(FUN_dat$dev)>=3&as.numeric(FUN_dat$dev)<=5,]
  valid_3_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))==34&as.numeric(FUN_dat$dev)>=3&as.numeric(FUN_dat$dev)<=4,]
  valid_4_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))>=35&as.numeric(as.character(FUN_dat$origin))<=37&as.numeric(FUN_dat$dev)>=2&as.numeric(FUN_dat$dev)<=3,]
  valid_5_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))>=38&as.numeric(as.character(FUN_dat$origin))<=39&as.numeric(FUN_dat$dev)==2,]
  valid_6_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))==40&as.numeric(FUN_dat$dev)==1,]
  valid_FC<-rbind(valid_1_FC,valid_2_FC,valid_3_FC,valid_4_FC,valid_5_FC,valid_6_FC)
  valid_FC<-valid_FC[order(as.numeric(as.character(valid_FC$origin))),]
  
  # Construction of Test Data
  test_1_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))>=2&as.numeric(as.character(FUN_dat$origin))<=35&as.numeric(as.character(FUN_dat$Calendar))>=39&as.numeric(as.character(FUN_dat$Calendar))<=41,]
  test_2_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))==36&as.numeric(as.character(FUN_dat$dev))>=4&as.numeric(as.character(FUN_dat$dev))<=5,]
  test_3_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))==37&as.numeric(as.character(FUN_dat$dev))==4,]
  test_4_FC<-FUN_dat[as.numeric(as.character(FUN_dat$origin))==38&as.numeric(as.character(FUN_dat$dev))==3,]
  test_FC<-rbind(test_1_FC,test_2_FC,test_3_FC,test_4_FC)
  test_FC<-test_FC[order(as.numeric(as.character(test_FC$origin))),]
  aug_valid_FC<-rbind(valid_FC,test_FC)
  aug_valid_FC<-aug_valid_FC[order(as.numeric(as.character(aug_valid_FC$origin))),]
  #Alternatively, fit a ODP on finalization claims that depends on development periods
  odp_FC<-glm(FCount~factor(dev),data=train_FC,family=quasipoisson(link="log"))
  new_train$PPCF=new_train$value/train_F$value
  #Create a column for operation time
  new_train$OT=train_cumF$value/N_rep
  #Remove the cells with zero finalized claim count but with positive payments; 
  new_train_PPCF<-na.omit(new_train[new_train$PPCF!=Inf,])
  ODP_PPCF_train<-glm(PPCF~OT,family=quasipoisson(link="log"),data=new_train_PPCF)
  #Density for PPCI and PPCF
  dens_PPCI2<-cal_dens_PPCI(fit_ODP_ppci,N,aug_valid,2,40)
  aug_valid_FC$dev=as.factor(aug_valid_FC$dev)
  dens_PPCF2<-cal_PPCF_dens(odp_FC,ODP_PPCF_train,train_cumF,aug_valid_FC,aug_valid,N,2,40) 
  
  
  
  mu_PPCI2<-cal_mu_PPCI(fit_ODP_ppci,N,aug_valid,2,40)
  mu_PPCF2<-cal_PPCF_mu(odp_FC,ODP_PPCF_train,train_cumF,aug_valid_FC,aug_valid,N,2,40)
  
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
  glm_ODP_Ho_In<-glm(formula=value~factor(origin)+log(as.numeric(as.character(dev)))+as.numeric(as.character(dev)),family=quasipoisson(link="log"),data=in_sample,trace=FALSE)
  
  #Density for Hoerl curve
  out_sample_Ho<-out_sample
  out_sample_Ho$dev<-as.numeric(as.character(out_sample$dev))
  
  
  #Under Gamma Assumption
  in_sample_Ho<-in_sample
  in_sample_Ho$dev<-as.numeric(as.character(in_sample$dev))
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
  dens_ODPGLM_out_sample<-cal_dens_ODP(ODP_GLM_in,out_sample)
  dens_GAGLM_out_sample<-cal_dens_GA(Ga_optimTau,tau_Ga,out_sample)
  dens_LNGLM_out_sample<-cal_dens_LN(LN_optimTau,tau_LN,out_sample)
  dens_ZAGA_out_sample<-cal_dens_ZAGA(ZAGA_in,out_sample)
  dens_ZALN_out_sample<-cal_dens_ZALN(LN_in,out_sample)
  
  out_sample_numeric<-out_sample
  out_sample_numeric$origin=as.numeric(out_sample$origin)
  out_sample_numeric$dev=as.numeric(out_sample$dev)
  out_sample_numeric$Calendar=as.numeric(out_sample$Calendar)
  
  #Density for Hoerl curve
  out_sample_Ho<-out_sample
  out_sample_Ho$dev<-as.numeric(as.character(out_sample$dev))
  dens_ODPHo_out_sample<-cal_dens_ODP(glm_ODP_Ho_In,out_sample_Ho)
  dens_GaHo_out_sample<-cal_dens_GA(glm_Ga_Ho_In,tau_Ga,out_sample_Ho)
  dens_LNHo_out_sample<-cal_dens_LN(glm_LN_Ho_In,tau_LN,out_sample_Ho)
  
  #Density for Calendar Periods
  dens_ODPCal_out_sample<-cal_dens_ODP(glm_ODP_Cal_In,out_sample)
  dens_GaCal_out_sample<-cal_dens_GA(glm_Ga_Cal_In,tau_Ga,out_sample)
  dens_LNCal_out_sample<-cal_dens_LN(glm_LN_Cal_In,tau_LN,out_sample)
  
  #Density for Smoothing Spline
  dens_SpNormal_out_sample<-cal_dens_Normal(sp_Normal_In,out_sample_numeric)
  dens_SpGamma_out_sample<-cal_dens_GA(sp_Gamma_In,tau_Ga,out_sample_numeric)
  dens_SpLN_out_sample<-cal_dens_LN(sp_LN_In,tau_LN,out_sample_numeric)
  
  #Density for GAMLSS
  out_sample_mu<-out_sample
  out_sample_sigma<-out_sample_numeric
  #Density for GAMLSS
  dens_GaGAMLSS_out_sample<-cal_dens_GA_Gamlss(gamlss2_GA_In,tau_Ga,out_sample_mu,out_sample_sigma)
  dens_LNGAMLSS_out_sample<-cal_dens_LN_Gamlss(gamlss2_LN_In,tau_LN,out_sample_mu,out_sample_sigma)
  
  
  ##Mean for Basic Models
  mu_ODPGLM_out_sample<-cal_mu_ODP(ODP_GLM_in,out_sample)
  mu_GAGLM_out_sample<-ifelse(cal_mu_GA(Ga_optimTau,tau_Ga,out_sample)>tau_Ga,cal_mu_GA(Ga_optimTau,tau_Ga,out_sample)-tau_Ga,0)
  mu_LNGLM_out_sample<-ifelse(cal_mu_LN(LN_optimTau,tau_LN,out_sample)>tau_LN,cal_mu_LN(LN_optimTau,tau_LN,out_sample)-tau_LN,0)
  mu_ZAGA_out_sample<-cal_mu_ZAGA(ZAGA_in,out_sample)
  mu_ZALN_out_sample<-cal_mu_ZALN(LN_in,out_sample)
  
  
  #Mean for Hoerl curve
  out_sample_Ho<-out_sample
  out_sample_Ho$dev<-as.numeric(as.character(out_sample$dev))
  mu_ODPHo_out_sample<-cal_mu_ODP(glm_ODP_Ho_In,out_sample_Ho)
  mu_GaHo_out_sample<-ifelse(cal_mu_GA(glm_Ga_Ho_In,tau_Ga,out_sample_Ho)>tau_Ga,cal_mu_GA(glm_Ga_Ho_In,tau_Ga,out_sample_Ho)-tau_Ga,0)
  mu_LNHo_out_sample<-ifelse(cal_mu_LN(glm_LN_Ho_In,tau_LN,out_sample_Ho)>tau_LN,cal_mu_LN(glm_LN_Ho_In,tau_LN,out_sample_Ho)-tau_LN,0)
  
  #Mean for Calendar Periods
  mu_ODPCal_out_sample<-cal_mu_ODP(glm_ODP_Cal_In,out_sample)
  mu_GaCal_out_sample<-ifelse(cal_mu_GA(glm_Ga_Cal_In,tau_Ga,out_sample)>tau_Ga,cal_mu_GA(glm_Ga_Cal_In,tau_Ga,out_sample)-tau_Ga,0)
  mu_LNCal_out_sample<-ifelse(cal_mu_LN(glm_LN_Cal_In,tau_LN,out_sample)>tau_LN,cal_mu_LN(glm_LN_Cal_In,tau_LN,out_sample)-tau_LN,0)
  
  #Mean for Smoothing Spline
  mu_SpNormal_out_sample<-cal_mu_Normal(sp_Normal_In,out_sample_numeric)
  mu_SpGamma_out_sample<-ifelse(cal_mu_GA(sp_Gamma_In,tau_Ga,out_sample_numeric)>tau_Ga,cal_mu_GA(sp_Gamma_In,tau_Ga,out_sample_numeric)-tau_Ga,0)
  mu_SpLN_out_sample<-ifelse(cal_mu_LN(sp_LN_In,tau_LN,out_sample_numeric)>tau_LN,cal_mu_LN(sp_LN_In,tau_LN,out_sample_numeric)-tau_LN,0)
  
  #Mean for GAMLSS
  out_sample_mu<-out_sample
  out_sample_sigma<-out_sample_numeric
  #Mean for GAMLSS
  mu_GaGAMLSS_out_sample<-ifelse(cal_mu_GA_Gamlss(gamlss2_GA_In,tau_Ga,out_sample_mu,out_sample_sigma)>tau_Ga,cal_mu_GA_Gamlss(gamlss2_GA_In,tau_Ga,out_sample_mu,out_sample_sigma)-tau_Ga,0)
  mu_LNGAMLSS_out_sample<-ifelse(cal_mu_LN_Gamlss(gamlss2_LN_In,tau_LN,out_sample_mu,out_sample_sigma)>tau_LN,cal_mu_LN_Gamlss(gamlss2_LN_In,tau_LN,out_sample_mu,out_sample_sigma)-tau_LN,0)
  
  
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
  addrow<-data.frame(origin=as.factor(40),dev=1,UplusN=UplusN[length(UplusN)],FCount=SettlCount_dt[SettlCount_dt$origin==40,]$value,Calendar=41)
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
  
  
  #Density for PPCI and PPCF
  dens_PPCI_out_sample<-cal_dens_PPCI(fit_ODP_ppci_In,N_in,out_sample,2,40)
  FUN_dat_Out$dev=as.factor(FUN_dat_Out$dev)
  dens_PPCF_out_sample<-cal_PPCF_dens(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out,out_sample,N_in,2,40)
  
  #Calculate Predictive Mean for Out-of-Sample Data
  
  
  #Mean for PPCI and PPCF
  mu_PPCI_out_sample<-cal_mu_PPCI(fit_ODP_ppci_In,N_in,out_sample,2,40)
  FUN_dat_Out$dev=as.factor(FUN_dat_Out$dev)
  mu_PPCF_out_sample<-cal_PPCF_mu(odp_FC_In,ODP_PPCF_In,cum_F_dt,FUN_dat_Out,out_sample,N_in,2,40)
  
 
  
  #Generate a data.frame to store  predictive mean in out-of-sample data
  meta_mu_out_0_aug_new1<-data.frame(origin=out_sample$origin,dev=out_sample$dev,mu_ODP_GLM=mu_ODPGLM_out_sample,mu_GAGLM=mu_GAGLM_out_sample,mu_LNGLM=mu_LNGLM_out_sample,mu_ZAGA=mu_ZAGA_out_sample,mu_ZALN=mu_ZALN_out_sample,mu_ODPHo=mu_ODPHo_out_sample,mu_GaHo=mu_GaHo_out_sample,mu_LNHo=mu_LNHo_out_sample,mu_ODPCal=mu_ODPCal_out_sample,mu_GaCal=mu_GaCal_out_sample,mu_LNCal=mu_LNCal_out_sample,mu_SpNormal=mu_SpNormal_out_sample,mu_SpGamma=mu_SpGamma_out_sample,mu_SpLN=mu_SpLN_out_sample,mu_GaGAMLSS=mu_GaGAMLSS_out_sample,mu_LNGAMLSS=mu_LNGAMLSS_out_sample,mu_PPCI=mu_PPCI_out_sample,mu_PPCF=mu_PPCF_out_sample)
  
  #Generate a Level 1 Meta data
  dens_PPCF2[dens_PPCF2==0]=min(dens_PPCF2[dens_PPCF2!=0])
  
  meta_dens_augvalid<-data.frame(origin=aug_valid$origin,dev=aug_valid$dev,dens_ODP_GLM=dens_ODPGLM2,dens_GAGLM=dens_GAGLM2,dens_LNGLM=dens_LNGLM2,dens_ZAGA=dens_ZAGA2,dens_ZALN=dens_ZALN2,dens_ODPHo=dens_ODPHo2,dens_GaHo=dens_GaHo2,dens_LNHo=dens_LNHo2,dens_ODPCal=dens_ODPCal2,dens_GaCal=dens_GaCal2,dens_LNCal=dens_LNCal2,dens_SpNormal=dens_SpNormal2,dens_SpGamma=dens_SpGamma2,dens_SpLN=dens_SpLN2,dens_GaGAMLSS=dens_GaGAMLSS2,dens_LNGAMLSS=dens_LNGAMLSS2,dens_PPCI=dens_PPCI2,dens_PPCF=dens_PPCF2)
  #Construct 2 Validation sets to train model weights
  #Valid1(contains predicted density from AC2 TO AC10),Valid2(AC11-30),Valid3(AC31-40)
  meta_dens_out<-data.frame(origin=out_sample$origin,dev=out_sample$dev,dens_ODP_GLM=dens_ODPGLM_out_sample,dens_GAGLM=dens_GAGLM_out_sample,dens_LNGLM=dens_LNGLM_out_sample,dens_ZAGA=dens_ZAGA_out_sample,dens_ZALN=dens_ZALN_out_sample,dens_ODPHo=dens_ODPHo_out_sample,dens_GaHo=dens_GaHo_out_sample,dens_LNHo=dens_LNHo_out_sample,dens_ODPCal=dens_ODPCal_out_sample,dens_GaCal=dens_GaCal_out_sample,dens_LNCal=dens_LNCal_out_sample,dens_SpNormal=dens_SpNormal_out_sample,dens_SpGamma=dens_SpGamma_out_sample,dens_SpLN=dens_SpLN_out_sample,dens_GaGAMLSS=dens_GaGAMLSS_out_sample,dens_LNGAMLSS=dens_LNGAMLSS_out_sample,dens_PPCI=dens_PPCI_out_sample,dens_PPCF=dens_PPCF_out_sample)
  

  LogS_valid_IndEn2_new1<-log(meta_dens_augvalid[,-c(1,2)])
  ###Record the Best Performance model in the validation set
  index_best<-which.max(apply(LogS_valid_IndEn2_new1,MARGIN=2,FUN=mean))
  
  w_init_1_3<-rep(1/18,18)
  #Calculate the density for equally weighted ensemble
  dens_Eneq<-as.matrix(meta_dens_out[,-c(1,2)])%*%w_init_1_3
  meta_dens_out_1_aug_new1<-data.frame(origin=out_sample$origin,dev=out_sample$dev,dens_ODP_GLM=dens_ODPGLM_out_sample,dens_GAGLM=dens_GAGLM_out_sample,dens_LNGLM=dens_LNGLM_out_sample,dens_ZAGA=dens_ZAGA_out_sample,dens_ZALN=dens_ZALN_out_sample,dens_ODPHo=dens_ODPHo_out_sample,dens_GaHo=dens_GaHo_out_sample,dens_LNHo=dens_LNHo_out_sample,dens_ODPCal=dens_ODPCal_out_sample,dens_GaCal=dens_GaCal_out_sample,dens_LNCal=dens_LNCal_out_sample,dens_SpNormal=dens_SpNormal_out_sample,dens_SpGamma=dens_SpGamma_out_sample,dens_SpLN=dens_SpLN_out_sample,dens_GaGAMLSS=dens_GaGAMLSS_out_sample,dens_LNGAMLSS=dens_LNGAMLSS_out_sample,dens_PPCI=dens_PPCI_out_sample,dens_PPCF=dens_PPCF_out_sample,dens_EqEn=dens_Eneq)
  
  
  meta_mu_out_0_aug_new1_original<-meta_mu_out_0_aug_new1[,-c(1,2)]*10000
  out_mu_EW<-as.matrix(meta_mu_out_0_aug_new1_original)%*%w_init_1_3
  meta_mu_out_1_aug_new1_original<-cbind(meta_mu_out_0_aug_new1_original,out_mu_EW)
  out_mu_PB<-meta_mu_out_1_aug_new1_original[,index_best]
  meta_mu_augvalid<-data.frame(origin=aug_valid$origin,dev=aug_valid$dev,mu_ODP_GLM=mu_ODPGLM2*10000,mu_GAGLM=mu_GAGLM2*10000,mu_LNGLM=mu_LNGLM2*10000,mu_ZAGA=mu_ZAGA2*10000,mu_ZALN=mu_ZALN2*10000,mu_ODPHo=mu_ODPHo2*10000,mu_GaHo=mu_GaHo2*10000,mu_LNHo=mu_LNHo2*10000,mu_ODPCal=mu_ODPCal2*10000,mu_GaCal=mu_GaCal2*10000,mu_LNCal=mu_LNCal2*10000,mu_SpNormal=mu_SpNormal2*10000,mu_SpGamma=mu_SpGamma2*10000,mu_SpLN=mu_SpLN2*10000,mu_GaGAMLSS=mu_GaGAMLSS2*10000,mu_LNGAMLSS=mu_LNGAMLSS2*10000,mu_PPCI=mu_PPCI2*10000,mu_PPCF=mu_PPCF2*10000,observed=aug_valid$value*10000)
  
 
  future_claims[,D]<-out_sample$value
  
  index_best_validmod2[D]<-index_best
  ##Record the Predictive Density at Validation Set
  valid_dens_list[[D]]<- meta_dens_augvalid
  
  ##Record the Predictive Density at Out-of-Sample
  out_dens_list[[D]]<-meta_dens_out_1_aug_new1
  
  #Record the predicted mean in validation data
  valid_mu_list[[D]]<-meta_mu_augvalid
  
  #Record the predicted mean in out-of-sample data
  out_mu_list[[D]]<-meta_mu_out_1_aug_new1_original

  
}  