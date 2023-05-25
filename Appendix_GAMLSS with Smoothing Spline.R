dens_GaGAMLSSSP_out_sample<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
dens_LNGAMLSSSP_out_sample<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
mu_GaGAMLSSSP_out_sample<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
mu_LNGAMLSSSP_out_sample<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)


dens_LNGAMLSSSP_in_sample<-matrix(NA,nrow=nrow(in_sample),ncol=ntri)
dens_LNSP_in_sample<-matrix(NA,nrow=nrow(in_sample),ncol=ntri)

#ntri<-100
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
  
  in_sample<-in_sample[order(as.numeric(as.character(in_sample$origin))),]
  
  in_sample_numeric<-in_sample
  in_sample_numeric$origin=as.numeric(as.character(in_sample$origin))
  in_sample_numeric$dev=as.numeric(as.character(in_sample$dev))
  in_sample_numeric$Calendar=as.numeric(as.character(in_sample$Calendar))
  
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
  
  ###GAMLSS
  #gamlss2_GA_SP<-gamlss(formula=(value+tau_Ga)~cs(as.numeric(as.character(dev)))+cs(as.numeric(as.character(origin))),data=in_sample,sigma.formula=~cs(as.numeric(as.character(dev))),family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
  gamlss2_LN_SP<-gamlss(formula=(value+tau_LN)~scs(as.numeric(dev))+scs(as.numeric(origin)),data=in_sample_numeric,sigma.formula=~scs(as.numeric(as.character(dev))),family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  sp_LN_In<-gamlss(formula=(value+tau_LN)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=in_sample_numeric,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
  
  ###Predictive density for in-sample data
  
  in_sample_mu<-in_sample
  in_sample_sigma<-in_sample_numeric
  
  dens_LNGAMLSSSP_in_sample[,D]<-cal_dens_LN_Gamlss(gamlss2_LN_SP,tau_LN,in_sample_numeric,in_sample_numeric)
  dens_LNSP_in_sample[,D]<-cal_dens_LN(sp_LN_In,tau_LN,in_sample_numeric)
  
  ###Predictive density for out-sample data
  out_sample<-out_sample[order(as.numeric(as.character(out_sample$origin))),]
  
  
  out_sample_numeric<-out_sample
  out_sample_numeric$origin=as.numeric(out_sample$origin)
  out_sample_numeric$dev=as.numeric(out_sample$dev)
  out_sample_numeric$Calendar=as.numeric(out_sample$Calendar)
  
  #Density for GAMLSS
  out_sample_mu<-out_sample
  out_sample_sigma<-out_sample_numeric
  
  
  #Density for GAMLSS
  
  #dens_GaGAMLSSSP_out_sample[,D]<-cal_dens_GA_Gamlss(gamlss2_GA_SP,tau_Ga,out_sample_mu,out_sample_sigma)
  #dens_LNGAMLSSSP_out_sample[,D]<-cal_dens_LN_Gamlss(gamlss2_LN_SP,tau_LN,out_sample_mu,out_sample_sigma)
 
 
  #Mean for GAMLSS
  
  #out_sample_mu<-out_sample
  #out_sample_sigma<-out_sample_numeric
  #Mean for GAMLSS
  #mu_GaGAMLSSSP_out_sample[,D]<-ifelse(cal_mu_GA_Gamlss(gamlss2_GA_SP,tau_Ga,out_sample_mu,out_sample_sigma)>tau_Ga,cal_mu_GA_Gamlss(gamlss2_GA_SP,tau_Ga,out_sample_mu,out_sample_sigma)-tau_Ga,0)
  #mu_LNGAMLSSSP_out_sample[,D]<-ifelse(cal_mu_LN_Gamlss(gamlss2_LN_SP,tau_LN,out_sample_mu,out_sample_sigma)>tau_LN,cal_mu_LN_Gamlss(gamlss2_LN_SP,tau_LN,out_sample_mu,out_sample_sigma)-tau_LN,0)
}



# save(dens_LNGAMLSSSP_in_sample,file="dens_LNGAMLSSSP_in_sample")
# save(dens_LNSP_in_sample,file='dens_LNSP_in_sample')
# save(dens_LNGAMLSSSP_out_sample,file="LNGAMLSSP_out_dens")
# save(mu_LNGAMLSSSP_out_sample,file="LNGAMLSSP_out_mu")