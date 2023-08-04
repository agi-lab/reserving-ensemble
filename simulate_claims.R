################################################################################
## Required variables to run module
################################################################################

data.dir = sprintf("./simulation/triangle_%s-data", tri.size)

if (!exists('n.sims')) {stop("Number of simulations 'n.sims' should be defined")}
if (!exists('tri.size')) {stop("Size of triangles 'tri.size' should be defined")}
if (!dir.exists(file.path(data.dir))){
  stop(sprintf("Please create subdirectory '%s' to store simulated data", data.dir))
}

print(sprintf('... Simulating %s size %s triangles...', n.sims, tri.size))

################################################################################
## Code for simulating claim data from SynthETIC
################################################################################

library("SynthETIC")
library("ChainLadder")

################################################################################
### Helper Functions
################################################################################

############## Claim size: Default power normal
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


notidel_param <- function(claim_size, occurrence_period){
    
    
    mean <- 2.47 * tri.size / 40
    cv <- 1.54
    shape <- get_Weibull_parameters(mean, cv)[1, ]
    scale <- get_Weibull_parameters(mean, cv)[2, ]
    c(shape = shape, scale = scale)
}


# Claim closure: Default Weibull

setldel_param <- function(claim_size, occurrence_period) {
    
    target_mean <- 11.74 * tri.size / 40
    
    # specify the target Weibull coefficient of variation
    target_cv <- 0.61
    
    c(shape = get_Weibull_parameters(target_mean, target_cv)[1, ],
      scale = get_Weibull_parameters(target_mean, target_cv)[2, ])
}


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

# Remove the Superimposed inflation
# 2) With respect to payment "time" (continuous scale)
# -> compounding by user-defined time unit

SI_payment <-  function(payment_time, claim_size){1}
SI_occurrence <- function(occurrence_time, claim_size){1}

data.to.claims.df <- function(data, na.rm = FALSE) {
    data<-as.data.frame(as.triangle(data), na.rm = na.rm)
    data$calendar <- as.numeric(data$origin)+as.numeric(data$dev)
    data[order(as.numeric(data$dev), as.numeric(data$origin)),]
    return(data)
}

################################################################################
### Main Simulation Code
################################################################################

infl <- 0.02
    
# SPLICE assumptions
ref_claim <- 200000
if (!tri.size %in% c(10, 20, 40)){stop("Only triangles of size 10, 20 and 40 are currently supported")}
time_unit <- 10/tri.size
years <- 10
I <- years / time_unit
E <- c(rep(12000, I)) # effective annual exposure rates
lambda <- c(rep(0.03, I))

# Package-wise global parameters
set_parameters(ref_claim = ref_claim, time_unit = time_unit)

simulate_splice <- function(seed) {
    # Creating claims simulations
    
    set.seed(seed)
    
    # Claim occurrence
    # Example 1.1: Claim occurrence: Constant exposure and frequency
    # Number of claims occurring for each period i
    n_vector <- claim_frequency(I = I, E = E, freq = lambda)
    
    # Occurrence time of each claim r, for each period i
    occurrence_times <- claim_occurrence(frequency_vector = n_vector)
    
    claim_sizes <- claim_size(frequency_vector = n_vector, simfun = S_df, type = "p", range = c(0, 1e24))
    
    ## Output
    # Notification Delay
    # By changing claim notification to transformed Gamma, the first development periods contain no zero value claims 
    notidel <- claim_notification(n_vector, claim_size(n_vector), 
                                  paramfun = notidel_param)
    
    # simulate the settlement delays from the Weibull with parameters above
    setldel <- claim_closure(n_vector, claim_sizes, rfun = rweibull, paramfun = setldel_param)
    
    # Claim partial payment: Default mixture distribution
    no_payments <- claim_payment_no(n_vector, claim_sizes)
    
    # Claim Partial payment (Without Inflation): Default Distribution
    payment_sizes <- claim_payment_size(n_vector, claim_sizes, no_payments)
    
    payment_delays <- claim_payment_delay(
        n_vector, claim_sizes, no_payments, setldel,
        rfun = r_pmtdel, paramfun = param_pmtdel,
        occurrence_period = rep(1:I, times = n_vector)
    )
    
    # payment times on a continuous time scale
    payment_times <- claim_payment_time(n_vector, occurrence_times, notidel, payment_delays)
    # payment times in periods
    payment_periods <- claim_payment_time(n_vector, occurrence_times, notidel, payment_delays,discrete = TRUE)
    
    # Claim inflation
    # Base inflation: a vector of quarterly rates
    # In this data set we set base inflation to be at 2% p.a. constant for both past and future
    demo_rate <- (1 + infl)^(1/4) - 1
    base_inflation_past <- 
        base_inflation_future <- rep(demo_rate, times = 40)
    base_inflation_vector <- c(base_inflation_past, base_inflation_future)
    
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
    
   
    
    agg_sqrt_withTail <- claim_output(
        n_vector,
        payment_times, 
        payment_inflated, 
        incremental = TRUE,
        future = TRUE,
        adjust = FALSE
    )
    
  
    agg_sqrt<-agg_sqrt_withTail[,-ncol(agg_sqrt_withTail)]
    colnames(agg_sqrt)<-c(1:tri.size)
    rownames(agg_sqrt)<-c(1:tri.size)
    
    agg_tri_withInflation <- claim_output(
        n_vector, 
        payment_times, 
        payment_inflated,
        incremental = TRUE,
        future = FALSE, 
        # adjust = FALSE  
    )
  
    colnames(agg_tri_withInflation)<-c(1:tri.size)
    rownames(agg_tri_withInflation)<-c(1:tri.size)
    
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
    transaction_dataset$calendar=transaction_dataset$ReportingPeriod+transaction_dataset$occurrence_period
    transaction_dataset$SettlementPeriod=ceiling(transaction_dataset$setldel)
    
    return (list(
        agg_sqrt = agg_sqrt, 
        agg_tri_withInflation = agg_tri_withInflation,
        all_claims = all_claims,
        transaction_dataset = transaction_dataset
    ))
}

################################################################################
### Main Simulation Code
################################################################################

for (D in 1:n.sims){
    
    claims.sim <- simulate_splice(20200130+D)
    
    agg_sqrt <- claims.sim$agg_sqrt
    agg_tri_withInflation <- claims.sim$agg_tri_withInflation
    all_claims <- claims.sim$all_claims
    transaction_dataset <- claims.sim$transaction_dataset
    
    #################################################
    ### Generate Aggregate Claims Data
    ### Full Triangle
    #################################################
    
    full_dat <- data.to.claims.df(agg_sqrt)
    full_dat$value=full_dat$value/10000
    full_dat <- df.to.factor(full_dat)
  
    #################################################
    ### Generate Aggregate Claims Data
    ### Past Triangle
    #################################################
    
    # convert the aggregate data generated by SynthETIC into triangle format: value, development year and accident year
    # only contains the past observations (the upper triangle) as the lower triangle are coded as N.A
    dt_withInflation_past<-data.to.claims.df(agg_tri_withInflation, na.rm=TRUE)
    dt_withInflation_past$value=dt_withInflation_past$value/10000
    dt_withInflation_past <- df.to.factor(dt_withInflation_past)
    
    #################################################
    ### Generate Notification Count Data
    #################################################
    
    # Past Triangle
    matri_NotiCount<-matrix(NA,nrow=tri.size,ncol=tri.size)
    # Full triangle
    matri_NotiCount_sqrt<-matrix(NA,nrow=tri.size,ncol=tri.size)
    
    for (i in 1:tri.size){
        for (j in 1:tri.size){
            if (i+j<=tri.size+1){matri_NotiCount[i,j]=nrow(transaction_dataset[transaction_dataset$occurrence_period==i&transaction_dataset$ReportingPeriod==j&transaction_dataset$calendar<=tri.size+1,])}
            else{matri_NotiCount[i,j]=NA}
            matri_NotiCount_sqrt[i,j]=nrow(transaction_dataset[transaction_dataset$occurrence_period==i&transaction_dataset$ReportingPeriod==j,])
        }
    }
    
    # Past Triangle
    NotiCount_dt<-data.to.claims.df(matri_NotiCount, na.rm = TRUE)
    
    # Full Triangle
    NotiCount_full_dt<-data.to.claims.df(matri_NotiCount_sqrt)
    
    ### Generate Cumulative Results
    # Past Triangle
    cum_Count_dt<-data.to.claims.df(incr2cum(as.triangle(matri_NotiCount)), na.rm = TRUE)
    
    # Full Triangle
    cum_Count_dt_full<-data.to.claims.df(incr2cum(as.triangle(matri_NotiCount_sqrt)))
    
    #################################################
    ### Generate Settlement Count Data
    #################################################
    
    # Past Triangle
    matri_SettlementCount<-matrix(NA,nrow=tri.size,ncol=tri.size)
    # Full Triangle
    sqrt_SettlementCount<-matrix(NA,nrow=tri.size,ncol=tri.size)
    
    for (i in 1:tri.size){
        for (j in 1:tri.size){
            if (i+j<=tri.size+1){matri_SettlementCount[i,j]=nrow(transaction_dataset[transaction_dataset$occurrence_period==i&transaction_dataset$SettlementPeriod==j&transaction_dataset$calendar<=tri.size+1,])}
            else{matri_SettlementCount[i,j]=NA}
          
            sqrt_SettlementCount[i,j]=nrow(transaction_dataset[transaction_dataset$occurrence_period==i&transaction_dataset$SettlementPeriod==j,])
        }
    }
    
    
    # Past Triangle
    SettlCount_dt <- data.to.claims.df(matri_SettlementCount, na.rm=TRUE)
    
    # Full Triangle
    SettlCount_full<-data.to.claims.df(sqrt_SettlementCount)
    
    ### Generate Cumulative Results
    # Past Triangle
    cum_F_dt<-data.to.claims.df(incr2cum(as.triangle(matri_SettlementCount)), na.rm=TRUE)
    
    # Full Triangle
    cum_F_dt_full<-data.to.claims.df(incr2cum(as.triangle(sqrt_SettlementCount)))
    
    #################################################
    ### Generate Unfinalised Claims Count Data
    #################################################
    
    # Past Triangle
    cum_U_dt<-data.frame(origin=cum_F_dt$origin,dev=cum_F_dt$dev,value=cum_Count_dt$value-cum_F_dt$value)
    
    # Full Triangle
    cum_U_dt_full<-data.frame(origin=cum_F_dt_full$origin,dev=cum_F_dt_full$dev,calendar=cum_F_dt_full$calendar,value=cum_Count_dt_full$value-cum_F_dt_full$value)
    
    
    UplusN<-matrix(NA,nrow=tri.size,ncol=tri.size)
    UplusN_full<-matrix(NA,nrow=tri.size,ncol=tri.size)
    for (i in 1:tri.size){
        for (j in 1:tri.size){
            if (i+j<=tri.size+1){UplusN[i, j] <- cum_U_dt[cum_U_dt$origin==i&cum_U_dt$dev==j, 'value'] + NotiCount_dt[NotiCount_dt$origin==i&NotiCount_dt$dev==j, 'value']}
            else{UplusN[i, j] <- NA}
          
            UplusN_full[i, j] <- cum_U_dt_full[cum_U_dt_full$origin==i&cum_U_dt_full$dev==j, 'value'] + NotiCount_full_dt[NotiCount_full_dt$origin==i&NotiCount_full_dt$dev==j, 'value']
        }
    }
    
    # Past Triangle
    UplusN <- data.to.claims.df(UplusN, na.rm=TRUE)
    
    # Full Triangle
    UplusN_full<-data.to.claims.df(UplusN_full)
    
    # Write CSV
    past_tri_data <- data.frame(
        origin=dt_withInflation_past$origin,
        dev=dt_withInflation_past$dev,
        calendar=dt_withInflation_past$calendar,
        aggregate_claims=dt_withInflation_past$value,
        notif_count=NotiCount_dt$value, 
        cum_notif_count=cum_Count_dt$value, 
        settle_count=SettlCount_dt$value, 
        cum_settle_count=cum_F_dt$value, 
        unfinal_count=UplusN$value
    )
    full_tri_data <- data.frame(
        origin=full_dat$origin,
        dev=full_dat$dev,
        calendar=full_dat$calendar,
        aggregate_claims=full_dat$value,
        notif_count=NotiCount_full_dt$value, 
        cum_notif_count=cum_Count_dt_full$value, 
        settle_count=SettlCount_full$value, 
        cum_settle_count=cum_F_dt_full$value, 
        unfinal_count=UplusN_full$value
      )
    
    write.csv(
        past_tri_data,
        sprintf("%s/sim%s-past-data.csv", data.dir, D),
        row.names = FALSE
    )
    
    
    write.csv(
        full_tri_data,
        sprintf("%s/sim%s-full-data.csv", data.dir, D),
        row.names = FALSE
    )
    
}
