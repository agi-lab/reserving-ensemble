
if (!exists('tri.size')) {stop("Size of triangles 'tri.size' should be defined")}

################################################################################
### Fitting Model on training sets
### For each model, we:
### 1) Fit the model
### 2) Calculate Densities
### 3) Calculate CDF
################################################################################
# Note that the training model uses:
# training set for fitting, and valid+test set for dens/cdf/mu
#
# The full data model uses:
# in-sample (train+valid), and test set for dens/cdf/mu
# 
# The columns required for analysis are:
# - origin
# - dev
# - calendar
# - aggregate_claims
# - notif_count
# - cum_notif_count
# - settle_count
# - cum_settle_count
# - unfinal_count
################################################################################
# There are 18 models considered in this file for triangle size 40
################################################################################

fit_all_component_models_40 <- function(train.data, test.data) {
    sink(nullfile())
    
    train.data.numeric <- df.to.numeric(train.data)
    test.data.numeric <- df.to.numeric(test.data)
    y.test <- test.data.numeric$aggregate_claims
    
    ###Fitting of all the models in the training sets
    ####Basic Models 
    #ODP
    ODP_GLM<-glm(formula=aggregate_claims~factor(origin)+factor(dev),family=quasipoisson(link="log"),data=train.data.numeric)
    param_ODPGLM<-fit_param_ODP(ODP_GLM,test.data)
    dens_ODPGLM<-cal_dens_ODP(y.test, param_ODPGLM)
    CDF_ODPGLM<-cal_CDF_ODP(y.test, param_ODPGLM)
    mu_ODPGLM<-cal_mu_ODP(param_ODPGLM)
    
    #Gamma 
    #tau_Ga<-5
    Ga_optimTau<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(origin)+factor(dev),data=train.data.numeric,family=GA(mu.link="log", sigma.link ="log"))
    param_GAGLM<-fit_param_GA(Ga_optimTau,train.data,test.data,tau_Ga)
    dens_GAGLM<-cal_dens_GA(y.test, param_GAGLM)
    CDF_GAGLM<-cal_CDF_GA(y.test, param_GAGLM)
    mu_GAGLM<-pmax(cal_mu_GA(param_GAGLM)-tau_Ga,0)
    
    #Log-Normal
    #tau_LN<-5
    LN_optimTau<-gamlss(formula=(aggregate_claims+tau_LN)~factor(origin)+factor(dev),data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"))
    param_LNGLM<-fit_param_LN(LN_optimTau,train.data,test.data,tau_LN)
    dens_LNGLM<-cal_dens_LN(y.test, param_LNGLM)
    CDF_LNGLM<-cal_CDF_LN(y.test, param_LNGLM)
    mu_LNGLM<-pmax(cal_mu_LN(param_LNGLM)-tau_LN,0)
    
    #ZAGA
    gamma_1<-gamlss(formula=aggregate_claims~factor(origin)+factor(dev),nu.formula=~as.numeric(as.character(dev)),data=train.data.numeric,family=ZAGA(mu.link="log",sigma.link = "log", nu.link = "logit"))
    param_ZAGA<-fit_param_ZAGA(gamma_1,train.data,test.data)
    dens_ZAGA<-cal_dens_ZAGA(y.test, param_ZAGA)
    CDF_ZAGA<-cal_CDF_ZAGA(y.test, param_ZAGA)
    mu_ZAGA<-cal_mu_ZAGA(param_ZAGA)
    
    #ZALN
    LN_1<-gamlssZadj(y=aggregate_claims,mu.formula = ~factor(origin)+factor(dev),xi0.formula=~as.numeric(as.character(dev)),data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"))
    param_ZALN<-fit_param_ZALN(LN_1,train.data,test.data)
    dens_ZALN<-cal_dens_ZALN(y.test, param_ZALN)
    CDF_ZALN<-cal_CDF_ZALN(y.test, param_ZALN)
    mu_ZALN<-cal_mu_ZALN(param_ZALN)
    
    #### Models With Hoerl Curve
    #Under ODP Assumption
    #R return warning for rank deficient fit for Hoerl Curve ODP
    glm_ODP_Ho_tr1<-glm(formula=aggregate_claims~factor(origin)+log(dev)+dev,family=quasipoisson(link="log"),data=train.data.numeric)
    param_ODPHo<-fit_param_ODP(glm_ODP_Ho_tr1,test.data.numeric)
    dens_ODPHo<-cal_dens_ODP(y.test, param_ODPHo)
    CDF_ODPHo<-cal_CDF_ODP(y.test, param_ODPHo)
    mu_ODPHo<-cal_mu_ODP(param_ODPHo)
    
    #Under Gamma Assumption
    glm_Ga_Ho_tr1<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(origin)+log(dev)+dev,data=train.data.numeric,family=GA(mu.link="log", sigma.link ="log"))
    param_GaHo<-fit_param_GA(glm_Ga_Ho_tr1,train.data.numeric,test.data.numeric,tau_Ga)
    dens_GaHo<-cal_dens_GA(y.test, param_GaHo)
    CDF_GaHo<-cal_CDF_GA(y.test, param_GaHo)
    mu_GaHo<-pmax(cal_mu_GA(param_GaHo)-tau_Ga,0)
    
    #Under Log-Normal Assumption
    glm_LN_Ho_tr1<-gamlss(formula=(aggregate_claims+tau_LN)~factor(origin)+log(dev)+dev,data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"))
    param_LNHo<-fit_param_LN(glm_LN_Ho_tr1,train.data.numeric,test.data.numeric,tau_LN)
    dens_LNHo<-cal_dens_LN(y.test, param_LNHo)
    CDF_LNHo<-cal_CDF_LN(y.test, param_LNHo)
    mu_LNHo<-pmax(cal_mu_LN(param_LNHo)-tau_LN,0)
    
    #### Models With Calendar Periods
    #Under ODP Assumption
    glm_ODP_Cal_tr1<-glm(formula=aggregate_claims~factor(dev)+calendar,family=quasipoisson(link="log"),data=train.data.numeric)
    param_ODPCal<-fit_param_ODP(glm_ODP_Cal_tr1,test.data.numeric)
    dens_ODPCal<-cal_dens_ODP(y.test, param_ODPCal)
    CDF_ODPCal<-cal_CDF_ODP(y.test, param_ODPCal)
    mu_ODPCal<-cal_mu_ODP(param_ODPCal)
    
    #Under Gamma Assumption
    glm_Ga_Cal_tr1<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(dev)+calendar,data=train.data.numeric,family=GA(mu.link="log", sigma.link ="log"))
    param_GaCal<-fit_param_GA(glm_Ga_Cal_tr1,train.data.numeric,test.data.numeric,tau_Ga)
    dens_GaCal<-cal_dens_GA(y.test, param_GaCal)
    CDF_GaCal<-cal_CDF_GA(y.test, param_GaCal)
    mu_GaCal<-pmax(cal_mu_GA(param_GaCal)-tau_Ga,0)
    
    #Under Log-Normal Assumption
    glm_LN_Cal_tr1<-gamlss(formula=(aggregate_claims+tau_LN)~factor(dev)+calendar,data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"))
    param_LNCal<-fit_param_LN(glm_LN_Cal_tr1,train.data.numeric,test.data.numeric,tau_LN)
    dens_LNCal<-cal_dens_LN(y.test, param_LNCal)
    CDF_LNCal<-cal_CDF_LN(y.test, param_LNCal)
    mu_LNCal<-pmax(cal_mu_LN(param_LNCal)-tau_LN,0)
    
    ####Smoothing Spline
    #Under Normal Assumption
    sp_Normal<-gamlss(formula=aggregate_claims~scs(origin)+scs(dev),data=train.data.numeric,family=NO(),trace=FALSE)
    sp_Normal <- refit(sp_Normal)
    param_SpNormal<-fit_param_NO(sp_Normal,train.data.numeric,test.data.numeric)
    dens_SpNormal<-cal_dens_Normal(y.test, param_SpNormal)
    CDF_SpNormal<-cal_CDF_Normal(y.test, param_SpNormal)
    mu_SpNormal<-cal_mu_Normal(param_SpNormal)
    
    #Under Gamma Assumption
    sp_Gamma<-gamlss(formula=(aggregate_claims+tau_Ga)~scs(origin)+scs(dev),data=train.data.numeric,family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
    sp_Gamma <- refit(sp_Gamma)
    param_SpGamma<-fit_param_GA(sp_Gamma,train.data.numeric,test.data.numeric,tau_Ga)
    dens_SpGamma<-cal_dens_GA(y.test, param_SpGamma)
    CDF_SpGamma<-cal_CDF_GA(y.test, param_SpGamma)
    mu_SpGamma<-pmax(cal_mu_GA(param_SpGamma)-tau_Ga,0)
    
    #Under Log-Normal Assumption
    sp_LN<-gamlss(formula=(aggregate_claims+tau_LN)~scs(origin)+scs(dev),data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
    sp_LN <- refit(sp_LN)
    param_SpLN<-fit_param_LN(sp_LN,train.data.numeric,test.data.numeric,tau_LN)
    dens_SpLN<-cal_dens_LN(y.test, param_SpLN)
    CDF_SpLN<-cal_CDF_LN(y.test, param_SpLN)
    mu_SpLN<-pmax(cal_mu_LN(param_SpLN)-tau_LN,0)
    
    ####GAMLSS
    #GAMLSS 2: Smooth Effects on the predictor for sigma term
    #Under Gamma Assumption
    gamlss_GA<-gamlss(formula=(aggregate_claims+tau_Ga)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=train.data.numeric,sigma.formula=~cs(as.numeric(as.character(dev))),family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
    gamlss_GA <- refit(gamlss_GA)
    param_GaGAMLSS<-fit_param_GAGamlss(gamlss_GA,train.data.numeric,test.data.numeric,tau_Ga)
    dens_GaGAMLSS<-cal_dens_GA_Gamlss(y.test, param_GaGAMLSS)
    CDF_GaGAMLSS<-cal_CDF_GA_Gamlss(y.test, param_GaGAMLSS)
    mu_GaGAMLSS<-pmax(cal_mu_GA_Gamlss(param_GaGAMLSS)-tau_Ga,0)
    
    #Under Log-Normal Assumption
    gamlss_LN<-gamlss(formula=(aggregate_claims+tau_LN)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=train.data.numeric,sigma.formula=~cs(as.numeric(as.character(dev))),family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
    gamlss_LN <- refit(gamlss_LN)
    param_LNGAMLSS<-fit_param_LNGamlss(gamlss_LN,train.data.numeric,test.data.numeric,tau_LN)
    dens_LNGAMLSS<-cal_dens_LN_Gamlss(y.test, param_LNGAMLSS)
    CDF_LNGAMLSS<-cal_CDF_LN_Gamlss(y.test, param_LNGAMLSS)
    mu_LNGAMLSS<-pmax(cal_mu_LN_Gamlss(param_LNGAMLSS)-tau_LN,0)
    
    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train.data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train.data[train.data$origin==1,]$notif_count)
    for (i in 2:tri.size){
      N[i]<- sum(train.data[train.data$origin==i,]$notif_count) + 
             sum(round(predict(fit_nc,newdata=test.data[test.data$origin==i,],type="response"),0)) 
    }
    
    ###PPCI
    ##Fit a ODP/Chain-Ladder model to predict total claim notification count at each accident period 
    PPCI <- calc_PPCI(N, train.data)
    ##Fitting a ODP model on PPCI
    fit_ODP_ppci<-glm(PPCI~factor(dev),family=quasipoisson(link="log"),data=train.data.numeric)
    param_PPCI<-fit_param_PPCI(fit_ODP_ppci,N,test.data)
    dens_PPCI<-cal_dens_PPCI(y.test, param_PPCI)
    CDF_PPCI<-cal_CDF_PPCI(y.test, param_PPCI)
    mu_PPCI<-cal_mu_PPCI(param_PPCI)
    
    ###PPCF
    #Alternatively, fit a ODP on finalization claims that depends on development periods
    odp_FC<-glm(settle_count~factor(dev),data=train.data,family=quasipoisson(link="log"))
    PPCF_data <- calc_PPCF(N, train.data) 
    ODP_PPCF<-glm(PPCF~OT,family=quasipoisson(link="log"),data=PPCF_data)
    param_PPCF <- fit_param_PPCF(odp_FC,ODP_PPCF, N, train.data,test.data) 
    dens_PPCF <- cal_dens_PPCF(y.test, param_PPCF)
    CDF_PPCF <- cal_CDF_PPCF(y.test, param_PPCF)
    mu_PPCF<-cal_mu_PPCF(param_PPCF)
    
    ################################################################################
    ################################################################################
    ################################################################################
    
    component_models <- list(
        fit_ODP_GLM = ODP_GLM,
        fit_GAGLM = Ga_optimTau,
        fit_LNGLM = LN_optimTau,
        fit_ZAGA = gamma_1,
        fit_ZALN = LN_1,
        fit_ODPHo = glm_ODP_Ho_tr1,
        fit_GaHo = glm_Ga_Ho_tr1,
        fit_LNHo = glm_LN_Ho_tr1,
        fit_ODPCal = glm_ODP_Cal_tr1,
        fit_GaCal = glm_Ga_Cal_tr1,
        fit_LNCal = glm_LN_Cal_tr1,
        fit_SpNormal = sp_Normal,
        fit_SpGamma = sp_Gamma,
        fit_SpLN = sp_LN,
        fit_GaGAMLSS = gamlss_GA,
        fit_LNGAMLSS = gamlss_LN,
        fit_PPCI = list(model = fit_ODP_ppci, N = N),
        fit_PPCF = list(model_subCount = odp_FC, model_subPayments = ODP_PPCF, N=N)
    )
    
    # Predictive Mean  
    meta_mu<-data.frame(
      origin=test.data$origin,
      dev=test.data$dev,
      mu_ODP_GLM=mu_ODPGLM,
      mu_GAGLM=mu_GAGLM,
      mu_LNGLM=mu_LNGLM,
      mu_ZAGA=mu_ZAGA,
      mu_ZALN=mu_ZALN,
      mu_ODPHo=mu_ODPHo,
      mu_GaHo=mu_GaHo,
      mu_LNHo=mu_LNHo,
      mu_ODPCal=mu_ODPCal,
      mu_GaCal=mu_GaCal,
      mu_LNCal=mu_LNCal,
      mu_SpNormal=mu_SpNormal,
      mu_SpGamma=mu_SpGamma,
      mu_SpLN=mu_SpLN,
      mu_GaGAMLSS=mu_GaGAMLSS,
      mu_LNGAMLSS=mu_LNGAMLSS,
      mu_PPCI=mu_PPCI,
      mu_PPCF=mu_PPCF
    )
    
    # Density
    meta_dens<-data.frame(
      origin=test.data$origin,
      dev=test.data$dev,
      dens_ODP_GLM=dens_ODPGLM,
      dens_GAGLM=dens_GAGLM,
      dens_LNGLM=dens_LNGLM,
      dens_ZAGA=dens_ZAGA,
      dens_ZALN=dens_ZALN,
      dens_ODPHo=dens_ODPHo,
      dens_GaHo=dens_GaHo,
      dens_LNHo=dens_LNHo,
      dens_ODPCal=dens_ODPCal,
      dens_GaCal=dens_GaCal,
      dens_LNCal=dens_LNCal,
      dens_SpNormal=dens_SpNormal,
      dens_SpGamma=dens_SpGamma,
      dens_SpLN=dens_SpLN,
      dens_GaGAMLSS=dens_GaGAMLSS,
      dens_LNGAMLSS=dens_LNGAMLSS,
      dens_PPCI=dens_PPCI,
      dens_PPCF=dens_PPCF
    )
    
    #CDF
    meta_CDF <-data.frame(
      origin=test.data$origin,
      dev=test.data$dev,
      CDF_ODP_GLM=CDF_ODPGLM,
      CDF_GAGLM=CDF_GAGLM,
      CDF_LNGLM=CDF_LNGLM,
      CDF_ZAGA=CDF_ZAGA,
      CDF_ZALN=CDF_ZALN,
      CDF_ODPHo=CDF_ODPHo,
      CDF_GaHo=CDF_GaHo,
      CDF_LNHo=CDF_LNHo,
      CDF_ODPCal=CDF_ODPCal,
      CDF_GaCal=CDF_GaCal,
      CDF_LNCal=CDF_LNCal,
      CDF_SpNormal=CDF_SpNormal,
      CDF_SpGamma=CDF_SpGamma,
      CDF_SpLN=CDF_SpLN,
      CDF_GaGAMLSS=CDF_GaGAMLSS,
      CDF_LNGAMLSS=CDF_LNGAMLSS,
      CDF_PPCI=CDF_PPCI,
      CDF_PPCF=CDF_PPCF
    )
    
    sink()
    return (list(component_models=component_models, meta_mu=meta_mu, meta_dens=meta_dens, meta_CDF=meta_CDF))
}

## Define a fit_all_component_models function for smaller triangles: 
### We exclude ZALN and ZAGA as there are no zero incremental claims for 10x10 and 20x20 

fit_all_component_models_20 <- function(train.data, test.data) {
    sink(nullfile())
    
    train.data.numeric <- df.to.numeric(train.data)
    test.data.numeric <- df.to.numeric(test.data)
    y.test <- test.data.numeric$aggregate_claims
    
    ###Fitting of all the models in the training sets
    ####Basic Models 
    #ODP
    ODP_GLM<-glm(formula=aggregate_claims~factor(origin)+factor(dev),family=quasipoisson(link="log"),data=train.data.numeric)
    param_ODPGLM<-fit_param_ODP(ODP_GLM,test.data)
    dens_ODPGLM<-cal_dens_ODP(y.test, param_ODPGLM)
    CDF_ODPGLM<-cal_CDF_ODP(y.test, param_ODPGLM)
    mu_ODPGLM<-cal_mu_ODP(param_ODPGLM)
    
    #Gamma 
    #tau_Ga<-5
    Ga_optimTau<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(origin)+factor(dev),data=train.data.numeric,family=GA(mu.link="log", sigma.link ="log"))
    param_GAGLM<-fit_param_GA(Ga_optimTau,train.data,test.data,tau_Ga)
    dens_GAGLM<-cal_dens_GA(y.test, param_GAGLM)
    CDF_GAGLM<-cal_CDF_GA(y.test, param_GAGLM)
    mu_GAGLM<-pmax(cal_mu_GA(param_GAGLM)-tau_Ga,0)
    
    #Log-Normal
    #tau_LN<-5
    LN_optimTau<-gamlss(formula=(aggregate_claims+tau_LN)~factor(origin)+factor(dev),data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"))
    param_LNGLM<-fit_param_LN(LN_optimTau,train.data,test.data,tau_LN)
    dens_LNGLM<-cal_dens_LN(y.test, param_LNGLM)
    CDF_LNGLM<-cal_CDF_LN(y.test, param_LNGLM)
    mu_LNGLM<-pmax(cal_mu_LN(param_LNGLM)-tau_LN,0)
    
    #### Models With Hoerl Curve
    #Under ODP Assumption
    #R return warning for rank deficient fit for Hoerl Curve ODP
    glm_ODP_Ho_tr1<-glm(formula=aggregate_claims~factor(origin)+log(dev)+dev,family=quasipoisson(link="log"),data=train.data.numeric)
    param_ODPHo<-fit_param_ODP(glm_ODP_Ho_tr1,test.data.numeric)
    dens_ODPHo<-cal_dens_ODP(y.test, param_ODPHo)
    CDF_ODPHo<-cal_CDF_ODP(y.test, param_ODPHo)
    mu_ODPHo<-cal_mu_ODP(param_ODPHo)
    
    #Under Gamma Assumption
    glm_Ga_Ho_tr1<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(origin)+log(dev)+dev,data=train.data.numeric,family=GA(mu.link="log", sigma.link ="log"))
    param_GaHo<-fit_param_GA(glm_Ga_Ho_tr1,train.data.numeric,test.data.numeric,tau_Ga)
    dens_GaHo<-cal_dens_GA(y.test, param_GaHo)
    CDF_GaHo<-cal_CDF_GA(y.test, param_GaHo)
    mu_GaHo<-pmax(cal_mu_GA(param_GaHo)-tau_Ga,0)
    
    #Under Log-Normal Assumption
    glm_LN_Ho_tr1<-gamlss(formula=(aggregate_claims+tau_LN)~factor(origin)+log(dev)+dev,data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"))
    param_LNHo<-fit_param_LN(glm_LN_Ho_tr1,train.data.numeric,test.data.numeric,tau_LN)
    dens_LNHo<-cal_dens_LN(y.test, param_LNHo)
    CDF_LNHo<-cal_CDF_LN(y.test, param_LNHo)
    mu_LNHo<-pmax(cal_mu_LN(param_LNHo)-tau_LN,0)
    
    #### Models With Calendar Periods
    #Under ODP Assumption
    glm_ODP_Cal_tr1<-glm(formula=aggregate_claims~factor(dev)+calendar,family=quasipoisson(link="log"),data=train.data.numeric)
    param_ODPCal<-fit_param_ODP(glm_ODP_Cal_tr1,test.data.numeric)
    dens_ODPCal<-cal_dens_ODP(y.test, param_ODPCal)
    CDF_ODPCal<-cal_CDF_ODP(y.test, param_ODPCal)
    mu_ODPCal<-cal_mu_ODP(param_ODPCal)
    
    #Under Gamma Assumption
    glm_Ga_Cal_tr1<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(dev)+calendar,data=train.data.numeric,family=GA(mu.link="log", sigma.link ="log"))
    param_GaCal<-fit_param_GA(glm_Ga_Cal_tr1,train.data.numeric,test.data.numeric,tau_Ga)
    dens_GaCal<-cal_dens_GA(y.test, param_GaCal)
    CDF_GaCal<-cal_CDF_GA(y.test, param_GaCal)
    mu_GaCal<-pmax(cal_mu_GA(param_GaCal)-tau_Ga,0)
    
    #Under Log-Normal Assumption
    glm_LN_Cal_tr1<-gamlss(formula=(aggregate_claims+tau_LN)~factor(dev)+calendar,data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"))
    param_LNCal<-fit_param_LN(glm_LN_Cal_tr1,train.data.numeric,test.data.numeric,tau_LN)
    dens_LNCal<-cal_dens_LN(y.test, param_LNCal)
    CDF_LNCal<-cal_CDF_LN(y.test, param_LNCal)
    mu_LNCal<-pmax(cal_mu_LN(param_LNCal)-tau_LN,0)
    
    ####Smoothing Spline
    #Under Normal Assumption
    sp_Normal<-gamlss(formula=aggregate_claims~scs(origin)+scs(dev),data=train.data.numeric,family=NO(),trace=FALSE)
    sp_Normal <- refit(sp_Normal)
    param_SpNormal<-fit_param_NO(sp_Normal,train.data.numeric,test.data.numeric)
    dens_SpNormal<-cal_dens_Normal(y.test, param_SpNormal)
    CDF_SpNormal<-cal_CDF_Normal(y.test, param_SpNormal)
    mu_SpNormal<-cal_mu_Normal(param_SpNormal)
    
    #Under Gamma Assumption
    sp_Gamma<-gamlss(formula=(aggregate_claims+tau_Ga)~scs(origin)+scs(dev),data=train.data.numeric,family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
    sp_Gamma <- refit(sp_Gamma)
    param_SpGamma<-fit_param_GA(sp_Gamma,train.data.numeric,test.data.numeric,tau_Ga)
    dens_SpGamma<-cal_dens_GA(y.test, param_SpGamma)
    CDF_SpGamma<-cal_CDF_GA(y.test, param_SpGamma)
    mu_SpGamma<-pmax(cal_mu_GA(param_SpGamma)-tau_Ga,0)
    
    #Under Log-Normal Assumption
    sp_LN<-gamlss(formula=(aggregate_claims+tau_LN)~scs(origin)+scs(dev),data=train.data.numeric,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
    sp_LN <- refit(sp_LN)
    param_SpLN<-fit_param_LN(sp_LN,train.data.numeric,test.data.numeric,tau_LN)
    dens_SpLN<-cal_dens_LN(y.test, param_SpLN)
    CDF_SpLN<-cal_CDF_LN(y.test, param_SpLN)
    mu_SpLN<-pmax(cal_mu_LN(param_SpLN)-tau_LN,0)
    
    ####GAMLSS
    #GAMLSS 2: Smooth Effects on the predictor for sigma term
    #Under Gamma Assumption
    gamlss_GA<-gamlss(formula=(aggregate_claims+tau_Ga)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=train.data.numeric,sigma.formula=~as.numeric(as.character(dev)),family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
    gamlss_GA <- refit(gamlss_GA)
    param_GaGAMLSS<-fit_param_GAGamlss(gamlss_GA,train.data.numeric,test.data.numeric,tau_Ga)
    dens_GaGAMLSS<-cal_dens_GA_Gamlss(y.test, param_GaGAMLSS)
    CDF_GaGAMLSS<-cal_CDF_GA_Gamlss(y.test, param_GaGAMLSS)
    mu_GaGAMLSS<-pmax(cal_mu_GA_Gamlss(param_GaGAMLSS)-tau_Ga,0)
    
    #Under Log-Normal Assumption
    gamlss_LN<-gamlss(formula=(aggregate_claims+tau_LN)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=train.data.numeric,sigma.formula=~as.numeric(as.character(dev)),family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
    gamlss_LN <- refit(gamlss_LN)
    param_LNGAMLSS<-fit_param_LNGamlss(gamlss_LN,train.data.numeric,test.data.numeric,tau_LN)
    dens_LNGAMLSS<-cal_dens_LN_Gamlss(y.test, param_LNGAMLSS)
    CDF_LNGAMLSS<-cal_CDF_LN_Gamlss(y.test, param_LNGAMLSS)
    mu_LNGAMLSS<-pmax(cal_mu_LN_Gamlss(param_LNGAMLSS)-tau_LN,0)
    
    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train.data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train.data[train.data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train.data[train.data$origin==i,]$notif_count) + 
            sum(round(predict(fit_nc,newdata=test.data[test.data$origin==i,],type="response"),0)) 
    }
    
    ###PPCI
    ##Fit a ODP/Chain-Ladder model to predict total claim notification count at each accident period 
    PPCI <- calc_PPCI(N, train.data)
    ##Fitting a ODP model on PPCI
    fit_ODP_ppci<-glm(PPCI~factor(dev),family=quasipoisson(link="log"),data=train.data.numeric)
    param_PPCI<-fit_param_PPCI(fit_ODP_ppci,N,test.data)
    dens_PPCI<-cal_dens_PPCI(y.test, param_PPCI)
    CDF_PPCI<-cal_CDF_PPCI(y.test, param_PPCI)
    mu_PPCI<-cal_mu_PPCI(param_PPCI)
    
    ###PPCF
    #Alternatively, fit a ODP on finalization claims that depends on development periods
    odp_FC<-glm(settle_count~factor(dev),data=train.data,family=quasipoisson(link="log"))
    PPCF_data <- calc_PPCF(N, train.data) 
    ODP_PPCF<-glm(PPCF~OT,family=quasipoisson(link="log"),data=PPCF_data)
    param_PPCF <- fit_param_PPCF(odp_FC,ODP_PPCF, N, train.data,test.data) 
    dens_PPCF <- cal_dens_PPCF(y.test, param_PPCF)
    CDF_PPCF <- cal_CDF_PPCF(y.test, param_PPCF)
    mu_PPCF<-cal_mu_PPCF(param_PPCF)
    
    ################################################################################
    ################################################################################
    ################################################################################
    
    component_models <- list(
        fit_ODP_GLM = ODP_GLM,
        fit_GAGLM = Ga_optimTau,
        fit_LNGLM = LN_optimTau,
        fit_ODPHo = glm_ODP_Ho_tr1,
        fit_GaHo = glm_Ga_Ho_tr1,
        fit_LNHo = glm_LN_Ho_tr1,
        fit_ODPCal = glm_ODP_Cal_tr1,
        fit_GaCal = glm_Ga_Cal_tr1,
        fit_LNCal = glm_LN_Cal_tr1,
        fit_SpNormal = sp_Normal,
        fit_SpGamma = sp_Gamma,
        fit_SpLN = sp_LN,
        fit_GaGAMLSS = gamlss_GA,
        fit_LNGAMLSS = gamlss_LN,
        fit_PPCI = list(model = fit_ODP_ppci, N = N),
        fit_PPCF = list(model_subCount = odp_FC, model_subPayments = ODP_PPCF, N=N)
    )
    
    # Predictive Mean  
    meta_mu<-data.frame(
        origin=test.data$origin,
        dev=test.data$dev,
        mu_ODP_GLM=mu_ODPGLM,
        mu_GAGLM=mu_GAGLM,
        mu_LNGLM=mu_LNGLM,
        mu_ODPHo=mu_ODPHo,
        mu_GaHo=mu_GaHo,
        mu_LNHo=mu_LNHo,
        mu_ODPCal=mu_ODPCal,
        mu_GaCal=mu_GaCal,
        mu_LNCal=mu_LNCal,
        mu_SpNormal=mu_SpNormal,
        mu_SpGamma=mu_SpGamma,
        mu_SpLN=mu_SpLN,
        mu_GaGAMLSS=mu_GaGAMLSS,
        mu_LNGAMLSS=mu_LNGAMLSS,
        mu_PPCI=mu_PPCI,
        mu_PPCF=mu_PPCF
    )
    
    # Density
    meta_dens<-data.frame(
        origin=test.data$origin,
        dev=test.data$dev,
        dens_ODP_GLM=dens_ODPGLM,
        dens_GAGLM=dens_GAGLM,
        dens_LNGLM=dens_LNGLM,
        dens_ODPHo=dens_ODPHo,
        dens_GaHo=dens_GaHo,
        dens_LNHo=dens_LNHo,
        dens_ODPCal=dens_ODPCal,
        dens_GaCal=dens_GaCal,
        dens_LNCal=dens_LNCal,
        dens_SpNormal=dens_SpNormal,
        dens_SpGamma=dens_SpGamma,
        dens_SpLN=dens_SpLN,
        dens_GaGAMLSS=dens_GaGAMLSS,
        dens_LNGAMLSS=dens_LNGAMLSS,
        dens_PPCI=dens_PPCI,
        dens_PPCF=dens_PPCF
    )
    
    #CDF
    meta_CDF <-data.frame(
        origin=test.data$origin,
        dev=test.data$dev,
        CDF_ODP_GLM=CDF_ODPGLM,
        CDF_GAGLM=CDF_GAGLM,
        CDF_LNGLM=CDF_LNGLM,
        CDF_ODPHo=CDF_ODPHo,
        CDF_GaHo=CDF_GaHo,
        CDF_LNHo=CDF_LNHo,
        CDF_ODPCal=CDF_ODPCal,
        CDF_GaCal=CDF_GaCal,
        CDF_LNCal=CDF_LNCal,
        CDF_SpNormal=CDF_SpNormal,
        CDF_SpGamma=CDF_SpGamma,
        CDF_SpLN=CDF_SpLN,
        CDF_GaGAMLSS=CDF_GaGAMLSS,
        CDF_LNGAMLSS=CDF_LNGAMLSS,
        CDF_PPCI=CDF_PPCI,
        CDF_PPCF=CDF_PPCF
    )
    
    sink()
    return (list(component_models=component_models, meta_mu=meta_mu, meta_dens=meta_dens, meta_CDF=meta_CDF))
}



