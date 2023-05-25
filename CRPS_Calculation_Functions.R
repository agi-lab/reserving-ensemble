z_u <- 2*round(max(out_sample$value),0)
z_l<-0
z<-z_l:z_u
N<-length(z_u:z_l)
I<-function(y,z) ifelse(y<=z,1,0)

# construct an indicator function within the CRPS

######Write a CRPS function for LN
cal_crps_LN<-function(z,model,newdata){
  y<- newdata$value
  mu<-predict(model,what="mu", newdata=newdata,type="response")
  sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  diff<-mapply(FUN=function(y,mu,sigma)sum((pLOGNO(q=z,mu=mu,sigma=sigma)-I(y,z))^2),y,mu,sigma)
  #for each data point, calculate its CRPS over the range of interest
  return(diff)
}





####Construct a CRPS function for ensemble

cal_CRPS_ensemble<-function(w,z,newdata,ODPGLM,GAGLM,LNGLM,ZAGA,ZALN,ODPHo,GaHo,LNHo,ODPCal,GaCAL,LNCAL,NoSp,GaSp,LNSp,GaGAMLSS,LNGAMLSS,PPCI,NO_PPCI,index_start,index_end,odp_FC,model_subPayments,train_cumF,newdataFC,newdata_Pay,NO_PPCF)
  {
  ##ODP GLM
  pred_ODPGLM<-predict.glm(ODPGLM,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu<- pred_ODPGLM$fit
  pred_ODP_phi<-(pred_ODPGLM$residual.scale)^2
  ##Gamma GLM
  pred_mu_GAGLM<-predict(GAGLM,what="mu",newdata=newdata,type="response")
  pred_sigma_GAGLM<-predict(GAGLM,what="sigma",newdata=newdata,type="response")
  ##LN GLM
  pred_mu_LNGLM<-predict(LNGLM,what="mu",newdata=newdata,type="response")
  pred_sigma_LNGLM<-predict(LNGLM,what="sigma",newdata=newdata,type="response")
  ##ZAGA
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu_ZAGA<-predict(ZAGA,what="mu",newdata=newdata,type="response")
  pred_sigma_ZAGA<-predict(ZAGA,what="sigma",newdata=newdata,type="response")
  pred_nu_ZAGA<-predict(ZAGA,what="nu",newdata=newdata_nu,type="response")
  ##ZALN
  pred_mu_ZALN<-predict(ZALN,parameter="mu",newdata=newdata,type="response")
  pred_sigma_ZALN<-predict(ZALN,parameter="sigma",newdata=newdata,type="response")
  pred_nu_ZALN<-predict(ZALN,parameter="xi0",newdata=newdata_nu,type="response")
  
  ##ODP Hoerl
  newdata_Ho<-newdata
  newdata_Ho$dev<-as.numeric(as.character(newdata$dev))
  pred_ODPHo<-predict.glm(ODPHo,newdata=newdata_Ho,type="response",se.fit=TRUE)
  pred_mu_ODPHo<-pred_ODPHo$fit
  pred_phi_ODPHo<-(pred_ODPHo$residual.scale)^2
  
  #GaHo
  pred_mu_GaHo<-predict(GaHo,what="mu",newdata=newdata_Ho,type="response")
  pred_sigma_GaHo<-predict(GaHo,what="sigma",newdata=newdata_Ho,type="response")
  
  
  #LNHo
  pred_mu_LNHo<-predict(LNHo,what="mu",newdata=newdata_Ho,type="response")
  pred_sigma_LNHo<-predict(LNHo,what="sigma",newdata=newdata_Ho,type="response")
  
  #ODPCal
  pred_ODPCal<-predict.glm(ODPCal,newdata=newdata,type="response",se.fit=TRUE)
  pred_mu_ODPCal<- pred_ODPCal$fit
  pred_phi_ODPCal<-(pred_ODPCal$residual.scale)^2
 
  #GaCAL
  pred_mu_GaCAL<-predict(GaCAL,what="mu",newdata=newdata,type="response")
  pred_sigma_GaCAL<-predict(GaCAL,what="sigma",newdata=newdata,type="response")
  
  #LNCAL
  pred_mu_LNCAL<-predict(LNCAL,what="mu",newdata=newdata,type="response")
  pred_sigma_LNCAL<-predict(LNCAL,what="sigma",newdata=newdata,type="response")
  
  #NoSp
  newdata_numeric<-newdata
  newdata_numeric$origin=as.numeric(as.character(newdata$origin))
  newdata_numeric$dev=as.numeric(as.character(newdata$dev))
  newdata_numeric$Calendar=as.numeric(as.character(newdata$Calendar))
  
  pred_mu_NoSp<-predict(NoSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_NoSp<-predict(NoSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #GaSp
  pred_mu_GaSp<-predict(GaSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_GaSp<-predict(GaSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #LNSP
  pred_mu_LNSp<-predict(LNSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_LNSp<-predict(LNSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #GaGAMLSS
  pred_mu_GaGAMLSS<-predict(GaGAMLSS,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_GaGAMLSS<-predict(GaGAMLSS,what="sigma",newdata=newdata_numeric,type="response")
  
  
  #LNGAMLSS
  pred_mu_LNGAMLSS<-predict(LNGAMLSS,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_LNGAMLSS<-predict(LNGAMLSS,what="sigma",newdata=newdata_numeric,type="response")
  
 
  #PPCI
  N_rep_valid<-rep(NO_PPCI[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_PPCI<-predict.glm(PPCI,newdata=newdata,type="response",se.fit=TRUE)
  pred_mu_PPCI<-pred_PPCI$fit*N_rep_valid
  pred_phi_PPCI<-(pred_PPCI$residual.scale)^2*(N_rep_valid)^2
  
  #PPCF
  N_rep_valid<-rep(NO_PPCF[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
  
  pred_F<-predict(odp_FC,newdata=newdataFC,type="response")
  pred_F_cum<-list()
  for(i in index_start:index_end){
    pf<-predict(odp_FC,newdata=newdataFC[newdataFC$origin==i,],type="response")
    pf_cum<-round(cumsum(c(train_cumF[train_cumF$origin==i,]$value[nrow(train_cumF[train_cumF$origin==i,])],pf))[-1],0)
    pred_F_cum[[i]]<-pf_cum
  }
  pred_F_cum<-as.vector(unlist(pred_F_cum[index_start:index_end]))
  
  pred_OT<-pred_F_cum/N_rep_valid
  pred_payment<-predict.glm(model_subPayments,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
  pred_mu_PPCF<-pred_payment$fit*pred_F
  pred_phi_PPCF<-mean((pred_payment$residual.scale)^2*(pred_F)^2)
  
  y<-newdata$value
  

  
  
  diff<-mapply(FUN=function(y,pred_ODP_mu,pred_ODP_phi,pred_mu_GAGLM,pred_sigma_GAGLM,pred_mu_LNGLM,pred_sigma_LNGLM,pred_mu_ZAGA,pred_sigma_ZAGA,pred_nu_ZAGA,pred_mu_ZALN,pred_sigma_ZALN,pred_nu_ZALN,pred_mu_ODPHo,pred_phi_ODPHo,pred_mu_GaHo,pred_sigma_GaHo,pred_mu_LNHo,pred_sigma_LNHo,pred_mu_ODPCal,pred_phi_ODPCal,pred_mu_GaCAL,pred_sigma_GaCAL,pred_mu_LNCAL,pred_sigma_LNCAL,pred_mu_NoSp,pred_sigma_NoSp,pred_mu_GaSp,pred_sigma_GaSp,pred_mu_LNSp,pred_sigma_LNSp,pred_mu_GaGAMLSS,pred_sigma_GaGAMLSS,pred_mu_LNGAMLSS,pred_sigma_LNGAMLSS,pred_mu_PPCI,pred_phi_PPCI,pred_mu_PPCF,pred_phi_PPC)sum((cbind(ptweedie(q=z,mu=pred_ODP_mu,phi=pred_ODP_phi,power=1),pGA(q=z,pred_mu_GAGLM,pred_sigma_GAGLM),pLNO(q=z,pred_mu_LNGLM,pred_sigma_LNGLM),pZAGA(q=z,pred_mu_ZAGA,pred_sigma_ZAGA,pred_nu_ZAGA),pZALN(q=z,pred_mu_ZALN,pred_sigma_ZALN,pred_nu_ZALN),ptweedie(q=z,mu=pred_mu_ODPHo,phi=pred_phi_ODPHo,power=1),pGA(q=z,pred_mu_GaHo,pred_sigma_GaHo),pLNO(q=z,pred_mu_LNHo,pred_sigma_LNHo),ptweedie(q=z,mu=pred_mu_ODPCal,phi=pred_phi_ODPCal,power=1),pGA(q=z,pred_mu_GaCAL,pred_sigma_GaCAL),pLNO(q=z,pred_mu_LNCAL,pred_sigma_LNCAL),pNO(q=z,pred_mu_NoSp,pred_sigma_NoSp),pGA(q=z,pred_mu_GaSp,pred_sigma_GaSp),pLNO(q=z,pred_mu_LNSp,pred_sigma_LNSp),pGA(q=z,pred_mu_GaGAMLSS,pred_sigma_GaGAMLSS),pLNO(q=z,pred_mu_LNGAMLSS,pred_sigma_LNGAMLSS),ptweedie(q=z,mu=pred_mu_PPCI,phi=pred_phi_PPCI,power=1),ptweedie(q=z,mu=pred_mu_PPCF,phi=pred_phi_PPCF,power=1))%*%w-I(y,z))^2),y,pred_ODP_mu,pred_ODP_phi,pred_mu_GAGLM,pred_sigma_GAGLM,pred_mu_LNGLM,pred_sigma_LNGLM,pred_mu_ZAGA,pred_sigma_ZAGA,pred_nu_ZAGA,pred_mu_ZALN,pred_sigma_ZALN,pred_nu_ZALN,pred_mu_ODPHo,pred_phi_ODPHo,pred_mu_GaHo,pred_sigma_GaHo,pred_mu_LNHo,pred_sigma_LNHo,pred_mu_ODPCal,pred_phi_ODPCal,pred_mu_GaCAL,pred_sigma_GaCAL,pred_mu_LNCAL,pred_sigma_LNCAL,pred_mu_NoSp,pred_sigma_NoSp,pred_mu_GaSp,pred_sigma_GaSp,pred_mu_LNSp,pred_sigma_LNSp,pred_mu_GaGAMLSS,pred_sigma_GaGAMLSS,pred_mu_LNGAMLSS,pred_sigma_LNGAMLSS,pred_mu_PPCI,pred_phi_PPCI,pred_mu_PPCF,pred_phi_PPCF)

  return (diff)
}





#####CRPS function for ADLP ensembles
cal_CRPS_ADLPensemble<-function(w1,w2,z,index_subset1,index_subset2,newdata,tau_Ga,tau_LN,ODPGLM,GAGLM,LNGLM,ZAGA,ZALN,ODPHo,GaHo,LNHo,ODPCal,GaCAL,LNCAL,NoSp,GaSp,LNSp,GaGAMLSS,LNGAMLSS,PPCI,NO_PPCI,index_start,index_end,odp_FC,model_subPayments,train_cumF,newdataFC,newdata_Pay,NO_PPCF)
{
  ##ODP GLM
  pred_ODPGLM<-predict.glm(ODPGLM,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu<- pred_ODPGLM$fit
  pred_ODP_phi<-(pred_ODPGLM$residual.scale)^2
  ##Gamma GLM
  pred_mu_GAGLM<-predict(GAGLM,what="mu",newdata=newdata,type="response")
  pred_sigma_GAGLM<-predict(GAGLM,what="sigma",newdata=newdata,type="response")
  ##LN GLM
  pred_mu_LNGLM<-predict(LNGLM,what="mu",newdata=newdata,type="response")
  pred_sigma_LNGLM<-predict(LNGLM,what="sigma",newdata=newdata,type="response")
  ##ZAGA
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu_ZAGA<-predict(ZAGA,what="mu",newdata=newdata,type="response")
  pred_sigma_ZAGA<-predict(ZAGA,what="sigma",newdata=newdata,type="response")
  pred_nu_ZAGA<-predict(ZAGA,what="nu",newdata=newdata_nu,type="response")
  ##ZALN
  pred_mu_ZALN<-predict(ZALN,parameter="mu",newdata=newdata,type="response")
  pred_sigma_ZALN<-predict(ZALN,parameter="sigma",newdata=newdata,type="response")
  pred_nu_ZALN<-predict(ZALN,parameter="xi0",newdata=newdata_nu,type="response")
  
  ##ODP Hoerl
  newdata_Ho<-newdata
  newdata_Ho$dev<-as.numeric(as.character(newdata$dev))
  pred_ODPHo<-predict.glm(ODPHo,newdata=newdata_Ho,type="response",se.fit=TRUE)
  pred_mu_ODPHo<-pred_ODPHo$fit
  pred_phi_ODPHo<-(pred_ODPHo$residual.scale)^2
  
  #GaHo
  pred_mu_GaHo<-predict(GaHo,what="mu",newdata=newdata_Ho,type="response")
  pred_sigma_GaHo<-predict(GaHo,what="sigma",newdata=newdata_Ho,type="response")
  
  
  #LNHo
  pred_mu_LNHo<-predict(LNHo,what="mu",newdata=newdata_Ho,type="response")
  pred_sigma_LNHo<-predict(LNHo,what="sigma",newdata=newdata_Ho,type="response")
  
  #ODPCal
  pred_ODPCal<-predict.glm(ODPCal,newdata=newdata,type="response",se.fit=TRUE)
  pred_mu_ODPCal<- pred_ODPCal$fit
  pred_phi_ODPCal<-(pred_ODPCal$residual.scale)^2
  
  #GaCAL
  pred_mu_GaCAL<-predict(GaCAL,what="mu",newdata=newdata,type="response")
  pred_sigma_GaCAL<-predict(GaCAL,what="sigma",newdata=newdata,type="response")
  
  #LNCAL
  pred_mu_LNCAL<-predict(LNCAL,what="mu",newdata=newdata,type="response")
  pred_sigma_LNCAL<-predict(LNCAL,what="sigma",newdata=newdata,type="response")
  
  #NoSp
  newdata_numeric<-newdata
  newdata_numeric$origin=as.numeric(as.character(newdata$origin))
  newdata_numeric$dev=as.numeric(as.character(newdata$dev))
  newdata_numeric$Calendar=as.numeric(as.character(newdata$Calendar))
  
  pred_mu_NoSp<-predict(NoSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_NoSp<-predict(NoSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #GaSp
  pred_mu_GaSp<-predict(GaSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_GaSp<-predict(GaSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #LNSP
  pred_mu_LNSp<-predict(LNSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_LNSp<-predict(LNSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #GaGAMLSS
 
  pred_mu_GaGAMLSS<-predict(GaGAMLSS,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_GaGAMLSS<-predict(GaGAMLSS,what="sigma",newdata=newdata_numeric,type="response")
   
  
  
  #LNGAMLSS
  #pred_mu_LNGAMLSS<-predict(LNGAMLSS,what="mu",newdata=newdata_mu,type="response")
  #pred_sigma_LNGAMLSS<-predict(LNGAMLSS,what="sigma",newdata=newdata_sigma,type="response")
  
  pred_mu_LNGAMLSS<-predict(LNGAMLSS,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_LNGAMLSS<-predict(LNGAMLSS,what="sigma",newdata=newdata_numeric,type="response")
  
  
  #PPCI
  N_rep_valid<-rep(NO_PPCI[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_PPCI<-predict.glm(PPCI,newdata=newdata,type="response",se.fit=TRUE)
  pred_mu_PPCI<-pred_PPCI$fit*N_rep_valid
  pred_phi_PPCI<-(pred_PPCI$residual.scale)^2*(N_rep_valid)^2
  
  #PPCF
  N_rep_valid<-rep(NO_PPCF[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
  
  pred_F<-predict(odp_FC,newdata=newdataFC,type="response")
  pred_F_cum<-list()
  for(i in index_start:index_end){
    pf<-predict(odp_FC,newdata=newdataFC[newdataFC$origin==i,],type="response")
    pf_cum<-round(cumsum(c(train_cumF[train_cumF$origin==i,]$value[nrow(train_cumF[train_cumF$origin==i,])],pf))[-1],0)
    pred_F_cum[[i]]<-pf_cum
  }
  pred_F_cum<-as.vector(unlist(pred_F_cum[index_start:index_end]))
  
  pred_OT<-pred_F_cum/N_rep_valid
  pred_payment<-predict.glm(model_subPayments,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
  pred_mu_PPCF<-pred_payment$fit*pred_F
  pred_phi_PPCF<-mean((pred_payment$residual.scale)^2*(pred_F)^2)
  
  y<-newdata$value
  ###do not apply index number to phi as it is a single number
  
  #####Inputs in the first subset
  # y1=y[index_subset1]
  # pred1_ODP_mu<-pred_ODP_mu[index_subset1]
  # pred1_ODP_phi<-pred_ODP_phi
  # 
  # pred1_mu_GAGLM<-ifelse(pred_mu_GAGLM[index_subset1]-tau_Ga>0,pred_mu_GAGLM[index_subset1]-tau_Ga,0.01)
  # pred1_sigma_GAGLM<-pred_sigma_GAGLM[index_subset1]
  # 
  # pred1_mu_LNGLM<-ifelse(pred_mu_LNGLM[index_subset1]-tau_LN>0,pred_mu_LNGLM[index_subset1]-tau_LN,0.01)
  # pred1_sigma_LNGLM<-pred_sigma_LNGLM[index_subset1]
  # 
  # pred1_mu_ZAGA<-pred_mu_ZAGA[index_subset1]
  # pred1_sigma_ZAGA<-pred_sigma_ZAGA[index_subset1]
  # pred1_nu_ZAGA<-pred_nu_ZAGA[index_subset1]
  # pred1_mu_ZALN<-pred_mu_ZALN[index_subset1]
  # pred1_sigma_ZALN<-pred_sigma_ZALN[index_subset1]
  # pred1_nu_ZALN<-pred_nu_ZALN[index_subset1]
  # pred1_mu_ODPHo<-pred_mu_ODPHo[index_subset1]
  # pred1_phi_ODPHo<-pred_phi_ODPHo
  # 
  # pred1_mu_GaHo<-ifelse(pred_mu_GaHo[index_subset1]-tau_Ga>0,pred_mu_GaHo[index_subset1]-tau_Ga,0.01)
  # pred1_sigma_GaHo<-pred_sigma_GaHo[index_subset1]
  # 
  # pred1_mu_LNHo<-ifelse(pred_mu_LNHo[index_subset1]-tau_LN>0,pred_mu_LNHo[index_subset1]-tau_LN,0.01)
  # pred1_sigma_LNHo<-pred_sigma_LNHo[index_subset1]
  # 
  # pred1_mu_ODPCal<-pred_mu_ODPCal[index_subset1]
  # pred1_phi_ODPCal<-pred_phi_ODPCal
  # 
  # pred1_mu_GaCAL<-ifelse(pred_mu_GaCAL[index_subset1]-tau_Ga>0,pred_mu_GaCAL[index_subset1]-tau_Ga,0.01)
  # pred1_sigma_GaCAL<-pred_sigma_GaCAL[index_subset1]
  # 
  # pred1_mu_LNCAL<-ifelse(pred_mu_LNCAL[index_subset1]-tau_LN>0,pred_mu_LNCAL[index_subset1]-tau_LN,0.01)
  # pred1_sigma_LNCAL<-pred_sigma_LNCAL[index_subset1]
  # 
  # pred1_mu_NoSp<-pred_mu_NoSp[index_subset1]
  # pred1_sigma_NoSp<-pred_sigma_NoSp[index_subset1]
  # 
  # pred1_mu_GaSp<-ifelse(pred_mu_GaSp[index_subset1]-tau_Ga>0,pred_mu_GaSp[index_subset1]-tau_Ga,0.01)
  # pred1_sigma_GaSp<-pred_sigma_GaSp[index_subset1]
  # 
  # pred1_mu_LNSp<-ifelse(pred_mu_LNSp[index_subset1]-tau_LN>0,pred_mu_LNSp[index_subset1]-tau_LN,0.01)
  # pred1_sigma_LNSp<-pred_sigma_LNSp[index_subset1]
  # 
  # pred1_mu_GaGAMLSS<-ifelse(pred_mu_GaGAMLSS[index_subset1]-tau_Ga>0,pred_mu_GaGAMLSS[index_subset1]-tau_Ga,0.01)
  # pred1_sigma_GaGAMLSS<-pred_sigma_GaGAMLSS[index_subset1]
  # 
  # pred1_mu_LNGAMLSS<-ifelse(pred_mu_LNGAMLSS[index_subset1]-tau_LN>0,pred_mu_LNGAMLSS[index_subset1]-tau_LN,0.01)
  # pred1_sigma_LNGAMLSS<-pred_sigma_LNGAMLSS[index_subset1]
  # 
  # pred1_mu_PPCI<-pred_mu_PPCI[index_subset1]
  # pred1_phi_PPCI<-pred_phi_PPCI[index_subset1]
  # pred1_mu_PPCF<-pred_mu_PPCF[index_subset1]
  # pred1_phi_PPCF<-pred_phi_PPCF
  y1=y[index_subset1]
  pred1_ODP_mu<-pred_ODP_mu[index_subset1]
  pred1_ODP_phi<-pred_ODP_phi
  
  pred1_mu_GAGLM<-pred_mu_GAGLM[index_subset1]
  pred1_sigma_GAGLM<-pred_sigma_GAGLM[index_subset1]
  
  pred1_mu_LNGLM<-pred_mu_LNGLM[index_subset1]
  pred1_sigma_LNGLM<-pred_sigma_LNGLM[index_subset1]
  
  pred1_mu_ZAGA<-pred_mu_ZAGA[index_subset1]
  pred1_sigma_ZAGA<-pred_sigma_ZAGA[index_subset1]
  pred1_nu_ZAGA<-pred_nu_ZAGA[index_subset1]
  pred1_mu_ZALN<-pred_mu_ZALN[index_subset1]
  pred1_sigma_ZALN<-pred_sigma_ZALN[index_subset1]
  pred1_nu_ZALN<-pred_nu_ZALN[index_subset1]
  pred1_mu_ODPHo<-pred_mu_ODPHo[index_subset1]
  pred1_phi_ODPHo<-pred_phi_ODPHo
  
  pred1_mu_GaHo<-pred_mu_GaHo[index_subset1]
  pred1_sigma_GaHo<-pred_sigma_GaHo[index_subset1]
  
  pred1_mu_LNHo<-pred_mu_LNHo[index_subset1]
  pred1_sigma_LNHo<-pred_sigma_LNHo[index_subset1]
  
  pred1_mu_ODPCal<-pred_mu_ODPCal[index_subset1]
  pred1_phi_ODPCal<-pred_phi_ODPCal
  
  pred1_mu_GaCAL<-pred_mu_GaCAL[index_subset1]
  pred1_sigma_GaCAL<-pred_sigma_GaCAL[index_subset1]
  
  pred1_mu_LNCAL<-pred_mu_LNCAL[index_subset1]
  pred1_sigma_LNCAL<-pred_sigma_LNCAL[index_subset1]
  
  pred1_mu_NoSp<-pred_mu_NoSp[index_subset1]
  pred1_sigma_NoSp<-pred_sigma_NoSp[index_subset1]
  
  pred1_mu_GaSp<-pred_mu_GaSp[index_subset1]
  pred1_sigma_GaSp<-pred_sigma_GaSp[index_subset1]
  
  pred1_mu_LNSp<-pred_mu_LNSp[index_subset1]
  pred1_sigma_LNSp<-pred_sigma_LNSp[index_subset1]
  
  pred1_mu_GaGAMLSS<-pred_mu_GaGAMLSS[index_subset1]
  pred1_sigma_GaGAMLSS<-pred_sigma_GaGAMLSS[index_subset1]
  
  pred1_mu_LNGAMLSS<-pred_mu_LNGAMLSS[index_subset1]
  pred1_sigma_LNGAMLSS<-pred_sigma_LNGAMLSS[index_subset1]
  
  pred1_mu_PPCI<-pred_mu_PPCI[index_subset1]
  pred1_phi_PPCI<-pred_phi_PPCI[index_subset1]
  pred1_mu_PPCF<-pred_mu_PPCF[index_subset1]
  pred1_phi_PPCF<-pred_phi_PPCF
  
  
  
  ####Inputs in the second subset
  y2=y[index_subset2]
  pred2_ODP_mu<-pred_ODP_mu[index_subset2]
  pred2_ODP_phi<-pred_ODP_phi
  
  pred2_mu_GAGLM<-pred_mu_GAGLM[index_subset2]
  pred2_sigma_GAGLM<-pred_sigma_GAGLM[index_subset2]
  
  pred2_mu_LNGLM<-pred_mu_LNGLM[index_subset2]
  pred2_sigma_LNGLM<-pred_sigma_LNGLM[index_subset2]
  
  pred2_mu_ZAGA<-pred_mu_ZAGA[index_subset2]
  pred2_sigma_ZAGA<-pred_sigma_ZAGA[index_subset2]
  pred2_nu_ZAGA<-pred_nu_ZAGA[index_subset2]
  pred2_mu_ZALN<-pred_mu_ZALN[index_subset2]
  pred2_sigma_ZALN<-pred_sigma_ZALN[index_subset2]
  pred2_nu_ZALN<-pred_nu_ZALN[index_subset2]
  pred2_mu_ODPHo<-pred_mu_ODPHo[index_subset2]
  pred2_phi_ODPHo<-pred_phi_ODPHo
  
  pred2_mu_GaHo<-pred_mu_GaHo[index_subset2]
  pred2_sigma_GaHo<-pred_sigma_GaHo[index_subset2]
  
  pred2_mu_LNHo<-pred_mu_LNHo[index_subset2]
  pred2_sigma_LNHo<-pred_sigma_LNHo[index_subset2]
  
  pred2_mu_ODPCal<-pred_mu_ODPCal[index_subset2]
  pred2_phi_ODPCal<-pred_phi_ODPCal
  
  pred2_mu_GaCAL<-pred_mu_GaCAL[index_subset2]
  pred2_sigma_GaCAL<-pred_sigma_GaCAL[index_subset2]
  
  pred2_mu_LNCAL<-pred_mu_LNCAL[index_subset2]
  pred2_sigma_LNCAL<-pred_sigma_LNCAL[index_subset2]
  
  pred2_mu_NoSp<-pred_mu_NoSp[index_subset2]
  pred2_sigma_NoSp<-pred_sigma_NoSp[index_subset2]
  
  pred2_mu_GaSp<-pred_mu_GaSp[index_subset2]
  pred2_sigma_GaSp<-pred_sigma_GaSp[index_subset2]
  
  pred2_mu_LNSp<-pred_mu_LNSp[index_subset2]
  pred2_sigma_LNSp<-pred_sigma_LNSp[index_subset2]
  
  pred2_mu_GaGAMLSS<-pred_mu_GaGAMLSS[index_subset2]
  pred2_sigma_GaGAMLSS<-pred_sigma_GaGAMLSS[index_subset2]
  
  pred2_mu_LNGAMLSS<-pred_mu_LNGAMLSS[index_subset2]
  pred2_sigma_LNGAMLSS<-pred_sigma_LNGAMLSS[index_subset2]
  
  pred2_mu_PPCI<-pred_mu_PPCI[index_subset2]
  pred2_phi_PPCI<-pred_phi_PPCI[index_subset2]
  pred2_mu_PPCF<-pred_mu_PPCF[index_subset2]
  pred2_phi_PPCF<-pred_phi_PPCF
  
  
  diff1<-mapply(FUN=function(y1,pred1_ODP_mu,pred1_ODP_phi,pred1_mu_GAGLM,pred1_sigma_GAGLM,pred1_mu_LNGLM,pred1_sigma_LNGLM,pred1_mu_ZAGA,pred1_sigma_ZAGA,pred1_nu_ZAGA,pred1_mu_ZALN,pred1_sigma_ZALN,pred1_nu_ZALN,pred1_mu_ODPHo,pred1_phi_ODPHo,pred1_mu_GaHo,pred1_sigma_GaHo,pred1_mu_LNHo,pred1_sigma_LNHo,pred1_mu_ODPCal,pred1_phi_ODPCal,pred1_mu_GaCAL,pred1_sigma_GaCAL,pred1_mu_LNCAL,pred1_sigma_LNCAL,pred1_mu_NoSp,pred1_sigma_NoSp,pred1_mu_GaSp,pred1_sigma_GaSp,pred1_mu_LNSp,pred1_sigma_LNSp,pred1_mu_GaGAMLSS,pred1_sigma_GaGAMLSS,pred1_mu_LNGAMLSS,pred1_sigma_LNGAMLSS,pred1_mu_PPCI,pred1_phi_PPCI,pred1_mu_PPCF,pred1_phi_PPCF)sum((cbind(ptweedie(q=z,mu=pred1_ODP_mu,phi=pred1_ODP_phi,power=1),pGA(q=z,pred1_mu_GAGLM,pred1_sigma_GAGLM),pLNO(q=z,pred1_mu_LNGLM,pred1_sigma_LNGLM),pZAGA(q=z,pred1_mu_ZAGA,pred1_sigma_ZAGA,pred1_nu_ZAGA),pZALN(q=z,pred1_mu_ZALN,pred1_sigma_ZALN,pred1_nu_ZALN),ptweedie(q=z,mu=pred1_mu_ODPHo,phi=pred1_phi_ODPHo,power=1),pGA(q=z,pred1_mu_GaHo,pred1_sigma_GaHo),pLNO(q=z,pred1_mu_LNHo,pred1_sigma_LNHo),ptweedie(q=z,mu=pred1_mu_ODPCal,phi=pred1_phi_ODPCal,power=1),pGA(q=z,pred1_mu_GaCAL,pred1_sigma_GaCAL),pLNO(q=z,pred1_mu_LNCAL,pred1_sigma_LNCAL),pNO(q=z,pred1_mu_NoSp,pred1_sigma_NoSp),pGA(q=z,pred1_mu_GaSp,pred1_sigma_GaSp),pLNO(q=z,pred1_mu_LNSp,pred1_sigma_LNSp),pGA(q=z,pred1_mu_GaGAMLSS,pred1_sigma_GaGAMLSS),pLNO(q=z,pred1_mu_LNGAMLSS,pred1_sigma_LNGAMLSS),ptweedie(q=z,mu=pred1_mu_PPCI,phi=pred1_phi_PPCI,power=1),ptweedie(q=z,mu=pred1_mu_PPCF,phi=pred1_phi_PPCF,power=1))%*%w1-I(y1,z))^2),y1=y1,pred1_ODP_mu,pred1_ODP_phi,pred1_mu_GAGLM,pred1_sigma_GAGLM,pred1_mu_LNGLM,pred1_sigma_LNGLM,pred1_mu_ZAGA,pred1_sigma_ZAGA,pred1_nu_ZAGA,pred1_mu_ZALN,pred1_sigma_ZALN,pred1_nu_ZALN,pred1_mu_ODPHo,pred1_phi_ODPHo,pred1_mu_GaHo,pred1_sigma_GaHo,pred1_mu_LNHo,pred1_sigma_LNHo,pred1_mu_ODPCal,pred1_phi_ODPCal,pred1_mu_GaCAL,pred1_sigma_GaCAL,pred1_mu_LNCAL,pred1_sigma_LNCAL,pred1_mu_NoSp,pred1_sigma_NoSp,pred1_mu_GaSp,pred1_sigma_GaSp,pred1_mu_LNSp,pred1_sigma_LNSp,pred1_mu_GaGAMLSS,pred1_sigma_GaGAMLSS,pred1_mu_LNGAMLSS,pred1_sigma_LNGAMLSS,pred1_mu_PPCI,pred1_phi_PPCI,pred1_mu_PPCF,pred1_phi_PPCF)
  
  diff2<-mapply(FUN=function(y2,pred2_ODP_mu,pred2_ODP_phi,pred2_mu_GAGLM,pred2_sigma_GAGLM,pred2_mu_LNGLM,pred2_sigma_LNGLM,pred2_mu_ZAGA,pred2_sigma_ZAGA,pred2_nu_ZAGA,pred2_mu_ZALN,pred2_sigma_ZALN,pred2_nu_ZALN,pred2_mu_ODPHo,pred2_phi_ODPHo,pred2_mu_GaHo,pred2_sigma_GaHo,pred2_mu_LNHo,pred2_sigma_LNHo,pred2_mu_ODPCal,pred2_phi_ODPCal,pred2_mu_GaCAL,pred2_sigma_GaCAL,pred2_mu_LNCAL,pred2_sigma_LNCAL,pred2_mu_NoSp,pred2_sigma_NoSp,pred2_mu_GaSp,pred2_sigma_GaSp,pred2_mu_LNSp,pred2_sigma_LNSp,pred2_mu_GaGAMLSS,pred2_sigma_GaGAMLSS,pred2_mu_LNGAMLSS,pred2_sigma_LNGAMLSS,pred2_mu_PPCI,pred2_phi_PPCI,pred2_mu_PPCF,pred2_phi_PPCF)sum((cbind(ptweedie(q=z,mu=pred2_ODP_mu,phi=pred2_ODP_phi,power=1),pGA(q=z,pred2_mu_GAGLM,pred2_sigma_GAGLM),pLNO(q=z,pred2_mu_LNGLM,pred2_sigma_LNGLM),pZAGA(q=z,pred2_mu_ZAGA,pred2_sigma_ZAGA,pred2_nu_ZAGA),pZALN(q=z,pred2_mu_ZALN,pred2_sigma_ZALN,pred2_nu_ZALN),ptweedie(q=z,mu=pred2_mu_ODPHo,phi=pred2_phi_ODPHo,power=1),pGA(q=z,pred2_mu_GaHo,pred2_sigma_GaHo),pLNO(q=z,pred2_mu_LNHo,pred2_sigma_LNHo),ptweedie(q=z,mu=pred2_mu_ODPCal,phi=pred2_phi_ODPCal,power=1),pGA(q=z,pred2_mu_GaCAL,pred2_sigma_GaCAL),pLNO(q=z,pred2_mu_LNCAL,pred2_sigma_LNCAL),pNO(q=z,pred2_mu_NoSp,pred2_sigma_NoSp),pGA(q=z,pred2_mu_GaSp,pred2_sigma_GaSp),pLNO(q=z,pred2_mu_LNSp,pred2_sigma_LNSp),pGA(q=z,pred2_mu_GaGAMLSS,pred2_sigma_GaGAMLSS),pLNO(q=z,pred2_mu_LNGAMLSS,pred2_sigma_LNGAMLSS),ptweedie(q=z,mu=pred2_mu_PPCI,phi=pred2_phi_PPCI,power=1),ptweedie(q=z,mu=pred2_mu_PPCF,phi=pred2_phi_PPCF,power=1))%*%w2-I(y2,z))^2),y2=y2,pred2_ODP_mu,pred2_ODP_phi,pred2_mu_GAGLM,pred2_sigma_GAGLM,pred2_mu_LNGLM,pred2_sigma_LNGLM,pred2_mu_ZAGA,pred2_sigma_ZAGA,pred2_nu_ZAGA,pred2_mu_ZALN,pred2_sigma_ZALN,pred2_nu_ZALN,pred2_mu_ODPHo,pred2_phi_ODPHo,pred2_mu_GaHo,pred2_sigma_GaHo,pred2_mu_LNHo,pred2_sigma_LNHo,pred2_mu_ODPCal,pred2_phi_ODPCal,pred2_mu_GaCAL,pred2_sigma_GaCAL,pred2_mu_LNCAL,pred2_sigma_LNCAL,pred2_mu_NoSp,pred2_sigma_NoSp,pred2_mu_GaSp,pred2_sigma_GaSp,pred2_mu_LNSp,pred2_sigma_LNSp,pred2_mu_GaGAMLSS,pred2_sigma_GaGAMLSS,pred2_mu_LNGAMLSS,pred2_sigma_LNGAMLSS,pred2_mu_PPCI,pred2_phi_PPCI,pred2_mu_PPCF,pred2_phi_PPCF)
  
  diff<-c(diff1,diff2)
  return (diff)
}



####################Some examples of CRPS calculation

# crps_LN<-cal_crps_LN(z,sp_LN_In,out_sample_numeric)
# 
# #Mean CRPS for Equally Weighted ensemble
# mean(cal_CRPS_ensemble(w=rep(1/18,18),z=z,newdata=out_sample,ODPGLM=ODP_GLM_in,GAGLM=Ga_optimTau,LNGLM=LN_optimTau,ZAGA=ZAGA_in,ZALN=LN_in,ODPHo=glm_ODP_Ho_In,GaHo=glm_Ga_Ho_In,LNHo=glm_LN_Ho_In,ODPCal=glm_ODP_Cal_In,GaCAL=glm_Ga_Cal_In,LNCAL=glm_LN_Cal_In,NoSp=sp_Normal_In,GaSp=sp_Gamma_In,LNSp=sp_LN_In,GaGAMLSS=gamlss2_GA_In,LNGAMLSS=gamlss2_LN_In,PPCI=fit_ODP_ppci_In,NO_PPCI=N_in,index_start=2,index_end=40,odp_FC=odp_FC_In,model_subPayments=ODP_PPCF_In,train_cumF=cum_F_dt,newdataFC=FUN_dat_Out,newdata_Pay=out_sample,NO_PPCF=N_in))
# 
# #Mean CRPS for SLP
# mean(cal_CRPS_ensemble(w=model_weights_simul_par0[[100]],z=z,newdata=out_sample,ODPGLM=ODP_GLM_in,GAGLM=Ga_optimTau,LNGLM=LN_optimTau,ZAGA=ZAGA_in,ZALN=LN_in,ODPHo=glm_ODP_Ho_In,GaHo=glm_Ga_Ho_In,LNHo=glm_LN_Ho_In,ODPCal=glm_ODP_Cal_In,GaCAL=glm_Ga_Cal_In,LNCAL=glm_LN_Cal_In,NoSp=sp_Normal_In,GaSp=sp_Gamma_In,LNSp=sp_LN_In,GaGAMLSS=gamlss2_GA_In,LNGAMLSS=gamlss2_LN_In,PPCI=fit_ODP_ppci_In,NO_PPCI=N_in,index_start=2,index_end=40,odp_FC=odp_FC_In,model_subPayments=ODP_PPCF_In,train_cumF=cum_F_dt,newdataFC=FUN_dat_Out,newdata_Pay=out_sample,NO_PPCF=N_in))
# 
# #Mean CRPS for ADLP8
# ind_subset1_par7_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=15
# ind_subset2_par7_new<-as.numeric(as.character(out_sample$origin))>=16&as.numeric(as.character(out_sample$origin))<=40
# crps_ADLP8<-cal_CRPS_ADLPensemble(w1=as.vector(unlist(model_weights_simul_par7_new[[100]][1])),w2=as.vector(unlist(model_weights_simul_par7_new[[100]][2])),tau_Ga=tau_Ga,tau_LN=tau_LN,index_subset1=ind_subset1_par7_new,index_subset2=ind_subset2_par7_new,z=z,newdata=out_sample,ODPGLM=ODP_GLM_in,GAGLM=Ga_optimTau,LNGLM=LN_optimTau,ZAGA=ZAGA_in,ZALN=LN_in,ODPHo=glm_ODP_Ho_In,GaHo=glm_Ga_Ho_In,LNHo=glm_LN_Ho_In,ODPCal=glm_ODP_Cal_In,GaCAL=glm_Ga_Cal_In,LNCAL=glm_LN_Cal_In,NoSp=sp_Normal_In,GaSp=sp_Gamma_In,LNSp=sp_LN_In,GaGAMLSS=gamlss2_GA_In,LNGAMLSS=gamlss2_LN_In,PPCI=fit_ODP_ppci_In,NO_PPCI=N_in,index_start=2,index_end=40,odp_FC=odp_FC_In,model_subPayments=ODP_PPCF_In,train_cumF=cum_F_dt,newdataFC=FUN_dat_Out,newdata_Pay=out_sample,NO_PPCF=N_in)
# mean(crps_ADLP8)

#plot(x=out_sample$value,y=SpLN_crps[,2])
#points(x=out_sample$value,y=ADLP8_crps[,2],col="blue")

##############Additional Function
#Return the CDF at a specified point

####CDF funtion of ensemble

cal_CDF_ensemble<-function(w,z,newdata,index,ODPGLM,GAGLM,LNGLM,ZAGA,ZALN,ODPHo,GaHo,LNHo,ODPCal,GaCAL,LNCAL,NoSp,GaSp,LNSp,GaGAMLSS,LNGAMLSS,PPCI,NO_PPCI,index_start,index_end,odp_FC,model_subPayments,train_cumF,newdataFC,newdata_Pay,NO_PPCF)
{
  
  ##ODP GLM
  pred_ODPGLM<-predict.glm(ODPGLM,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu<- pred_ODPGLM$fit
  pred_ODP_phi<-(pred_ODPGLM$residual.scale)^2
  ##Gamma GLM
  pred_mu_GAGLM<-predict(GAGLM,what="mu",newdata=newdata,type="response")
  pred_sigma_GAGLM<-predict(GAGLM,what="sigma",newdata=newdata,type="response")
  ##LN GLM
  pred_mu_LNGLM<-predict(LNGLM,what="mu",newdata=newdata,type="response")
  pred_sigma_LNGLM<-predict(LNGLM,what="sigma",newdata=newdata,type="response")
  ##ZAGA
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu_ZAGA<-predict(ZAGA,what="mu",newdata=newdata,type="response")
  pred_sigma_ZAGA<-predict(ZAGA,what="sigma",newdata=newdata,type="response")
  pred_nu_ZAGA<-predict(ZAGA,what="nu",newdata=newdata_nu,type="response")
  ##ZALN
  pred_mu_ZALN<-predict(ZALN,parameter="mu",newdata=newdata,type="response")
  pred_sigma_ZALN<-predict(ZALN,parameter="sigma",newdata=newdata,type="response")
  pred_nu_ZALN<-predict(ZALN,parameter="xi0",newdata=newdata_nu,type="response")
  
  ##ODP Hoerl
  newdata_Ho<-newdata
  newdata_Ho$dev<-as.numeric(as.character(newdata$dev))
  pred_ODPHo<-predict.glm(ODPHo,newdata=newdata_Ho,type="response",se.fit=TRUE)
  pred_mu_ODPHo<-pred_ODPHo$fit
  pred_phi_ODPHo<-(pred_ODPHo$residual.scale)^2
  
  #GaHo
  pred_mu_GaHo<-predict(GaHo,what="mu",newdata=newdata_Ho,type="response")
  pred_sigma_GaHo<-predict(GaHo,what="sigma",newdata=newdata_Ho,type="response")
  
  
  #LNHo
  pred_mu_LNHo<-predict(LNHo,what="mu",newdata=newdata_Ho,type="response")
  pred_sigma_LNHo<-predict(LNHo,what="sigma",newdata=newdata_Ho,type="response")
  
  #ODPCal
  pred_ODPCal<-predict.glm(ODPCal,newdata=newdata,type="response",se.fit=TRUE)
  pred_mu_ODPCal<- pred_ODPCal$fit
  pred_phi_ODPCal<-(pred_ODPCal$residual.scale)^2
  
  #GaCAL
  pred_mu_GaCAL<-predict(GaCAL,what="mu",newdata=newdata,type="response")
  pred_sigma_GaCAL<-predict(GaCAL,what="sigma",newdata=newdata,type="response")
  
  #LNCAL
  pred_mu_LNCAL<-predict(LNCAL,what="mu",newdata=newdata,type="response")
  pred_sigma_LNCAL<-predict(LNCAL,what="sigma",newdata=newdata,type="response")
  
  #NoSp
  newdata_numeric<-newdata
  newdata_numeric$origin=as.numeric(as.character(newdata$origin))
  newdata_numeric$dev=as.numeric(as.character(newdata$dev))
  newdata_numeric$Calendar=as.numeric(as.character(newdata$Calendar))
  
  pred_mu_NoSp<-predict(NoSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_NoSp<-predict(NoSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #GaSp
  pred_mu_GaSp<-predict(GaSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_GaSp<-predict(GaSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #LNSP
  pred_mu_LNSp<-predict(LNSp,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_LNSp<-predict(LNSp,what="sigma",newdata=newdata_numeric,type="response")
  
  #GaGAMLSS
  pred_mu_GaGAMLSS<-predict(GaGAMLSS,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_GaGAMLSS<-predict(GaGAMLSS,what="sigma",newdata=newdata_numeric,type="response")
  
  
  #LNGAMLSS
  pred_mu_LNGAMLSS<-predict(LNGAMLSS,what="mu",newdata=newdata_numeric,type="response")
  pred_sigma_LNGAMLSS<-predict(LNGAMLSS,what="sigma",newdata=newdata_numeric,type="response")
  
  
  #PPCI
  N_rep_valid<-rep(NO_PPCI[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_PPCI<-predict.glm(PPCI,newdata=newdata,type="response",se.fit=TRUE)
  pred_mu_PPCI<-pred_PPCI$fit*N_rep_valid
  pred_phi_PPCI<-(pred_PPCI$residual.scale)^2*(N_rep_valid)^2
  
  #PPCF
  N_rep_valid<-rep(NO_PPCF[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
  
  pred_F<-predict(odp_FC,newdata=newdataFC,type="response")
  pred_F_cum<-list()
  for(i in index_start:index_end){
    pf<-predict(odp_FC,newdata=newdataFC[newdataFC$origin==i,],type="response")
    pf_cum<-round(cumsum(c(train_cumF[train_cumF$origin==i,]$value[nrow(train_cumF[train_cumF$origin==i,])],pf))[-1],0)
    pred_F_cum[[i]]<-pf_cum
  }
  pred_F_cum<-as.vector(unlist(pred_F_cum[index_start:index_end]))
  
  pred_OT<-pred_F_cum/N_rep_valid
  pred_payment<-predict.glm(model_subPayments,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
  pred_mu_PPCF<-pred_payment$fit*pred_F
  pred_phi_PPCF<-mean((pred_payment$residual.scale)^2*(pred_F)^2)
  
  
  y=newdata$value[index]
  CDF<-mapply(FUN=function(y,pred_ODP_mu,pred_ODP_phi,pred_mu_GAGLM,pred_sigma_GAGLM,pred_mu_LNGLM,pred_sigma_LNGLM,pred_mu_ZAGA,pred_sigma_ZAGA,pred_nu_ZAGA,pred_mu_ZALN,pred_sigma_ZALN,pred_nu_ZALN,pred_mu_ODPHo,pred_phi_ODPHo,pred_mu_GaHo,pred_sigma_GaHo,pred_mu_LNHo,pred_sigma_LNHo,pred_mu_ODPCal,pred_phi_ODPCal,pred_mu_GaCAL,pred_sigma_GaCAL,pred_mu_LNCAL,pred_sigma_LNCAL,pred_mu_NoSp,pred_sigma_NoSp,pred_mu_GaSp,pred_sigma_GaSp,pred_mu_LNSp,pred_sigma_LNSp,pred_mu_GaGAMLSS,pred_sigma_GaGAMLSS,pred_mu_LNGAMLSS,pred_sigma_LNGAMLSS,pred_mu_PPCI,pred_phi_PPCI,pred_mu_PPCF,pred_phi_PPC) (cbind(ptweedie(q=z,mu=pred_ODP_mu,phi=pred_ODP_phi,power=1),pGA(q=z,pred_mu_GAGLM,pred_sigma_GAGLM),pLNO(q=z,pred_mu_LNGLM,pred_sigma_LNGLM),pZAGA(q=z,pred_mu_ZAGA,pred_sigma_ZAGA,pred_nu_ZAGA),pZALN(q=z,pred_mu_ZALN,pred_sigma_ZALN,pred_nu_ZALN),ptweedie(q=z,mu=pred_mu_ODPHo,phi=pred_phi_ODPHo,power=1),pGA(q=z,pred_mu_GaHo,pred_sigma_GaHo),pLNO(q=z,pred_mu_LNHo,pred_sigma_LNHo),ptweedie(q=z,mu=pred_mu_ODPCal,phi=pred_phi_ODPCal,power=1),pGA(q=z,pred_mu_GaCAL,pred_sigma_GaCAL),pLNO(q=z,pred_mu_LNCAL,pred_sigma_LNCAL),pNO(q=z,pred_mu_NoSp,pred_sigma_NoSp),pGA(q=z,pred_mu_GaSp,pred_sigma_GaSp),pLNO(q=z,pred_mu_LNSp,pred_sigma_LNSp),pGA(q=z,pred_mu_GaGAMLSS,pred_sigma_GaGAMLSS),pLNO(q=z,pred_mu_LNGAMLSS,pred_sigma_LNGAMLSS),ptweedie(q=z,mu=pred_mu_PPCI,phi=pred_phi_PPCI,power=1),ptweedie(q=z,mu=pred_mu_PPCF,phi=pred_phi_PPCF,power=1))%*%w),y,pred_ODP_mu[index],pred_ODP_phi,pred_mu_GAGLM[index],pred_sigma_GAGLM[index],pred_mu_LNGLM[index],pred_sigma_LNGLM[index],pred_mu_ZAGA[index],pred_sigma_ZAGA[index],pred_nu_ZAGA[index],pred_mu_ZALN[index],pred_sigma_ZALN[index],pred_nu_ZALN[index],pred_mu_ODPHo[index],pred_phi_ODPHo,pred_mu_GaHo[index],pred_sigma_GaHo[index],pred_mu_LNHo[index],pred_sigma_LNHo[index],pred_mu_ODPCal[index],pred_phi_ODPCal,pred_mu_GaCAL[index],pred_sigma_GaCAL[index],pred_mu_LNCAL[index],pred_sigma_LNCAL[index],pred_mu_NoSp[index],pred_sigma_NoSp[index],pred_mu_GaSp[index],pred_sigma_GaSp[index],pred_mu_LNSp[index],pred_sigma_LNSp[index],pred_mu_GaGAMLSS[index],pred_sigma_GaGAMLSS[index],pred_mu_LNGAMLSS[index],pred_sigma_LNGAMLSS[index],pred_mu_PPCI[index],pred_phi_PPCI[index],pred_mu_PPCF[index],pred_phi_PPCF)
  
  return (CDF)
}

cal_CDF_LN<-function(z,model,newdata,index){
  y<- newdata$value[index]
  mu<-predict(model,what="mu", newdata=newdata[index,],type="response")
  sigma<-predict(model,what="sigma",newdata=newdata[index,],type="response")
  CDF<-pLOGNO(q=z,mu=mu,sigma=sigma)
  #for each data point, calculate its CRPS over the range of interest
  return(CDF)
}






