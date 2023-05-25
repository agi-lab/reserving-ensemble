cal_CRPS_FMatrix<-function(z,newdata,ODPGLM,GAGLM,LNGLM,ZAGA,ZALN,ODPHo,GaHo,LNHo,ODPCal,GaCAL,LNCAL,NoSp,GaSp,LNSp,GaGAMLSS,LNGAMLSS,PPCI,NO_PPCI,index_start,index_end,odp_FC,model_subPayments,train_cumF,newdataFC,newdata_Pay,NO_PPCF)
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
  
####Create a matrix to store the predictive cummulative probabilities  
F_mat_list<-list()  
for (i in 1:length(y)){
F_mat_list[[i]]<-cbind(ptweedie(q=z,mu=pred_ODP_mu[i],phi=pred_ODP_phi,power=1),
      pGA(q=z,pred_mu_GAGLM[i],pred_sigma_GAGLM[i]),
      pLNO(q=z,pred_mu_LNGLM[i],pred_sigma_LNGLM)[i],
      pZAGA(q=z,pred_mu_ZAGA[i],pred_sigma_ZAGA[i],pred_nu_ZAGA[i]),
      pZALN(q=z,pred_mu_ZALN[i],pred_sigma_ZALN[i],pred_nu_ZALN[i]),
      ptweedie(q=z,mu=pred_mu_ODPHo[i],phi=pred_phi_ODPHo,power=1),
      pGA(q=z,pred_mu_GaHo[i],pred_sigma_GaHo[i]),
      pLNO(q=z,pred_mu_LNHo[i],pred_sigma_LNHo[i]),
      ptweedie(q=z,mu=pred_mu_ODPCal[i],phi=pred_phi_ODPCal,power=1),
      pGA(q=z,pred_mu_GaCAL[i],pred_sigma_GaCAL[i]),
      pLNO(q=z,pred_mu_LNCAL[i],pred_sigma_LNCAL[i]),
      pNO(q=z,pred_mu_NoSp[i],pred_sigma_NoSp[i]),
      pGA(q=z,pred_mu_GaSp[i],pred_sigma_GaSp[i]),
      pLNO(q=z,pred_mu_LNSp[i],pred_sigma_LNSp[i]),
      pGA(q=z,pred_mu_GaGAMLSS[i],pred_sigma_GaGAMLSS[i]),
      pLNO(q=z,pred_mu_LNGAMLSS[i],pred_sigma_LNGAMLSS[i]),
      ptweedie(q=z,mu=pred_mu_PPCI[i],phi=pred_phi_PPCI[i],power=1),
      ptweedie(q=z,mu=pred_mu_PPCF[i],phi=pred_phi_PPCF,power=1))
} 
F_matrix<-do.call(rbind, F_mat_list)


####Create a vector to store the indicator values
I_list<-list()
for (i in 1:length(y)){
  I_list[[i]]<-ifelse(z>=y[i],1,0)
}

I_ECDF<-as.vector(unlist(I_list))

  return (list(F_matrix=F_matrix,I_ECDF=I_ECDF))
}





####################################################
#------------------------------------------------------------------------------------------------
##Some testing of codes

####################################################

test<-cal_CRPS_FMatrix(z=z,
                      newdata=out_sample,
                      ODPGLM=ODP_GLM_in,
                      GAGLM=Ga_optimTau,
                      LNGLM=LN_optimTau,
                      ZAGA=ZAGA_in,
                      ZALN=LN_in,
                      ODPHo=glm_ODP_Ho_In,
                      GaHo=glm_Ga_Ho_In,
                      LNHo=glm_LN_Ho_In,
                      ODPCal=glm_ODP_Cal_In,
                      GaCAL=glm_Ga_Cal_In,
                      LNCAL=glm_LN_Cal_In,
                      NoSp=sp_Normal_In,
                      GaSp=sp_Gamma_In,
                      LNSp=sp_LN_In,
                      GaGAMLSS=gamlss2_GA_In,
                      LNGAMLSS=gamlss2_LN_In,
                      PPCI=fit_ODP_ppci_In,
                      NO_PPCI=N_in,
                      index_start=2,
                      index_end=40,odp_FC=odp_FC_In,model_subPayments=ODP_PPCF_In,train_cumF=cum_F_dt,newdataFC=FUN_dat_Out,newdata_Pay=out_sample,NO_PPCF=N_in)

F_mat<-test$F_matrix
I<-test$I_ECDF
w_init<-rep(1/18,18)

###Calculate the crps using the above function
crps_ind<-as.vector(F_mat%*%w_init-I)
crps_test<-(crps_ind%*%crps_ind)/nrow(out_sample)
crps_test