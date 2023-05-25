##Specify the number of simulations: can be altered by the users
ntri<-100


count_by_origin <- function(data,start,end){
  cbo<-c()
  for (i in start:end){
    cbo[i]<-nrow(data[data$origin==i,]) 
  }
  return(cbo)
}


#Writing the Minorization-Maximization Algorithm


MM_func<-function(w,dat){
  wd<-dat%*%w
  dens<-matrix(NA,nrow=nrow(dat),ncol=ncol(dat))
  for (i in 1:nrow(dat)){
    dens[i,]<-dat[i,]/wd[i]
  }
  mean<-apply(dens,MARGIN=2,FUN=mean)
  return(mean) 
}

MM_optim<-function(w_init,dat,testdat,nIters){
  w<-matrix(NA,nrow=nIters+1,ncol=length(w_init))
  TL<-c()
  TestL<-c()
  w[1,]<-w_init
  TL[1]<--mean(log(dat%*%w_init))
  TestL[1]<--mean(log(testdat%*%w_init))
  for (i in 2:(nIters+1)){
    w[i,]<-w[i-1,]*MM_func(w[i-1,],dat)
    TL[i]<--mean(log(dat%*%w[i,]))
    TestL[i]<--mean(log(testdat%*%w[i,]))
  }
  return(list(params=w,NLL=TL,TestNLL=TestL,finalparams=w[nIters+1,],finalNLL=TL[nIters+1]))
}



#create a ODP density function
dODP<-function(y,lambda,phi) {
  x<-floor(y/phi)
  lambda2<-floor(lambda/phi)
  dpois(x,lambda2)/phi
}
#Create density function for ZALN
dZALN<-Zadj.d(family="LOGNO")
pZALN<-Zadj.p(family="LOGNO")




#Create functions that calculate the predictive density for each component model

##Calculate ODP density from ODP GLM
cal_dens_ODP<-function(model,newdata){
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu<-pred_ODP$fit
  pred_ODP_phi<-(pred_ODP$residual.scale)^2
  y<-newdata$value
  return(dODP(y,lambda=pred_ODP_mu,phi=pred_ODP_phi))
}

##Calculate Gamma density from Gamma GLM
cal_dens_GA<-function(model,tau,newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(dGA(x=newdata$value+tau,mu=pred_mu,sigma=pred_sigma))
}
##Calculate Log-Normal density from Log-Normal GLM
cal_dens_LN<-function(model,tau,newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(dLNO(x=newdata$value+tau,mu=pred_mu,sigma=pred_sigma))
}
##Calculate ZAGA density from ZAGA GLM
cal_dens_ZAGA<-function(model,newdata){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  pred_nu<-predict(model,what="nu",newdata=newdata_nu,type="response")
  return(dZAGA(x=newdata$value,mu=pred_mu,sigma=pred_sigma,nu=pred_nu))
}

##Calculate ZALN density from ZALN GLM
cal_dens_ZALN<-function(model,newdata){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu<-predict(model,parameter="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,parameter="sigma",newdata=newdata,type="response")
  pred_nu<-predict(model,parameter="xi0",newdata=newdata_nu,type="response")
  return(dZALN(x=newdata$value,mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))
}
##Calculate Normal density from Normal GLM
cal_dens_Normal<-function(model, newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(dNO(x=newdata$value,mu=pred_mu,sigma=pred_sigma))
}
##Calculate ODP density from PPCI model
cal_dens_PPCI<-function(model,N,newdata,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  mu<-pred_ODP$fit*N_rep_valid
  y<-newdata$value
  dens<-c()
  for (i in 1:nrow(newdata)){
    phi<-(pred_ODP$residual.scale)^2*(N_rep_valid[i])^2
    dens[i]<-dtweedie(y=y[i],mu=mu[i],phi=phi,power=1)/phi
  }
  return(dens)
}
##Calculate ODP density from PPCF GLM
cal_PPCF_dens<-function(model_subCount,model_subPayments,train_cumF,newdataFC,newdata_Pay,N,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
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
  mu<-pred_payment$fit*pred_F
  phi<-(pred_payment$residual.scale)^2*(pred_F)^2
  y<-newdata_Pay$value
  
  dens_PPCF<-c()
  for (i in 1:nrow(newdata_Pay)){
    ds<-dtweedie(y=y[i],mu=mu[i],phi=phi[i],power=1)/phi[i]
    dens_PPCF[i]<-ds
  }
  dens_PPCF[dens_PPCF==0]=min(dens_PPCF[dens_PPCF!=0])
  #assign a small number to the zero density in order to prevent -Inf Log Score
  return(dens_PPCF)
}

##Calculate Log-Normal density for Log-Normal GAMLSS
cal_dens_LN_Gamlss<-function(model,tau,newdata_mu,newdata_sigma){
  pred_mu<-predict(model,what="mu",newdata=newdata_mu,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  return(dLNO(x=newdata_mu$value+tau,mu=pred_mu,sigma=pred_sigma))
}

##Calculate Gamma density for Gamma GAMLSS
cal_dens_GA_Gamlss<-function(model,tau,newdata_mu,newdata_sigma){
  pred_mu<-predict(model,what="mu",newdata=newdata_mu,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  return(dGA(x=newdata_mu$value+tau,mu=pred_mu,sigma=pred_sigma))
}





#Create functions that calculate the predictive mean for each component model

##Calculate ODP mean for ODP GLM
cal_mu_ODP<-function(model,newdata){
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu<-pred_ODP$fit
  return(pred_ODP_mu)
}
##Calculate GA mean for GA GLM
cal_mu_GA<-function(model,tau,newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  return(pred_mu)
}

##Calculate LN mean for LN GLM
cal_mu_LN<-function(model,tau,newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  mean<-exp(pred_mu+pred_sigma^2/2)
  return(mean)
}

##Calculate ZAGA mean for ZAGA GLM
cal_mu_ZAGA<-function(model,newdata){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_nu<-predict(model,what="nu",newdata=newdata_nu,type="response")
  mean<-(1-pred_nu)*pred_mu
  return(mean)
}

##Calculate ZALN mean for ZALN GLM
cal_mu_ZALN<-function(model,newdata){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu<-predict(model,parameter="mu",newdata=newdata,type="response")
  pred_nu<-predict(model,parameter="xi0",newdata=newdata_nu,type="response")
  pred_sigma<-predict(model,parameter="sigma",newdata=newdata,type="response")
  mean<-(1-pred_nu)*exp(pred_mu+pred_sigma^2/2)
  return(mean)
}

##Calculate Normal mean for Normal GLM
cal_mu_Normal<-function(model, newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  return(pred_mu)
}

##Calculate ODP mean for PPCI model
cal_mu_PPCI<-function(model,N,newdata,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  mu<-pred_ODP$fit*N_rep_valid
  return(mu)
}

##Calculate ODP mean for PPCF model
cal_PPCF_mu<-function(model_subCount,model_subPayments,train_cumF,newdataFC,newdata_Pay,N,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
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
  mu<-pred_payment$fit*pred_F
  return(mu)
}

##Calculate GA mean for GA GAMLSS

cal_mu_GA_Gamlss<-function(model,tau,newdata_mu,newdata_sigma){
  pred_mu<-predict(model,what="mu",newdata=newdata_mu,type="response")
  return(pred_mu)
}

##Calculate LN mean for LN GAMLSS
cal_mu_LN_Gamlss<-function(model,tau,newdata_mu,newdata_sigma){
  pred_mu<-predict(model,what="mu",newdata=newdata_mu,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  mean<-exp(pred_mu+1/2*pred_sigma^2)
  return(mean)
}









# Create function that calculate the CDF values -------------------------------------------------------------------------

#create a ODP CDF function
pODP<-function(y,lambda,phi) {
  x<-floor(y/phi)
  lambda2<-floor(lambda/phi)
  ppois(x,lambda2)/phi
}

pODP_pre <-function(y,lambda,phi) {
  x <- floor(y/phi)-1
  lambda2<-floor(lambda/phi)
  ppois(x,lambda2)/phi
}



#Create CDF function for ZALN
dZALN<-Zadj.d(family="LOGNO")
pZALN<-Zadj.p(family="LOGNO")




#Create functions that calculate the predictive CDF for each component model

##Calculate ODP CDF from ODP GLM
cal_CDF_ODP<-function(model,newdata){
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu<-pred_ODP$fit
  pred_ODP_phi<-(pred_ODP$residual.scale)^2
  y<-newdata$value
  return(pODP(y,lambda=pred_ODP_mu,phi=pred_ODP_phi))
}

cal_CDF_ODP_pre <- function(model,newdata){
  pred_ODP <- predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu <- pred_ODP$fit
  pred_ODP_phi <- (pred_ODP$residual.scale)^2
  y<-newdata$value
  return(pODP_pre(y,lambda=pred_ODP_mu,phi=pred_ODP_phi))
}



##Calculate Gamma CDF from Gamma GLM
cal_CDF_GA<-function(model,tau,newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(pGA(q=newdata$value+tau,mu=pred_mu,sigma=pred_sigma))
}

##Calculate Log-Normal CDF from Log-Normal GLM
cal_CDF_LN<-function(model,tau,newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(pLNO(q=newdata$value+tau,mu=pred_mu,sigma=pred_sigma))
}
##Calculate ZAGA CDF from ZAGA GLM
cal_CDF_ZAGA<-function(model,newdata){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  pred_nu<-predict(model,what="nu",newdata=newdata_nu,type="response")
  return(pZAGA(q=newdata$value,mu=pred_mu,sigma=pred_sigma,nu=pred_nu))
}

##Calculate ZALN CDF from ZALN GLM
cal_CDF_ZALN<-function(model,newdata){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_mu<-predict(model,parameter="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,parameter="sigma",newdata=newdata,type="response")
  pred_nu<-predict(model,parameter="xi0",newdata=newdata_nu,type="response")
  return(pZALN(q=newdata$value,mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))
}
##Calculate Normal CDF from Normal GLM
cal_CDF_Normal<-function(model, newdata){
  pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(pNO(q=newdata$value,mu=pred_mu,sigma=pred_sigma))
}
##Calculate ODP CDF from PPCI model

cal_CDF_PPCI_pre<-function(model,N,newdata,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  mu<-pred_ODP$fit*N_rep_valid
  y<-newdata$value
  CDF<-c()
  for (i in 1:nrow(newdata)){
    phi<-(pred_ODP$residual.scale)^2*(N_rep_valid[i])^2
    CDF[i]<-ptweedie(q=y[i]-1,mu=mu[i],phi=phi,power=1)/phi
  }
  return(CDF)
}


cal_CDF_PPCI<-function(model,N,newdata,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  mu<-pred_ODP$fit*N_rep_valid
  y<-newdata$value
  CDF<-c()
  for (i in 1:nrow(newdata)){
    phi<-(pred_ODP$residual.scale)^2*(N_rep_valid[i])^2
    CDF[i]<-ptweedie(q=y[i],mu=mu[i],phi=phi,power=1)/phi
  }
  return(CDF)
}
##Calculate ODP CDF from PPCF GLM
cal_PPCF_CDF_pre<-function(model_subCount,model_subPayments,train_cumF,newdataFC,newdata_Pay,N,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
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
  mu<-pred_payment$fit*pred_F
  phi<-(pred_payment$residual.scale)^2*(pred_F)^2
  y<-newdata_Pay$value
  
  CDF_PPCF<-c()
  for (i in 1:nrow(newdata_Pay)){
    ds<-ptweedie(q=y[i]-1,mu=mu[i],phi=phi[i],power=1)/phi[i]
    CDF_PPCF[i]<-ds
  }
  CDF_PPCF[CDF_PPCF==0]=min(CDF_PPCF[CDF_PPCF!=0])
  #assign a small number to the zero CDF in order to prevent -Inf Log Score
  return(CDF_PPCF)
}


cal_PPCF_CDF<-function(model_subCount,model_subPayments,train_cumF,newdataFC,newdata_Pay,N,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
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
  mu<-pred_payment$fit*pred_F
  phi<-(pred_payment$residual.scale)^2*(pred_F)^2
  y<-newdata_Pay$value
  
  CDF_PPCF<-c()
  for (i in 1:nrow(newdata_Pay)){
    ds<-ptweedie(q=y[i],mu=mu[i],phi=phi[i],power=1)/phi[i]
    CDF_PPCF[i]<-ds
  }
  CDF_PPCF[CDF_PPCF==0]=min(CDF_PPCF[CDF_PPCF!=0])
  #assign a small number to the zero CDF in order to prevent -Inf Log Score
  return(CDF_PPCF)
}

##Calculate Log-Normal CDF for Log-Normal GAMLSS
cal_CDF_LN_Gamlss<-function(model,tau,newdata_mu,newdata_sigma){
  pred_mu<-predict(model,what="mu",newdata=newdata_mu,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  return(pLNO(q=newdata_mu$value+tau,mu=pred_mu,sigma=pred_sigma))
}

##Calculate Gamma CDF for Gamma GAMLSS
cal_CDF_GA_Gamlss<-function(model,tau,newdata_mu,newdata_sigma){
  pred_mu<-predict(model,what="mu",newdata=newdata_mu,type="response")
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  return(pGA(q=newdata_mu$value+tau,mu=pred_mu,sigma=pred_sigma))
}



# Create function that calcualte the predicted variance:

#Create functions that calculate the predictive mean for each component model

#Create functions that calculate the predictive mean for each component model

##Calculate ODP mean for ODP GLM
cal_sigma_ODP<-function(model,newdata){
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  pred_ODP_mu<-pred_ODP$fit
  pred_ODP_phi<-(pred_ODP$residual.scale)^2
  pred_ODP_sigma<-pred_ODP_phi*pred_ODP_mu
  return(pred_ODP_sigma)
}
##Calculate GA mean for GA GLM
cal_sigma_GA<-function(model,tau,newdata){
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(pred_sigma)
}

##Calculate LN mean for LN GLM
cal_sigma_LN<-function(model,tau,newdata){
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(pred_sigma)
}

##Calculate ZAGA mean for ZAGA GLM
cal_sigma_ZAGA<-function(model,newdata){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(pred_sigma)
}

##Calculate ZALN mean for ZALN GLM
cal_sigma_ZALN<-function(model,newdata){
  newdata_nu<-newdata
  newdata_nu$dev=as.numeric(as.character(newdata$dev))
  pred_sigma <- predict(model,parameter="sigma",newdata=newdata,type="response")
  return(pred_sigma)
}

##Calculate Normal mean for Normal GLM
cal_sigma_Normal<-function(model, newdata){
  pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
  return(pred_sigma)
}

##Calculate ODP mean for PPCI model
cal_sigma_PPCI<-function(model,N,newdata,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata,index_start,index_end)[index_start:index_end])
  pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
  mu<-pred_ODP$fit*N_rep_valid
  y<-newdata$value
  sigma<-c()
  for (i in 1:nrow(newdata)){
    sigma[i]<-(pred_ODP$residual.scale)^2*(N_rep_valid[i])^2*mu[i]
    
  }
  return(sigma)
}

##Calculate ODP mean for PPCF model
##Calculate ODP CDF from PPCF GLM
cal_PPCF_sigma<-function(model_subCount,model_subPayments,train_cumF,newdataFC,newdata_Pay,N,index_start,index_end){
  N_rep_valid<-rep(N[index_start:index_end],times=count_by_origin(newdata_Pay,index_start,index_end)[index_start:index_end])
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
  mu<-pred_payment$fit*pred_F
  phi<-(pred_payment$residual.scale)^2*(pred_F)^2
  return(phi*mu)
}

##Calculate GA mean for GA GAMLSS

cal_sigma_GA_Gamlss<-function(model,tau,newdata_sigma){
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  return(pred_sigma)
}

##Calculate LN mean for LN GAMLSS
cal_sigma_LN_Gamlss<-function(model,tau,newdata_sigma){
  pred_sigma<-predict(model,what="sigma",newdata=newdata_sigma,type="response")
  return(pred_sigma)
}
