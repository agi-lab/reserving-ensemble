# load("index_bestvalid")
# load("out_dens_list_par1")

#out_dens_list<-out_dens_list_par1
#ntri<-100

EW_LS<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
for(i in 1:ntri){
  EW_LS[,i]<-log(out_dens_list[[i]][,21])
}
BMV_LS<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
for(i in 1:ntri){
  BMV_LS[,i]<-log(out_dens_list[[i]][,index_best_validmod2[i]+2])
}



#DM Test for SLP v.s EW and BMV

###DM Test statistics for SLP v.s EW
S_bar_0<-apply(log(OW_out_dens_par0),MARGIN=2,FUN=mean)
G_bar_0<-c()
for (i in 1:ntri){
  G_bar_0[i]<-mean(EW_LS[,i])
}
tn_0_DM_EW<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(log(OW_out_dens_par0))*sum((log(OW_out_dens_par0)[,i]-EW_LS[,i])^2))
  tn_0_DM_EW[i]<-sqrt(nrow(log(OW_out_dens_par0)))*(S_bar_0[i]-G_bar_0[i])/sigma2_1[i]
}


###DMTest statistics for SLP v.s BMV
S_bar_0<-apply(log(OW_out_dens_par0),MARGIN=2,FUN=mean)
G_bar_0<-apply(BMV_LS,MARGIN=2,FUN=mean)
tn_0_DM_BMV<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(log(OW_out_dens_par0))*sum((log(OW_out_dens_par0)[,i]-BMV_LS[,i])^2))
  tn_0_DM_BMV[i]<-sqrt(nrow(log(OW_out_dens_par0)))*(S_bar_0[i]-G_bar_0[i])/sigma2_1[i]
}


#Adjusted DM Test for SLP v.s EW and BMV

## Adjusted DM tests for SLP v.s.EW
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(log(OW_out_dens_par0[,i])-EW_LS[,i])
}
##Caculate the Spectral Density at Zero
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-log(OW_out_dens_par0[,i])-EW_LS[,i]
  fd_0[i]<-spectrum0.ar(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_0_spec_EW<-c()
for(i in 1:ntri){
  tn_0_spec_EW[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(log(OW_out_dens_par0)))
  ###fd_0[i]/nrow(log(OW_out_dens_par0)) provides an estimation of the variance of mean log score differential
}

## Adjusted DM tests for SLP v.s.BMV
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(log(OW_out_dens_par0[,i])-BMV_LS[,i])
}
##Caculate the Spectral Density at Zero
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-log(OW_out_dens_par0[,i])-BMV_LS[,i]
  fd_0[i]<-spectrum0.ar(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_0_spec_BMV<-c()
for(i in 1:ntri){
  tn_0_spec_BMV[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(log(OW_out_dens_par0)))
  ###fd_0[i]/nrow(log(OW_out_dens_par0)) provides an estimation of the variance of mean log score differential
}


#DM Test for Partition Strategy 8 v.s EW and BMV

##DM Test for Partition Strategy 8 v.s EW 

S_bar_8<-apply(log(OW_out_dens_par8_new),MARGIN=2,FUN=mean)
G_bar_8<-c()
for (i in 1:ntri){
  G_bar_8[i]<-mean(EW_LS[,i])
}
tn_8_DM_EW<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(log(OW_out_dens_par8_new))*sum((log(OW_out_dens_par8_new)[,i]-EW_LS[,i])^2))
  tn_8_DM_EW[i]<-sqrt(nrow(log(OW_out_dens_par8_new)))*(S_bar_8[i]-G_bar_8[i])/sigma2_1[i]
}

##DM Test for Partition Strategy 8 v.s BMV
S_bar_8<-apply(log(OW_out_dens_par8_new),MARGIN=2,FUN=mean)
G_bar_8<-apply(BMV_LS,MARGIN=2,FUN=mean)
tn_8_DM_BMV<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(log(OW_out_dens_par8_new))*sum((log(OW_out_dens_par8_new)[,i]-BMV_LS[,i])^2))
  tn_8_DM_BMV[i]<-sqrt(nrow(log(OW_out_dens_par8_new)))*(S_bar_8[i]-G_bar_8[i])/sigma2_1[i]
}

#Adjusted DM Test for Partition Strategy 8 v.s EW and BMV

##Adjusted DM Test for Partition Strategy 8 v.s EW
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(log(OW_out_dens_par8_new[,i])-EW_LS[,i])
}
##Caculate the Spectral Density at Zero
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-log(OW_out_dens_par8_new[,i])-EW_LS[,i]
  fd_0[i]<-spectrum0.ar(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_8_spec_EW<-c()
for(i in 1:ntri){
  tn_8_spec_EW[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(log(OW_out_dens_par8_new)))
  ###fd_0[i]/nrow(log(OW_out_dens_par8_new)) provides an estimation of the variance of mean log score differential
}

##Adjusted DM Test for Partition Strategy 8 v.s BMV
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(log(OW_out_dens_par8_new[,i])-BMV_LS[,i])
}
##Calculate the Spectral Density at Zero
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-log(OW_out_dens_par8_new[,i])-BMV_LS[,i]
  fd_0[i]<-spectrum0.ar(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_8_spec_BMV<-c()
for(i in 1:ntri){
  tn_8_spec_BMV[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(log(OW_out_dens_par8_new)))
  ###fd_0[i]/nrow(log(OW_out_dens_par8_new)) provides an estimation of the variance of mean log score differential
}



#DM Test for Partition Strategy 8 v.s SLP
S_bar_8v0<-apply(log(OW_out_dens_par8_new),MARGIN=2,FUN=mean)
G_bar_8v0<-apply(log(OW_out_dens_par0),MARGIN=2,FUN=mean)
tn_8_DM_vs0<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(log(OW_out_dens_par8_new))*sum((log(OW_out_dens_par8_new)[,i]-log(OW_out_dens_par0)[,i])^2))
  tn_8_DM_vs0[i]<-sqrt(nrow(log(OW_out_dens_par8_new)))*(S_bar_8v0[i]-G_bar_8v0[i])/sigma2_1[i]
}

#Adjusted DM Test for Partition Strategy 8 v.s SLP
####Best Model in Validation Set
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(log(OW_out_dens_par8_new[,i])-log(OW_out_dens_par0)[,i])
}
##Calculate the Spectral Density at Zero
library("coda")
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-log(OW_out_dens_par8_new[,i])-log(OW_out_dens_par0)[,i]
  fd_0[i]<-spectrum0.ar(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_8_spec_8v0<-c()
for(i in 1:ntri){
  tn_8_spec_8v0[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(log(OW_out_dens_par8_new)))
  ###fd_0[i]/nrow(log(OW_out_dens_par8_new)) provides an estimation of the variance of mean log score differential
}






##################################################################
######################Statistical Tests for CRPS######################
#######################################################################


##DM Test for Partition Strategy 8 v.s BMV
S_bar_8<-apply(ADLP8_crps,MARGIN=2,FUN=mean)
G_bar_8<-apply(SpLN_crps,MARGIN=2,FUN=mean)
tn_8_DM_BMV_crps<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(ADLP8_crps)*sum((ADLP8_crps[,i]-SpLN_crps[,i])^2))
  tn_8_DM_BMV_crps[i]<-sqrt(nrow(ADLP8_crps))*(G_bar_8[i]-S_bar_8[i])/sigma2_1[i]
}
mean(tn_8_DM_BMV_crps>qnorm(0.975))
#Adjusted DM Test for Partition Strategy 8 v.s EW and BMV


##Adjusted DM Test for Partition Strategy 8 v.s BMV
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(SpLN_crps[,i]-ADLP8_crps[,i])
}
##Calculate the Spectral Density at Zero
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-SpLN_crps[,i]-ADLP8_crps[,i]
  fd_0[i]<-spectrum0.ar(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_8_spec_BMV_crps<-c()
for(i in 1:ntri){
  tn_8_spec_BMV_crps[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(ADLP8_crps))
  ###fd_0[i]/nrow(log(OW_out_dens_par8_new)) provides an estimation of the variance of mean log score differential
}
mean(tn_8_spec_BMV_crps>qnorm(0.975))

####ADLP5 and ADLP8
S_bar<-apply(ADLP5_crps,MARGIN=2,FUN=mean)
G_bar<-apply(ADLP8_crps,MARGIN=2,FUN=mean)
tn_8vs5_DM_BMV_crps<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(ADLP8_crps)*sum((ADLP8_crps[,i]-ADLP5_crps[,i])^2))
  tn_8vs5_DM_BMV_crps[i]<-sqrt(nrow(ADLP8_crps))*(G_bar[i]-S_bar[i])/sigma2_1[i]
}
mean(tn_8vs5_DM_BMV_crps>qnorm(0.975))


#####ADLP5 and ADLP8 (Adjusted DM)
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(ADLP8_crps[,i]-ADLP5_crps[,i])
}
##Calculate the Spectral Density at Zero
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-ADLP8_crps[,i]-ADLP5_crps[,i]
  fd_0[i]<-spectrum0.ar(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_8_spec_5_crps<-c()
for(i in 1:ntri){
  tn_8_spec_5_crps[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(ADLP8_crps))
  ###fd_0[i]/nrow(log(OW_out_dens_par8_new)) provides an estimation of the variance of mean log score differential
}
mean(tn_8_spec_5_crps>qnorm(0.975))




####SLP and ADLP8
S_bar<-apply(ADLP8_crps,MARGIN=2,FUN=mean)
G_bar<-apply(SLP_crps,MARGIN=2,FUN=mean)
tn_8vsSLP_DM_BMV_crps<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(ADLP8_crps)*sum((ADLP8_crps[,i]-SLP_crps[,i])^2))
  tn_8vsSLP_DM_BMV_crps[i]<-sqrt(nrow(ADLP8_crps))*(G_bar[i]-S_bar[i])/sigma2_1[i]
}
mean(tn_8vsSLP_DM_BMV_crps>qnorm(0.975))


#####SLP and ADLP8 (Adjusted DM)
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(SLP_crps[,i]-ADLP8_crps[,i])
}
##Calculate the Spectral Density at Zero
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-SLP_crps[,i]-ADLP8_crps[,i]
  fd_0[i]<-spectrum0(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_8_spec_SLP_crps<-c()
for(i in 1:ntri){
  tn_8_spec_SLP_crps[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(ADLP8_crps))
  ###fd_0[i]/nrow(log(OW_out_dens_par8_new)) provides an estimation of the variance of mean log score differential
}
mean(tn_8_spec_SLP_crps>qnorm(0.975))





####EW and ADLP8
S_bar<-apply(ADLP8_crps,MARGIN=2,FUN=mean)
G_bar<-apply(EqEns_crps,MARGIN=2,FUN=mean)
tn_8vsEW_DM_crps<-c()
sigma2_1<-c()

for (i in 1:ntri){
  sigma2_1[i]<-sqrt(1/nrow(ADLP8_crps)*sum((ADLP8_crps[,i]-EqEns_crps[,i])^2))
  tn_8vsEW_DM_crps[i]<-sqrt(nrow(ADLP8_crps))*(G_bar[i]-S_bar[i])/sigma2_1[i]
}
mean(tn_8vsEW_DM_crps[tn_8vsEW_DM_crps>=0]>qnorm(0.975))


#####EW and ADLP8 (Adjusted DM)
d_bar<-c()
for (i in 1:ntri){
  d_bar[i]<-mean(EqEns_crps[,i]-ADLP8_crps[,i])
}
##Calculate the Spectral Density at Zero
fd_0<-c()
for (i in 1:ntri){
  loss_diff<-EqEns_crps[,i]-ADLP8_crps[,i]
  fd_0[i]<-spectrum0(loss_diff)$spec
  ##estimate the spectral density at zero by fitting an autoregressive model;
}
##Calculate the test statistic
tn_8_spec_EqEns_crps<-c()
for(i in 1:ntri){
  tn_8_spec_EqEns_crps[i]<-d_bar[i]/sqrt(fd_0[i]/nrow(ADLP8_crps))
  ###fd_0[i]/nrow(log(OW_out_dens_par8_new)) provides an estimation of the variance of mean log score differential
}
mean(tn_8_spec_EqEns_crps[tn_8_spec_EqEns_crps>=0]>qnorm(0.975))


