
OW_out_dens_par0<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
OW_out_CDF_par0 <- matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par0<-list()
w_init_1_3<-rep(1/18,18)

##Caculate the Density of the Optimized Ensemble 
for (D in 1:ntri){
  #D <- 1
  meta_dens_augvalid<-valid_dens_list[[D]][,-21]
  meta_dens_out <- out_dens_list[[D]][,-c(21,22)]
  meta_CDF_out <- out_CDF_list[[D]][,-c(21,22)]
  
  #Train the model weights using the MM Algorithm
  mat_valid<-as.matrix(meta_dens_augvalid[,3:ncol(meta_dens_augvalid)])
  mat_test<-as.matrix(meta_dens_out[,3:ncol(meta_dens_out)])
  optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
  w_MM_aug_new1<-optim_MM$params
  finalw_MM_aug_new1<-optim_MM$finalparams
  TrainL_MM_aug_new1<-optim_MM$NLL
  FinalTrainL_MM_aug_new1<-optim_MM$finalNLL
  TestL_MM_aug_new1<--log(mat_test%*%finalw_MM_aug_new1)
  
  
  out_meta<-as.matrix(meta_dens_out[,3:ncol(meta_dens_out)])
  out_meta_CDF <- as.matrix(meta_CDF_out[,3:ncol(meta_CDF_out)])
  Out_sample_EnDens_aug_new1<-out_meta%*%as.vector(finalw_MM_aug_new1)
  Out_sample_EnCDF_aug_new1 <- out_meta_CDF%*%as.vector(finalw_MM_aug_new1)
  
  #Calculate the predictive density by the ensemble
  OW_out_dens_par0[,D]<-as.vector(Out_sample_EnDens_aug_new1)
  OW_out_CDF_par0[,D]<-as.vector(Out_sample_EnCDF_aug_new1)
  model_weights_simul_par0[[D]]<-finalw_MM_aug_new1
  
}

### Calculate the Log-score associated with Partion strategy 0:

origin<-as.numeric(as.character(out_sample$origin))
dev<-as.numeric(as.character(out_sample$dev))

OW_out_LS_par0_dat <- as.data.frame(cbind(origin,dev,OW_out_dens_par0))

OW_out_LS_par0_dat <- as.data.frame(cbind(origin,dev,OW_out_dens_par0))
LS_acc_OW_par0 <- matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par0[i-1,]<-apply(log(OW_out_LS_par0_dat[OW_out_LS_par0_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}

LS_acc_OW_par0_avg <- apply(LS_acc_OW_par0,MARGIN=1,FUN=mean)


#save(model_weights_simul_par0, file = "model_weights_simul_par0")
#save(OW_out_dens_par0, file = "OW_out_dens_par0")
