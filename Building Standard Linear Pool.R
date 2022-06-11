
OW_out_dens_par0<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par0<-list()
w_init_1_3<-rep(1/18,18)
##Caculate the Density of the Optimized Ensemble 
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-21]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
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
  Out_sample_EnDens_aug_new1<-out_meta%*%as.vector(finalw_MM_aug_new1)
  
  #Calculate the predictive density by the ensemble
  OW_out_dens_par0[,D]<-as.vector(Out_sample_EnDens_aug_new1)
  model_weights_simul_par0[[D]]<-finalw_MM_aug_new1
  
}