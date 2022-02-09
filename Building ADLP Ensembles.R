#Define the index for subset
ind_subset1_par1_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=28
ind_subset2_par1_new<-as.numeric(as.character(out_sample$origin))>=29&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par2_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=31
ind_subset2_par2_new<-as.numeric(as.character(out_sample$origin))>=32&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par3_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=33
ind_subset2_par3_new<-as.numeric(as.character(out_sample$origin))>=34&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par4_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=26
ind_subset2_par4_new<-as.numeric(as.character(out_sample$origin))>=27&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par5_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=23
ind_subset2_par5_new<-as.numeric(as.character(out_sample$origin))>=24&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par6_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=19
ind_subset2_par6_new<-as.numeric(as.character(out_sample$origin))>=20&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par7_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=15
ind_subset2_par7_new<-as.numeric(as.character(out_sample$origin))>=16&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par8_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=13
ind_subset2_par8_new<-as.numeric(as.character(out_sample$origin))>=14&as.numeric(as.character(out_sample$origin))<=40





#Partition Strategy 1

OW_out_dens_par1_new<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par1_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=28,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=28,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=29&as.numeric(as.character(meta_dens_out$origin))<=40,]
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:2){
    mat_valid<-as.matrix(meta_valid_list_aug_new1[[i]][,3:ncol(meta_valid_list_aug_new1[[i]])])
    mat_test<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
    w_MM_aug_new1[[i]]<-optim_MM$params
    finalw_MM_aug_new1[[i]]<-optim_MM$finalparams
    TrainL_MM_aug_new1[[i]]<-optim_MM$NLL
    FinalTrainL_MM_aug_new1[[i]]<-optim_MM$finalNLL
    TestL_MM_aug_new1[[i]]<--log(mat_test%*%finalw_MM_aug_new1[[i]])
  }
  Out_sample_EnDens_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_par1_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par1_new[[D]]<-finalw_MM_aug_new1
  
}


#Partion Strategy 2

OW_out_dens_par2_new<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par2_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=31,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=31,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=32&as.numeric(as.character(meta_dens_out$origin))<=40,]
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:2){
    mat_valid<-as.matrix(meta_valid_list_aug_new1[[i]][,3:ncol(meta_valid_list_aug_new1[[i]])])
    mat_test<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
    w_MM_aug_new1[[i]]<-optim_MM$params
    finalw_MM_aug_new1[[i]]<-optim_MM$finalparams
    TrainL_MM_aug_new1[[i]]<-optim_MM$NLL
    FinalTrainL_MM_aug_new1[[i]]<-optim_MM$finalNLL
    TestL_MM_aug_new1[[i]]<--log(mat_test%*%finalw_MM_aug_new1[[i]])
  }
  Out_sample_EnDens_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_par2_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par2_new[[D]]<-finalw_MM_aug_new1
  
}

#Partition Strategy 3


OW_out_dens_par3_new<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par3_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=33,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=33,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=34&as.numeric(as.character(meta_dens_out$origin))<=40,]
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:2){
    mat_valid<-as.matrix(meta_valid_list_aug_new1[[i]][,3:ncol(meta_valid_list_aug_new1[[i]])])
    mat_test<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
    w_MM_aug_new1[[i]]<-optim_MM$params
    finalw_MM_aug_new1[[i]]<-optim_MM$finalparams
    TrainL_MM_aug_new1[[i]]<-optim_MM$NLL
    FinalTrainL_MM_aug_new1[[i]]<-optim_MM$finalNLL
    TestL_MM_aug_new1[[i]]<--log(mat_test%*%finalw_MM_aug_new1[[i]])
  }
  Out_sample_EnDens_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_par3_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par3_new[[D]]<-finalw_MM_aug_new1
  
}


#Partition Strategy 4

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par4_new<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par4_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=26,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=26,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=27&as.numeric(as.character(meta_dens_out$origin))<=40,]
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:2){
    mat_valid<-as.matrix(meta_valid_list_aug_new1[[i]][,3:ncol(meta_valid_list_aug_new1[[i]])])
    mat_test<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
    w_MM_aug_new1[[i]]<-optim_MM$params
    finalw_MM_aug_new1[[i]]<-optim_MM$finalparams
    TrainL_MM_aug_new1[[i]]<-optim_MM$NLL
    FinalTrainL_MM_aug_new1[[i]]<-optim_MM$finalNLL
    TestL_MM_aug_new1[[i]]<--log(mat_test%*%finalw_MM_aug_new1[[i]])
  }
  Out_sample_EnDens_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_par4_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par4_new[[D]]<-finalw_MM_aug_new1
  
}


#Partition Strategy 5

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par5_new<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par5_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=23,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=23,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=24&as.numeric(as.character(meta_dens_out$origin))<=40,]
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:2){
    mat_valid<-as.matrix(meta_valid_list_aug_new1[[i]][,3:ncol(meta_valid_list_aug_new1[[i]])])
    mat_test<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
    w_MM_aug_new1[[i]]<-optim_MM$params
    finalw_MM_aug_new1[[i]]<-optim_MM$finalparams
    TrainL_MM_aug_new1[[i]]<-optim_MM$NLL
    FinalTrainL_MM_aug_new1[[i]]<-optim_MM$finalNLL
    TestL_MM_aug_new1[[i]]<--log(mat_test%*%finalw_MM_aug_new1[[i]])
  }
  Out_sample_EnDens_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_par5_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par5_new[[D]]<-finalw_MM_aug_new1
  
}




#Partition Strategy 6

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par6_new<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par6_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=19,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=19,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=20&as.numeric(as.character(meta_dens_out$origin))<=40,]
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:2){
    mat_valid<-as.matrix(meta_valid_list_aug_new1[[i]][,3:ncol(meta_valid_list_aug_new1[[i]])])
    mat_test<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
    w_MM_aug_new1[[i]]<-optim_MM$params
    finalw_MM_aug_new1[[i]]<-optim_MM$finalparams
    TrainL_MM_aug_new1[[i]]<-optim_MM$NLL
    FinalTrainL_MM_aug_new1[[i]]<-optim_MM$finalNLL
    TestL_MM_aug_new1[[i]]<--log(mat_test%*%finalw_MM_aug_new1[[i]])
  }
  Out_sample_EnDens_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_par6_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par6_new[[D]]<-finalw_MM_aug_new1
  
}



#Partition Strategy 7

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par7_new<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par7_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=15,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=15,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=16&as.numeric(as.character(meta_dens_out$origin))<=40,]
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:2){
    mat_valid<-as.matrix(meta_valid_list_aug_new1[[i]][,3:ncol(meta_valid_list_aug_new1[[i]])])
    mat_test<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
    w_MM_aug_new1[[i]]<-optim_MM$params
    finalw_MM_aug_new1[[i]]<-optim_MM$finalparams
    TrainL_MM_aug_new1[[i]]<-optim_MM$NLL
    FinalTrainL_MM_aug_new1[[i]]<-optim_MM$finalNLL
    TestL_MM_aug_new1[[i]]<--log(mat_test%*%finalw_MM_aug_new1[[i]])
  }
  Out_sample_EnDens_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_par7_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par7_new[[D]]<-finalw_MM_aug_new1
  
}




#Partition Strategy 8

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par8_new<-matrix(NA,nrow=780,ncol=ntri)
model_weights_simul_par8_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=13,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=13,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=14&as.numeric(as.character(meta_dens_out$origin))<=40,]
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:2){
    mat_valid<-as.matrix(meta_valid_list_aug_new1[[i]][,3:ncol(meta_valid_list_aug_new1[[i]])])
    mat_test<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    optim_MM<-MM_optim(w_init_1_3,dat=mat_valid,testdat=mat_test,nIters=500)
    w_MM_aug_new1[[i]]<-optim_MM$params
    finalw_MM_aug_new1[[i]]<-optim_MM$finalparams
    TrainL_MM_aug_new1[[i]]<-optim_MM$NLL
    FinalTrainL_MM_aug_new1[[i]]<-optim_MM$finalNLL
    TestL_MM_aug_new1[[i]]<--log(mat_test%*%finalw_MM_aug_new1[[i]])
  }
  Out_sample_EnDens_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_par8_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par8_new[[D]]<-finalw_MM_aug_new1
  
}