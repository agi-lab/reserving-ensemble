##Caculate the Density of the Optimized Ensemble 
#ntri<-100


##Partition Strategy 1
OW_out_dens_3par_PS1<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_3par_PS1<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=5,]
  augvalid_meta_2_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=15,]
  augvalid_meta_3_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=40,]
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=5,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=6&as.numeric(as.character(meta_dens_out$origin))<=14,]
  out_meta_3_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=15&as.numeric(as.character(meta_dens_out$origin))<=40,]
  
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1,augvalid_meta_3_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1,out_meta_3_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:3){
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
  for (i in 1:3){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_3par_PS1[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_3par_PS1[[D]]<-finalw_MM_aug_new1
}




##Partition Strategy 2
OW_out_dens_3par_PS2<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_3par_PS2<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=15,]
  augvalid_meta_2_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=29,]
  augvalid_meta_3_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=40,]
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=15,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=16&as.numeric(as.character(meta_dens_out$origin))<=29,]
  out_meta_3_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=30&as.numeric(as.character(meta_dens_out$origin))<=40,]
  
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1,augvalid_meta_3_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1,out_meta_3_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:3){
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
  for (i in 1:3){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_3par_PS2[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_3par_PS2[[D]]<-finalw_MM_aug_new1
}


##Partition Strategy 3
OW_out_dens_3par_PS3<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_3par_PS3<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=17,]
  augvalid_meta_2_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=31,]
  augvalid_meta_3_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=40,]
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=17,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=18&as.numeric(as.character(meta_dens_out$origin))<=31,]
  out_meta_3_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=32&as.numeric(as.character(meta_dens_out$origin))<=40,]
  
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1,augvalid_meta_3_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1,out_meta_3_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:3){
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
  for (i in 1:3){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_3par_PS3[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_3par_PS3[[D]]<-finalw_MM_aug_new1
}


##Partition Strategy 4
OW_out_dens_3par_PS4<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_3par_PS4<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=23,]
  augvalid_meta_2_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=33,]
  augvalid_meta_3_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=40,]
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=23,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=24&as.numeric(as.character(meta_dens_out$origin))<=33,]
  out_meta_3_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=34&as.numeric(as.character(meta_dens_out$origin))<=40,]
  
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1,augvalid_meta_3_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1,out_meta_3_aug_new1)
  TrainL_MM_aug_new1<-list()
  FinalTrainL_MM_aug_new1<-list()
  TestL_MM_aug_new1<-list()
  w_MM_aug_new1<-list()
  finalw_MM_aug_new1<-list()
  for(i in 1:3){
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
  for (i in 1:3){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
  }
  OW_out_dens_3par_PS4[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_3par_PS4[[D]]<-finalw_MM_aug_new1
}

