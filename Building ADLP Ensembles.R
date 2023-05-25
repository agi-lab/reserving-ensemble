###########################
##Note: the numbering of the ADLP ensemble is not the same as the numbering of ADLP ensemble in the paper due to better presentation in paper
##########################
####please refer to the following index section for the actual partition of subsets

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


ind_subset1_par9_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=11
ind_subset2_par9_new<-as.numeric(as.character(out_sample$origin))>=12&as.numeric(as.character(out_sample$origin))<=40


ind_subset1_par10_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=9
ind_subset2_par10_new<-as.numeric(as.character(out_sample$origin))>=10&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par11_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=7
ind_subset2_par11_new<-as.numeric(as.character(out_sample$origin))>=8&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par12_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=5
ind_subset2_par12_new<-as.numeric(as.character(out_sample$origin))>=6&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par13_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=3
ind_subset2_par13_new<-as.numeric(as.character(out_sample$origin))>=4&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par14_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=16
ind_subset2_par14_new<-as.numeric(as.character(out_sample$origin))>=17&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par15_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=17
ind_subset2_par15_new<-as.numeric(as.character(out_sample$origin))>=18&as.numeric(as.character(out_sample$origin))<=40

ind_subset1_par16_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=18
ind_subset2_par16_new<-as.numeric(as.character(out_sample$origin))>=19&as.numeric(as.character(out_sample$origin))<=40


ind_subset1_par17_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=4
ind_subset2_par17_new<-as.numeric(as.character(out_sample$origin))>=5&as.numeric(as.character(out_sample$origin))<=40


ind_subset1_par18_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=14
ind_subset2_par18_new<-as.numeric(as.character(out_sample$origin))>=15&as.numeric(as.character(out_sample$origin))<=40




#Partition Strategy 1

OW_out_dens_par1_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
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

OW_out_dens_par2_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
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


OW_out_dens_par3_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
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

OW_out_dens_par4_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
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

OW_out_dens_par5_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
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

OW_out_dens_par6_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
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

OW_out_dens_par7_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par7_new<-list()
w_init_1_3 <- rep(1/18, 18)
# # ntri <- 100
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




#Partition Strategy 8 (the optimal strategy identified) ------------------------------------------------------------------------------------------------------------

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par8_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
OW_out_CDF_par8_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)

model_weights_simul_par8_new<-list()
for (D in 1:ntri){
  #D <- 1
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=13,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=13,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=14&as.numeric(as.character(meta_dens_out$origin))<=40,]
  
  out_meta_1_aug_new1_CDF<-meta_CDF_out[as.numeric(as.character(meta_CDF_out$origin))>=2&as.numeric(as.character(meta_CDF_out$origin))<=13,]
  out_meta_2_aug_new1_CDF<-meta_CDF_out[as.numeric(as.character(meta_CDF_out$origin))>=14&as.numeric(as.character(meta_CDF_out$origin))<=40,]
  
  
  #Train the model weights
  meta_valid_list_aug_new1<-list(augvalid_meta_1_new1,augvalid_meta_2_new1)
  meta_test_list_aug_new1<-list(out_meta_1_aug_new1,out_meta_2_aug_new1)
  meta_test_list_aug_new1_CDF<-list(out_meta_1_aug_new1_CDF,out_meta_2_aug_new1_CDF)
  
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
  Out_sample_EnCDF_aug_new1<-list()
  for (i in 1:2){
    out_meta<-as.matrix(meta_test_list_aug_new1[[i]][,3:ncol(meta_test_list_aug_new1[[i]])])
    Out_sample_EnDens_aug_new1[[i]]<-out_meta%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
    
    out_meta_CDF<-as.matrix(meta_test_list_aug_new1_CDF[[i]][,3:ncol(meta_test_list_aug_new1_CDF[[i]])])
    Out_sample_EnCDF_aug_new1[[i]]<-out_meta_CDF%*%as.vector(unlist(finalw_MM_aug_new1[[i]]))
    
    
  }
  OW_out_dens_par8_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  OW_out_CDF_par8_new[,D]<-as.vector(unlist(Out_sample_EnCDF_aug_new1))
  model_weights_simul_par8_new[[D]]<-finalw_MM_aug_new1
  
}


#Partition Strategy 9

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par9_new <- matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par9_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=11,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=11,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=12&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par9_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par9_new[[D]]<-finalw_MM_aug_new1
  
}


#Partition Strategy 10

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par10_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par10_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=9,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=9,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=10&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par10_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par10_new[[D]]<-finalw_MM_aug_new1
  
}



#Partition Strategy 11

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par11_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par11_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=7,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=7,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=8&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par11_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par11_new[[D]]<-finalw_MM_aug_new1
  
}


#Partition Strategy 12

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par12_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par12_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=5,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=5,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=6&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par12_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par12_new[[D]]<-finalw_MM_aug_new1
  
}




#Partition Strategy 13

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par13_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par13_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=3,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=3,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=4&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par13_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par13_new[[D]]<-finalw_MM_aug_new1
  
}




#Partition Strategy 14

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par14_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par14_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=16,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=16,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=17&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par14_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par14_new[[D]]<-finalw_MM_aug_new1
  
}

#Partition Strategy 15

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par15_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par15_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=17,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=17,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=18&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par15_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par15_new[[D]]<-finalw_MM_aug_new1
  
}


#Partition Strategy 16

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par16_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par16_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=18,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=18,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=19&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par16_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par16_new[[D]]<-finalw_MM_aug_new1
  
}


#Partition Strategy 17

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par17_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par17_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=4,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=4,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=5&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par17_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par17_new[[D]]<-finalw_MM_aug_new1
  
}


#Partition Strategy 18

##Redo the Caculation of the Density of the Optimized Ensemble 

OW_out_dens_par18_new<-matrix(NA,nrow=nrow(out_sample),ncol=ntri)
model_weights_simul_par18_new<-list()
for (D in 1:ntri){
  meta_dens_augvalid<-valid_dens_list[[D]][,-c(21,22)]
  meta_dens_out<-out_dens_list[[D]][,-c(21,22)]
  
  augvalid_meta_1_new1<-meta_dens_augvalid[as.numeric(as.character(meta_dens_augvalid$origin))>=2&as.numeric(as.character(meta_dens_augvalid$origin))<=14,]
  augvalid_meta_2_new1<-meta_dens_augvalid
  
  
  #Construct the Corresponding Out-of-Sample Test sets
  out_meta_1_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=2&as.numeric(as.character(meta_dens_out$origin))<=14,]
  out_meta_2_aug_new1<-meta_dens_out[as.numeric(as.character(meta_dens_out$origin))>=15&as.numeric(as.character(meta_dens_out$origin))<=40,]
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
  OW_out_dens_par18_new[,D]<-as.vector(unlist(Out_sample_EnDens_aug_new1))
  model_weights_simul_par18_new[[D]]<-finalw_MM_aug_new1
  
}


#Calculate (out-of-sample) Log Score by Accident periods

origin<-as.numeric(as.character(out_dens_list[[1]]$origin))
dev<-as.numeric(as.character(out_dens_list[[1]]$dev))


OW_out_LS_par0_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par0))
LS_acc_OW_par0<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par0[i-1,]<-apply(log(OW_out_LS_par0_dat[OW_out_LS_par0_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par0_avg<-apply(LS_acc_OW_par0,MARGIN=1,FUN=mean)

OW_out_LS_par1_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par1_new))
LS_acc_OW_par1_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par1_new[i-1,]<-apply(log(OW_out_LS_par1_new_dat[OW_out_LS_par1_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par1_new_avg<-apply(LS_acc_OW_par1_new,MARGIN=1,FUN=mean)

OW_out_LS_par2_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par2_new))
LS_acc_OW_par2_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par2_new[i-1,]<-apply(log(OW_out_LS_par2_new_dat[OW_out_LS_par2_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par2_new_avg<-apply(LS_acc_OW_par2_new,MARGIN=1,FUN=mean)

OW_out_LS_par3_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par3_new))
LS_acc_OW_par3_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par3_new[i-1,]<-apply(log(OW_out_LS_par3_new_dat[OW_out_LS_par3_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par3_new_avg<-apply(LS_acc_OW_par3_new,MARGIN=1,FUN=mean)

OW_out_LS_par4_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par4_new))
LS_acc_OW_par4_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par4_new[i-1,]<-apply(log(OW_out_LS_par4_new_dat[OW_out_LS_par4_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par4_new_avg<-apply(LS_acc_OW_par4_new,MARGIN=1,FUN=mean)



OW_out_LS_par5_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par5_new))
LS_acc_OW_par5_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par5_new[i-1,]<-apply(log(OW_out_LS_par5_new_dat[OW_out_LS_par5_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par5_new_avg<-apply(LS_acc_OW_par5_new,MARGIN=1,FUN=mean)


OW_out_LS_par6_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par6_new))
LS_acc_OW_par6_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par6_new[i-1,]<-apply(log(OW_out_LS_par6_new_dat[OW_out_LS_par6_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par6_new_avg<-apply(LS_acc_OW_par6_new,MARGIN=1,FUN=mean)


OW_out_LS_par7_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par7_new))
LS_acc_OW_par7_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par7_new[i-1,]<-apply(log(OW_out_LS_par7_new_dat[OW_out_LS_par7_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par7_new_avg<-apply(LS_acc_OW_par7_new,MARGIN=1,FUN=mean)


OW_out_LS_par8_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par8_new))
LS_acc_OW_par8_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par8_new[i-1,]<-apply(log(OW_out_LS_par8_new_dat[OW_out_LS_par8_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par8_new_avg<-apply(LS_acc_OW_par8_new,MARGIN=1,FUN=mean)

OW_out_LS_par9_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par9_new))

LS_acc_OW_par9_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par9_new[i-1,]<-apply(log(OW_out_LS_par9_new_dat[OW_out_LS_par9_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par9_new_avg<-apply(LS_acc_OW_par9_new,MARGIN=1,FUN=mean)


OW_out_LS_par10_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par10_new))
LS_acc_OW_par10_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par10_new[i-1,]<-apply(log(OW_out_LS_par10_new_dat[OW_out_LS_par10_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par10_new_avg<-apply(LS_acc_OW_par10_new,MARGIN=1,FUN=mean)



OW_out_LS_par11_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par11_new))
LS_acc_OW_par11_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par11_new[i-1,]<-apply(log(OW_out_LS_par11_new_dat[OW_out_LS_par11_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par11_new_avg<-apply(LS_acc_OW_par11_new,MARGIN=1,FUN=mean)


OW_out_LS_par12_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par12_new))
LS_acc_OW_par12_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par12_new[i-1,]<-apply(log(OW_out_LS_par12_new_dat[OW_out_LS_par12_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par12_new_avg<-apply(LS_acc_OW_par12_new,MARGIN=1,FUN=mean)

OW_out_LS_par13_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par13_new))
LS_acc_OW_par13_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par13_new[i-1,]<-apply(log(OW_out_LS_par13_new_dat[OW_out_LS_par13_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par13_new_avg<-apply(LS_acc_OW_par13_new,MARGIN=1,FUN=mean)


OW_out_LS_par14_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par14_new))
LS_acc_OW_par14_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par14_new[i-1,]<-apply(log(OW_out_LS_par14_new_dat[OW_out_LS_par14_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par14_new_avg<-apply(LS_acc_OW_par14_new,MARGIN=1,FUN=mean)


OW_out_LS_par15_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par15_new))
LS_acc_OW_par15_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par15_new[i-1,]<-apply(log(OW_out_LS_par15_new_dat[OW_out_LS_par15_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par15_new_avg<-apply(LS_acc_OW_par15_new,MARGIN=1,FUN=mean)

OW_out_LS_par16_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par16_new))
LS_acc_OW_par16_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par16_new[i-1,]<-apply(log(OW_out_LS_par16_new_dat[OW_out_LS_par16_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par16_new_avg<-apply(LS_acc_OW_par16_new,MARGIN=1,FUN=mean)


OW_out_LS_par17_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par17_new))
LS_acc_OW_par17_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par17_new[i-1,]<-apply(log(OW_out_LS_par17_new_dat[OW_out_LS_par17_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par17_new_avg<-apply(LS_acc_OW_par17_new,MARGIN=1,FUN=mean)



OW_out_LS_par18_new_dat<-as.data.frame(cbind(origin,dev,OW_out_dens_par18_new))
LS_acc_OW_par18_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  LS_acc_OW_par18_new[i-1,]<-apply(log(OW_out_LS_par18_new_dat[OW_out_LS_par18_new_dat$origin==i,][,-c(1,2)]),MARGIN=2,FUN=mean)
}
LS_acc_OW_par18_new_avg<-apply(LS_acc_OW_par18_new,MARGIN=1,FUN=mean)


