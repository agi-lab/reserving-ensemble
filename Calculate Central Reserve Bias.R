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

ind_subset1_par16_new<-as.numeric(as.character(out_sample$origin))>=2&as.numeric(as.character(out_sample$origin))<=13
ind_subset2_par16_new<-as.numeric(as.character(out_sample$origin))>=14&as.numeric(as.character(out_sample$origin))<=40


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

load("out_mu_list")
load('true_reserve')
out_mu_list<-out_mu_list_par1
# SLP
ntri<-70
mu_OW_par0<-matrix(NA,nrow=780,ncol=ntri)


for (i in 1:ntri){
  mu_OW_par0[,i]<-as.matrix(out_mu_list[[i]][,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par0[[i]]))
}


res_OW_par0<-apply(mu_OW_par0,MARGIN=2,FUN=sum)
res_bias_par0<-abs((res_OW_par0-true_reserve[1:ntri])/true_reserve[1:ntri])

#Equally Weighted

mu_OW_EW<-matrix(NA,nrow=780,ncol=ntri)
for (i in 1:ntri){
  mu_OW_EW[,i]<-as.matrix(out_mu_list[[i]][,-c(19,20)])%*%c(rep(1/18,18))
}


res_OW_EW<-apply(mu_OW_EW,MARGIN=2,FUN=sum)
res_bias_EW<-abs((res_OW_EW-true_reserve[1:ntri])/true_reserve[1:ntri])

#BMV
mu_OW_BMV<-matrix(NA,nrow=780,ncol=ntri)
for (i in 1:ntri){
  mu_OW_BMV[,i]<-out_mu_list[[i]][,14]
}

res_OW_BMV<-apply(mu_OW_BMV,MARGIN=2,FUN=sum)
res_bias_BMV<-abs((res_OW_BMV-true_reserve[1:ntri])/true_reserve[1:ntri])
mean(res_bias_BMV[28:37])
#Partition Strategy 1
mu_OW_out1_par1_new<-matrix(NA,nrow=sum(ind_subset1_par1_new),ncol=ntri)
mu_OW_out2_par1_new<-matrix(NA,nrow=sum(ind_subset2_par1_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par1_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par1_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par1_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par1_new[,i]<-as.matrix(out_mu_list_par1[[i]][ind_subset2_par1_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par1_new[[i]][2]))
}
res_OW_par1_new<-apply(rbind(mu_OW_out1_par1_new,mu_OW_out2_par1_new),MARGIN=2,FUN=sum)
res_bias_par1<-abs((res_OW_par1_new-true_reserve[1:ntri])/true_reserve[1:ntri])

#Predicted Mean for Optimized Ensemble: par2_new


mu_OW_out1_par2_new<-matrix(NA,nrow=sum(ind_subset1_par2_new),ncol=ntri)
mu_OW_out2_par2_new<-matrix(NA,nrow=sum(ind_subset2_par2_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par2_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par2_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par2_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par2_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par2_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par2_new[[i]][2]))
}
res_OW_par2_new<-apply(rbind(mu_OW_out1_par2_new,mu_OW_out2_par2_new),MARGIN=2,FUN=sum)
res_bias_par2<-abs((res_OW_par2_new-true_reserve[1:ntri])/true_reserve[1:ntri])


#Predicted Mean for Optimized Ensemble: par3_new


mu_OW_out1_par3_new<-matrix(NA,nrow=sum(ind_subset1_par3_new),ncol=ntri)
mu_OW_out2_par3_new<-matrix(NA,nrow=sum(ind_subset2_par3_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par3_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par3_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par3_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par3_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par3_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par3_new[[i]][2]))
}
res_OW_par3_new<-apply(rbind(mu_OW_out1_par3_new,mu_OW_out2_par3_new),MARGIN=2,FUN=sum)
res_bias_par3<-abs((res_OW_par3_new-true_reserve[1:ntri])/true_reserve[1:ntri])

#Predicted Mean for Optimized Ensemble: par4_new

mu_OW_out1_par4_new<-matrix(NA,nrow=sum(ind_subset1_par4_new),ncol=ntri)
mu_OW_out2_par4_new<-matrix(NA,nrow=sum(ind_subset2_par4_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par4_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par4_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par4_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par4_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par4_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par4_new[[i]][2]))
}
res_OW_par4_new<-apply(rbind(mu_OW_out1_par4_new,mu_OW_out2_par4_new),MARGIN=2,FUN=sum)
res_bias_par4<-abs((res_OW_par4_new-true_reserve[1:ntri])/true_reserve[1:ntri])


#Predicted Mean for Optimized Ensemble: par5_new


mu_OW_out1_par5_new<-matrix(NA,nrow=sum(ind_subset1_par5_new),ncol=ntri)
mu_OW_out2_par5_new<-matrix(NA,nrow=sum(ind_subset2_par5_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par5_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par5_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par5_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par5_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par5_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par5_new[[i]][2]))
}
res_OW_par5_new<-apply(rbind(mu_OW_out1_par5_new,mu_OW_out2_par5_new),MARGIN=2,FUN=sum)
res_bias_par5<-abs((res_OW_par5_new-true_reserve[1:ntri])/true_reserve[1:ntri])


#Predicted Mean for Optimized Ensemble: par6_new


mu_OW_out1_par6_new<-matrix(NA,nrow=sum(ind_subset1_par6_new),ncol=ntri)
mu_OW_out2_par6_new<-matrix(NA,nrow=sum(ind_subset2_par6_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par6_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par6_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par6_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par6_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par6_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par6_new[[i]][2]))
}
res_OW_par6_new<-apply(rbind(mu_OW_out1_par6_new,mu_OW_out2_par6_new),MARGIN=2,FUN=sum)
res_bias_par6<-abs((res_OW_par6_new-true_reserve[1:ntri])/true_reserve[1:ntri])

#Predicted Mean for Optimized Ensemble: par7_new


mu_OW_out1_par7_new<-matrix(NA,nrow=sum(ind_subset1_par7_new),ncol=ntri)
mu_OW_out2_par7_new<-matrix(NA,nrow=sum(ind_subset2_par7_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par7_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par7_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par7_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par7_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par7_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par7_new[[i]][2]))
}
res_OW_par7_new<-apply(rbind(mu_OW_out1_par7_new,mu_OW_out2_par7_new),MARGIN=2,FUN=sum)
res_bias_par7<-abs((res_OW_par7_new-true_reserve[1:ntri])/true_reserve[1:ntri])

#Predicted Mean for Optimized Ensemble: par8_new


mu_OW_out1_par8_new<-matrix(NA,nrow=sum(ind_subset1_par8_new),ncol=ntri)
mu_OW_out2_par8_new<-matrix(NA,nrow=sum(ind_subset2_par8_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par8_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par8_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par8_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par8_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par8_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par8_new[[i]][2]))
}
res_OW_par8_new<-apply(rbind(mu_OW_out1_par8_new,mu_OW_out2_par8_new),MARGIN=2,FUN=sum)
res_bias_par8<-abs((res_OW_par8_new-true_reserve[1:ntri])/true_reserve[1:ntri])

#Predicted Mean for Optimized Ensemble: par9_new


mu_OW_out1_par9_new<-matrix(NA,nrow=sum(ind_subset1_par9_new),ncol=ntri)
mu_OW_out2_par9_new<-matrix(NA,nrow=sum(ind_subset2_par9_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par9_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par9_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par9_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par9_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par9_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par9_new[[i]][2]))
}
res_OW_par9_new<-apply(rbind(mu_OW_out1_par9_new,mu_OW_out2_par9_new),MARGIN=2,FUN=sum)
res_bias_par9<-abs((res_OW_par9_new-true_reserve[1:ntri])/true_reserve[1:ntri])

#Predicted Mean for Optimized Ensemble: par10_new


mu_OW_out1_par10_new<-matrix(NA,nrow=sum(ind_subset1_par10_new),ncol=ntri)
mu_OW_out2_par10_new<-matrix(NA,nrow=sum(ind_subset2_par10_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par10_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par10_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par10_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par10_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par10_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par10_new[[i]][2]))
}
res_OW_par10_new<-apply(rbind(mu_OW_out1_par10_new,mu_OW_out2_par10_new),MARGIN=2,FUN=sum)
res_bias_par10<-abs((res_OW_par10_new-true_reserve[1:ntri])/true_reserve[1:ntri])

#Predicted Mean for Optimized Ensemble: par11_new


mu_OW_out1_par11_new<-matrix(NA,nrow=sum(ind_subset1_par11_new),ncol=ntri)
mu_OW_out2_par11_new<-matrix(NA,nrow=sum(ind_subset2_par11_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par11_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par11_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par11_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par11_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par11_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par11_new[[i]][2]))
}
res_OW_par11_new<-apply(rbind(mu_OW_out1_par11_new,mu_OW_out2_par11_new),MARGIN=2,FUN=sum)
res_bias_par11<-abs((res_OW_par11_new-true_reserve[1:ntri])/true_reserve[1:ntri])



#Predicted Mean for Optimized Ensemble: par12_new


mu_OW_out1_par12_new<-matrix(NA,nrow=sum(ind_subset1_par12_new),ncol=ntri)
mu_OW_out2_par12_new<-matrix(NA,nrow=sum(ind_subset2_par12_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par12_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par12_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par12_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par12_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par12_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par12_new[[i]][2]))
}
res_OW_par12_new<-apply(rbind(mu_OW_out1_par12_new,mu_OW_out2_par12_new),MARGIN=2,FUN=sum)
res_bias_par12<-abs((res_OW_par12_new-true_reserve[1:ntri])/true_reserve[1:ntri])



#Predicted Mean for Optimized Ensemble: par13_new


mu_OW_out1_par13_new<-matrix(NA,nrow=sum(ind_subset1_par13_new),ncol=ntri)
mu_OW_out2_par13_new<-matrix(NA,nrow=sum(ind_subset2_par13_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par13_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par13_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par13_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par13_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par13_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par13_new[[i]][2]))
}
res_OW_par13_new<-apply(rbind(mu_OW_out1_par13_new,mu_OW_out2_par13_new),MARGIN=2,FUN=sum)
res_bias_par13<-abs((res_OW_par13_new-true_reserve[1:ntri])/true_reserve[1:ntri])


#Predicted Mean for Optimized Ensemble: par14_new


mu_OW_out1_par14_new<-matrix(NA,nrow=sum(ind_subset1_par14_new),ncol=ntri)
mu_OW_out2_par14_new<-matrix(NA,nrow=sum(ind_subset2_par14_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par14_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par14_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par14_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par14_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par14_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par14_new[[i]][2]))
}
res_OW_par14_new<-apply(rbind(mu_OW_out1_par14_new,mu_OW_out2_par14_new),MARGIN=2,FUN=sum)
res_bias_par14<-abs((res_OW_par14_new-true_reserve[1:ntri])/true_reserve[1:ntri])



#Predicted Mean for Optimized Ensemble: par15_new


mu_OW_out1_par15_new<-matrix(NA,nrow=sum(ind_subset1_par15_new),ncol=ntri)
mu_OW_out2_par15_new<-matrix(NA,nrow=sum(ind_subset2_par15_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par15_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par15_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par15_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par15_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par15_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par15_new[[i]][2]))
}
res_OW_par15_new<-apply(rbind(mu_OW_out1_par15_new,mu_OW_out2_par15_new),MARGIN=2,FUN=sum)
res_bias_par15<-abs((res_OW_par15_new-true_reserve[1:ntri])/true_reserve[1:ntri])




#Predicted Mean for Optimized Ensemble: par16_new


mu_OW_out1_par16_new<-matrix(NA,nrow=sum(ind_subset1_par16_new),ncol=ntri)
mu_OW_out2_par16_new<-matrix(NA,nrow=sum(ind_subset2_par16_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par16_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par16_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par16_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par16_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par16_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par16_new[[i]][2]))
}
res_OW_par16_new<-apply(rbind(mu_OW_out1_par16_new,mu_OW_out2_par16_new),MARGIN=2,FUN=sum)
res_bias_par16<-abs((res_OW_par16_new-true_reserve[1:ntri])/true_reserve[1:ntri])


#Predicted Mean for Optimized Ensemble: par17_new


mu_OW_out1_par17_new<-matrix(NA,nrow=sum(ind_subset1_par17_new),ncol=ntri)
mu_OW_out2_par17_new<-matrix(NA,nrow=sum(ind_subset2_par17_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par17_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par17_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par17_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par17_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par17_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par17_new[[i]][2]))
}
res_OW_par17_new<-apply(rbind(mu_OW_out1_par17_new,mu_OW_out2_par17_new),MARGIN=2,FUN=sum)
res_bias_par17<-abs((res_OW_par17_new-true_reserve[1:ntri])/true_reserve[1:ntri])

#par18
mu_OW_out1_par18_new<-matrix(NA,nrow=sum(ind_subset1_par18_new),ncol=ntri)
mu_OW_out2_par18_new<-matrix(NA,nrow=sum(ind_subset2_par18_new),ncol=ntri)

for (i in 1:ntri){
  mu_OW_out1_par18_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset1_par18_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par18_new[[i]][1]))
}

for (i in 1:ntri){
  mu_OW_out2_par18_new[,i]<-as.matrix(out_mu_list[[i]][ind_subset2_par18_new,-c(19,20)])%*%as.vector(unlist(model_weights_simul_par18_new[[i]][2]))
}
res_OW_par18_new<-apply(rbind(mu_OW_out1_par18_new,mu_OW_out2_par18_new),MARGIN=2,FUN=sum)
res_bias_par18<-abs((res_OW_par18_new-true_reserve[1:ntri])/true_reserve[1:ntri])


#Create a vector that contains all reserve bias

rb_central_comparision_ADLP<-cbind(res_bias_par0,res_bias_par1,res_bias_par2,res_bias_par3,res_bias_par4,res_bias_par5,res_bias_par6,res_bias_par7,res_bias_par8,res_bias_par9,res_bias_par10,res_bias_par11,res_bias_par12,res_bias_par13,res_bias_par14,res_bias_par15,res_bias_par16,res_bias_par18)

save(rb_central_comparision_ADLP,file="Central Relative Bias")