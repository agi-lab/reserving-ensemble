
####### Log Score plottings ###########################


### Calculate the distribution of Log Score:
LS_dist_SLP<-apply(OW_out_dens_par0,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP1<-apply(OW_out_dens_par1_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP2<-apply(OW_out_dens_par2_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP3<-apply(OW_out_dens_par3_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP4<-apply(OW_out_dens_par4_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP5<-apply(OW_out_dens_par5_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP6<-apply(OW_out_dens_par6_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP7<-apply(OW_out_dens_par7_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP8<-apply(OW_out_dens_par8_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP9<-apply(OW_out_dens_par9_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP10<-apply(OW_out_dens_par10_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP11<-apply(OW_out_dens_par11_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP12<-apply(OW_out_dens_par12_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP13<-apply(OW_out_dens_par13_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP14<-apply(OW_out_dens_par14_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP15<-apply(OW_out_dens_par15_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP16<-apply(OW_out_dens_par16_new,MARGIN=2,FUN=function(x) mean(log(x)))
LS_dist_ADLP18<-apply(OW_out_dens_par18_new,MARGIN=2,FUN=function(x) mean(log(x)))


### Plot the mean Log Score with different split points
pdf("Mean_Log_Score_withDifferentSplitPoints.pdf")
corres_splitpoints_SLPincluded<-c(0,28,31,33,26,23,19,13,15,11,9,7,5,3,16,17,18,14)
par(mar = c(3.8, 3.8, 3.8, 2))
LS_comparision_ADLP <-cbind(LS_dist_SLP,LS_dist_ADLP1,LS_dist_ADLP2,LS_dist_ADLP3,LS_dist_ADLP4,LS_dist_ADLP5,LS_dist_ADLP6,LS_dist_ADLP7,LS_dist_ADLP8,LS_dist_ADLP9,LS_dist_ADLP10,LS_dist_ADLP11,LS_dist_ADLP12,LS_dist_ADLP13,LS_dist_ADLP14,LS_dist_ADLP15,LS_dist_ADLP16,LS_dist_ADLP18)
plot(x=corres_splitpoints_SLPincluded,y= apply(LS_comparision_ADLP,FUN=mean,MARGIN=2),type="p",xlab="Accident Periods",ylab="Mean Log Score",main="Mean Log Score v.s. Different split points")
points(x=corres_splitpoints_SLPincluded[which.max(apply(LS_comparision_ADLP,FUN=mean,MARGIN=2))],y=max(apply(LS_comparision_ADLP,FUN=mean,MARGIN=2)),type="p",col="black",pch=16,cex=1.6)
points(x=0,y=apply(LS_comparision_ADLP,FUN=mean,MARGIN=2)[1],type="p",col="grey",pch=16,cex=1.6)
legend(x=15,y=-3.88,legend=c("SLP","Optimal ADLP"),col=c("grey","black"),pch=c(16,16),cex=1.5)
dev.off()

### Plot the Log Score by accident periods:
pdf("LogScore_by_AccidentPeriods.pdf")
par(mar = c(3.8, 3.8, 3.8, 0))
plot(x=2:40,y=avg_Log_Score_acc_par1[,1],col="grey",ylim=c(-8,0),type="l",main="Log-Score by Accident Periods ",ylab="Mean Log Score",xlab="Accident Periods")
for (i in 2:18){
  points(x=2:40,y=avg_Log_Score_acc_par1[,i],col="grey",type="l")
}
points(x=2:40,y=LS_acc_OW_par0_avg,col="yellow",type="l",lwd=2)
points(x=2:40,y=LS_acc_OW_par3_new_avg,col="#2417DA",type="l",lwd=1)
points(x=2:40,y=LS_acc_OW_par8_new_avg,col="#6D4691",type="l",lwd=1)
points(x=2:40,y=LS_acc_OW_par13_new_avg,col="#FFA500",type="l",lwd=1)
points(x=2:40,y=avg_Log_Score_acc_par1[,20],col="red",type="l",lwd=1)
points(x=2:40,y=LS_bestvalid_out_acc_mean_par1,col="green",type="l",lwd=1)
legend(x=20,y=0,legend=c("Component Models","Best Models in Validation Set","SLP","ADLP 1(AP 2-3)","ADLP 8(AP2-15)","ADLP 17(AP2-33)","Equally Weighted Ensemble"),col=c("grey","green","yellow","#FFA500","#6D4691","#2417DA","red"),lty=1,cex=0.6)
dev.off()



####### CRPS plottings ###########################

### Calculate the distribution of CRPS: 

MeanCRPS_BMV<-apply(SpLN_crps,MARGIN=2,FUN=mean)
MeanCRPS_EW<-apply(EqEns_crps,MARGIN=2,FUN=mean)
MeanCRPS_SLP<-apply(SLP_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP1<-apply(ADLP1_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP2<-apply(ADLP2_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP3<-apply(ADLP3_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP4<-apply(ADLP4_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP5<-apply(ADLP5_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP6<-apply(ADLP6_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP7<-apply(ADLP7_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP8<-apply(ADLP8_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP9<-apply(ADLP9_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP10<-apply(ADLP10_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP11<-apply(ADLP11_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP12<-apply(ADLP12_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP13<-apply(ADLP13_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP14<-apply(ADLP14_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP15<-apply(ADLP15_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP16<-apply(ADLP16_crps,MARGIN=2,FUN=mean)
MeanCRPS_ADLP17<-apply(ADLP17_crps,MARGIN=2,FUN=mean)

### Calculate the CRPS by accident periods:

origin<-as.numeric(as.character(out_dens_list[[1]]$origin))
dev<-as.numeric(as.character(out_dens_list[[1]]$dev))

OW_out_CRPS_ADLP5_new_dat<-as.data.frame(cbind(origin,dev,ADLP5_crps))
CRPS_acc_OW_ADLP5_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  CRPS_acc_OW_ADLP5_new[i-1,]<-apply(OW_out_CRPS_ADLP5_new_dat[OW_out_CRPS_ADLP5_new_dat$origin==i,][,-c(1,2)],MARGIN=2,FUN=mean)
}
CRPS_acc_OW_ADLP5_new_avg<-apply(CRPS_acc_OW_ADLP5_new,MARGIN=1,FUN=mean)


OW_out_CRPS_ADLP8_new_dat<-as.data.frame(cbind(origin,dev,ADLP8_crps))
CRPS_acc_OW_ADLP8_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  CRPS_acc_OW_ADLP8_new[i-1,]<-apply(OW_out_CRPS_ADLP8_new_dat[OW_out_CRPS_ADLP8_new_dat$origin==i,][,-c(1,2)],MARGIN=2,FUN=mean)
}
CRPS_acc_OW_ADLP8_new_avg<-apply(CRPS_acc_OW_ADLP8_new,MARGIN=1,FUN=mean)

OW_out_CRPS_SpLN_new_dat<-as.data.frame(cbind(origin,dev,SpLN_crps))
CRPS_acc_OW_SpLN_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  CRPS_acc_OW_SpLN_new[i-1,]<-apply(OW_out_CRPS_SpLN_new_dat[OW_out_CRPS_SpLN_new_dat$origin==i,][,-c(1,2)],MARGIN=2,FUN=mean)
}
CRPS_acc_OW_SpLN_new_avg<-apply(CRPS_acc_OW_SpLN_new,MARGIN=1,FUN=mean)


OW_out_CRPS_SLP_new_dat<-as.data.frame(cbind(origin,dev,SLP_crps))
CRPS_acc_OW_SLP_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  CRPS_acc_OW_SLP_new[i-1,]<-apply(OW_out_CRPS_SLP_new_dat[OW_out_CRPS_SLP_new_dat$origin==i,][,-c(1,2)],MARGIN=2,FUN=mean)
}
CRPS_acc_OW_SLP_new_avg<-apply(CRPS_acc_OW_SLP_new,MARGIN=1,FUN=mean)

OW_out_CRPS_EqEns_new_dat<-as.data.frame(cbind(origin,dev,EqEns_crps))
CRPS_acc_OW_EqEns_new<-matrix(NA,nrow=39,ncol=ntri)
for (i in 2:40){
  CRPS_acc_OW_EqEns_new[i-1,]<-apply(OW_out_CRPS_EqEns_new_dat[OW_out_CRPS_EqEns_new_dat$origin==i,][,-c(1,2)],MARGIN=2,FUN=mean)
}
CRPS_acc_OW_EqEns_new_avg<-apply(CRPS_acc_OW_EqEns_new,MARGIN=1,FUN=mean)


### Plot the mean CRPS with different split points

pdf("Mean_CRPS_withDifferentSplitPoints.pdf")
CRPS_comparision_ADLP<-cbind(MeanCRPS_SLP,MeanCRPS_ADLP1,MeanCRPS_ADLP2,MeanCRPS_ADLP3,MeanCRPS_ADLP4,MeanCRPS_ADLP5,MeanCRPS_ADLP6,MeanCRPS_ADLP7,MeanCRPS_ADLP8,MeanCRPS_ADLP9,MeanCRPS_ADLP10,MeanCRPS_ADLP11,MeanCRPS_ADLP12,MeanCRPS_ADLP13,MeanCRPS_ADLP14,MeanCRPS_ADLP15,MeanCRPS_ADLP16,MeanCRPS_ADLP17)
corres_splitpoints_SLPincluded<-c(0,28,31,33,26,23,19,13,15,11,9,7,5,3,16,17,18,4)
par(mar = c(3.8, 3.8, 3.8, 2))
plot(x=corres_splitpoints_SLPincluded,y=apply(CRPS_comparision_ADLP,FUN=mean,MARGIN=2),type="p",xlab="Accident Periods",ylab="Mean CRPS",main="Mean CRPS v.s. Different split points")
points(x=corres_splitpoints_SLPincluded[which.min(apply(CRPS_comparision_ADLP,FUN=mean,MARGIN=2))],y=min(apply(CRPS_comparision_ADLP,FUN=mean,MARGIN=2)),type="p",pch=16,cex=1.6,col="black")
points(x=0,y=apply(CRPS_comparision_ADLP,FUN=mean,MARGIN=2)[1],type="p",pch=16,cex=1.6,col="grey")
points(x=15,y=apply(CRPS_comparision_ADLP,FUN=mean,MARGIN=2)[9],type="p",pch=16,cex=1.6,col="brown")
legend("topright",legend=c("SLP","ADLP8","ADLP13"),col=c("grey","brown","black"),pch=c(16,16,16), cex = 0.8)
dev.off()

### Plot the mean CRPS by accident periods: 

pdf("Mean_CRPS_byAccidentPeriods.pdf")
par(mar = c(3.8, 3.8, 3.8, 0))
plot(x=2:40,y=CRPS_acc_OW_ADLP5_new_avg,col="blue",type="l",main="CRPS by Accident Periods ",ylab="Mean CRPS",xlab="Accident Periods",ylim=c(2,25))
points(x=2:40,y=CRPS_acc_OW_ADLP8_new_avg,col="#6D4691",type="l")
points(x=2:40,y=CRPS_acc_OW_SpLN_new_avg,col="green",type="l")
points(x=2:40,y=CRPS_acc_OW_SLP_new_avg,col="orange",type="l")
legend("bottomright",legend=c("ADLP13","ADLP8","BMV","SLP"),col=c("blue","#6D4691","green","orange"),lty=1, cex = 0.7)
dev.off()


#corres_splitpoints <- c(0,28,31,33,26,23,19,15,13,11,9,7,5,3,16,17,18,4)
#plot(x = corres_splitpoints, y = reserve_bias75_comparison, xlab = "Accident periods", ylab = "Relative Reserve bias at 75th quantile")
