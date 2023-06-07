rm(list=ls())
#Define the reference DoE of size 18(54)
xrf1<-c(0.015,0.0375,0.06,0.0825,0.105,0.1275,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.25,1.75,2.25,3.5) 
xrf2<-rep(xrf1,each=3)

#Define all the LOD and BOD
xd10<-rep(c(0.0675,0.4   ,1.6025,3.5),c(4,5,4,5))
xd11<-rep(c(0.065 ,0.3875,1.57  ,3.5),c(5,4,4,5))
xd12<-rep(c(0.0525,0.19  ,0.3825,1.51,3.5),c(4,1,4,4,5))
xa10<-rep(c(0.045 ,0.3925,1.665 ,3.5),c(5,4,5,4))
xa11<-rep(c(0.055 ,0.4125,1.7275,3.5),c(5,5,5,3))
xa12<-rep(c(0.0525,0.24  ,0.325 ,0.565,1.72,3.5),c(5,2,1,3,5,2))
xd20<-rep(c(0.0675,0.4   ,1.6025,3.5),c(13,14,13,14))
xd21<-rep(c(0.065 ,0.3875,1.57  ,3.5),c(14,13,13,14))
xd22<-rep(c(0.055 ,0.3555,1.4925,3.5),c(14,13,13,14))
xa20<-rep(c(0.045 ,0.3975,1.6475,3.5),c(14,13,16,11))
xa21<-rep(c(0.055 ,0.4175,1.7375,3.5),c(15,14,16,9))
xa22<-rep(c(0.05  ,0.25  ,0.52  ,1.73 ,3.5),c(17,7,9,14,7))

#choose if we assume covariance matrix for the parameter prior distribution
Method=2
#draw a SAMPLE of size 300000
R=300000
GammaParams <- function(mu, var) {
  s <- var/mu
  a <- mu/s
  return(Params = list(a = a, s = s))
}

Prior<-c(0.72,1.68,0.13,3.10)
Prior[4]=Prior[4]-Prior[3]
SD1<-0.35
SD2<-0.6
nparams=length(Prior)
Prior_1=matrix(NA,nrow=R,ncol=nparams)
Prior_2=Prior_1

if        (Method==1){ #no covariances
for (i in 1:nparams){
GParams<-GammaParams(Prior[i],(Prior[i]*SD1)^2)
Prior_1[,i]<-rgamma(R,GParams$a,scale=GParams$s)
GParams<-GammaParams(Prior[i],(Prior[i]*SD2)^2)
Prior_2[,i]<-rgamma(R,GParams$a,scale=GParams$s)
}}else if (Method==2){ #covariances and correlations
Cor=matrix(0.8,nparams,nparams)+diag(nparams)*0.2
Cor[2,1]=0.6; Cor[1,2]=0.6;
Cor[2,3]=0.6; Cor[3,2]=0.6;
Cor[4,3]=0.9; Cor[3,4]=0.9;
require(MASS)
TEMP<-pnorm(mvrnorm(R,rep(0,nparams),Cor))
for (i in 1:nparams){
GParams<-GammaParams(Prior[i],(Prior[i]*SD1)^2)
Prior_1[,i]<-qgamma(TEMP[,i],GParams$a,scale=GParams$s)
GParams<-GammaParams(Prior[i],(Prior[i]*SD2)^2)
Prior_2[,i]<-qgamma(TEMP[,i],GParams$a,scale=GParams$s)
}}
Prior_1[,4]=Prior_1[,3]+Prior_1[,4]
Prior_2[,4]=Prior_2[,3]+Prior_2[,4]

#simulate initial rates and define NLS for regression
NLSfit0<-function(Pvec,x,n,warn){
  Y<-Pvec[1]*x/(Pvec[3]+x)+Pvec[2]*x/(Pvec[4]+x)+rnorm(n,sd=(Pvec[1]*1/(Pvec[3]+1)+Pvec[2]*1/(Pvec[4]+1))*0.02)
  Y[Y<0]<-0
  return(coef(nls(Y~v1*x/(k1+x)+v2*x/(k2+x),start=list(v1=0.72, v2=1.68, k1=0.13, k2=3.10),control=warn)))
}

NLSfit1<-function(Pvec,x,n,warn){
  Y<-Pvec[1]*x/(Pvec[3]+x)+Pvec[2]*x/(Pvec[4]+x)+rnorm(n,sd=(Pvec[1]*1/(Pvec[3]+1)+Pvec[2]*1/(Pvec[4]+1))*0.02)
  Y[Y<0]<-0
  return(coef(nls(Y~v1*x/(k1+x)+v2*x/(k2+x),start=list(v1=0.72, v2=1.68, k1=0.13, k2=3.10),control=warn,algorithm="port",lower=0.01)))
}

NLSfit2<-function(Pvec,x,n,warn){
  Y<-Pvec[1]*x/(Pvec[3]+x)+Pvec[2]*x/(Pvec[4]+x)*(1+rnorm(n,sd=0.02))
  Y[Y<0]<-0
  return(coef(nls(Y~v1*x/(k1+x)+v2*x/(k2+x),start=list(v1=0.72, v2=1.68, k1=0.13, k2=3.10),control=warn,algorithm="port",lower=0.01)))
}

NLSfit3<-function(Pvec,x,n,warn){
  Y<-Pvec[1]*x/(Pvec[3]+x)+Pvec[2]*x/(Pvec[4]+x)+rnorm(n,sd=(Pvec[1]*1/(Pvec[3]+1)+Pvec[2]*1/(Pvec[4]+1))*0.02)
  Y[Y<0]<-0
  return(coef(nls(Y~v1*x/(k1+x)+v2*x/(k2+x),start=list(v1=0.72, v2=1.68, k1=0.13, k2=3.10),control=warn,algorithm="port",lower=c(0.5,0.5,0.1,0.1))))
}

#simulation function
simu<-function(index,NLSfit){
      if (index== 1){x=xrf1; Pvec=Prior_1; n=18
}else if (index== 2){x=xd10; Pvec=Prior_1; n=18
}else if (index== 3){x=xd11; Pvec=Prior_1; n=18
}else if (index== 4){x=xa10; Pvec=Prior_1; n=18
}else if (index== 5){x=xa11; Pvec=Prior_1; n=18
}else if (index== 6){x=xrf2; Pvec=Prior_1; n=54
}else if (index== 7){x=xd20; Pvec=Prior_1; n=54
}else if (index== 8){x=xd21; Pvec=Prior_1; n=54
}else if (index== 9){x=xa20; Pvec=Prior_1; n=54
}else if (index==10){x=xa21; Pvec=Prior_1; n=54
}else if (index==11){x=xrf1; Pvec=Prior_2; n=18
}else if (index==12){x=xd10; Pvec=Prior_2; n=18
}else if (index==13){x=xd12; Pvec=Prior_2; n=18
}else if (index==14){x=xa10; Pvec=Prior_2; n=18
}else if (index==15){x=xa12; Pvec=Prior_2; n=18
}else if (index==16){x=xrf2; Pvec=Prior_2; n=54
}else if (index==17){x=xd20; Pvec=Prior_2; n=54
}else if (index==18){x=xd22; Pvec=Prior_2; n=54
}else if (index==19){x=xa20; Pvec=Prior_2; n=54
}else if (index==20){x=xa22; Pvec=Prior_2; n=54
}
Hat=apply(Pvec,1,NLSfit,x=x,n=n,warn=nls.control(minFactor=1/32768,warnO=TRUE))
return(Hat)
}

#use Parallel Computing for parametric bootstrapping: one hour for each
#save all the parameter estimates
require(parallel) 

cl <- makeCluster(detectCores())
clusterExport(cl,c("Prior_1","Prior_2","xrf1","xrf2","xd10","xd11","xd12","xd20","xd21","xd22","xa10","xa11","xa12","xa20","xa21","xa22"))
Hatset0<-parLapply(cl,1:20,simu,NLSfit=NLSfit0)
stopCluster(cl)

cl <- makeCluster(detectCores())
clusterExport(cl,c("Prior_1","Prior_2","xrf1","xrf2","xd10","xd11","xd12","xd20","xd21","xd22","xa10","xa11","xa12","xa20","xa21","xa22"))
Hatset1<-parLapply(cl,1:20,simu,NLSfit=NLSfit1)
stopCluster(cl)

cl <- makeCluster(detectCores())
clusterExport(cl,c("Prior_1","Prior_2","xrf1","xrf2","xd10","xd11","xd12","xd20","xd21","xd22","xa10","xa11","xa12","xa20","xa21","xa22"))
Hatset2<-parLapply(cl,1:20,simu,NLSfit=NLSfit2)
stopCluster(cl)

cl <- makeCluster(detectCores())
clusterExport(cl,c("Prior_1","Prior_2","xrf1","xrf2","xd10","xd11","xd12","xd20","xd21","xd22","xa10","xa11","xa12","xa20","xa21","xa22"))
Hatset3<-parLapply(cl,1:20,simu,NLSfit=NLSfit3)
stopCluster(cl)