# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

source("simu_Kw-GAS-COV.R")
source("GASKwFit-COV.R")
source("simu_Beta-GAS-COV.R")
source("GASBetaFit-COV.R")
library(forecast)
############################
#### initial quantities ####
############################
start_time <- Sys.time()
vn<-c(500)
R<-1000
bug<-0
true<-"Beta"
# true<-"Kw"
scenarios<-c("bbonita","capivara","avermelha")
scenario=scenarios[2]
if(scenario==scenarios[1]){
  scen=1
  # bbonita
  w<-0.5;A=0.1;B=0.7;phi=10.1;beta=c(0.6,-0.2) 
  cov<-"sin&cos"
  X<-sin(2*pi*t/12) # in sample
  X_hat<-sin(2*pi*t_hat/12) # out of sample
  theta=c(w,A,B,phi,beta)
}else{
  if(scenario==scenarios[2]){
    scen=2
    w=1;A=0.1;B=.3;phi=2.5;beta=-0.4
    cov<-"sin"
    theta=c(w,A,B,phi,beta)
  }else{
    if(scenario==scenarios[3])
      scen=3
    # avermelha without covariates
    w<--.1;A=0.1;B=0.7;phi=2.7;beta=0
    cov<-NA
    X=NA
    X_hat=NA
    theta=c(w,A,B,phi)
  }
}
estimbeta<-estim<-array(NA,c(R,length(theta),length(vn)))
pv_LB<-AIC<-RMSER<-MAPE<-array(NA,c(R,2,length(vn)))
mean_MLEs<-RMSE<-RB<-matrix(NA,length(vn),length(theta))
pv_LB05<-pv_LB010<-pv_LB001<-contAIC<-contRMSER<-contMAPE<-matrix(NA,length(vn),2)
##########################
#### simulation study ####
##########################
for(j in 1:length(vn)){
  n<-vn[j]
  t=1:n
  t_hat=n+1
  i<-1
  if(scenario==scenarios[1]){
    # bbonita
    X<-cbind(sin(2*pi*t/12),cos(2*pi*t/12))
    X_hat<-cbind(sin(2*pi*t_hat/12),cos(2*pi*t_hat/12)) 
  }else{
    if(scenario==scenarios[2]){
      # capivara without covariate
      X<-sin(2*pi*t/12) # in sample
      X_hat<-sin(2*pi*t_hat/12) # out of sample
    }else{
      if(scenario==scenarios[3])
        # avermelha
        X=NA
        X_hat=NA
    }
  }
  set.seed(10)
  while(i<=R){
    if(true=="Kw"){
      y <- simu.GASKw(n,A=A,B=B,w=w,phi=phi,beta=beta,X=cov)
    }else{
      y <- simu.GASbeta(n,A=A,B=B,w=w,alpha=phi,beta=beta,X=cov)
    }
    result <- try(GASKw.fit(y,1,1,X=X,X_hat=X_hat,h=1),silent = T)
    resultbeta <- try(GASbeta.fit(y,1,1,X=X,X_hat=X_hat,h=1),silent = T)
    if(class(result)=="try-error" ||  result$conv != 0 ||
       class(resultbeta)=="try-error" ||  resultbeta$conv != 0){
      bug<-bug+1
      print(bug)
    }else{
      estim[i,,j]<-result$model[,1]
      estimbeta[i,,j]<-resultbeta$model[,1]
      AIC[i,,j]<-c(result$aic,resultbeta$aic)
      RMSER[i,,j]<-c(accuracy(result$fitted, y)[2],
                     accuracy(resultbeta$fitted, y)[2])
      MAPE[i,,j]<-c(accuracy(result$fitted, y)[5],
                    accuracy(resultbeta$fitted, y)[5])
      pv_LB[i,,j]<-c(Box.test(result$residuals, lag = 20, type = "Ljung")$p.value,
                     Box.test(resultbeta$residuals, lag = 20, type = "Ljung")$p.value)
      i<-i+1
    }
  }
  mean_MLEs[j,] <- apply(estim[,,j], 2, mean)
  RB[j,] <- (mean_MLEs[j,]-theta)/theta*100 
  RMSE[j,] <- sqrt(apply(estim[,,j],2,var)+(mean_MLEs[j,]-theta)^2)
  contAIC[j,]<-apply(apply(AIC[,,j],1,rank)==1,1,sum)
  contRMSER[j,]<-apply(apply(RMSER[,,j],1,rank)==1,1,sum)
  contMAPE[j,]<-apply(apply(MAPE[,,j],1,rank)==1,1,sum)
  pv_LB010<-apply(pv_LB[,,j]>0.1,2,sum)
  pv_LB05<-apply(pv_LB[,,j]>0.05,2,sum)
  pv_LB001<-apply(pv_LB[,,j]>0.01,2,sum)
}

saving<-paste0("results_simu/simu_true_",true,"_scen",scen,"vn",vn,"R",R,".RData")
save.image(saving)
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)
