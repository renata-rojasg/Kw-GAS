# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

source("simu_Kw-GAS-COV.R")
source("GASKwFit-COV.R")
############################
#### initial quantities ####
############################
start_time <- Sys.time()
vn<-c(100,200,500,1000)
R<-5000
bug<-0
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
    # "capivara"
    w=-0.1;A=0.1;B=.9;phi=4.1;beta=0.5
    cov<-"sin"
    theta=c(w,A,B,phi,beta)
  }else{
    if(scenario==scenarios[3])
      scen=3
    # avermelha without covariates
    w<--.1;A=0.1;B=0.7;phi=2.7
    cov<-NA
    X=NA
    X_hat=NA
    theta=c(w,A,B,phi)
  }
}
estim<-array(NA,c(R,length(theta),length(vn)))
RMSE<-mean_MLEs<-RB<-matrix(NA,length(vn),length(theta))
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
      # "capivara"
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
    y <- simu.GASKw(n,A=A,B=B,w=w,phi=phi,beta=beta,X=cov)
    result <- try(GASKw.fit(y,1,1,X=X,X_hat=X_hat,h=1),silent = T)
    if(class(result)=="try-error" || 
       result$conv != 0){
      bug<-bug+1
    }else{
      estim[i,,j]<-result$model[,1]
      i<-i+1
    }
  }
  mean_MLEs[j,] <- apply(estim[,,j], 2, mean)
  RB[j,] <- (mean_MLEs[j,]-theta)/theta*100 
  RMSE[j,] <- sqrt(apply(estim[,,j],2,var)+(mean_MLEs[j,]-theta)^2)
}

saving<-paste0("results_simu/simu_scen",scen,"R",R,".RData")
save.image(saving)
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)
