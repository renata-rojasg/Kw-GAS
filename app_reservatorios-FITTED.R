# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

library(forecast)
library(readr)
library(xtable)
source("GASKwFit-COV.R")
source("GASBetaFit-COV.R")
############################
#### initial quantities ####
############################
w1<-4.5 # width for plots 
h11<-4.5 # height for plots
## Data
southeast_reservoir <- as.data.frame(read_csv("reservoirs_data.csv")[,-1])

R<-1:dim(southeast_reservoir)[2]
acuracia1<-array(NA,c(2,length(R),4))
pv_LB_Kw<-pv_LB_beta<-
  nome<-modelo<-modelo_beta<-c()
final<-final_beta<-list()
start_time <- Sys.time()
j<-0
for(i in R){
  ################################
  #### defining the variables ####
  ################################
  y<-na.omit(southeast_reservoir[,i])
  if(i==2){
    y<-ts(y[1:(length(y)-6)],
          frequency = 12,start = c(2013,1))
  }else{
    y<-ts(y[1:(length(y)-6)],
          frequency = 12,start = c(2011,1))
  }
  t<-1:length(y)
  t_hat<-length(y)+1
  Xsin<-sin(2*pi*t/12)
  Xsin_hat<-sin(2*pi*t_hat/12)
  Xcos<-cos(2*pi*t/12)
  Xcos_hat<-cos(2*pi*t_hat/12)
  X<-cbind(sin(2*pi*t/12),cos(2*pi*t/12))
  X_hat<-cbind(sin(2*pi*t_hat/12),cos(2*pi*t_hat/12))
  ########################
  #### fitting Kw-GAS ####
  ########################
  ressin<-try(GASKw.fit(y,1,1,X=Xsin,X_hat=Xsin_hat,h1=1),silent=T)
  if(class(ressin)=="try-error" || ressin$conv != 0){ressin$aic<-.Machine$double.xmax}
  rescos<-GASKw.fit(y,1,1,X=Xcos,X_hat=Xcos_hat,h1=1)
  ressincos<-try(GASKw.fit(y,1,1,X=X,X_hat = X_hat,h1=1),silent=T)
  if(class(ressincos)=="try-error" || 
     ressincos$conv != 0 ||
     sum(is.nan(ressincos$model[,4]))>0
  ){ressincos$aic<-.Machine$double.xmax}
  res<-GASKw.fit(y,1,1,h1=1)
  minKw<-min(ressin$aic,rescos$bic,ressincos$aic,res$aic)
  if(minKw==ressin$aic){final1<-ressin; modelo1<-"Kw with sin covariate"} else{
    if(minKw==rescos$aic){final1<-rescos; modelo1<-"Kw with cos covariate"} else{
      if(minKw==ressincos$aic){final1<-ressincos; modelo1<-"Kw with sin and cos covariate"} else{
        if(minKw==res$aic){final1<-res; modelo1<-"Kw without covariate"} 
      }}}
  ##########################
  #### fitting Beta-GAS ####  
  ##########################
  resbetasin<-GASbeta.fit(y,1,1,X=Xsin,X_hat=Xsin_hat,h1=1)
  resbetacos<-GASbeta.fit(y,1,1,X=Xcos,X_hat=Xcos_hat,h1=1)
  resbetasincos<-GASbeta.fit(y,1,1,X=X,X_hat = X_hat,h1=1)
  resbeta<-GASbeta.fit(y,1,1,h1=1)
  minbeta<-min(resbetasin$aic,resbetacos$bic,resbetasincos$aic,resbeta$aic)
  if(minbeta==resbetasin$aic){final_beta1<-resbetasin; modelo_beta1<-"beta with sin covariate"} else{
    if(minbeta==resbetacos$aic){final_beta1<-resbetacos; modelo_beta1<-"beta with cos covariate"} else{
      if(minbeta==resbetasincos$aic){final_beta1<-resbetasincos; modelo_beta1<-"beta with sin and cos covariate"} else{
        if(minbeta==resbeta$aic){final_beta1<-resbeta; modelo_beta1<-"beta without covariate"} 
      }}}        
  ###########################
  #### residual analysis ####  
  ###########################
  pv_LB_Kw1<-Box.test(final1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_beta1<-Box.test(final_beta1$residuals, lag = 20, type = "Ljung")$p.value
  #############################
  #### selecting resevoirs ####  
  #############################
  if(((pv_LB_Kw1>0.05)&
      (pv_LB_beta1>0.05))>0){
    j<-j+1
    final[[j]]<-final1
    final_beta[[j]]<-final_beta1
    modelo[j]<-modelo1
    modelo_beta[j]<-modelo_beta1
    nome[j]<-substring(colnames(southeast_reservoir)[i],1,nchar(colnames(southeast_reservoir)[i])-13)
    print(nome[j])
    pv_LB_Kw[j]<-Box.test(final1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_beta[j]<-Box.test(final_beta1$residuals, lag = 20, type = "Ljung")$p.value
    acuracia1[1,j,]<-c(accuracy(final[[j]]$fitted, y)[,c(2,3,5)],final[[j]]$aic)
    acuracia1[2,j,]<-c(accuracy(final_beta[[j]]$fitted, y)[,c(2,3,5)],final_beta[[j]]$aic)
    ###############
    #### plots ####  
    ###############
    #---------------------------------
    # time series
    time_series<-paste0("plots/time_series",i,".eps")
    postscript(time_series,width = w1, height = h11,family = "Times")
    par(mfrow=c(1,1))
    plot(y,main=nome[j])
    lines(final[[j]]$fitted, col=2,lty=2,lwd=2)
    lines(final_beta[[j]]$fitted, col=4,lty=3,lwd=1.5)
    legend("top", 
           c("Original","Kw-GAS","Beta-GAS"),
           col = c(1,2,4),
           lty= c(1,2,3),
           lwd = c(1,2,1.5), bty="n", cex = 1)
    dev.off()
    #---------------------------------
    # seasonality
    months<-paste0("plots/months",i,".eps")
    postscript(months,width = w1, height = h11,family = "Times")
    par(mfrow=c(1,1))
    monthplot(y,main=nome[j])
    dev.off()
    #---------------------------------
    resid_plot<-paste0("plots/resid_plot",i,".eps")
    postscript(resid_plot,width = w1*2, height = h11*2,family = "Times")
    # acf Kw residuals 
    par(mfrow=c(1,2))
    acf(final[[j]]$residuals)
    # acf Beta residuals 
    acf(final_beta[[j]]$residuals)
    dev.off()
  }
}
########################
#### final measures ####  
########################
acuracia1<-acuracia1[,1:j,]
rankRMSE<-apply((apply(acuracia1[,,1],2,rank)==1),1,sum)
rankMAE<-apply((apply(acuracia1[,,2],2,rank)==1),1,sum)
rankMAPE<-apply((apply(acuracia1[,,3],2,rank)==1),1,sum)
rankAIC<-apply((apply(acuracia1[,,4],2,rank)==1),1,sum)

result<-rbind(rankAIC,
              rankMAPE,
              rankRMSE
              )
prop.table(result,1)

type_model_Kw<-table(modelo)
type_model_beta<-table(modelo_beta)

count_LB_Kw<-sum(round(pv_LB_Kw,4)>0.05)
count_LB_beta<-sum(round(pv_LB_beta,4)>0.05)

results<-data.frame(
  AIC_Kw=acuracia1[1,,4],
  AIC_beta=acuracia1[2,,4],
  MAPE_Kw=acuracia1[1,,3],
  MAPE_beta=acuracia1[2,,3],
  RMSE_Kw=acuracia1[1,,1],
  RMSE_beta=acuracia1[2,,1],
  pv_LB_Kw=round(pv_LB_Kw,4),
  pv_LB_beta=round(pv_LB_beta,4))
rownames(results)<-nome

saving<-paste0("app_reservoirs.RData")
save.image(saving)
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)