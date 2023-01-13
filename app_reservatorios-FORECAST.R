# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

####################
#### R packages ####
####################
library(forecast)
source("GASKwFit-COV.R")
source("GASBetaFit-COV.R")
library(readr)
############################
#### initial quantities ####
############################
start_time <- Sys.time()
## Data
southeast_reservoir <- as.data.frame(read_csv("reservoirs_data.csv")[,-1])

R<-dim(southeast_reservoir)[2]
h<-6
out_forecast<-array(NA,c(h,2,R))
ac<-array(NA,c(2,R,2))
nome<-c()
model<-c("Kw with sin and cos covariate", "Kw with sin and cos covariate",
         "Kw with sin and cos covariate", "Kw with sin covariate",        
         "Kw with sin and cos covariate", "Kw with sin covariate",        
         "Kw with sin covariate"        , "Kw with sin covariate",        
         "Kw with sin and cos covariate", "Kw with sin and cos covariate",
         "Kw with sin and cos covariate")
model_beta<-c("beta with sin and cos covariate", "beta with sin and cos covariate",
              "beta with sin and cos covariate", "beta with sin and cos covariate",
              "beta without covariate"         , "beta with sin covariate",
              "beta with sin covariate"        , "beta with sin covariate",        
              "beta with sin and cos covariate", "beta with sin covariate",        
              "beta with sin covariate"    )

for(i in 1:R){
  nome[i]<-substring(colnames(southeast_reservoir)[i],1,nchar(colnames(southeast_reservoir)[i])-13)
  y1<-na.omit(southeast_reservoir[,i])
  y_prev<-y1[(length(y1)[1]-h+1):length(y1)[1]]
  print(i)
  for(j in 1:h){
    y<-ts(y1[1:(length(y1)[1]-h+j-1)],frequency = 12)
    t<-1:length(y)
    t_hat<-length(y)+1
    if(model[i]=="Kw with sin covariate"){
      X<-sin(2*pi*t/12)
      X_hat<-sin(2*pi*t_hat/12)
    }else{
      if(model[i]=="Kw with sin and cos covariate"){
        X<-cbind(sin(2*pi*t/12),cos(2*pi*t/12))
        X_hat<-cbind(sin(2*pi*t_hat/12),cos(2*pi*t_hat/12))
      }else{
        if(model[i]=="Kw without covariate")
          X<-X_hat<-NA
      }
    }
    if(model_beta[i]=="beta with sin covariate"){
      Xbeta<-sin(2*pi*t/12)
      Xbeta_hat<-sin(2*pi*t_hat/12)
    }else{
      if(model_beta[i]=="beta with sin and cos covariate"){
        Xbeta<-cbind(sin(2*pi*t/12),cos(2*pi*t/12))
        Xbeta_hat<-cbind(sin(2*pi*t_hat/12),cos(2*pi*t_hat/12))
      }else{
        if(model_beta[i]=="beta without covariate")
          Xbeta<-X_beta<-NA
      }
    }
    res<-GASKw.fit(y,1,1,X=X,X_hat=X_hat,h1=1)
    resbeta<-GASbeta.fit(y,1,1,X=Xbeta,X_hat=Xbeta_hat,h1=1,d=0)
    out_forecast[j,1,i]<-res$forecast
    out_forecast[j,2,i]<-resbeta$forecast
  }
  assign(paste0("y_prev",i),y_prev)
  ac[1,i,]<- accuracy(out_forecast[,1,i], y_prev)[,c(2,5)]
  ac[2,i,]<- accuracy(out_forecast[,2,i], y_prev)[,c(2,5)]
}

rankRMSE_out<-apply((apply(ac[,,1],2,rank)==1),1,sum)
rankMAPE_out<-apply((apply(ac[,,2],2,rank)==1),1,sum)
result_out<-rbind(rankRMSE_out,
                  rankMAPE_out)
print(result_out)
print(prop.table(result_out,1))
xtable(t(matrix(paste0(round(result_out,3)," (",
                       round(prop.table(result_out,1),4)*100,")"),2,2)))

results_out<-data.frame(
  MAPE_Kw_out=ac[1,,2],
  MAPE_beta_out=ac[2,,2],
  RMSE_Kw_out=ac[1,,1],
  RMSE_beta_out=ac[2,,1]
)

rownames(results_out)<-nome
print(results_out)

saving<-paste0("app_reservoirs-FORECAST1step.RData")
save.image(saving)

end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)
