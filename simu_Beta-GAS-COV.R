# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

simu.GASbeta <- function(n,A=NA,B=NA,w=1,alpha=2, beta=0, X=NA,freq=12,d=0)
{
  ar<-NA
  ma<-NA
  beta<-as.matrix(beta,1,2)
  
  ##### X definitions
  if(is.na(X)){
    X<-matrix(0, c(n,1))
  }else{
    if(X=="cos"){
      X=as.matrix(cos(2*pi*(1:n)/12))
      if(beta==0) stop("Inform the value of beta")
    }else
      if(X=="sin"){
        X=as.matrix(sin(2*pi*(1:n)/12))
        if(beta==0) stop("Inform the value of beta")
      }else
        if(X=="sin&cos"){
          X=cbind(sin(2*pi*(1:n)/12),cos(2*pi*(1:n)/12))
          if(dim(beta)[1]!=2) stop("Inform the value of beta 2")
        }
  }
  
  ###### GASp model
  if(any(is.na(A)==F))
  {
    ar <- 1:length(A)
  }
  
  ###### GASq model
  if(any(is.na(B)==F))
  {
    ma <- 1:length(B)
  }
  
  ###### GASpq model
  if(any(is.na(A)==F) && any(is.na(B)==F))
  {
    p <- max(ar)
    q <- max(ma)
    m <- max(p,q,na.rm=T)
    p1 <- length(ar)
    q1 <- length(ma)
    
    st <- f <- rep(0,n+m) # E(error)=0 
    b <- y <- rep(NA,n+m)
    
    for(i in (m+1):(n+m))
    {
      f[i]  <- w + as.numeric(A%*%st[i-ar]) + 
        as.numeric(B%*%f[i-ma]) + X[i-m,]%*%beta
      b[i]   <- exp(f[i])
      y[i]    <- rbeta(1, b[i],alpha)
      st[i]<- (log(y[i])-digamma(b[i])+digamma(b[i]+alpha))/
        (b[i]^(1-2*d)*(trigamma(b[i])-trigamma(b[i]+alpha))^(1-d))
    }
    
  }
  return( ts(y[(m+1):(n+m)],frequency=freq) )
}



