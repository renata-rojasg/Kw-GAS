# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

source("kum-q-phi.R")
simu.GASKw <- function(n, A=NA, B=NA, w=1 , phi=0.5, beta=0, X=NA, tau=0.5, link = "logit"){
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog"))){
    stats <- make.link(linktemp)
  }  else {
    stop(paste(linktemp, "link not available, available links are \"logit\", ","\"probit\" and \"cloglog\""))
  }
  
  link = linktemp 
  linkfun = stats$linkfun
  linkinv = stats$linkinv 
  mu.eta = stats$mu.eta 
  diflink = mu.eta
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
  if (any(is.na(A) ==F ))
  {
    ar <- 1:length(A)
  }
  # ##### GASq model
  if(any(is.na(B) ==F ))
  {
    ma <-1 : length(B)
  }
  # ##### GASpq model
  if(any(is.na(A) ==F ) && any (is.na(B) ==F ))
  {
    p <- max ( ar )
    q <- max ( ma )
    m <- max ( p, q, na.rm=T )
    p1 <- length( ar )
    q1 <- length( ma )
    # inicializa os vetores dos dados
    st <- f <- rep(0, n+m) # E(erro) =0
    mu <- y <- rep(NA, n+m)
    
    for ( i in (m+1) : (n+m)){
      f[i] <- w + as.numeric(A%*%st[i-ar]) + as.numeric(B%*%f[i-ma]) + X[i-m,]%*%beta
      mu[i] <- linkinv(f[i])
      y[i] <- rkum(1,mu[i],phi)
      st[i] <- ((1-mu[i]^phi)*log(1-mu[i]^phi))/(phi*mu[i]^(phi-1))*
        (log(1-tau)/log(1-mu[i]^phi)*log(1-y[i]^phi)+1)*(diflink(mu[i]))^(-3)
    }
  }
  
  return(ts(y[-(1:m)],frequency = 12))
}
