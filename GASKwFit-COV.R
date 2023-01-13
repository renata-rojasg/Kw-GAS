# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

source("kum-q-phi.R")
GASKw.fit <- function (y, ar, ma,X=NA, X_hat=NA, tau=0.5 ,link = "logit", h1=1)
{
  if (min(y) <= 0 || max(y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  
  maxit1<-10000
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
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
  
  ###########################################################
  names_A <- c(paste("A", ar, sep = ""))
  names_B <- c(paste("B", ma, sep = ""))
  
  if(any(is.na(X))==FALSE){
    if(any(is.na(X_hat))==TRUE) 
      stop("You need to inform X_hat")
    X<-as.matrix(X)
    X_hat<-as.matrix(X_hat)
    k = ncol(X)
    names_beta <- c(paste("beta", 1:k, sep = ""))
  }else{
    X <- matrix(0, c(n,1))
    X_hat<- as.matrix(rep(0,h1+1))
    k=0
    names_beta <- NULL
  }
  ###########################################################  
  ######### GAS pq model
  { 
    # print("GAS model",quote=F)
    loglik <- function(z) 
    {
      w <- z[1]
      A <- z[2:(p1+1)]
      B <- z[(p1+2):(p1+q1+1)]
      phi <- z[p1+q1+2]
      if(k==0)  beta <- as.matrix(0) else beta <- as.matrix(z[(p1+q1+3):(p1+q1+2+k)])
      st<-f<-rep(0,n+m) 
      mu<- rep(0,n+m)
      
      for(i in (m+1):n)
      {
        f[i]  <- w + as.numeric(A%*%st[i-p]) + 
          as.numeric(B%*%f[i-q]) + X[i,]%*%beta
        mu[i]   <- linkinv(f[i])
        st[i] <- ((1-mu[i]^phi)*log(1-mu[i]^phi))/(phi*mu[i]^(phi-1))*
        (log(1-tau)/log(1-mu[i]^phi)*log(1-y[i]^phi)+1)*(diflink(mu[i]))^(-3)
        
      }
      mu   <- linkinv(f[(m+1):n])
      y1<-y[(m+1):n]

      ll <- suppressWarnings(log(dkum(y1, mu, phi,tau=tau)))
      sum(ll)
    } 
    
    opt <- try(
      optim(c(rep(.1,(p1+q1+k+2))), loglik, method = "BFGS", hessian = T,
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    ,silent = T)
    
    if(class(opt)=="try-error"){
      opt <- suppressWarnings(
        optim(c(rep(.2,(p1+q1+k+2))), loglik, method = "BFGS", hessian = T,
              control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
      )
    }
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    ########################################################### 
    #### RETURN
    z <- c()
    z$conv <- opt$conv
    coef <- opt$par
    w <- coef[1]
    A <- coef[2:(p1+1)]
    B <- coef[(p1+2):(p1+q1+1)]
    phi <- coef[p1+q1+2]
    if(k==0)  beta <- as.matrix(0) else beta <- coef[(p1+q1+3):(p1+q1+2+k)]
    
    
    names(coef)<-c("omega",names_A,names_B,"phi",names_beta)
    z$coeff <- coef
    J_inv <- solve(-(opt$hessian))
    
    z$w <- w
    z$A <- A
    z$B <- B
    z$phi <- phi
    z$beta <- beta
    
    z$stderror<-sqrt(diag(J_inv))
    z$zstat <- abs(z$coef/z$stderror)
    z$pvalues <- 2*(1 - pnorm(z$zstat) )
    
    z$loglik<-opt$value
    z$k<- (p1+q1+2+k)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
    
    model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
    colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
    z$model <- model_presentation
    
    ########################################################### 
    #### FORECAST
    fhat <- sthat <-rep(0,n) 
    muhat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      fhat[i] <- w + as.numeric(A%*%sthat[i-p]) + 
        as.numeric(B%*%fhat[i-q]) + X[i,]%*%beta
      muhat[i]   <- linkinv(fhat[i])
      sthat[i] <- ((1-muhat[i]^phi)*log(1-muhat[i]^phi))/(phi*muhat[i]^(phi-1))*
        (log(1-tau)/log(1-muhat[i]^phi)*log(1-y[i]^phi)+1)*(diflink(muhat[i]))^(-3)
          }
    y1<-y[(m+1):n]
    
    z$fitted <- ts(muhat,start=start(y),frequency=frequency(y))
    z$fhat <- fhat
    z$sthat <- sthat
    z$serie <- y
    
    
    ynew_prev <- c(y[n],rep(NA,h1))
    
    fhatf <- c(fhat[n],rep(0,h1))
    sthatf <-c(sthat[n],rep(0,h1))
    muhatf<- c(muhat[n],rep(0,h1))
    
    for(i in 2:(h1+1))
    {
      fhatf[i] <-  w + as.numeric(A%*%sthatf[i-p]) + 
        as.numeric(B%*%fhatf[i-q]) +  X_hat[i-1,]%*%beta
      muhatf[i]   <- linkinv(fhatf[i])
      ynew_prev[i] <- muhatf[i]
      nabla <-phi*muhatf[i]^(phi-1)/((1-muhatf[i]^phi)*log(1-muhatf[i]^phi))*(log(.5)/log(1-muhatf[i]^phi)*log(1-ynew_prev[i]^phi)+1)*diflink(muhatf[i])^(-1)
      St <- (phi^2*muhatf[i]^(2*phi-2)/((1-muhatf[i]^phi)^2*log(1-muhatf[i]^phi)^2))^(-1)*(diflink(muhatf[i]))^(-2)
      sthatf[i] <-nabla*St
      
    }
    
    z$forecast <- ynew_prev[2:(h1+1)]
    
    ##############################################
    # residuals
    z$residuals<-rep(0,n) 
      for(i in (m+1):n){
        z$residuals[i] <- qnorm(pkum(y[i],z$fitted[i],phi,tau=tau))
      }
  }
  return(z)
}