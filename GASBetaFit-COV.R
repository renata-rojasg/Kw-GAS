# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

GASbeta.fit <- function (y, ar, ma, X=NA, X_hat=NA, h1=6,d=0)#, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid)
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
    loglik <- function(z) 
    {
      w <- z[1]
      A <- z[2:(p1+1)]
      B <- z[(p1+2):(p1+q1+1)]
      alpha = z[p1+q1+2]
      if(k==0)  beta <- as.matrix(0) else beta <- as.matrix(z[(p1+q1+3):(p1+q1+2+k)])
      st<-f<-rep(0,n+m) 
      b<- rep(0,n+m)
      
      for(i in (m+1):n)
      {
        f[i]  <- w + as.numeric(A%*%st[i-ar]) + 
          as.numeric(B%*%f[i-ma]) + X[i,]%*%beta
        b[i]   <- exp(f[i])
        st[i]<- (log(y[i])-digamma(b[i])+digamma(b[i]+alpha))/
          (b[i]^(1-2*d)*(trigamma(b[i])-trigamma(b[i]+alpha))^(1-d))
      }
      b<-exp(f[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, b, alpha, log=TRUE))
      sum(ll)
    } 
    opt <- suppressWarnings(
      optim(c(rep(.1,(p1+q1+k+2))), 
            loglik, method = "BFGS", hessian = T,
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    )
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
   
    ########################################################### 
    #### return Z
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)
    z$coeff <- coef
    w <- coef[1]
    A <- coef[2:(p1+1)]
    B <- coef[(p1+2):(p1+q1+1)]
    alpha <- coef[p1+q1+2]
    if(k==0)  beta <- as.matrix(0) else beta <- coef[(p1+q1+3):(p1+q1+2+k)]
    
    names(coef)<-c("omega",names_A,names_B,"alpha",names_beta)
    z$coeff <- coef
    J_inv <- solve(-(opt$hessian))
    
    z$w <- w
    z$A <- A
    z$B <- B
    z$alpha <- alpha
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
    bhat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      fhat[i] <- w + as.numeric(A%*%sthat[i-ar]) + 
        as.numeric(B%*%fhat[i-ma]) + X[i,]%*%beta
      bhat[i]   <- exp(fhat[i])  
      sthat[i] <- (log(y[i])-digamma(bhat[i])+digamma(bhat[i]+alpha))/
        (bhat[i]^(1-2*d)*(trigamma(bhat[i])-trigamma(bhat[i]+alpha))^(1-d))
    }
    y1<-y[(m+1):n]
    
    esperanca <- bhat/(bhat+alpha)
    z$fitted <- ts(c(esperanca),start=start(y),frequency=frequency(y))
    z$fhat <- fhat
    z$sthat <- sthat
    z$serie <- y
    
    ynew_prev <- c(y[n],rep(NA,h1))

    fhatf <- c(fhat[n],rep(0,h1))
    sthatf <-c(sthat[n],rep(0,h1))
    bhatf<- c(bhat[n],rep(0,h1))

    for(i in 2:(h1+1))
    {
      fhatf[i] <-  w + as.numeric(A%*%sthatf[i-ar]) 
      + as.numeric(B%*%fhatf[i-ma]) +  X_hat[i-1,]%*%beta
      bhatf[i]   <- exp(fhatf[i])
      ynew_prev[i] <- bhatf[i]/(bhatf[i]+alpha)
      sthatf[i] <- (log(ynew_prev[i])-digamma(bhatf[i])+digamma(bhatf[i]+alpha))/
        (bhatf[i]^(1-2*d)*(trigamma(bhatf[i])-trigamma(bhatf[i]+alpha))^(1-d))
    }

    z$forecast <- ynew_prev[2:(h1+1)]
    
    ##############################################
    # residuals
    z$residuals<-rep(0,n) 
    for(i in (m+1):n){
      z$residuals[i] <- qnorm(pbeta(y[i], bhat[i], alpha))
    }
  }
  return(z)
}
