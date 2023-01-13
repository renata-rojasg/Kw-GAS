# Reference: Kw-GAS
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

# density function
dkum<-function(y,q,phi,tau=.5)
{
  d<-log(1-tau)/log(1-q^phi)*
    phi*y^(phi-1)*(1-y^phi)^(log(1-tau)/log(1-q^phi)-1)
  d
}

# cumulative distribution function
pkum<-function(y,q,phi,tau=.5)
{
  p<- 1-(1-y^phi)^(log(1-tau)/log(1-q^phi))
  p
}

# quantile function
qkum<-function(u,q,phi,tau=.5)
{
  q<- (1-(1-u)^(log(1-q^phi)/log(1-tau)))^(1/phi)
  q
}

# inversion method for random generation
rkum<-function(n,q,phi,tau=.5)
{
  u<- runif(n)
  y<-(1-(1-u)^(log(1-q^phi)/log(1-tau)))^(1/phi)
  y
}


