
###Composite endpoints

samplesize<-function(mean.lat,mean.bin,var.lat,var.bin,alpha,beta){
  
  n.lat=(var.lat*(qnorm(1-beta)-(qnorm(alpha)))^2)/(mean.lat)^2
  n.bin=(var.bin*(qnorm(1-beta)-(qnorm(alpha)))^2)/(mean.bin)^2
  
  return(c(n.lat,n.bin))
}

powerfunc<-function(mean,var,alpha,n){
  
  power.lat=pnorm((mean/(sqrt(var/n)))-qnorm(1-alpha))
  
  return(power.lat)
}

###Co-primary endpoints 

zfunc<-function(zalph,mean,sigma,n){
  z=(mean/sigma)*sqrt(n/2)-zalph
  return(z)
}

powercoprim<-function(z,Sigma,n){
  
  powerco=pmvnorm(-z,mean=0,sigma = Sigma)
  
  return(powerco)
  
}

##Multiple primary 

multfunc2<-function(z,sigmat){
  pmult = pnorm(z[1])+pnorm(z[2])-pmvnorm(z,mean=c(0,0),lower=c(-Inf,-Inf),sigma=sigmat)
  
  return(pmult)
}

multfunc3<-function(z,sigmat){
  pmult = pnorm(z[1])+pnorm(z[2])+pnorm(z[3])-pmvnorm(c(z[1],z[2]),mean=c(0,0),lower=c(-Inf,-Inf),sigma=matrix(c(sigmat[1,1],sigmat[1,2],sigmat[2,1],sigmat[2,2]),nrow=2))-
    pmvnorm(c(z[1],z[3]),mean=c(0,0),lower=c(-Inf,-Inf),sigma=matrix(c(sigmat[1,1],sigmat[1,3],sigmat[3,1],sigmat[3,3]),nrow=2))-
    pmvnorm(c(z[2],z[3]),mean=c(0,0),lower=c(-Inf,-Inf),sigma=matrix(c(sigmat[2,2],sigmat[2,3],sigmat[3,2],sigmat[3,3]),nrow=2))+
  pmvnorm(z,mean=c(0,0,0),lower=c(-Inf,-Inf,-Inf),sigma=sigmat)
  
  return(pmult)
}



  
