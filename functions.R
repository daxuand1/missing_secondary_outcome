# function used to calculate variance, mean, 
# and derivative of random variables followed by 
# binomial, poisson, and gaussian
v<-function(x,theta,dist) 
{ 
  n<-length(x[,1])
  reg<-x%*%theta

  if (dist=="gaussian")
  {
    v<-rep(1,n)
    der<-x
    mu<-reg
  }
  else if (dist=="binomial")
  {
    pi<-exp(reg)/(1+exp(reg))
    v<-pi*(1-pi)
    der<-x*as.vector(v)
    mu<-as.vector(pi)
  }
  else if (dist=="poisson")
  {
    v<-exp(reg)
    der<-x*as.vector(exp(reg))
    mu<-as.vector(exp(reg))
  }
  
  list(v=v,der=der,mu=mu)
}

##based on adjusted empirical likelihood (h function in the paper)
wgeef<-function(theta,adata,r=r,id=id,dist=dist,time=time)
{ #full wgee
  y<-adata[,1]
  x<-adata[,-1]
  n<-length(unique(id))
  A<-diag(1,time,time)
  R1<-R2<-R3<-R4<-diag(0,time,time)
  R1[1,1]<-1
  R2[2,2]<-1
  R3[3,3]<-1
  R4[4,4]<-1
  W<-diag(1,time,time)
  V<-v(x,theta,dist)
  wgeef<-rep()
  dwgeef<-rep()
  sum_dwgee<-0
  
  y[which(is.na(y))]<-0
  x[which(is.na(x))]<-0
  for (i in 1:n)
  {
    index_obs<-which(r[((i-1)*time+1):(i*time)]==1)
    AA<-(A*V$v[((i-1)*time+1):(i*time)]^(-0.5))[index_obs,index_obs]
    WW<-W*r[((i-1)*time+1):(i*time)]
    WW_obs<-WW[index_obs,index_obs]
    R_obs1<-R1[index_obs,index_obs]
    R_obs2<-R2[index_obs,index_obs]
    R_obs3<-R3[index_obs,index_obs]
    R_obs4<-R4[index_obs,index_obs]
    e<-(y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)])
    e_obs<-e[index_obs]
    
    if(length(e_obs)==0){
      wgeei1<-wgeei2<-wgeei3<-wgeei4<-matrix(0, nrow = ncol(x), ncol = 1)
      dwgeei1<-dwgeei2<-dwgeei3<-dwgeei4<-matrix(0, nrow = ncol(x), ncol = ncol(x))
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
    }
    else if (length(e_obs)==1)
    {
      wgeei1<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*e_obs)
      wgeei2<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*e_obs)
      wgeei3<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*e_obs)
      wgeei4<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*e_obs)
      dwgeei1<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei2<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei3<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei4<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
      #wgeei134<-rbind(wgeei1,wgeei2)
      #dwgeei134<-rbind(dwgeei1,dwgeei2)
    }else
    {
      wgeei1<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%e_obs
      wgeei2<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%e_obs
      wgeei3<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%e_obs
      wgeei4<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%e_obs
      
      dwgeei1<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei2<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei3<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei4<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
      #wgeei134<-rbind(wgeei1,wgeei2)
      #dwgeei134<-rbind(dwgeei1,dwgeei2)
    }
    wgeef<-cbind(wgeef,wgeei134)
    dwgeef<-cbind(dwgeef,dwgeei134)
    sum_dwgee<-sum_dwgee+dwgeei134
  }
  wgeef_adjusted<-cbind(wgeef,-max(log(n)/2,1)*apply(wgeef,1,mean))
  #l_index<-rep(1:time,time=n)
  dwgeef_adjusted<-cbind(dwgeef,-max(log(n)/2,1)*sum_dwgee/n)
  return(cbind(wgeef_adjusted,dwgeef_adjusted))
}

##based on empirical likelihood (h function in the paper)
wgeef_oracle<-function(theta,adata,r=r,id=id,dist=dist,time=time)
{ #full wgee
  y<-adata[,1]
  x<-adata[,-1]
  n<-length(unique(id))
  A<-diag(1,time,time)
  R1<-R2<-R3<-R4<-diag(0,time,time)
  R1[1,1]<-1
  R2[2,2]<-1
  R3[3,3]<-1
  R4[4,4]<-1
  W<-diag(1,time,time)
  V<-v(x,theta,dist)
  wgeef<-rep()
  dwgeef<-rep()
  sum_dwgee<-0
  
  #z.col<-ncol(z)
  
  y[which(is.na(y))]<-0
  x[which(is.na(x))]<-0
  for (i in 1:n)
  {
    index_obs<-which(r[((i-1)*time+1):(i*time)]==1)
    AA<-(A*V$v[((i-1)*time+1):(i*time)]^(-0.5))[index_obs,index_obs]
    WW<-W*r[((i-1)*time+1):(i*time)]
    WW_obs<-WW[index_obs,index_obs]
    R_obs1<-R1[index_obs,index_obs]
    R_obs2<-R2[index_obs,index_obs]
    R_obs3<-R3[index_obs,index_obs]
    R_obs4<-R4[index_obs,index_obs]
    e<-(y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)])
    e_obs<-e[index_obs]
    
    if(length(e_obs)==0){
      wgeei1<-wgeei2<-wgeei3<-wgeei4<-matrix(0, nrow = ncol(x), ncol = 1)
      dwgeei1<-dwgeei2<-dwgeei3<-dwgeei4<-matrix(0, nrow = ncol(x), ncol = ncol(x))
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
    }
    else if (length(e_obs)==1)
    {
      wgeei1<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*e_obs)
      wgeei2<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*e_obs)
      wgeei3<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*e_obs)
      wgeei4<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*e_obs)
      dwgeei1<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei2<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei3<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei4<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
      #wgeei134<-rbind(wgeei1,wgeei2)
      #dwgeei134<-rbind(dwgeei1,dwgeei2)
    }else
    {
      wgeei1<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%e_obs
      wgeei2<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%e_obs
      wgeei3<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%e_obs
      wgeei4<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%e_obs
      
      dwgeei1<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei2<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei3<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei4<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
      #wgeei134<-rbind(wgeei1,wgeei2)
      #dwgeei134<-rbind(dwgeei1,dwgeei2)
    }
    wgeef<-cbind(wgeef,wgeei134)
    dwgeef<-cbind(dwgeef,dwgeei134)
    #sum_dwgee<-sum_dwgee+dwgeei12
  }
  #wgeef_adjusted<-cbind(wgeef,-max(log(n)/2,1)*apply(wgeef,1,mean))
  #l_index<-rep(1:time,time=n)
  #dwgeef_adjusted<-cbind(dwgeef,-max(log(n)/2,1)*sum_dwgee/n)
  return(cbind(wgeef,dwgeef))
}

###basic functions to solve empirical likelihood

##first derivative of -log EL
R1der<-function(lambda,ZZ)
{
  apply(ZZ,2,function(xx)
  {as.matrix(xx,ncol=1) / 
      as.vector((1+t(lambda) %*% 
                   as.matrix(xx,ncol=1)))}) %*% 
    rep(1,ncol(ZZ))
}

##second derivative of -log EL
R2der<-function(lambda,ZZ)
{
  r2der<-0
  for(i in 1:ncol(ZZ))
  {
    r2der_i<--as.matrix(ZZ[,i],ncol=1) %*%
      t(as.matrix(ZZ[,i],ncol=1)) /
      as.vector(1+t(lambda) %*% 
                  as.matrix(ZZ[,i],ncol=1))^2
    r2der<-r2der+r2der_i
  }
  r2der
}

##-log EL
R0der<-function(lambda,ZZ)
{
  apply(ZZ,2, function (xx){
    log(as.vector(1+t(lambda)%*%as.matrix(xx,ncol=1)))}) %*% rep(1,ncol(ZZ))
}

#function to find lambda, given theta (calculate the tuning parameter)
lambda_find<-function(theta, maxit = 100)
{
  ZZ<-wgeef(theta=theta,adata,r=r,id=id,dist=dist,time=time)[,1:(n+1)]
  
  gamma<-1
  k<-0
  lambda<-rep(0,nrow(ZZ))
  tol<-1e-7
  
  repeat{
    # step 3
    rl<-R1der(lambda,ZZ)
    rll<-R2der(lambda,ZZ) 
    Delta<--ginv(rll)%*%rl
    
    if(mean(abs(Delta))<tol | k > 100)
      {break}
    else{
      # step 4
      repeat{
        delta = gamma * Delta
        
        index_1<-apply(ZZ,2,function(xx){
          ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)})
        index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
        
        if (sum(index_1)>0 | index_2>0){
          gamma<-gamma/2
        }
        else{
          break
        }
      }
      
      # step 5
      lambda<-lambda+delta
      k<-k+1
      gamma<-(k)^(-0.5)
    }
  }
  lambda
}

#### multiroot_method (calculate the parameter in the secondary model)
theta_ee<-function(theta)
{
  total<-wgeef(theta=theta,adata,r=r,id=id,dist=dist,time=time)
  ZZ<-total[,1:(n+1)]
  ZZ_d<-total[,(n+2):((ncol(x)+1)*n+1+ncol(x))]
  theta_ee<-0
  lambda=lambda_find(theta)
  for (i in 1:(n+1))
  {
    scaler<-(1/(1+t(matrix(lambda,ncol=1))%*%ZZ[,i]))
    ee_i<-matrix((t(ZZ_d[,((i-1)*ncol(x)+1):(i*ncol(x))])%*%
                    matrix(lambda,ncol=1)), nrow=ncol(x))*as.vector(scaler)
    theta_ee<-theta_ee+ee_i
  }
  theta_ee
}

# EL estimation with wgee

##based on adjusted empirical likelihood (h function in the paper)
wgeef_wgee<-function(theta,adata,r=r,id=id,dist=dist,time=time)
{ #full wgee
  y<-adata[,1]
  x<-adata[,-1]
  n<-length(unique(id))
  A<-diag(1,time,time)
  R1<-R2<-R3<-R4<-diag(0,time,time)
  R1[1,1]<-1
  R2[2,2]<-1
  R3[3,3]<-1
  R4[4,4]<-1
  W<-diag(1,time,time)
  V<-v(x,theta,dist)
  wgeef<-rep()
  dwgeef<-rep()
  sum_dwgee<-0
  
  #z.col<-ncol(z)
  
  y[which(is.na(y))]<-0
  x[which(is.na(x))]<-0
  for (i in 1:n)
  {
    index_obs<-which(r[((i-1)*time+1):(i*time)]==1)
    AA<-(A*V$v[((i-1)*time+1):(i*time)]^(-0.5))[index_obs,index_obs]
    WW<-diag(weight_0[((i-1)*time+1):(i*time)])
    WW_obs<-WW[index_obs,index_obs]
    R_obs1<-R1[index_obs,index_obs]
    R_obs2<-R2[index_obs,index_obs]
    R_obs3<-R3[index_obs,index_obs]
    R_obs4<-R4[index_obs,index_obs]
    e<-(y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)])
    e_obs<-e[index_obs]
    
    if(length(e_obs)==0){
      wgeei1<-wgeei2<-wgeei3<-wgeei4<-matrix(0, nrow = ncol(x), ncol = 1)
      dwgeei1<-dwgeei2<-dwgeei3<-dwgeei4<-matrix(0, nrow = ncol(x), ncol = ncol(x))
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
    }
    else if (length(e_obs)==1)
    {
      wgeei1<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*e_obs)
      wgeei2<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*e_obs)
      wgeei3<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*e_obs)
      wgeei4<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*e_obs)
      dwgeei1<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei2<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei3<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei4<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
      #wgeei134<-rbind(wgeei1,wgeei2)
      #dwgeei134<-rbind(dwgeei1,dwgeei2)
    }else
    {
      wgeei1<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%e_obs
      wgeei2<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%e_obs
      wgeei3<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%e_obs
      wgeei4<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%e_obs
      
      dwgeei1<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei2<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei3<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei4<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
      #wgeei134<-rbind(wgeei1,wgeei2)
      #dwgeei134<-rbind(dwgeei1,dwgeei2)
    }
    wgeef<-cbind(wgeef,wgeei134)
    dwgeef<-cbind(dwgeef,dwgeei134)
    sum_dwgee<-sum_dwgee+dwgeei134
  }
  wgeef_adjusted<-cbind(wgeef,-max(log(n)/2,1)*apply(wgeef,1,mean))
  #l_index<-rep(1:time,time=n)
  dwgeef_adjusted<-cbind(dwgeef,-max(log(n)/2,1)*sum_dwgee/n)
  return(cbind(wgeef_adjusted,dwgeef_adjusted))
}

##based on empirical likelihood (h function in the paper)
wgeef_oracle_wgee<-function(theta,adata,r=r,id=id,dist=dist,time=time)
{ #full wgee
  y<-adata[,1]
  x<-adata[,-1]
  n<-length(unique(id))
  A<-diag(1,time,time)
  R1<-R2<-R3<-R4<-diag(0,time,time)
  R1[1,1]<-1
  R2[2,2]<-1
  R3[3,3]<-1
  R4[4,4]<-1
  W<-diag(1,time,time)
  V<-v(x,theta,dist)
  wgeef<-rep()
  dwgeef<-rep()
  sum_dwgee<-0
  
  #z.col<-ncol(z)
  
  y[which(is.na(y))]<-0
  x[which(is.na(x))]<-0
  for (i in 1:n)
  {
    index_obs<-which(r[((i-1)*time+1):(i*time)]==1)
    AA<-(A*V$v[((i-1)*time+1):(i*time)]^(-0.5))[index_obs,index_obs]
    WW<-diag(weight_0[((i-1)*time+1):(i*time)])
    WW_obs<-WW[index_obs,index_obs]
    R_obs1<-R1[index_obs,index_obs]
    R_obs2<-R2[index_obs,index_obs]
    R_obs3<-R3[index_obs,index_obs]
    R_obs4<-R4[index_obs,index_obs]
    e<-(y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)])
    e_obs<-e[index_obs]
    
    if(length(e_obs)==0){
      wgeei1<-wgeei2<-wgeei3<-wgeei4<-matrix(0, nrow = ncol(x), ncol = 1)
      dwgeei1<-dwgeei2<-dwgeei3<-dwgeei4<-matrix(0, nrow = ncol(x), ncol = ncol(x))
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
    }
    else if (length(e_obs)==1)
    {
      wgeei1<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*e_obs)
      wgeei2<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*e_obs)
      wgeei3<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*e_obs)
      wgeei4<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*e_obs)
      dwgeei1<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei2<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei3<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei4<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
      #wgeei134<-rbind(wgeei1,wgeei2)
      #dwgeei134<-rbind(dwgeei1,dwgeei2)
    }else
    {
      wgeei1<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%e_obs
      wgeei2<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%e_obs
      wgeei3<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%e_obs
      wgeei4<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%e_obs
      
      dwgeei1<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei2<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei3<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei4<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      
      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
      #wgeei134<-rbind(wgeei1,wgeei2)
      #dwgeei134<-rbind(dwgeei1,dwgeei2)
    }
    wgeef<-cbind(wgeef,wgeei134)
    dwgeef<-cbind(dwgeef,dwgeei134)
    #sum_dwgee<-sum_dwgee+dwgeei12
  }
  #wgeef_adjusted<-cbind(wgeef,-max(log(n)/2,1)*apply(wgeef,1,mean))
  #l_index<-rep(1:time,time=n)
  #dwgeef_adjusted<-cbind(dwgeef,-max(log(n)/2,1)*sum_dwgee/n)
  return(cbind(wgeef,dwgeef))
}

#function to find lambda, given theta and an (calculate the tuning parameter)
lambda_find_wgee<-function(theta)
{
  ZZ<-wgeef_wgee(theta=theta,adata,r=r,id=id,dist=dist,time=time)[,1:(n+1)]
  
  gamma<-1
  c<-0
  lambda<-rep(0,nrow(ZZ))
  tol<-10e-8
  
  repeat{
    rl<-R1der(lambda,ZZ)
    rll<-R2der(lambda,ZZ) 
    theta<--ginv(rll)%*%rl
    if(mean(abs(theta))<tol)
    {break}
    else{
      repeat{
        mm<-0
        repeat{
          theta<-gamma*theta
          index_1<-apply(ZZ,2,function(xx)
          {ifelse(1+t(lambda+theta)%*%as.matrix(xx,ncol=1)<=0,1,0)})
          if (sum(index_1)>0)
          {gamma<-gamma/2
          mm<-mm+1}
          else
          {break}
        }
        index_2<-ifelse(R0der(lambda+theta,ZZ)-R0der(lambda,ZZ)<0,1,0)
        if (index_2==1)
        {gamma<-gamma/2}
        else
        {break}
      }
    }
    lambda<-lambda+theta
    c<-c+1
    gamma<-(c)^(-0.5)
  }
  lambda
}

#### multiroot_method (calculate the parameter in the secondary model)
theta_ee_wgee<-function(theta)
{
  total<-wgeef_wgee(theta=theta,adata,r=r,id=id,dist=dist,time=time)
  ZZ<-total[,1:(n+1)]
  ZZ_d<-total[,(n+2):((ncol(x)+1)*n+1+ncol(x))]
  theta_ee<-0
  lambda=lambda_find_wgee(theta)
  for (i in 1:(n+1))
  {
    scaler<-(1/(1+t(matrix(lambda,ncol=1))%*%ZZ[,i]))
    ee_i<-matrix((t(ZZ_d[,((i-1)*ncol(x)+1):(i*ncol(x))])%*%matrix(lambda,ncol=1)),nrow=ncol(x))*as.vector(scaler)
    theta_ee<-theta_ee+ee_i
  }
  theta_ee
}

# solve weighted estimating equation without ipw
score_function_en1<-function (beta)
{
  sf=rep()
  for (i in 1:n)
  {
    x_i=as.matrix(x_reg[i,],ncol=1)
    mu_i=as.vector(1/(1+exp(-t(x_i)%*%beta)))
    sf_i=x_i*(yy[i]-mu_i)*Prop_scores[i]*no_mis_ind[i]
    sf=cbind(sf,sf_i)
  }
  sf
}
nf_en1<-function (beta)
{
  score_function_en1(beta)%*%rep(1,n)
}

# solve unweighted estimating equation without ipw
score_function<-function (beta)
{
  sf=rep()
  for (i in 1:n)
  {
    x_i=as.matrix(x_reg[i,],ncol=1)
    mu_i=as.vector(1/(1+exp(-t(x_i)%*%beta)))
    sf_i=x_i*(yy[i]-mu_i)*no_mis_ind[i]
    sf=cbind(sf,sf_i)
  }
  sf
}
nf<-function (beta)
{
  score_function(beta)%*%rep(1,n)
}

# solve weighted estimating equation with ipw using FJ method.
score_function_ipw_fj<-function (beta, eta)
{
  sf=rep()
  
  for (i in 1:n)
  {
    if(no_mis_ind[i] == 0){
      sf_i = rep(0, length(beta))
    }
    else{
      x_i=as.matrix(x_reg[i,],ncol=1)
      mu_i=as.vector(1/(1+exp(-t(x_i)%*%beta)))
      predict_i = prod(1/(1+exp(-x_mis[(4*i-2):(4*i),] %*% eta)))
      
      sf_i=x_i*(yy[i]-mu_i)/predict_i
    }
    
    sf_j = rep()
    for(j in (4*i-2):(4*i)){
      if(j %in% adjusted_idx){
        x_j=as.matrix(x_mis[j,],ncol=1)
        mu_j=as.vector(1/(1+exp(-t(x_j)%*%eta)))
        sf_j=cbind(sf_j, x_j*(no_mis_ind_0[j]-mu_j))
      }
      else{
        sf_j = cbind(sf_j, rep(0, length(eta)))
      }
    }
    sf_j = rowSums(sf_j)
    
    sf=cbind(sf,c(sf_i, sf_j))
  }
  
  sf
}

nf_en1_ipw_fj<-function (alpha)
{
  beta = alpha[1:length(betaT)]
  eta = alpha[-c(1:length(betaT))]
  score_function_ipw_fj(beta, eta)%*%Prop_scores
}

# solve unweighted estimating equation with ipw using FJ method.
nf_ipw_fj<-function(alpha)
{
  beta = alpha[1:length(betaT)]
  eta = alpha[-c(1:length(betaT))]
  score_function_ipw_fj(beta, eta)%*%rep(1,n)
}

# solve weighted estimating equation with ipw using PP method.
score_function_en1_ipw_pp<-function (beta)
{
  sf=rep()
  for (i in 1:n)
  {
    x_i=as.matrix(x_reg[i,],ncol=1)
    mu_i=as.vector(1/(1+exp(-t(x_i)%*%beta)))
    sf_i=x_i*(yy[i]-mu_i)*Prop_scores[i]*weight[i]
    sf=cbind(sf,sf_i)
  }
  sf
}
nf_en1_ipw_pp<-function (beta)
{
  score_function_en1_ipw_pp(beta)%*%rep(1,n)
}

# solve unweighted estimating equation with ipw using PP.
score_function_ipw_pp<-function(beta)
{
  sf=rep()
  for (i in 1:n)
  {
    x_i=as.matrix(x_reg[i,],ncol=1)
    mu_i=as.vector(1/(1+exp(-t(x_i)%*%beta)))
    sf_i=x_i*(yy[i]-mu_i)*weight[i]
    sf=cbind(sf,sf_i)
  }
  sf
}
nf_ipw_pp<-function(beta)
{
  score_function_ipw_pp(beta)%*%rep(1,n)
}


# calculate SE of unweighted estimator with ipw using FJ method
partial_der_A = function(beta, eta, w=rep(1,n)){
  der = 0
  
  for (i in 1:n)
  {
    if(no_mis_ind[i] == 0){
      der_i = matrix(0, length(beta), length(beta))
    }
    else{
      x_i=as.matrix(x_reg[i,],ncol=1)
      mu_i=as.vector(1/(1+exp(-t(x_i)%*%beta)))
      predict_i = prod(1/(1+exp(-x_mis[(4*i-2):(4*i),] %*% eta)))
      
      der_i=-x_i%*%t(x_i)*mu_i*(1-mu_i)/predict_i * w[i]
    }
    
    der = der + der_i
  }
  
  der/n
}

partial_der_B = function(beta, eta, w=rep(1,n)){
  der = 0
  
  for (i in 1:n)
  {
    if(no_mis_ind[i] == 0){
      der_i = matrix(0, length(beta), length(eta))
    }
    else{
      x_i=as.matrix(x_reg[i,],ncol=1)
      mu_i=as.vector(1/(1+exp(-t(x_i)%*%beta)))
      
      x = x_mis[(4*i-2):(4*i),]
      x = rbind(x, x[1,]+x[2,], x[1,]+x[3,],
                x[2,]+x[3,], x[1,]+x[2,]+x[3,])
      x = -colSums(x * as.vector(exp(-x%*%eta))) 

      der_i = (yy[i]-mu_i)*x_i%*%x * w[i]
    }
    
    der = der + der_i
  }
  
  der/n
}

partial_der_C = function(beta, eta, w=rep(1,n)){
  der = 0
  
  for(i in adjusted_idx){
    x_i=as.matrix(x_mis[i,],ncol=1)
    mu_i=as.vector(1/(1+exp(-t(x_i)%*%eta)))
    der_i=-x_i%*%t(x_i)*mu_i*(1-mu_i) * w[ceiling(i/4)]
    
    der = der + der_i
  }
  
  der/n
}

partial_der = function(beta, eta, w=rep(1,n)){
  A = partial_der_A(beta, eta, w=w)
  B = partial_der_B(beta, eta, w=w)
  C = partial_der_C(beta, eta, w=w)
  
  cbind(rbind(A, matrix(0, length(eta), length(beta))),
        rbind(B, C))
}

Var_ols_ipw_fj<-function(beta, eta)
{
  Gamma = partial_der(beta, eta)
  sf = score_function_ipw_fj(beta, eta)
  Sigma = sf %*% t(sf)/n
  
  ginv(Gamma) %*% Sigma %*% t(ginv(Gamma))
}

# evaluate weighted estimator with ipw using FJ method
infer_en1_ipw_fj = function(beta, eta, rho=1){
  sf = score_function_ipw_fj(beta, eta)
  Gamma<-partial_der(beta, eta)
  Sigma<-sf %*% t(sf)/n
  Lambda<-sf %*% t(ZZ)/n
  
  s11<-0
  s12<-0
  for (i in 1:n)
  {
    s11_i<-(ZZ[,i]%*%t(ZZ[,i]))
    s11<-s11+s11_i
    s12_i<-ZZ_der[,((i-1)*ncol(x)+1):((i)*ncol(x))]
    s12<-s12+s12_i
  }
  s11<-s11/n
  s12<-s12/n
  omega<-t(s12)%*%ginv(s11)%*%s12
  S<-ginv(s11)-ginv(s11)%*%s12%*%ginv(omega)%*%t(s12)%*%t(ginv(s11))
  
  Ven1_ipw_fj = ginv(Gamma)%*%(Sigma-rho*Lambda%*%S%*%t(Lambda))%*%
    t(ginv(Gamma))/n
  
  criteria_ipw_fj = diag(diag(Vols_ipw_fj*n)^{-0.5})%*%
    ginv(Gamma)%*%Lambda%*%S%*%t(Lambda)%*%
    t(ginv(Gamma))%*% 
    diag(diag(Vols_ipw_fj*n)^{-0.5})
  
  criteria_ipw_fj<-rho*tr(criteria_ipw_fj[1:length(betaT), 1:length(betaT)])
  
  index_beta_ipw_fj<-rep()
  for(i in 1:length(betaT))
  {
    index_beta_ipw_fj_i<-ifelse(betaT[i] <= beta[i]+
                              1.96*sqrt(Ven1_ipw_fj[i,i]) &
                              betaT[i] >= beta[i]-
                              1.96*sqrt(Ven1_ipw_fj[i,i]),1,0)
    index_beta_ipw_fj<-c(index_beta_ipw_fj,index_beta_ipw_fj_i)
  }
  
  return(list(v = Ven1_ipw_fj, c = criteria_ipw_fj, 
              p = index_beta_ipw_fj,
              g = Gamma,
              l = Lambda,
              s = S))
}

# calculate SE of unweighted estimator with ipw using PP method
Var_ols_ipw_pp<-function (beta)
{
  Gamma = partial_der(beta, eta_initial)
  Gamma_11 = Gamma[1:length(beta), 1:length(beta)]
  Gamma_12 = Gamma[1:length(beta), -c(1:length(beta))]
  
  sf = score_function_ipw_fj(beta, eta_initial)
  Sigma = sf %*% t(sf)/n
  Sigma_11 = Sigma[1:length(beta), 1:length(beta)]
  Sigma_22 = Sigma[-c(1:length(beta)), -c(1:length(beta))]
  
  ginv(Gamma_11) %*% 
    (Sigma_11 - Gamma_12 %*% ginv(Sigma_22) %*% t(Gamma_12)) %*% 
    t(ginv(Gamma_11))
}

# evaluate weighted estimator with ipw using FJ method
infer_en1_ipw_pp = function(beta, rho=1){
  Gamma = partial_der(beta, eta_initial)
  Gamma_11 = Gamma[1:length(beta), 1:length(beta)]
  Gamma_12 = Gamma[1:length(beta), -c(1:length(beta))]
  
  sf = score_function_ipw_fj(beta, eta_initial)
  Sigma = sf %*% t(sf)/n
  Sigma_11 = Sigma[1:length(beta), 1:length(beta)]
  Sigma_22 = Sigma[-c(1:length(beta)), -c(1:length(beta))]
  
  Lambda<-sf %*% t(ZZ)/n
  Lambda_1 = Lambda[1:length(beta), ]  
  Lambda_2 = Lambda[-c(1:length(beta)), ]
  
  s11<-0
  s12<-0
  for (i in 1:n)
  {
    s11_i<-(ZZ[,i]%*%t(ZZ[,i]))
    s11<-s11+s11_i
    s12_i<-ZZ_der[,((i-1)*ncol(x)+1):((i)*ncol(x))]
    s12<-s12+s12_i
  }
  s11<-s11/n
  s12<-s12/n
  omega<-t(s12)%*%ginv(s11)%*%s12
  S<-ginv(s11)-ginv(s11)%*%s12%*%ginv(omega)%*%t(s12)%*%t(ginv(s11))
  
  Phi = Gamma_12 %*% ginv(Sigma_22) %*% Lambda_2 %*% S %*% t(Lambda_1)
  
  Ven1_ipw_pp = ginv(Gamma_11)%*%
    (Sigma_11 - Gamma_12 %*% ginv(Sigma_22) %*% t(Gamma_12) -
       rho*(Lambda_1%*%S%*%t(Lambda_1) + Phi + t(Phi)))%*%
    t(ginv(Gamma_11))/n
  
  Vols_ipw_pp = ginv(Gamma_11)%*%
    (Sigma_11 - Gamma_12 %*% ginv(Sigma_22) %*% t(Gamma_12))%*%
    t(ginv(Gamma_11))/n
  
  criteria_ipw_pp<-rho*tr(diag(diag(Vols_beta_ipw_pp*n)^{-0.5})%*%
                        ginv(Gamma_11)%*%
                        (Lambda_1%*%S%*%t(Lambda_1)+Phi+t(Phi))%*%
                        t(ginv(Gamma_11))%*% 
                        diag(diag(Vols_beta_ipw_pp*n)^{-0.5}))
  
  index_beta_ipw_pp<-rep()
  for(i in 1:length(betaT))
  {
    index_beta_ipw_pp_i<-ifelse(betaT[i] <= beta[i]+
                              1.96*sqrt(Ven1_ipw_pp[i,i]) &
                              betaT[i] >= beta[i]-
                              1.96*sqrt(Ven1_ipw_pp[i,i]),1,0)
    index_beta_ipw_pp<-c(index_beta_ipw_pp,index_beta_ipw_pp_i)
  }
  
  return(list(v = Ven1_ipw_pp, vols = Vols_ipw_pp,
              c = criteria_ipw_pp, 
              p = index_beta_ipw_pp,
              g = Gamma_11,
              l = Lambda_1,
              s = S,
              phi = Phi))
}


# other functions

comb_fun = function(list1, list2){
  mapply(cbind, list1, list2, SIMPLIFY=F)
}

data_sim2 <- function(data, id, x_mis, para, lag_level){
  
  complete=F
  x_final=list()
  prob_miss=0
  if(complete==F){
    data$ind <- 1
    #data$prey_mis <- 0
    data$y_first <- NULL;
    N=length(unique(id))
    for (i in 1:N){
      nsubj=max(table(id))
      for (k in 2:nsubj){
        lagy=NULL
        if(lag_level!=0){
          for(lag_i in 1:lag_level){
            ylagname=paste("ylag",lag_i,sep="")
            y_i=data$response[which(id==i)]
            y_new=c(rep(NA,lag_i),y_i[1:(length(y_i)-lag_i)])
            lagy=cbind(lagy,y_new)
            colnames(lagy)[lag_i]=ylagname
          }
        }
        
        x_mis_i=cbind(x_mis[which(id==i),],lagy=lagy)
        x_final[[i]]=x_mis_i
        
        p <- 1/(1+exp(sum(-x_mis_i[k,]*para,na.rm=T)))
        ###depends on the first and second observation;
        ind_p <- rbinom(1,1,p)
        data[data$id==i,]$ind[k] <- ind_p 
        if (data[data$id==i,]$ind[k-1]==0){
          data[data$id==i,]$ind[(k):nsubj] <- rep(0,(nsubj-(k)+1))
        }
      }  
      data$y_first[(i*nsubj-(nsubj-1)):(i*nsubj)]=rep(data[data$id==i,]$response[1],nsubj)
    }
    x_final=do.call("rbind",x_final)
    response_mis=ifelse(data$ind==0,NA,data$response)
    data_final=cbind(data[,c("id","response","ind")],response_mis,x,x_final)
    data<-as.data.frame(data_final)
  }
  return(data)
}

my.double.factor <- function(x) {as.numeric(levels(x))[x]}
