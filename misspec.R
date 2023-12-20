
setwd("~/work/data_borrow/")
#setwd("C:/Users/Daxuan Deng/Desktop/data borrow/sim/")
rm(list=ls(all=TRUE))
library(rootSolve)
library(geepack)
library(MASS)
library(wgeesel)
library(bindata)
library(psych)
library(tidyverse)
library(foreach)
library(doParallel)
source("functions.R")

### main data parameter config
n=500 # sample size, 2000/1000/500
r_0<-0.89 # control main and secondary outcome correlation

time=4 # four timepoints
correlation<-0.5 # longitudinal main outcome correlation across timepoints
variance<-1
id<-rep(1:n,each=time)
ind = (1:length(id))%%time==0 # index for different timepoints

betaT<-c(1, -1, 0.3, -0.5) # true parameter for the main outcome
thetaT<-c(2, -0.5, 1, -1) # true parameter for the secondary outcomes
etaT = c(2.65, -3, 3) # true parameter for drop-out missingness

# correlation for correlated covariates, x3check and x3acute
# 0.1 across timepoints
# correlated between x3 and x3tilde within each timepoint, 0.6
r1 = matrix(0.1, time, time)
diag(r1) = variance
r2 = matrix(0.1, time, time)
diag(r2) = 0.6
RR = rbind(cbind(r1, r2), cbind(r2, r1))

# parallel computation
plist = c("rootSolve", "geepack", "MASS", "wgeesel",
          "bindata", "psych", "tidyverse")
n.cores <- 4

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)

iteration<-1000
timestart<-Sys.time()
result = foreach(tt = 1:iteration, .combine = comb_fun,
                 .packages = plist) %dopar% 
  {
    set.seed(tt*10)                 
    source("functions.R")
    rslt = list() # results of each replication saved in this variable
    
    # (x_reg, yy) main data, (x,y) secondary data
    
    # generate main data with missingness
    x1 = rep(runif(n), each = time)
    x2 = rep(log10(1:time), n) + rnorm(length(id), sd=0.3)
    
    x3_0 = rmvbin(n,c(rep(0.4,time), rep(0.5,time)), bincorr=RR)
    x3_m = as.vector(t(x3_0[, 1:time]))
    x3_a = as.vector(t(x3_0[, -c(1:time)]))
    
    x_main = cbind(1, x1, x2, x3_m)
    x_mis_0 = cbind(1, x2)
    x = cbind(1, x1, x2, x3_a)
    x_reg = x_main[ind, ]
    
    # generate longitudinal main outcomes with missingness
    dt_main = data_sim(id = id, rho = correlation, phi = 1,
                       x = x_main, beta = betaT, 
                       x_mis = x_mis_0, para = etaT,
                       corstr = "exchangeable",
                       family = "binary", lag_level = 1)
    
    df_main = cbind(dt_main$data[,c("id", "response", "ind",
                                    "response_mis", "ylag1")],
                    x_main[, -1])
    yy = df_main$response[ind]
    yy_mis = df_main$response_mis[ind]
    x_mis = cbind(x_mis_0, df_main$ylag1)
    no_mis_ind_0 = df_main$ind
    
    # generate secondary data
    p_main = 1/(1+exp(-x_main %*% betaT))
    eps = (df_main$response - p_main)/sqrt(p_main*(1-p_main))
    eps_tilde = r_0 * eps + sqrt(1-r_0^2) * rnorm(length(id))
    
    y = x %*% thetaT + eps_tilde # secondary response
    
    
    # compute weight in IPW method
    adjusted_idx=lapply(1:length(unique(id)),function(x){
      idx=which(id==unique(id)[x])
      mis_idx=which(is.na(df_main[idx, "response_mis"])==T)
      if (length(mis_idx)>0) a=idx[2:min(mis_idx)]
      else a=idx[2:length(idx)]
      a
    })
    adjusted_idx=unlist(adjusted_idx)
    
    mis_fit=glm(ind~x2+ylag1,
                data=df_main[adjusted_idx,],family=binomial())
    rslt$eta_ipw_pp = mis_fit$coefficients
    rslt$Vols_eta_ipw_pp = sqrt(diag(vcov(mis_fit)))
    
    adjusted_w=lapply(1:length(unique(id)),function(x){
      idx=which(id==unique(id)[x])
      predict_d=df_main[idx,]
      predict_w=predict(mis_fit,newdata=df_main[idx,],type="response")
      predict_w[1]=1
      predict_w=unlist(lapply(1:length(idx),function(m){
        prod(predict_w[1:m])}))
      res=df_main[idx, "ind"]/predict_w
      res[which(is.na(res))]=0
      return(res)
    })
    weight=unlist(adjusted_w)[ind]
    no_mis_ind = (weight != 0)
    rslt$prob_no_mis = mean(no_mis_ind)
    rslt$weight = summary(weight)
    
    
    # data borrow part begins here
    r<-rep(1,nrow(x))
    dist="gaussian"
    x = x[, -2] # remove x1 for model misspecification
    # x = cbind(1, x1, x2, x3_m) # another misspecification
    adata<-cbind(y,x)
    
    ##get initial values for the algorithm
    sec_fit<-geese(y~x-1,id=id)
    theta_initial<-sec_fit$beta
    
    ##get estimated theta
    theta<-multiroot(f = theta_ee, start = as.vector(theta_initial))$root
    lambda<-lambda_find(theta)
    rslt$lambda = lambda
    
    #calculate informative scores (weights) based on oracle
    total<-wgeef_oracle(theta=theta,adata,r=r,id=id,dist=dist,time=time)
    ZZ<-total[,1:(n)]
    ZZ_der<-total[,(n+1):((ncol(x)+1)*n)]
    Prop_scores<-apply(ZZ,2,function(xx){
      1/(1+t(matrix(lambda,ncol=1))%*%xx)/n})
    # data borrow part ends here
    
    # get initial values for main regression analysis
    fit<-glm(yy_mis~x_reg-1, family = "binomial")
    beta_initial<-fit$coefficients # main model parameter
    
    eta_initial<-mis_fit$coefficients # missing model parameter
    
    # calculate enhanced (weighted) estimator without ipw
    
    beta_en1<-multiroot(f = nf_en1, start = as.vector(beta_initial))$root
    rslt$beta_en1_null = beta_en1
    
    # calculate unweighted estimator without ipw
    
    beta<-multiroot(f = nf, start = as.vector(beta_initial))$root
    rslt$beta_null= beta
    
    #calculate enhanced (weighted) estimator with ipw using FJ method
    alpha_en1_ipw_fj<-multiroot(f = nf_en1_ipw_fj,
                                start = c(beta_initial, eta_initial))$root
    beta_en1_ipw_fj = alpha_en1_ipw_fj[1:length(betaT)]
    eta_en1_ipw_fj = alpha_en1_ipw_fj[-c(1:length(betaT))]
    rslt$beta_en1_ipw_fj = beta_en1_ipw_fj
    rslt$eta_en1_ipw_fj = eta_en1_ipw_fj
    
    #calculate unweighted estimator with ipw using FJ method
    alpha_ipw_fj<-multiroot(f = nf_ipw_fj,
                            start = c(beta_initial, eta_initial))$root
    beta_ipw_fj = alpha_ipw_fj[1:length(betaT)]
    eta_ipw_fj = alpha_ipw_fj[-c(1:length(betaT))]
    rslt$beta_ipw_fj = beta_ipw_fj
    rslt$eta_ipw_fj = eta_ipw_fj
    
    #calculate enhanced (weighted) estimator with ipw using PP method
    
    beta_en1_ipw_pp<-multiroot(f = nf_en1_ipw_pp,
                               start = as.vector(beta_initial))$root
    rslt$beta_en1_ipw_pp = beta_en1_ipw_pp
    
    #calculate unweighted estimator with ipw using PP method
    
    beta_ipw_pp<-multiroot(f = nf_ipw_pp,
                           start = as.vector(beta_initial))$root
    rslt$beta_ipw_pp = beta_ipw_pp
    
    
    #calculate SE of unweighted estimator with ipw using FJ method
    Vols_ipw_fj<-Var_ols_ipw_fj(beta_ipw_fj, eta_ipw_fj)/n
    rslt$Vols_beta_ipw_fj = sqrt(diag(Vols_ipw_fj))[1:length(betaT)]
    rslt$Vols_eta_ipw_fj = sqrt(diag(Vols_ipw_fj))[-c(1:length(betaT))]
    
    #evaluate enhanced (weighted) estimator with ipw using FJ method
    res_ipw_fj = infer_en1_ipw_fj(beta_en1_ipw_fj, eta_en1_ipw_fj)
    V = sqrt(diag(res_ipw_fj$v))
    rslt$Ven1_beta_ipw_fj = V[1:length(betaT)]
    rslt$Ven1_eta_ipw_fj = V[-c(1:length(betaT))]
    rslt$criteria_ipw_fj = res_ipw_fj$c
    rslt$prop_ipw_fj = res_ipw_fj$p
    
    # calculate SE of unweighted estimator with ipw using PP method
    Vols_beta_ipw_pp<-Var_ols_ipw_pp(beta_ipw_pp)/n
    rslt$Vols_beta_ipw_pp = sqrt(diag(Vols_beta_ipw_pp))
    
    # evaluate enhanced (weighted) estimator with ipw using PP method
    res_ipw_pp = infer_en1_ipw_pp(beta_en1_ipw_pp)
    rslt$Ven1_beta_ipw_pp = sqrt(diag(res_ipw_pp$v))
    rslt$criteria_ipw_pp = res_ipw_pp$c
    rslt$prop_ipw_pp = res_ipw_pp$p
    
    return(rslt)
    }
timeend<-Sys.time()
timeend-timestart

parallel::stopCluster(cl = my.cluster)

save.image("result/result_misspec_500.RData")