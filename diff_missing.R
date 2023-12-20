
setwd("~/work/data_borrow/")
#setwd("C:/Users/Daxuan Deng/Desktop/data borrow/sim/")
rm(list=ls(all=TRUE))
library(jomo)
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
tauT = c(6.5, 1, -2) # true parameter for additional drop-out
n_mi = 40 # number of multiple imputation

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
n.cores <- 40

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK" # PSOCK for Windows platform and FORK for UNIX
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)

iteration<-1:1000
timestart<-Sys.time()
result = foreach(tt = iteration, .combine = comb_fun,
                 .packages = plist) %dopar% 
  {
    set.seed(tt*10)                 
    source("functions.R")
    rslt = list()                
    
    # (x_reg, yy) main data, (x,y) secondary data
    
    # generate main data with missingness
    x1 = rep(runif(n), each = time)
    x2 = rep(log10(1:time), n) + rnorm(length(id), sd=0.3)
    
    x3_0 = rmvbin(n,c(rep(0.4,time), rep(0.5,time)), bincorr=RR)
    x3_m = as.vector(t(x3_0[, 1:time]))
    x3_a = as.vector(t(x3_0[, -c(1:time)]))
    
    x_main = cbind(1, x1, x2, x3_m)
    x_mis_0 = cbind(1, x2)
    x_mis_1 = cbind(1, x3_a)
    x = cbind(1, x1, x2, x3_a)
    x_reg = x_main[ind, ]
    
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
    
    # add missing on the secondary outcome
    aux_data = as.data.frame(cbind(id, x, y))
    colnames(aux_data) = c("id", colnames(x), "response")
    
    dt_aux = data_sim2(data = aux_data, id = id, 
                       x_mis = x_mis_1, para = tauT, lag_level = 1)
    no_mis_ind_1 = dt_aux$ind
    no_mis_ind_2 = no_mis_ind_0 * no_mis_ind_1
    y = dt_aux$response
    y[no_mis_ind_2 == 0] = NA
    
    rslt$prob_no_mis_aux = mean(no_mis_ind_2[ind])
    
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
    weight_0 = unlist(adjusted_w)
    weight=weight_0[ind]
    no_mis_ind = (weight != 0)
    rslt$prob_no_mis = mean(no_mis_ind)
    rslt$weight = summary(weight)
    
    
    # get initial values for main regression analysis
    fit<-glm(yy_mis~x_reg-1, family = "binomial")
    beta_initial<-fit$coefficients # main model parameter
    
    eta_initial<-mis_fit$coefficients # missing model parameter
    
    
    # calculate unweighted estimator with ipw using PP method
    beta_ipw_pp<-multiroot(f = nf_ipw_pp,
                           start = as.vector(beta_initial))$root
    rslt$beta_ipw_pp = beta_ipw_pp
    
    
    # calculate SE of unweighted estimator with ipw using PP method
    Vols_beta_ipw_pp<-Var_ols_ipw_pp(beta_ipw_pp)/n
    rslt$Vols_beta_ipw_pp = sqrt(diag(Vols_beta_ipw_pp))
    
    # generate MI
    # as the outcomes are longitudinal, we impute timepoint-by-timepoint
    x_mi = cbind(1, x1[ind],
                 matrix(x2, nrow = n, byrow = T),
                 x3_0)
    colnames(x_mi) = c("inpt", "x1", paste0("x2",1:4),
                       paste0("x3m",1:4), paste0("x3a",1:4))
    
    dt_jm1 = data.frame(yy = as.factor(df_main$response_mis),
                    id = id) %>%
      mutate(timepoint = rep(c("yy1", "yy2", "yy3", "yy4"), n)) %>%
      spread(timepoint, yy) %>%
      select(-id)
    
    dt_jm2 = data.frame(y = y, id = id) %>%
      mutate(timepoint = rep(c("y1", "y2", "y3", "y4"), n)) %>%
      spread(timepoint, y) %>%
      select(-id)
    
    dt_jm = dt_jm1 %>% cbind(dt_jm2)
    
    # impute sequentially
    Y = dt_jm %>% select(y2, yy2) 
    X = dt_jm %>% select(y1, yy1) %>% cbind(x_mi)
    
    fit_jm2 = jomo(Y = Y, X = X, nimp = n_mi, output = 0)
    
    
    beta_ipw_pool = rep()
    ven1_ipw_pool = rep()
    dev_pool = rep()
    vdev_pool = rep()
    # beta_en1_mi_pool = rep()
    # ven1_mi_pool = rep()
    # beta_mi_pool = rep()
    # vols_mi_pool = rep()
    
    for(i in 1:n_mi){
      Y = dt_jm %>% select(y3, yy3) 
      X = fit_jm2 %>% filter(Imputation == i) %>%
        select(y1, yy1, y2, yy2) %>%
        cbind(x_mi)
      
      fit_jm3 =  jomo(Y = Y, X = X, nimp = 1, output = 0)
      
      Y = dt_jm %>% select(y4, yy4) 
      X = fit_jm3 %>% filter(Imputation == 1) %>%
        select(y1, yy1, y2, yy2, y3, yy3) %>%
        cbind(x_mi)
      
      fit_jm4 = jomo(Y = Y, X = X, nimp = 1, output = 0)
      
      y_imp = fit_jm4 %>% filter(Imputation == 1) %>%
        select(y1, y2, y3, y4)
      y_imp = as.vector(t(y_imp))
      
      
      # data borrow part begins here
      r = rep(1, nrow(x))
      dist="gaussian"
      adata<-cbind(y_imp,x)
      
      ##get initial values for the algorithm
      sec_fit<-geese(y_imp~x-1,id=id)
      theta_initial<-sec_fit$beta
      
      ##get estimated theta
      theta<-multiroot(f = theta_ee, start = as.vector(theta_initial))$root
      lambda<-lambda_find(theta)
      
      #calculate informative scores (weights) based on oracle
      total<-wgeef_oracle(theta=theta,adata,r=r,id=id,
                          dist=dist,time=time)
      ZZ<-total[,1:(n)]
      ZZ_der<-total[,(n+1):((ncol(x)+1)*n)]
      Prop_scores<-apply(ZZ,2,function(xx){
        1/(1+t(matrix(lambda,ncol=1))%*%xx)/n})
      # data borrow part ends here
      
      
      #calculate weighted estimator with ipw using PP method
      beta_en1_ipw_pp<-multiroot(f = nf_en1_ipw_pp,
                                 start = as.vector(beta_initial))$root
      beta_ipw_pool = cbind(beta_ipw_pool, beta_en1_ipw_pp)
      
      #calculate asymptotic variance etc of en1 with ipw using PP
      res_ipw_pp = infer_en1_ipw_pp(beta_en1_ipw_pp, rho=rho)
      ven1_ipw_pool = cbind(ven1_ipw_pool, diag(res_ipw_pp$v))
      
      gamma_11 = res_ipw_pp$g
      lambda_1 = res_ipw_pp$l
      S = res_ipw_pp$s
      
      dev = ginv(gamma_11) %*% lambda_1 %*% lambda
      dev_pool = cbind(dev_pool, dev)
      vdev = ginv(gamma_11) %*% lambda_1 %*% 
        (S/n) %*% t(lambda_1) %*% t(ginv(gamma_11))
      vdev_pool = cbind(vdev_pool, vdev)
    }
    
    # pooling MI results using Rubin's rule
    # ipw_pp
    beta_ipw = rowMeans(beta_ipw_pool)
    rslt$beta_en1_ipw_pp = beta_ipw
    
    ven1_ipw = rowMeans(ven1_ipw_pool) + 
      diag((1 + 1/n_mi) * var(t(as.matrix(beta_ipw_pool))))
    rslt$Ven1_beta_ipw_pp = sqrt(ven1_ipw)
    
    cbar_ipw = mean((diag(Vols_beta_ipw_pp) - ven1_ipw) / 
      diag(Vols_beta_ipw_pp))
    
    rslt$criteria_ipw_pp = cbar_ipw
    
    #calculate 95cp
    index_beta_ipw<-as.numeric(betaT <= beta_ipw + 1.96 * sqrt(ven1_ipw) &
                                 betaT >= beta_ipw - 1.96 * sqrt(ven1_ipw))
    rslt$prop_ipw_pp = index_beta_ipw

    return(rslt)
  }
timeend<-Sys.time()
timeend-timestart

parallel::stopCluster(cl = my.cluster)

save.image("result/result_diff_missing.RData")