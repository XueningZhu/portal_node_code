
#################################################################################
### 1. Screening

get.R2<-function(Ymat, xlist, Amat)
{
  # calculate N and Time
  Time = ncol(Ymat)
  N = nrow(Ymat)
  # generate diag(A, A, ..., A) (T times)
  A_list = rep(list(Amat), Time)
  AmatT = bdiag(A_list)
  # get (Y1', Y2', ..., Y_T')' (an N*T dimension vector)
  In = Diagonal(n = N)
  dY = Diagonal(x = as.vector(Ymat))
  # get A%*%Y_t and stack for t = 1, ..., T
  AY = AmatT%*%dY%*%kronecker(rep(1, Time), In)
  Y = as.vector(Ymat)
  X = do.call(rbind, Xlist)
  X = cbind(1, X) # add intercept to X
  Xvar = cbind(Y, X) # this Y tilde
  
  XY = t(AY)%*%Xvar # calculate t(Y tilde)%*%X (each column is tilde_Y'%*%X_j)
  R2_1 = rowSums((XY%*%ginv(crossprod(Xvar)))*(XY)) # the numerator
  R2_2 = colSums(AY^2) # the denominator
  R2 = R2_1/R2_2 # calculate R squares
  R2[!is.finite(R2)] = 0
  return(R2)
}

#################################################################################
### 2. Estimation

### 2.1 Penalized Functions
### LASSO for univariate optimization
lasso.d1<-function(z, lamb)
{
  zlam = abs(z) - lamb
  zlam[zlam<0] = 0
  theta = sign(z)*zlam # given z and lambda, calculate theta: Fan and Li (2001)
  return(theta)
}

### LLA for scad 
scad.L1<-function(z, theta, lamb)
{
  wei = scad(theta, lamb, a = 3.7) # the weight is first order derivative of scad function
  est = lasso.d1(z, wei) # p'(lambda)|lambda|, which equals to lasso; see Zou (2008)
  return(est)
}

### scad first-order derivative
scad <- function(theta, lamb, a = 3.7) # the value of a is suggested by Fan and Li (2001)
{
  theta = abs(theta)
  alam = a*lamb - theta
  alam[alam<0] = 0
  pp = lamb*((theta<=lamb) + alam/((a-1)*lamb)*(theta>lamb))
  return(pp)
}


### LLA for mcp
mcp.L1<-function(z, theta, lamb, gam)
{
  wei = mcp(theta, lamb, gam) # the weight is first order derivative of mcp function
  est = lasso.d1(z, wei) # p'(lambda)|lambda|, which equals to lasso; see Zou (2008)
  return(est)
}


### mcp first-order derivative
mcp <- function(theta, lamb, gam) # the value of a is suggested by Fan and Li (2001)
{
  theta = abs(theta)
  alam = lamb - theta/gam
  alam[alam<0] = 0
  return(alam)
}

### adaptive Lasso
ada<-function(theta, lamb, gam = 1)
{
  return(1/abs(theta)^gam*lamb)
}

### LLA for adaptive lasso
ada.L1<-function(z, theta, lamb, gam)
{
  wei = ada(theta, lamb, gam)
  est = lasso.d1(z, wei) # p'(lambda)|lambda|, which equals to lasso; see Zou (2008)
  return(est)
}

### 2.2 Estimation Methods
### estimate L1 problem
portal.L1<-function(Ymat, Xlist, W, lamb, Mind, est = NULL, method, verbose = T, verbose_pre)
{
  N = nrow(Ymat)
  
  # Add intercept
  if (is.null(Xlist[[1]]))
    Xlist = lapply(Xlist, function(x) matrix(1, nrow = N, ncol = 1))
  else
    Xlist = lapply(Xlist, function(x) cbind(1, x))
  
  ### set the initial values of paramters 
  set.initial<-function(est, Xlist, Mind, N, rand = F)
  {
    if (rand)
    {
      # if rand, then set the initial value to follow U[-0.2, 0.2]
      d1 = rep(0, N)
      d1[Mind] = runif(length(Mind), min = -0.2, max = 0.2)
      gamma = rep(0, ncol(Xlist[[1]]))
      return(list(d = d1, gamma = gamma))
    }
    if (is.null(est))
    {
      # if est = NULL, start from 0; this is the case for lasso estimation
      gamma = rep(0, ncol(Xlist[[1]]))
      d1 = rep(0, N)
    }
    else
    {
      # if this is not the lasso case, start from the pre-estimated (by lasso) initial estimators
      gamma = est$gamma
      d1 = rep(0, N)
      d1[Mind] = est$d
    }
    return(list(d = d1, gamma = gamma))
  }
  inis = set.initial(est, Xlist, Mind, N)
  d1 = inis$d; d2 = d1; gamma = inis$gamma
  
  ### use coordinate descent to conduct estimation
  delta = 1
  alpha = 0.5
  iter = 0
  while(delta>10^{-4}&iter<=80)
  {
    if (delta>1.5)
    {
      # if delta is too large, then set the initial estimators randomly
      alpha = alpha/2
      inis = set.initial(est, Xlist, Mind, N, rand = T)
      d1 = inis$d; gamma = inis$gamma
    }
    ### reserve the estimator in m'th step
    d0 = d1
    gamma0 = gamma
    ### set initial value of z v1 v2
    z = rep(0, length(beta))
    v1 = z; v2 = z
    
    for (j in Mind)
    {
      # calculate basic matrices and vectors
      beta = d1 # the length of d1 is N
      D = Diagonal(x = beta)
      WD = W%*%D
      In = Diagonal(N)
      S = In - WD
      S1 = solve(S)
      SY = S%*%Ymat
      Z = do.call(rbind, Xlist)
      
      ### given d, estimate gamma
      gamma = ginv(crossprod(Z))%*%crossprod(Z, as.vector(SY))
      Zgamma = do.call(cbind, lapply(Xlist, function(x) x%*%gamma))
      SYZ = SY - Zgamma
      sig2 = mean(SYZ^2)
      ### conduct Taylor's expansion of l with respect to dj (v1j and v2j are first- and second-order derivatives)
      ### delta_1j v1_1 v2_j are defined just the same as in the article
      delta1j = sig2^{-1}*sum((t(SYZ)%*%W[,j])*(Ymat[j,]))/(N*Time)
      eS1Wj = (S1%*%W[,j])[j]
      wwj = sum(W[,j]^2)
      v1[j] = -Time*eS1Wj + delta1j*N*Time
      v2[j] = Time*eS1Wj^2 + (sig2^{-1}*wwj*sum(Ymat[j,]^2)/(N*Time) - 2*delta1j^2)*N*Time
      
      ### update the new estimator
      z[j] = (v2[j]^{-1}*v1[j]) * alpha + d1[j]
      
      if (method=="lasso")
      {
        ### if using lasso estimation, conduct univariate L1 estimation
        d1[j] = lasso.d1(z[j], lamb)
      }
      if (method=="scad")
      {
        ### if using scad estimation, conduct LLA estimation
        d1[j] = scad.L1(z[j], d1[j], lamb)
      }
      if (method=="mcp")
      {
        ### if using mcp estimation, conduct LLA estimation
        d1[j] = mcp.L1(z[j], d1[j], lamb, gam = 1.5)
      }
      if (method=="alasso" & d2[j]!=0)
      {
        ### if using adaptive lasso estimation, conduct LLA estimation
        d1[j] = ada.L1(z[j], d2[j], lamb, gam = 1)
      }
      
    }
    iter = iter + 1
    
    delta = mean(c(abs(d1-d0), abs(gamma-gamma0)))
    
    if (verbose)
      cat(verbose_pre ,"\t",  delta, iter, "\r")
  }
  return(list(d = d1[Mind], gamma = gamma, iter = iter))
}


### the CCCP algorithm: for ALasso, MCP and SCAD
CCCP<-function(Ymat, Xlist, W, methods, lamb, Mind, verbose = T, verbose_pre, 
               lamb_range = list(c(0.0001, 0.008), c(0.14, 0.2), c(0.025, 0.1)))
{
  ### set appropriate lambda ranges for different methods (should be adjusted according to different dataset)
  lamb_range_all = do.call(rbind, lamb_range)
  lamb_range_mets = lamb_range_all[match(methods, c("alasso", "mcp", "scad")),, drop = F]
  
  ### LLA algorithm
  # Lasso estimation
  N = nrow(Ymat)
  est = portal.L1(Ymat, Xlist, W, lamb/log(N), Mind, est = NULL, method = "lasso", 
                  verbose, verbose_pre = paste(lamb, "CCCP Lasso: "))
  # LLA estimation for each method
  mm = length(methods)
  est_list = list()
  for (m in 1:mm)
  {
    met = methods[m]
    if (lamb_range_mets[m,1]<=lamb&lamb<=lamb_range_mets[m,2])
    {
      ### if lambda is within the pre-specified range, conduct LLA estimation
      est_list[[m]] = portal.L1(Ymat, Xlist, W, lamb = lamb, Mind, est = est, method = met, 
                                verbose, verbose_pre = paste(verbose_pre, met, ":"))
    }else{
      # return all zero vector
      est_list[[m]] = list(d = rep(0, length(Mind)), gamma = rep(0, ncol(Xlist[[1]])+1), iter = 0)
    }
  }
  return(est_list)
}

### 2.3 HBIC criterion
### Log likelihood (first part of hbic)
log.likelihood = function(Ymat, Xlist, W, bhat, Mind, gamma) {
  N = nrow(W)
  Time = ncol(Ymat)
  
  beta = rep(0, N)
  beta[Mind] = bhat
  D = Diagonal(x=beta)
  In = Diagonal(N)
  
  SM = In - W%*%D
  term_1 = log(det(SM))
  
  term_2 = 0
  for (t in 1:Time) {
    term_2 = term_2 + crossprod(SM %*% Ymat[,t] - Xlist[[t]]%*%gamma)
  }
  term_2 = N*Time*log(term_2/(N * Time))/2
  
  res = term_2[1,1]/N/Time*2 - 2/N * term_1
  return(res)
}

### HBIC modified from Wang et al (2013, AOS)
hbic.calc = function(Ymat, Xlist, W, bhat, Mind, gamma) {
  N = nrow(W)
  Time = ncol(Ymat)
  p = length(Mind)
  
  if (is.null(Xlist[[1]]))
    Xlist = lapply(Xlist, function(x) matrix(1, nrow = N, ncol = 1))
  else
    Xlist = lapply(Xlist, function(x) cbind(1, x))
  
  ind = which(bhat != 0)
  # N_ast is the effective sample size for d_j
  if (length(ind)==1)
    N_ast = sum(W[,Mind[ind]]^2)
  else
  {
    N_ast = mean(colSums(W[,Mind[ind]]^2))
    #N_ast = max(colSums(W[,Mind[ind]]^2))
  }
  sample_size = N_ast * Time
  
  # calculate the log-likelihood
  llh = log.likelihood(Ymat, Xlist, W, bhat, Mind, gamma)
  # calculate the penalty for model complexity
  penalty = sum(bhat != 0) * log(log(sample_size)) * log(p) / (sample_size)
  if (!is.finite(penalty))
    penalty = 0
  # penalty = 2*sum(bhat != 0) * log(sample_size)^2 * log(p) / (sample_size)
  # penalty = 2*sum(bhat != 0) * (log(N)+2*log(Time))^2 * log(p) / (sample_size)
  # penalty = sum(bhat != 0) * (log(sample_size) * log(p)) / (N * Time)
  # penalty = sum(bhat != 0) * (log(log(N*Time)) * log(p)) / (N * Time)
  res = llh + penalty
  return(c(llh, penalty, N_ast/N))
}



### HBIC selection for all the methods in simulation: TP FP TM, and MSE are compared
HBIC.select.lambda = function(Ymat, Xlist, W, Mind, methods, lambs, 
                              lamb_range = list(c(0.025, 0.1), c(0.0001, 0.008), c(0.025, 0.1), c(0.14, 0.2)),
                              true_beta, true_gam, 
                              verbose = T, verbose.details = F, verbose_pre)
{
  ### set appropriate lambda ranges for different methods
  lamb_all = lapply(1:length(lamb_range), function(i) lambs[lambs>=lamb_range[[i]][1]&lambs<=lamb_range[[i]][2]])
  lasso_lambs = lamb_all[[1]]
  lamb_mets = lamb_all[-1][match(methods, c("alasso", "mcp", "scad"))]
  
  n_met = length(methods)+1 # number of methods
  Time = ncol(Ymat) 
  
  ### estimation measurements
  est_iter = matrix(0, nrow = n_met, ncol = length(lambs))
  est_tp = matrix(0, nrow = n_met, ncol = length(lambs))
  est_fp = matrix(0, nrow = n_met, ncol = length(lambs))
  est_tm = matrix(0, nrow = n_met, ncol = length(lambs))
  
  est_beta_mse = matrix(0, nrow = n_met, ncol = length(lambs))
  est_gam_mse = matrix(0, nrow = n_met, ncol = length(lambs))
  est_hbic = matrix(0, nrow = n_met, ncol = length(lambs))
  
  llh_list = rep(list(matrix(0, nrow = 3, ncol = length(lambs))), n_met)
  lasso_est = list()
  mets_est = list() # for all the methods besides lasso
  for (i in 1:length(lambs))
  {
    lamb = lambs[i]
    if (is.element(lamb, lasso_lambs))
    {
      ### conduct lasso and scad estimation respectively
      lasso_est[[i]] = portal.L1(Ymat, Xlist, W, lamb, Mind, 
                                 est = NULL, method = "lasso", 
                                 verbose = T, verbose_pre = paste(verbose_pre, lamb, " Lasso: "))
      
      ### record simulation measurements
      est_iter[1,i] = lasso_est[[i]]$iter
      est_tp[1,i] = true.pos(true_beta[Mind], lasso_est[[i]]$d)
      est_fp[1,i] = false.pos(true_beta[Mind], lasso_est[[i]]$d)
      est_tm[1,i] = true.model(true_beta[Mind], lasso_est[[i]]$d)
      est_beta_mse[1,i] = mse(true_beta[Mind], lasso_est[[i]]$d)/length(Mind)
      est_gam_mse[1,i] = mse(true_gam, lasso_est[[i]]$gamma)/length(Mind)
      llh_list[[1]][,i] = hbic.calc(Ymat, Xlist, W, lasso_est[[i]]$d, Mind, lasso_est[[i]]$gamma)
      est_hbic[1,i] = llh_list[[1]][1,i] + llh_list[[1]][2,i]
    }else{
      est_hbic[1,i] = 10000
    }
    
    
    ### conduct estimation for different methods
    mets_est[[i]] = CCCP(Ymat, Xlist, W, methods, lamb, Mind, verbose = T, verbose_pre = paste(verbose_pre, lamb),
                         lamb_range = lamb_range[-1])
    
    for (m in 1:length(methods))
    {
      ### determine that whether the lambs in the range of corresponding methods
      if (is.element(lamb, lamb_mets[[m]]))
      {
        met_est = mets_est[[i]][[m]]
        est_iter[m+1,i] = met_est$iter
        
        ### record the estimation performance
        est_tp[m+1,i] = true.pos(true_beta[Mind], met_est$d)
        est_fp[m+1,i] = false.pos(true_beta[Mind], met_est$d)
        est_tm[m+1,i] = true.model(true_beta[Mind], met_est$d)
        est_beta_mse[m+1,i] = mse(true_beta[Mind], met_est$d)/length(Mind)
        est_gam_mse[m+1,i] = mse(true_gam, met_est$gamma)/length(Mind)
        llh_list[[m+1]][,i] = hbic.calc(Ymat, Xlist, W, met_est$d, Mind, met_est$gamma)
        est_hbic[m+1,i] = llh_list[[m+1]][1,i] + llh_list[[m+1]][2,i]
      }else{
        ### if lamb is not in the range of lambs of a method, just give hbic a very large number
        est_hbic[m+1,i] = 10000
      }
      
    }
    if (verbose.details)
    {
      ### print the estimation details
      cat(lamb, "\t", est_tp[,i], "\t", est_fp[,i], "\t", est_beta_mse[,i], "\t",
          est_gam_mse[,i], "\t", est_hbic[,i], 
          "\n")
    }
    
  }
  
  ### select lambda
  #ind_lasso = which.min(est_hbic[1,])
  ind_met = apply(est_hbic, 1, which.min)
  ind_arr = cbind(1:n_met, ind_met)
  est_list = list(lasso_est[[ind_met[1]]])
  for (m in 2:n_met)
  {
    est_list[[m]] = mets_est[[ind_met[m]]][[m-1]]
  }
  lamb = sapply(ind_met, function(x) lambs[x])
  
  return(list(est = est_list, 
              lamb = lamb,
              est_tp = est_tp[ind_arr],
              est_fp = est_fp[ind_arr],
              est_tm = est_tm[ind_arr],
              est_beta_mse = est_beta_mse[ind_arr],
              est_gam_mse = est_gam_mse[ind_arr],
              iter = est_iter[ind_arr],
              mets_est = mets_est, lasso_est = lasso_est
              ))
}



### 2.4 Simulation mesurements
### simulation evaluation measurements
# TP
true.pos = function(beta, bhat) {
  sum((bhat != 0) * (beta != 0))
}

# FP
false.pos = function(beta, bhat) {
  sum((bhat != 0) * (beta == 0))
}

# TM
true.model = function(beta, bhat) {
  tm = all((bhat != 0)==(beta != 0))
  return(tm)
}

# estimation MSE
mse = function(beta, bhat) {
  sum((beta-bhat)^2)
}

# estimation RMSE
RMSE = function(x){
  sqrt(mean(x^2))
}

# median for each row
rowMedian<-function(x)
{
  apply(x,1,median)
}


