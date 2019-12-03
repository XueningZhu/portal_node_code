
library(MASS)
library(Matrix)
library(poweRlaw)                                                                                     ### for power-law distribution
library(methods)


##############################################################################
### 1. generate X and Y

### generate X covariates
get.Xlist<-function(Time, p, rho = 0.5)
{
  Xsig = rho^abs(outer(1:p, 1:p,"-"))
  Xlist = lapply(1:Time, function(i){
    X = mvrnorm(n = N, mu = rep(0,nrow(Xsig)), Sigma = Xsig) # generate X matrix
  })
  return(Xlist)
}

### generate the reponses with autoregression effects
get.AR.Ymat<-function(N, Time, S1, Xlist, gamma, beta, burnT = 50, intercept = 1)
{
  ### S1 = solve(I - WD)
  lagY = length(beta)
  Ymat = matrix(0, nrow = N, ncol = Time + burnT)
  for (t in (lagY+1):(Time + burnT))
  {
    eps = rnorm(N)
    
    ### generate Y_t by S^{-1}(intercept+ Xt%*%gamma+ beta1*Y_{t-1} + beta2*Y_{t-2}+...+eps_t)
    Ymat[,t] = as.vector(S1%*%(intercept+Xlist[[t]]%*%gamma + 
                                 Ymat[,(t-1):(t-lagY), drop = F]%*%beta + eps))
    if (t>50)
    {
      ### update Xlist to include Y lag terms
      Xlist[[t]] = cbind(Xlist[[t]], Ymat[,(t-1):(t-lagY), drop = F])
    }
      
  }
  Xlist = Xlist[-(1:50)]
  Ymat = Ymat[,-(1:50)]
  
  return(list(Ymat = Ymat, Xlist = Xlist))
}


##############################################################################
### 1. generate W

### generate Dyad Independence Network
getDyadW<-function(N, N1 = 10, delta = 1.2, normalize = T)                                             ### simulate Dyad network W; N1: number of mutual pairs 2N1*N, delta: P((0,1)) = P((1,0)) = 0.5*N^{delta}
{
  A = matrix(0, nrow = N, ncol = N)                                                                    ### use A to store network structure
  
  ######################################### mutual follower ###########################################
  ind = which(upper.tri(A), arr.ind = T)                                                               ### record all the index of upper triangular matrix
  indM = ind[sample(1:nrow(ind), N*N1),]                                                               ### sample N*N1 as mutual pairs in the upper triangular matrix
  A[indM] = 1                                                                                          ### the corresponding links are set to be 1
  A[indM[,2:1]] = 1                                                                                    ### the following matrix is set to be symmetric, as a result, mutual pairs are 2N1*N
  
  ######################################### single relationship #######################################
  ind1 = which(A==0&upper.tri(A), arr.ind = T)                                                         ### record all the zero pairs in the upper triangular matrix
  indS = ind1[sample(1:nrow(ind1), N^delta),]                                                          ### choose N^delta index as single relations
  tmp = sample(1:nrow(indS), floor(N^delta/2))                                                         ### randomly choose 0.5*N^delta to be in the lower triangular matrix
  indS[tmp,] = indS[tmp, 2:1]                                                                          ### change the corresponding index to be inverse
  A[indS] = 1                                                                                          ### the single following relation is set to be 
  diag(A) = 0                                                                                          ### aii = 0
  if (!normalize)
    return(A)
  W = A/rowSums(A)                                                                                     ### W is row-normalized
  W = as(W, "dgCMatrix")
  return(W)
}

### generate Power-law Distribution Network
getPowerLawW<-function(N, alpha, normalize = T)                                                        ### get power-law network W
{
  Nfollowers = rpldis(N, 1, alpha)                                                                     ### generate N random numbers following power-law(1, alpha): k1-kN
  A = sapply(Nfollowers, function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  diag(A) = 0
  ind = which(rowSums(A)==0)                                                                           ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1                                                                ### for those node, randomly select 3 followees
  }
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  W = as(W, "dgCMatrix")
  return(W)
}


### generate Stochastic Block Network
getBlockW<-function(N, Nblock, normalize = T)                                                          ### get block network
{
  if (N%%Nblock==0){                                                                                   ### if N mod Nblock is integer
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)                        ### obtain the diagnal block list
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),                       ### generate following relations within the blocks
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)        ### if N mod Nblock is not integer
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-0.8}),                  ### generate following relations within the blocks
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-0.8}),                 ### generate following relations within the blocks
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)                                                                           ### combine the blocks in matrix
  offDiag = which(isDiag == 0, arr.ind = T)                                                            ### to calculate the index of the off digonal indexes
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                                          ### people between blocks have 0.3 prob to follow
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  ################ transform bA to be a symmetric matrix ##############################################
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  
  diag(bA) = 0
  ind = which(rowSums(bA)==0)                                                                          ### in case some row sums are zero
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                                                               ### for those node, randomly select 3 followees
  }
  
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                                                                   ### row normalize bA
  W = as(W, "dgCMatrix")
  return(W)
}



### generate W for stochastic block model in portal nodes paper
getBlockWS<-function(N, s, delta, alpha, Nblock, normalize = T, r = 2)
{
  A = matrix(0, N, N)
  ### generate portal nodes followers by N^delta
  Nfols = rep(floor(r*N^delta), s)
  ### generate network structure related to portal nodes (the first s columns for A)
  As = sapply(Nfols, function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  A[,1:s] = As
  ### generate the network structure between non-portal nodes by block model
  A[(s+1):N, (s+1):N] = getBlockW(N-s, Nblock, normalize = F)
  
  
  ### generate portal nodes with lower in-degrees
  nfol = colSums(A)
  n_deg = ceiling(median(nfol[-(1:s)]))
  A_portal_low = sapply(rep(n_deg, 2), function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  A[,5:6] = A_portal_low
  diag(A) = 0
  
  ind = which(rowSums(A)==0)                                                                           ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1                                                                ### for those node, randomly select 3 followees
  }
  
  
  diag(A) = 0
  
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  W = as(W, "dgCMatrix")
  return(W)
}

### generate W for dyad independence model in portal nodes paper
getDyadWS <- function(N, s, delta, alpha, N1 = 3, delta1 = 1.2, normalize = T, r = 1)
{
  A = matrix(0, N, N)
  ### generate portal nodes followers by N^delta
  Nfols = rep(floor(r*N^delta), s)
  ### generate network structure related to portal nodes (the first s columns for A)
  As = sapply(Nfols, function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  A[,1:s] = As
  ### generate the network structure between non-portal nodes by dyad independence model
  A[(s+1):N, (s+1):N] = getDyadW(N - s, N1 = N1, delta = delta1, normalize = F)
  
  ### generate portal nodes with lower in-degrees
  nfol = colSums(A)
  n_deg = ceiling(median(nfol[-(1:s)]))
  A_portal_low = sapply(rep(n_deg, 2), function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  A[,5:6] = A_portal_low
  diag(A) = 0
  ind = which(rowSums(A)==0)                                                                           ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1                                                                ### for those node, randomly select 3 followees
  }
  
  
  diag(A) = 0
  
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  W = as(W, "dgCMatrix")
  return(W)
}

### generate W for power-law distribution model in portal nodes paper
getPowerLawWs<-function(N, s, delta, alpha, normalize = T, r = 1)                                                        ### get power-law network W
{
  ### generate portal nodes followers by N^delta
  Nfols = rep(floor(r*N^delta), s)
  Nfollowers = c(Nfols, rpldis(N - s, 2, alpha))                                                                     ### generate N random numbers following power-law(1, alpha): k1-kN
  A = sapply(Nfollowers, function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  diag(A) = 0
  
  ### generate portal nodes with lower in-degrees
  nfol = colSums(A)
  n_deg = ceiling(median(nfol[-(1:s)]))
  A_portal_low = sapply(rep(n_deg, 2), function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  A[,5:6] = A_portal_low
  
  diag(A) = 0
  ind = which(rowSums(A)==0)                                                                           ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1                                                                ### for those node, randomly select 3 followees
  }
  
  diag(A) = 0
  
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  W = as(W, "dgCMatrix")
  return(W)
}

