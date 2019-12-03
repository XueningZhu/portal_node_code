
### This is the demo file of simulation
### we specify lagY as covariates here

source("simulator.R") # functions for generating responses and network structure
source("select_func.R") # functions for screening and post estimation

### simulation setup
N = 100 # N (network size)
Time = 100
K = 5 # number of blocks for block model

### true parameters
ns = 10  # number of portal nodes
ds = 0.2 + 0.05*(0:9)  # portal nodes coefficients (dj)
gamma = c(1, 2, 3, 4, 5) # exogenous covariates
p = length(gamma)
betaY = c(0.15, 0.05) # autoregression coefficients (lag 1 & 2)
delta_case = 1 # delta = 1/2


if (delta_case == 1)
{
  ### for delta = 1/2
  delta = 1/2
  lambs = sort(c(0.001, 0.003, 0.005, 0.008, 0.01,seq(0.025, 0.1, 0.03), 0.04, 0.07, 0.1, 0.12))
  lamb_range = list(c(0.01, 0.1), c(0.0001, 0.008), c(0.04, 0.2), c(0.025, 0.1))
}
if (delta_case == 2)
{
  ### for delta = 1/4
  delta = 1/4
  lambs = sort(c(0.0008, 0.001, 0.004, 0.005, 0.008, seq(0.01, 0.08, 0.02), 0.1))
  lamb_range = list(c(0.01, 0.07), c(0.0001, 0.008), c(0.01, 0.1), c(0.01, 0.1))
}

### estimation methods
Methods = c("alasso", "mcp", "scad")
n_met = length(Methods) + 1 # plus lasso
set.seed(1234)

Amat = getPowerLawWs(N, s = ns*2, delta = delta, alpha = 2.5, normalize = F) # generate adjacency matrix
W = Amat/rowSums(Amat) # row-normalized adjacency matrix
ss = 1:ns # portal nodes indexes

### true beta parameter
beta = rep(0, N)
beta[ss] = ds
n_sel = ceiling(N/log(N)) # number of selected nodes

### calculate D WD and (I - WD)^{-1}
D = Diagonal(x = beta)
WD = W%*%D
In = Diagonal(N)
S1 = solve(In - WD)

Xlist = get.Xlist(Time+50, p) # generate X variables with Time + 50 (burn time periods)
dat = get.AR.Ymat(N, Time, S1, Xlist, gamma, betaY, burnT = 50) # generate Y (the first 50 are burning period)
Ymat = dat$Ymat
Xlist = dat$Xlist

### conduct screening
r2 = get.R2(Ymat, Xlist, W)  # calculate R2
r2_order = order(r2, decreasing = T)
ss_hat = r2_order[1:n_sel] # obtain nodes with highest r2
(max(match(ss, r2_order)))
(sum(is.element(ss_hat, ss))) # TP for the portal nodes screening
min(r2[ss]) # minimal r2 for all true portal nodes

### post estimation
Mind = sort(ss_hat)
# use hbic to select lambdas for all methods (lamb_range gives different lamb ranges for different methods)
est = HBIC.select.lambda(Ymat, Xlist, W, Mind, 
                         methods = c("alasso", "mcp", "scad"),
                         lambs = lambs, lamb_range = lamb_range,
                         true_beta = beta, true_gam = c(1,gamma, betaY), 
                         verbose = T, verbose.details = F, verbose_pre = paste(""))

# estimation performance
est$est_tp
est$est_fp
est$est_tm
est$est_beta_mse
est$est_gam_mse
