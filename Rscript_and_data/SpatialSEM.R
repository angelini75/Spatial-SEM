############### Script of spatial SEM # ## ### #### ###########
# Purpose        : Include the gestatistical model into SEM
# Maintainer     : Marcos E. Angelini  (angelini75@gmail.com); 
# Contributions  : Gerard B. M. Heuvelink
# Status         : beta
# Note           : 
# sessionInfo(@RStudio desktop)  lenovo ThinkPad T430 (4 cores)
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.3 LTS

# pre-settings ####
rm(list=ls()[])
# Packages
require(lavaan)
require(sp)

# Set working directory
setwd("<YOUR WORKING DIRECTORY>")
# Some useful functions 
name <- function(x) { as.data.frame(names(x))}
# standardised data set 
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}

# Load data
data <- read.csv("data.csv")[,-1]
boxplot(data)
# load variable statistics
ST <- read.csv("stat.csv")
rownames(ST) <- ST$var
ST <- ST[,-1]
ST <- as.matrix(ST)

#### lavaan model ##############################################################
my.model.lv <- '
# Measurement model (lamda and epsilon)
#--------------------#
CEC.Ar =~ 1*CEC.A
CEC.Br =~ 1*CEC.B
CEC.Cr =~ 1*CEC.C
OC.Ar =~ 1*OC.A
OC.Br =~ 1*OC.B
OC.Cr =~ 1*OC.C
Clay.Ar =~ 1*Clay.A
Clay.Br =~ 1*Clay.B
Clay.Cr =~ 1*Clay.C
## Measurement error #
#CEC.A ~~ 0.05 * CEC.A
#CEC.B ~~ 0.05 * CEC.B
# CEC.C ~~ 0.05 * CEC.C
OC.A ~~ 0.05 * OC.A
OC.B ~~ 0.05 * OC.B
OC.C ~~ 0.05 * OC.C
Clay.A ~~ 0.05 * Clay.A
Clay.B ~~ 0.05 * Clay.B
Clay.C ~~ 0.05 *Clay.C

#--------------------#
# Structural model (gamma and betta matrices)
#--------------------#
Clay.Cr ~ dem + vdchn + X + lstm 
Clay.Ar ~ Clay.Cr + evisd + lstm + ndwi.b #+ Y 
Clay.Br ~ Clay.Ar + Clay.Cr + vdchn + twi + ndwi.b + X #+ Y

OC.Ar ~ Clay.Ar + evisd + lstm + ndwi.b 
OC.Br ~ OC.Ar + Clay.Br + evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + Clay.Ar 
CEC.Br ~ Clay.Br + 0*OC.Br
CEC.Cr ~ Clay.Cr + 0*OC.Cr

#------------------#
# Model error covariance (Psi)
#------------------#
CEC.Ar ~~ CEC.Br + 0*CEC.Cr
CEC.Cr ~~ CEC.Br
#OC.Cr ~~ 0*CEC.Br + 0*CEC.Cr + 0*CEC.Ar 

#------------------#
# lavaan suggestions
#------------------#
Clay.Br  ~   lstm
OC.Br  ~      dem
Clay.Ar  ~    twi
OC.Cr  ~      dem
OC.Ar  ~   ndwi.a

OC.Ar  ~  Clay.Br
OC.Br  ~~  Clay.Ar

CEC.Ar  ~     dem
CEC.Cr  ~   evisd
CEC.Br  ~     dem
CEC.Br  ~       X
CEC.Br  ~   evisd

OC.Cr ~~ Clay.Cr
CEC.Ar ~~ Clay.Br
#------------------#
'


# lavaan model calibration ####
fit <- sem(model = my.model.lv,data = data, meanstructure = FALSE, 
                    fixed.x = TRUE) # fixed.x means deterministic Xs

# lavaan model information ####
summary(fit)
# inspect() shows the lavaan matrices
par.list <- inspect(fit) # nonzero values are free parameters
inspect(fit, "est") # shows estimates
# partable is another function to extract parameter information
partable(fit)


#### Including geostatistical model in SEM #####################################
# first we create some objects and functions

# number of samples
N <- nrow(data)
# number of variables
p <- ncol(data)-1 # Y is not used as covariate

# MLIST: this object is a list with all the matrices of Eq. (3) 
MLIST <- inspect(fit, "est") # as follows:
# Beta and Gamma are included in
MLIST$beta
# lambda and kappa (Eq. 2) are included in 
MLIST$lambda
# Theta-delta and Theta-epsilon are included in
MLIST$theta
# System error (Psi and Phi) are included in
MLIST$psi
# we add two more elements: alpha and a (Eq. 10)
MLIST$alpha <- 0.5 #this is a starting value
MLIST$a <- 0.5 #another starting value

# list of free parameters. They will be starting values for the spatial model
par <- lavaan:::lav_model_get_parameters(lavmodel = fit@Model)
start.x <- c(par, MLIST$alpha, MLIST$a)
# Number of free parameters
free <- length(start.x)

# vector of soil property measures and covariate measures
z.all <- as.vector(as.matrix(data[,-18])) # Y is not used as a covariate

# Since we cannot longer use sem() function for calibration, we use the lavaan
# approach (PORT routine) to optimise the maximum likelihood function. For this,
# function to include free parameters into the list of matrices structure MLIST
x2MLIST <- function(x, MLIST) {
  lambda.x <- x[as.vector(par.list$lambda)[as.vector(par.list$lambda)!=0]]
  theta.x <- x[as.vector(par.list$theta)[as.vector(par.list$theta)!=0]]
  psi.x <- x[as.vector(par.list$psi)[as.vector(par.list$psi)!=0]]
  beta.x <- x[as.vector(par.list$beta)[as.vector(par.list$beta)!=0]]
  alpha.x <- x[free-2]
  a.x <- x[free-1]
  
  MLIST$lambda[which(as.vector(par.list$lambda)!=0)]     <- lambda.x
  MLIST$theta[which(as.vector(par.list$theta)!=0)]       <- theta.x
  MLIST$psi[which(as.vector(par.list$psi)!=0)]           <- psi.x
  MLIST$beta[which(as.vector(par.list$beta)!=0)]         <- beta.x
  MLIST$alpha                                            <- alpha.x
  MLIST$a                                                <- a.x
  MLIST
}

# get h: matrix of discances between samples. 
# warning: distances must be in planar units (meters). X and Y can be used 
# standardized but we need to use the same factor for standardization
xy <- data[, c("X","Y")] # coordinates of data
xy[,"Y"] <- xy[,"Y"] * ST["Y","std.dev"] / ST["X","std.dev"]
coordinates(xy) <- ~X+Y
h <- sp::spDists(xy) # h of cal data
# coordinate system was originally set as:
# NAD83.N <- CRS("+init=epsg:2796")

# get.RHO function: this function get the C matrix of Eq. (12), 
# wich depend on alpha and a
get.RHO <- function(MLIST = NULL, h = NULL) {
  a <- MLIST$a
  n <- nrow(h)
  alpha <- MLIST$alpha
  RHO <- matrix(rep(NA,n^2), nrow = n)
  for(i in seq_along(RHO)) {
    RHO[i] <- (1-alpha) * exp(-h[i]/a)
  }
  diag(RHO) <- 1
  RHO
}

# Function to compute Sigma-hat (Eq. 3) = SIGMA0. This is Sigma-hat at distance 
# zero, that has to be used in Eq. (12).
computeSigmaHat.LISREL <- lavaan:::computeSigmaHat.LISREL

# Objective function that has to be optimised
# Objective function:
objective_ML <- function(x, MLIST = NULL) { # x is a vector with the free param.
  MLIST <- x2MLIST(x = x, MLIST = MLIST) # we load x to the model structure 
  # compute Sigma.hat
  SIGMA0 <- computeSigmaHat.LISREL(MLIST = MLIST) 
  RHO <- get.RHO(MLIST,h) 
  # avoid non positive definite matrix by chacking that eigen values are positive
  if (all(eigen(SIGMA0)$values >0) & (all(eigen(RHO)$values >0))) {
    SIGMA0.inv <- chol2inv(chol(SIGMA0))
    RHO.inv <- chol2inv(chol(RHO))
    SIGMA.all.inv <- kronecker(SIGMA0.inv, RHO.inv) # Sigma-all Eq. (12)
    dL.S.R <- append((diag(chol(SIGMA0)))^N, (diag(chol(RHO)))^p)
    logdetSIGMA.all = 2*sum(log(dL.S.R))
    # Eq. (11):
    objective <- -1 * (-1/2*p*N*log(2*pi) - 1/2*logdetSIGMA.all -
                         1/2 * crossprod(z.all, SIGMA.all.inv) %*%z.all)
    cat("objective = ", objective, "\n")
    objective
  } else {
    objective <- Inf
    objective
  }
}

# Optimization of the objective function (calibration)
# warning: it may take several minutes. 
# Not run
out <- nlminb(start = start.x, objective = objective_ML,
              MLIST = MLIST, control = list(iter.max = 200))
# end(Not run)
load("out.RData") # load the outcome of the optimization function

# Extract parameters from out
MLIST.obs <- x2MLIST(out$par, MLIST) 
# Now MLIST.obs contain all model parameters for model prediction
# these coeficients are represented in Fig. 3.

#### Prediction ################################################################
# First we need a grid with locations of interest that contains the values of the
# covariates
s0.un <- read.csv("xy.csv")[,-1] # you have to uncompress the file (~23MB)
s0.un <- s0.un[,names(data)[10:18]]
s0 <- std(x = s0.un, st = ST[10:18,])
summary(s0)
xy.s0 <- s0[, c("X","Y")] # coordinates of data
xy.s0[,"Y"] <- xy.s0[,"Y"] * ST["Y","std.dev"] / ST["X","std.dev"]
# not run
plot(xy.s0) # to see the point distribution
# end(not run)

#### Functions for prediction: ####

# We need to compute a C matrix (Eq. 12) that include the prediction locations s0.
# We call this matrix RHO0. Function to get RHO0 including one prediction location
get.RHO0 <- function(MLIST = NULL, h0 = NULL) {
  a <- MLIST$a
  alpha <- MLIST$alpha
  n <- nrow(h0) # we compute h0 (distances to the s0 location) later on
  RHO0 <- matrix(rep(NA,n^2), nrow = nrow(h0), ncol = ncol(h0))
  for(i in seq_along(RHO0)) {
    RHO0[i] <- (1-alpha) * exp(-h0[i]/a)
  }
  #diag(RHO) <- 1
  RHO0
}

# Function to get prediction at s0 location from lavaan model. Based on Eq. (7).
get.pred <- function (MLIST = NULL, covar = NULL){
  var.names <- c("CEC.Ar","CEC.Br","CEC.Cr","OC.Ar","OC.Br","OC.Cr",
                 "Clay.Ar","Clay.Br","Clay.Cr","dem","vdchn","X",
                 "lstm","evisd","ndwi.b","twi","ndwi.a")
  m <- MLIST
  A <- m$beta[1:9,10:p] # this is Gamma matrix
  B <- m$beta[1:9,1:9]  # this is Beta matrix
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  k <- covar[,var.names[10:p]] # values of covariates at s0
  pr <- as.vector(as.matrix(k)) 
  pred <- t(IB.inv %*% A %*% pr) # this is Eq. (7)
  colnames(pred) <- var.names[1:9]
  pred 
}

# Function to get residuals. 
get.res <- function (m = NULL, z = NULL){
  A <- m$beta[1:9,10:p]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  sp <- z[,1:9]
  pr <- z[,10:p]
  res <- matrix(data = NA, nrow = N, ncol = 9)
  for(i in seq_along(pr[,1])){
    res[i,] <- t(sp[i,] - (IB.inv %*% A %*% pr[i,]))
  }
  colnames(res) <- colnames(sp)
  res
}

# Here I will extract only one location to make the prediction,
# but below I leave the code to predict all locations 
s01 <- s0[1,]

covar.st <- s0[1,c("dem","vdchn","X", "lstm","evisd","ndwi.b","twi","ndwi.a")]
# prediction based on lavaan model:
pred.lm <- get.pred(MLIST = MLIST.obs, covar = covar.st) 
# Now we have to add the kriging effect:
# Locations of calibration profiles (xy) and predition (ll)
ll <- s01[c("X","Y")] # location of s0
ll[,"Y"] <- ll[,"Y"] * ST["Y","std.dev"] / ST["X","std.dev"] # same factor for X and Y
xy <- as.data.frame(xy) # xy are the locations of calibration data
xy.ll <- rbind(xy,ll) # add the s0 location to the list of calibration locations
coordinates(xy.ll) <- ~X+Y # all coordinates (N+1 = 148)
h.all <- sp::spDists(xy.ll) # h.all is a distance matrix including s0
h0 <- matrix(h.all[1:N,N+1], ncol = 1, nrow = N) # h0 is a matrix of N x 1 
RHO0 <- get.RHO0(MLIST.obs, h0 = h0) # RHO0 has the same dimensions than h0

# obtain residuals for kriging
res <- get.res(m = MLIST.obs, z = as.matrix(data))
y.all <- as.vector(res) # vector of residuals

# Compute matrix SIGMA.yy (Eq. 16) (y=U, x=V)
SIGMA.yy <- computeSigmaHat.LISREL(MLIST = MLIST.obs)[1:9,1:9] # q x q
# Compute matrix SIGMA.xx (Eq. 17)
RHO <- get.RHO(MLIST = MLIST.obs, h = h)
SIGMA.xx <- kronecker(SIGMA.yy, RHO)   # qN x qN
# Compute SIGMA.xy (Eq. 18)
SIGMA.xy <- kronecker(SIGMA.yy, RHO0) # q x qN
# kriging of residuals
k.res <- t(crossprod(SIGMA.xy, chol2inv(chol(SIGMA.xx))) %*% y.all) 
colnames(k.res) <- paste0(colnames(res),".res")
# total prediction
predicted <- pred.lm + k.res # add kriging of residuals to lavaan prediction
predicted <- cbind(predicted, ll) # include coordinates to prediction
predicted
# Note that we can also estimate kriging variance

# Computing prediction error variance (Eq. 15)
PSI <- MLIST.obs$psi[1:9,1:9] # system error of SP
B <- MLIST.obs$beta[1:9,1:9]
I <- diag(nrow = 9, ncol = 9)
IB.inv <- solve(I - B)
theta <- MLIST.obs$theta[1:9,1:9] # measurement error
var.SIGMA.yy <- PSI + IB.inv %*% tcrossprod(theta, IB.inv) #
var.SIGMA.xx <- kronecker(var.SIGMA.yy, RHO)
var.SIGMA.xy <- kronecker(var.SIGMA.yy, RHO0)
# C_{UU} − C_{UV} C_{VV}^{−1} C_{VU} (Eq. 15)
var.zeta <- var.SIGMA.yy - 
  crossprod(var.SIGMA.xy,
            chol2inv(chol(var.SIGMA.xx))) %*% var.SIGMA.xy
# Vaiances for each soil property at each location can be obtain by:
variance <- matrix(diag(var.zeta), nrow = 1) 
colnames(variance) <- paste0(colnames(res),".var")
variance
# However, the off-diagonal contain covariances that might be interesting to 
# be analysed spatially



#### Loop to predict several locations usin doParallel ####
# Not run. This may take several hours
library(doParallel)
XX <- number of core processors
registerDoParallel(cores = XX) # with 48 cores it took 6.5 hours
result <-
  foreach(i = icount(nrow(s0)), .combine = rbind, # s0 = all prediction locations
          .packages = "sp") %dopar% {
            # covar for prediction
            covar.st <- s0[i,c("dem","vdchn","X", "lstm",
                               "evisd","ndwi.b","twi","ndwi.a")]
            # prediction from linear model and residuals
            pred.lm <- get.pred(MLIST = MLIST.obs, covar = covar.st) 
            res <- get.res(m = MLIST.obs, z = s)
            y.all <- as.vector(res) 
            ll <- s0[i,c("X","Y")]
            ll[,"Y"] <- ll[,"Y"] * ST["Y","std.dev"] / ST["X","std.dev"]
            xy <- as.data.frame(xy)
            xy.ll <- rbind(xy,ll)
            coordinates(xy.ll) <- ~X+Y
            h.all <- sp::spDists(xy.ll) 
            h0 <- matrix(h.all[1:N,N+1], ncol = 1, nrow = N)
            RHO0 <- get.RHO0(MLIST.obs, h0 = h0) 
            # get SIGMA.xy
            SIGMA.xy <- kronecker(SIGMA.yy, RHO0)
            # kriging of residuals
            k.res <- t(crossprod(SIGMA.xy, chol2inv(chol(SIGMA.xx))) %*% y.all) 
            colnames(k.res) <- paste0(colnames(res),".res")
            # total prediction
            predicted <- pred.lm + k.res
            predicted <- cbind(predicted, ll)
            # Computing prediction variance (could be in a function get.var())
            # PSI <- MLIST.obs$psi[1:9,1:9] # system error of SP
            # B <- MLIST.obs$beta[1:9,1:9] 
            # I <- diag(nrow = 9, ncol = 9)
            # IB.inv <- solve(I - B)
            # theta <- MLIST.obs$theta[1:9,1:9] # measurement error
            # var.SIGMA.yy <- PSI + IB.inv %*% tcrossprod(theta, IB.inv) #
            # var.SIGMA.xx <- kronecker(var.SIGMA.yy, RHO)
            # var.SIGMA.xy <- kronecker(var.SIGMA.yy, RHO0)
            # var.zeta <- var.SIGMA.yy -
            #   crossprod(var.SIGMA.xy,
            #             chol2inv(chol(var.SIGMA.xx))) %*% var.SIGMA.xy
            # variance <- matrix(diag(var.zeta), nrow = 1)
            # colnames(variance) <- paste0(colnames(res),".var")
            # #as.vector(var.zeta)
            result <- predicted #cbind(predicted, k.res, variance)
            result
          }
doParallel::stopImplicitCluster()
# end(not run)
