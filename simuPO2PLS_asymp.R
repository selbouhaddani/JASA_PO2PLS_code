#### ORGANIZATION: first all functions are defined, then all parameters and options are set
#### After approx line 138 the options are defined and can be changed
#### TIME: approx 15 minutes
#### HOW TO run this script:
#### Adjust parameters below related to speed, and source the file
#### This generates SIMUasymp***.RData files
#### Go to simuPO2PLS_asymp_aggr.R to generate the plots
#### NOTES: 
#### Some lines of code are modified for the sake of speed
#### Original lines of code are indicated with *original*
#### The printed output is the rejection proportion for 
####   asymptotic SE, parametric, non-parametric bootstrapped SE, and permutations
####   in the first two scenarios, the bootstrap and permutations are "switched off" 
####   and render non-interpretable output

#### Choose the number of replicates and cores in this script
nr_cores = 4
nr_replic = 16
nr_permrep = 10

library(parallel)
library(tidyverse)
library(magrittr)
library(OmicsPLS)
library(PO2PLS)

#### Non parametric bootstrapping
bootstrap.po2m <- function(fit, X, Y, perm.cores = 1, perm.K = 80, perm.reduce = FALSE, ...){
  if(perm.K == 0) return(list(fit$par,fit$par))
  if(perm.reduce) {
    perm.keepx <- with(fit$par, rowSums((W %*% B)^2) %>% order(decreasing = TRUE) %>% tail(19))
    X <- X[,which(1:ncol(X) %in% perm.keepx)]
  }
  r <- ncol(fit$par$W)
  rx <- ncol(fit$par$Wo)*sign(ssq(fit$par$Wo))
  ry <- ncol(fit$par$Co)*sign(ssq(fit$par$Co))

  boot.par <- mclapply(mc.cores = perm.cores, 1:perm.K,
                       FUN = function(niks) {
                         perm.indx <- sample(nrow(X), replace = TRUE)
                         fit_perm <- PO2PLS((X[perm.indx, ]), (Y[perm.indx, ]), r, rx, ry, ...)
                         return(fit_perm$parameters)
                       })
  return(boot.par)
}

#### Permutations
perm.po2m <- function(fit, X, Y, perm.cores = 1, perm.K = 80, perm.reduce = FALSE, ...){
  if(perm.K == 0) return(list(fit$par,fit$par))
  if(perm.reduce) {
    perm.keepx <- with(fit$par, rowSums((W %*% B)^2) %>% order(decreasing = TRUE) %>% tail(19))
    X <- X[,which(1:ncol(X) %in% perm.keepx)]
  }
  r <- ncol(fit$par$W)
  rx <- ncol(fit$par$Wo)*sign(ssq(fit$par$Wo))
  ry <- ncol(fit$par$Co)*sign(ssq(fit$par$Co))

  boot.par <- mclapply(mc.cores = perm.cores, 1:perm.K,
                       FUN = function(niks) {
                         perm.indx <- sample(nrow(X), replace = FALSE)
                         fit_perm <- PO2PLS((X), (Y[perm.indx, ]), r, rx, ry, ...)
                         return(fit_perm$parameters)
                       })
  return(boot.par)
}

#### Parametric bootstrapping
param_boot.po2m <- function(fit, X, Y, perm.cores = 1, perm.K = 80, perm.reduce = FALSE, ...){
  if(perm.K == 0) return(list(fit$par,fit$par))
  if(perm.reduce) {
    perm.keepx <- with(fit$par, rowSums((W %*% B)^2) %>% order(decreasing = TRUE) %>% tail(19))
    X <- X[,which(1:ncol(X) %in% perm.keepx)]
  }
  r <- ncol(fit$par$W)
  rx <- ncol(fit$par$Wo)*sign(ssq(fit$par$Wo))
  ry <- ncol(fit$par$Co)*sign(ssq(fit$par$Co))

  boot.par <- mclapply(mc.cores = perm.cores, 1:perm.K,
                       FUN = function(niks) {
                         dat <- generate_data(nrow(Y), fit$par)
                         fit_perm <- with(dat, PO2PLS(X, Y, r, rx, ry, ...))
                         return(fit_perm$parameters)
                       })
  return(boot.par)
}

#### Replication function
f <- function(boot.K=0, N=1e2, p=20, q=5, r=2, rx=1, ry=0, alpha=0.1){
  dat <- generate_data(N, parms)
  X <- dat$X
  Y <- dat$Y
  
  fit <- suppressMessages(invisible(PO2PLS(X, Y, r, rx, ry, steps=1e3, tol = 1e-7, init_param=parms, verbose=F)))
  
  asymp.se.B <- variances_inner.po2m(fit, X, Y)
  boot.fits <- suppressMessages((param_boot.po2m(fit, X, Y, perm.K = round(boot.K/2), steps = 100, tol=0.01, verbose=F)))
  if(r == 1){
  boot.se.B <- sapply(boot.fits, `[[`, "B") %>% sd
  } else {
    attr(boot.fits,"whichmax") <- sapply(boot.fits, function(e) crossprod(parms$W[,1],e$W)) %>% apply(2,which.max)
    attr(boot.fits,"ord") <- sapply(boot.fits, function(e) crossprod(parms$W[,1],e$W)) %>% apply(2,order, decreasing = TRUE)
    boot.se.B <- sapply(1:length(boot.fits), function(ii) with(boot.fits[[ii]],diag(B)[attr(boot.fits,"ord")[,ii]])) %>% apply(1,sd)
  }
  
  boot.non.fits <- suppressMessages((bootstrap.po2m(fit, X, Y, perm.K = round(boot.K/2), steps = 100, tol=0.01, verbose=F)))
  if(r==1){
    boot.non.se.B <- sapply(boot.non.fits, `[[`, "B") %>% sd
  } else {
    attr(boot.non.fits,"whichmax") <- sapply(boot.non.fits, function(e) crossprod(parms$W[,1],e$W)) %>% apply(2,which.max)
    attr(boot.non.fits,"ord") <- sapply(boot.non.fits, function(e) crossprod(parms$W[,1],e$W)) %>% apply(2,order, decreasing = TRUE)
    boot.non.se.B <- sapply(1:length(boot.non.fits), function(ii) with(boot.non.fits[[ii]],diag(B)[attr(boot.non.fits,"ord")[,ii]])) %>% apply(1,sd)
  }
    
  perm.fits <- suppressMessages(perm.po2m(fit, X, Y, perm.K = boot.K, steps = 100, tol=0.01, verbose=F))
  perm.B <- sapply(perm.fits, function(e) diag(e$B))
  if(r==1) perm.B %<>% t
  perm.pval <- mean(c(fit$par$B[1] <= perm.B[1,],TRUE))
  
  #boot.fits.se.W <- sapply(1:length(boot.fits), function(ii) with(boot.fits[[ii]],(W%*%B))[,attr(boot.fits,"whichmax")[ii]]) %>% apply(1,sd)
  #boot.fits.se.WBC <- sapply(1:length(boot.fits), function(ii) with(boot.fits[[ii]],(W%*%B%*%t(C)))) %>% apply(1,sd)
  
  #asymp.se.W <- NA
  #asymp.se.W <- variances.PO2PLS(fit, Reduce(cbind,dat), "complete")$Iobs %>% solve %>% diag %>% `[`(1:p) 
  #asymp.NA <- asymp.se.W < 0
  #asymp.se.W %<>% abs 
  #asymp.se.W <- (asymp.se.W %*% diag(asymp.se.B^2,r)) %>% sqrt
  
  #Ws.best <- fit$par$W[,which.max(crossprod(parms$W[,1],fit$par$W))]
  #Cs.best <- fit$par$C[,which.max(crossprod(parms$W[,1],fit$par$W))]  
  Bs.best <- diag(fit$par$B)[order(-crossprod(parms$W[,1],fit$par$W)^2)]
  if(is.null(fit$par)) warning("fit-par was null")
  return(list(boot.fits.B = sapply(boot.fits,function(e) {diag(e$B)}), se.B = list(asymp = asymp.se.B, boot = boot.se.B, nonboot = boot.non.se.B), 
    Bs = Bs.best, Bs.true = diag(parms$B), samp = N, dims = p, pwr = pwr, 
    pvals.B = list(asymp=pchisq(with(fit$par, (c(Bs.best[1]) / asymp.se.B[1] )^2), df=1, lower.tail=FALSE),
                    boot=pchisq(with(fit$par, (c(Bs.best[1]) / boot.se.B[1] )^2), df=1,  lower.tail=FALSE),
                    nonboot=pchisq(with(fit$par, (c(Bs.best[1]) / boot.non.se.B[1] )^2), df=1,  lower.tail=FALSE),
                    perm=perm.pval)))
}

tic = Sys.time()

r=2; rx=1; ry=0; 
############ Sample Size ########
cat("\n  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= \n")
cat("\n Increasing Sample Size\n\n")
for(samp_size in c(50, 500, 5000, 10000)[1:3]){
dim_p = 50
dim_q = 50
cat("Now at sample size ", samp_size, ", ... \n", sep="")
parms0 <- generate_params(dim_p, dim_q, r=r, rx=rx, ry=ry, alpha=0.5) 
parms <- parms0
parms$SigH <- with(parms0, SigH + SigT%*%B^2)
for(pwr in c(0, 0.1, 0.2, 0.5, 1)[1]){
  parms$B <- diag(pwr,ncol(parms$B))
  parms$SigH <- with(parms, SigH - SigT%*%B^2)
  cat("\n")
  # *original* # time_run <- system.time(outp.H0 <- mclapply(mc.cores=50,1:5e4, function(niks) {if(niks %% 10 == 0) cat(niks," "); f(r=r,rx=rx,ry=ry,q=dim_q,boot.K = 0, N = samp_size, p = dim_p)}))
  time_run <- system.time(outp.H0 <- mclapply(mc.cores=nr_cores, 1:nr_replic, function(niks) {if(niks %% 10 == 0) cat(niks," "); f(r=r,rx=rx,ry=ry,q=dim_q,boot.K = 0, N = samp_size, p = dim_p)}))
  cat("\n")
  save(outp.H0, parms, file = paste0("SIMUasympH0_N",samp_size,".RData"))
  print(time_run)
  cat("\n\n ----> Please only interpret the asymp column \n")
  sapply(outp.H0, `[[`, "pvals.B") %>% `<=`(0.05) %>% apply(1,mean) %>% print
}}


############ Dimensions ########
cat("\n  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= \n")
cat("\n Increasing Dimensions\n\n")
for(dim_p in c(20, 200, 2000)){
samp_size = 500
cat("Now at dimension ", dim_p, ", ... \n", sep="")
parms0 <- generate_params(dim_p, 5, r=r, rx=rx, ry=ry, alpha=0.5) 
parms <- parms0
parms$SigH <- with(parms, SigH + SigT%*%B^2)
for(pwr in c(0, 0.1, 0.2, 0.5, 1)[1]){
  parms$B <- diag(pwr,ncol(parms$B))
  parms$SigH <- with(parms, SigH - SigT%*%B^2)
  cat("\n")
  # *original* # time_run <- system.time(outp.H0 <- mclapply(mc.cores = 50, 1:2e3, function(niks) {if(niks %% 10 == 0) cat(niks," "); f(boot.K = 0, N = samp_size, p = dim_p)}))
  time_run <- system.time(outp.H0 <- mclapply(mc.cores=nr_cores, 1:nr_replic, function(niks) {if(niks %% 10 == 0) cat(niks," "); f(r=r,rx=rx,ry=ry,q=dim_q,boot.K = 0, N = samp_size, p = dim_p)}))
  cat("\n")
  save(outp.H0, parms, file = paste0("SIMUasympH0_p",dim_p,".RData"))
  print(time_run)
  cat("\n\n ----> Please only interpret the asymp column \n")
  sapply(outp.H0, `[[`, "pvals.B") %>% `<=`(0.05) %>% apply(1,mean) %>% print
}}


############ Power ########
cat("\n  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= \n")
cat("\n Increasing Power\n\n")
for(samp_size in c(50, 500)){
  dim_p = ifelse(samp_size==50, 200, 20)
  cat("Now at dimension ", dim_p, ", ... \n", sep="")
  cat("   Now at sample size ", samp_size, ", ... \n", sep="")
  parms0 <- generate_params(dim_p, 5, r=r, rx=rx, ry=ry, alpha=0.5) 
  parms <- parms0
  parms$SigH <- with(parms, SigH + SigT%*%B^2)
  for(pwr in c(0, 0.2, 0.5, 1)){
    cat("    Now at pwr ", pwr, ", ... \n", sep="")
    parms$B <- diag(pwr,ncol(parms$B))
    parms$SigH <- with(parms, SigH - SigT%*%B^2)
    cat("\n")
    # *original* # time_run <- system.time(outp.H0 <- mclapply(mc.cores = 50, 1:ifelse(pwr==0,500,500), function(niks) {if(niks %% 10 == 0) cat(niks," "); f(boot.K = 500, N = samp_size, p = dim_p)}))
    time_run <- system.time(outp.H0 <- mclapply(mc.cores = nr_cores, 1:nr_replic, function(niks) {if(niks %% 10 == 0) cat(niks," "); f(boot.K = nr_permrep, N = samp_size, p = dim_p)}))
    cat("\n")
    save(outp.H0, parms, file = paste0("SIMUasympAll_N",samp_size,"_p",dim_p,"_pwr",pwr,"_.RData"))
    print(time_run)
    sapply(outp.H0, `[[`, "pvals.B") %>% `<=`(0.05) %>% apply(1,mean) %>% print
}}

cat("\n\n  success: end \n\n")
print(Sys.time() - tic)

