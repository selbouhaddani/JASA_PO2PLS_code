cat("START")

#### ORGANIZATION: first all functions are defined, then all parameters and options are set
#### After approx line 300 the options are defined and can be changed
#### TIME: about 15 - 30 minutes per study 
#### HOW TO run this script:
#### Select one of four studies (normal, rank, distribution, or case-control outcome)
#### Set doOutcome, doDistr, etc. to FALSE or TRUE
#### Normal study: all do*** are FALSE
#### Choose the number of cores used by mclapply, see 'outp <- mclapply(...)' at the end
#### Also choose the number of replications
#### NOTE: the current #replications and cores are modified to obtain a reasonable runtime on a modest laptop
#### NOTE: the number of EM steps and test samples are lower than in the original simulation (see nr_steps and nr_testdata below)
#### Source the script per study, which generates SIMU***.RData files
#### Go to simuPO2PLS_aggr.R to generate the plots

## R/3.5.1
tic <- proc.time()

## Parameters given by the SLURM computing cluster
#arrID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

library(parallel)
library(tidyverse)
library(magrittr)
library(OmicsPLS)
library(PO2PLS)
library(pryr)

## List of distributions to consider
distr_list <- list(
  normal = rnorm, 
  student_t = function(e) {rt(e, df=2)}, 
  poiss = function(e) {rpois(e, lambda=1)},
  binom = function(e) {rbinom(e, size = 2, prob = 0.25)} 
)

#### Function, Convert O2PLS output to PO2PLS output
o2m_2_po2m <- function(fit){
  tmpsvd <- svd(fit$P_Yosc.,nu=0)
  To <- with(fit, T_Yosc %*% tmpsvd$v %*% diag(tmpsvd$d,ncol(tmpsvd$v)) %*% t(tmpsvd$v))
  tmpsvd <- svd(fit$P_Xosc.,nu=0)
  Uo <- with(fit, U_Xosc %*% tmpsvd$v %*% diag(tmpsvd$d,ncol(tmpsvd$v)) %*% t(tmpsvd$v))
  parameters <- list(
    W = fit$W.,
    Wo= suppressWarnings(orth(fit$P_Yosc.)),
    C = fit$C.,
    Co= suppressWarnings(orth(fit$P_Xosc.)),
    B = abs(fit$B_T. * diag(1,ncol(fit$B_T.))),
    SigT = cov(fit$Tt)* diag(1,ncol(fit$B_T.)),
    SigTo= cov(To)* diag(1,ncol(fit$P_Yosc.)),
    SigUo= cov(Uo)* diag(1,ncol(fit$P_Xosc.)),
    SigH = cov(fit$H_UT)* diag(1,ncol(fit$B_T.)) +0.1,
    sig2E= max(ssq(fit$E)/prod(dim(fit$E)),0.25),
    sig2F= max(0.25, ssq(fit$Ff)/prod(dim(fit$Ff)))
  )
  latent_vars <- list(
    mu_T = tibble(mu_T=fit$Tt),
    mu_U = tibble(mu_U=fit$U),
    mu_To= tibble(mu_To=To),
    mu_Uo= tibble(mu_Uo=Uo)
  )
  logl <- list(all_vals=NA,last_val=NA,df=NA)
  class(logl) <- "loglik.po2m"
  meta_data <- list(
    loglikelihood = logl,
    explained_vars = with(fit, c(R2X=R2X, R2Y=R2Y, R2Xj=R2Xcorr, R2Yj=R2Ycorr, R2Xhat=R2Xhat, R2Yhat=R2Yhat)),
    comps = with(fit, c(r=flags$n,rx=flags$nx,ry=flags$ny)),
    time = fit$flags$time.elapsed,
    call = fit$flags$call,
    convergence = TRUE, 
    ssqs = c(fit$flags$ssqX, fit$flags$ssqY)
  )
  outp <- list(parameters = parameters, latent_vars = latent_vars, meta_data=meta_data)
  class(outp) <- "po2m"
  return(outp)
}

#### Function, Generate parameters 
gen_par <- function(X, Y, r, rx, ry, alpha = 0.1, type = c("random", "o2m","unit"), homogen = FALSE){
  p = ifelse(is.matrix(X), ncol(X), X)
  q = ifelse(is.matrix(Y), ncol(Y), Y)
  
  if (length(alpha) == 1)
    alpha <- rep(alpha, 3)
  if (!(length(alpha) %in% c(1, 3)))
    stop("length alpha should be 1 or 3")
  outp <- list(
    W = orth(matrix(rnorm(p * r), p, r) + 1),
    Wo = suppressWarnings(sign(rx) * orth(matrix(rnorm(p * max(1, rx)), p, max(1, rx)) + seq(-p/2, p/2, length.out = p))), 
    C = orth(matrix(rnorm(q * r), q, r) + 1), 
    Co = suppressWarnings(sign(ry) * orth(matrix(rnorm(q * max(1, rx)), q, max(1, ry)) + seq(-q/2, q/2, length.out = q))), 
    B = B <- diag(1, r), #diag(sort(runif(r, 1, 4), decreasing = TRUE), r), 
    SigT = diag(sort(runif(r, 1, 2), decreasing = TRUE), r), 
    SigTo = sign(rx) * diag(sort(runif(max(1, rx), 1, 3), decreasing = TRUE), max(1, rx)), 
    SigUo = sign(ry) * diag(sort(runif(max(1, ry), 1, 3), decreasing = TRUE), max(1, ry)),
    beta = 3*(r:1)
  )
  if(homogen) {
    #outp$B <- diag(1, r)
    outp$SigH = 1e-4 #diag(alpha[3]/(1 - alpha[3]) * (mean(diag(outp$SigT %*% outp$B))), r)
    outp$SigTo <- alpha[3]/(1 - alpha[3]) * outp$SigTo / tr(outp$SigTo)*tr(outp$SigT)
    outp$SigUo / tr(outp$SigUo)*tr(outp$SigH)
  }
  outp$SigH = diag(alpha[3]/(1 - alpha[3]) * (mean(diag(outp$SigT %*% outp$B))), r)
  outp$SigU <- with(outp, SigT %*% B + SigH)
  return(with(outp, {
    c(outp, 
      sig2E = alpha[1]/(1 - alpha[1]) * (mean(diag(SigT)) + mean(diag(SigTo)))/p, 
      sig2F = alpha[2]/(1 - alpha[2]) * (mean(diag(SigT %*% B^2 + SigH)) + mean(diag(SigUo)))/q)
  }))
}

#### Function Generate data 
#### NOT FOR THE selected outcome design
gen_dat <- function (N, params, distr = rnorm)
{
  N2 <- N
  W = params$W
  C = params$C
  Wo = params$Wo
  Co = params$Co
  B = params$B
  beta = params$beta
  SigT = params$SigT
  SigTo = params$SigTo + 1e-06 * SigT[1] * (params$SigTo[1] == 0)
  SigH = params$SigH
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigU = SigT %*% B^2 + SigH
  SigUo = params$SigUo + 1e-06 * SigU[1] * (params$SigUo[1] == 0)
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  Gamma = rbind(cbind(W, matrix(0, p, r), Wo, matrix(0, p,
    ry)), cbind(matrix(0, q, r), C, matrix(0, q, rx), Co))
  VarZ = blockm(blockm(blockm(SigT, SigT %*% B, SigU), matrix(0,
    2 * r, rx), SigTo), matrix(0, 2 * r + rx, ry), SigUo)
  Z <- scale(matrix(distr(N2 * (2 * r + rx + ry)), N2))
  Z <- Z %*% chol(VarZ)
  Z[, 2 * r + 1:rx] <- sign(ssq(Wo)) * Z[, 2 * r + 1:rx]
  Z[, 2 * r + rx + 1:ry] <- sign(ssq(Co)) * Z[, 2 * r + rx +
                                                1:ry]
  EF <- scale(matrix(distr(N2 * p), N2)) * sqrt(sig2E)
  EF2 <- scale(matrix(distr(N2 * q), N2)) * sqrt(sig2F)
  EF <- cbind(EF, EF2)
  rm(EF2)
  dat <- Z %*% t(Gamma)
  dat <- dat + EF
  rm(EF)
  outc <- Z[,1:(2*r)] %*% c(beta,beta) #apply(Z[,1:(2*r)], 2, function(e) ifelse(e<0, mean(e[e<0]), mean(e[e>0]))) %*% c(beta,beta)
  outc <- outc + rnorm(length(outc), sd=.5*sd(outc))
  ii_keep <- order(-outc^2)[1:N] 
  outc <- scale(outc[ii_keep] > 0)
  X <- scale(dat[ii_keep, 1:p],scale=F)
  Y <- scale(dat[ii_keep, -(1:p)],scale=F)
  rm(dat)
  return(list(X = X, Y = Y, outc = outc))
}

#### Function, Generate data 
#### ONLY FOR THE selected outcome design
gen_dat2 <- function (N, params, distr = rnorm)
{
  N2 <- 2*N
  W = params$W
  C = params$C
  Wo = params$Wo
  Co = params$Co
  B = params$B
  beta = params$beta
  SigT = params$SigT
  SigTo = params$SigTo + 1e-06 * SigT[1] * (params$SigTo[1] == 0)
  SigH = params$SigH
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigU = SigT %*% B^2 + SigH
  SigUo = params$SigUo + 1e-06 * SigU[1] * (params$SigUo[1] == 0)
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  Gamma = rbind(cbind(W, matrix(0, p, r), Wo, matrix(0, p,
    ry)), cbind(matrix(0, q, r), C, matrix(0, q, rx), Co))
  VarZ = blockm(blockm(blockm(SigT, SigT %*% B, SigU), matrix(0,
    2 * r, rx), SigTo), matrix(0, 2 * r + rx, ry), SigUo)
  Z <- scale(matrix(distr(N2 * (2 * r + rx + ry)), N2))
  Z <- Z %*% chol(VarZ)
  Z[, 2 * r + 1:rx] <- sign(ssq(Wo)) * Z[, 2 * r + 1:rx]
  Z[, 2 * r + rx + 1:ry] <- sign(ssq(Co)) * Z[, 2 * r + rx + 1:ry]
  EF <- cbind(scale(matrix(distr(N2 * p), N2)) * sqrt(sig2E),
              scale(matrix(distr(N2 * q), N2)) * sqrt(sig2F))
  dat <- Z %*% t(Gamma) + EF
  outc <- Z[,1:(2*r)] %*% c(beta,beta) #apply(Z[,1:(2*r)], 2, function(e) ifelse(e<0, mean(e[e<0]), mean(e[e>0]))) %*% c(beta,beta)
  outc <- outc + rnorm(length(outc), sd=.5*sd(outc))
  ii_keep <- order(-outc^2)[1:N] 
  outc <- scale(outc[ii_keep] > 0)
  return(list(X = scale(dat[ii_keep, 1:p],scale=F), Y = scale(dat[ii_keep, -(1:p)],scale=F), outc = outc))
}

#### Replication function
f <- function(parms, N=1e2, distr_i = rnorm){
  #### Retrieve # components
  r <- ncol(parms$W)
  rx <- ncol(parms$Wo)*sign(ssq(parms$Wo))
  ry <- ncol(parms$Co)*sign(ssq(parms$Co))
  
  ## IN CASE OF *RANK MISSPECIFICATION* subtract one component
  if(doRankMisspec){  
    message("RANK MISSPECIFICATION TRUE")
    r <- ncol(parms$W) - 1
    rx <- (ncol(parms$Wo)-1)*sign(ssq(parms$Wo))
    ry <- (ncol(parms$Co)-1)*sign(ssq(parms$Co))
  }
  
  #### Generate data
  if(!doOutcome) dat <- gen_dat(N, parms, distr=distr_i) #### NOT FOR OUTCOME
  if(doOutcome) dat <- gen_dat2(N, parms) #### *ONLY FOR OUTCOME*
  X <- scale(dat$X, scale=F)
  Y <- scale(dat$Y, scale=F)
  
  #### Retrieve dimensions
  p <- ncol(X)
  q <- ncol(Y)
  
  #### Fit PLS, O2PLS from OmicsPLS package, PO2PLS and PPLS from PO2PLS package
  fit_pls <- suppressWarnings(invisible(o2m(X, Y, r, 0, 0) %>% o2m_2_po2m)); print("PLS")
  fit_o2pls <- suppressWarnings(invisible(o2m(X, Y, r, rx, ry) %>% o2m_2_po2m)); print("O2PLS")
  fit_po2pls <- suppressMessages(invisible(PO2PLS(X, Y, r, rx, ry, steps=nr_steps, tol = 1e-6, verbose=F,init_param=fit_o2pls$par))); print("PO2PLS")
  fit_ppls <- suppressMessages(invisible(PO2PLS(X, Y, r, 0, 0, steps=nr_steps, tol = 1e-6, verbose=F,init_param=fit_pls$par))); print("PPLS")
  
  ## "Disable" fit_sifa by replicating another fit object
  fit_sifa <- fit_po2pls
  
  #### PLEASE NOTE that the original SIFA is matlab code and requires Octave and the RcppOctave package
  #### The simulations are based on this matlab code, but was very slow
  #### For those that don't have Octave or Matlab, we added the u=t restriction to PO2PLS, giving an alternative SIFA fit based on our faster implementation
  #### To use this, uncomment the following line, and comment out the block of code after it
  # fit_sifa <- suppressMessages(invisible(PO2PLS(X, Y, r, rx, ry, steps=1e3, tol = 1e-6, verbose=F, homogen_joint = TRUE,init_param=fit_o2pls$par))); print("SIFA")
  # time_SIFA <- system.time(print(pryr::mem_change({
  #   o_load(list(X=X,Y=Y))
  #   o_eval("addpath('~/selbouhaddani/said/PO2PLS/SIFA_Matlab')")
  #   o_eval(paste0("[~,~,Vj,Vo, se2, Sf0, Sf]=SIFA_A(ones(size(X,1),1), {X, Y}, ",r,", [",rx,",",ry,"], struct('sparsity',0, 'Niter', ",90,"));"))
  #   fit_sifa <- o_get("Vj","Vo","se2","Sf0","Sf")
  # })))
  # fit_sifa$par$W <- orth(fit_sifa$Vj[1:p,])
  # fit_sifa$par$C <- orth(fit_sifa$Vj[-(1:p),])
  # fit_sifa$par$Wo <- orth(fit_sifa$Vo[[1]])
  # fit_sifa$par$Co <- orth(fit_s$Vo[[2]])
  # fit_sifa$par$SigT <- fit_sifa$Sf0
  # fit_sifa$par$SigU <- fit_sifa$Sf0
  # fit_sifa$par$SigTo <- as.matrix(fit_sifa$Sf[[1]])
  # fit_sifa$par$SigUo <- as.matrix(fit_sifa$Sf[[2]])
  # fit_sifa$par$SigH <- 0*fit_sifa$par$SigU
  # fit_sifa$par$B <- diag(1,r)
  # fit_s$par$sig2E <- fit_s$se2[1]
  # fit_s$par$sig2F <- fit_s$se2[2]

  #### Gather the fits and add true 'fit'
  fits <- list(PO2PLS=fit_po2pls, PPLS=fit_ppls, PLS=fit_pls, O2PLS=fit_o2pls, SIFA=fit_sifa)
  fits2 <- fits
  fits2$TRuE <- list(parameters = parms)
  
  #### Generate test data with N=1e4
  #### Then calculate training/test prediction error
  print("pred")
  dat.tst <- gen_dat(ifelse(doOutcome, nr_testdata/2, nr_testdata), parms)
  Xtst <- dat.tst$X %>% scale(scale=F)
  Ytst <- dat.tst$Y %>% scale(scale=F)
  preds <- sapply(simplify=F, fits2, function(e) ssq(with(e$par, Ytst - (Xtst-(Xtst%*%Wo)%*%t(Wo))%*%W%*%B%*%t(C) ))/ssq(Ytst) )
  preds2 <- sapply(simplify=F, fits2, function(e) ssq(with(e$par, Y - (X-(X%*%Wo)%*%t(Wo))%*%W%*%B%*%t(C) ))/ssq(Y) )
  
  #### Calculate TPR
  ordr <- sapply(fits, function(ee) crossprod(ee$par$W, parms$W) %>% abs %>% apply(2, which.max))
  estim.tpr <- lapply(fits, function(ee) apply(ee$par$W, 2, function(eee) order(-eee^2)[1:round(nrow(parms$W)/4)] ))
  estim.tpr <- sapply(estim.tpr, function(ee) sapply(1:ncol(fit_po2pls$par$W), function(ii) mean(true.tpr[,ordr[ii]] %in% ee[,ii]) ))
  print("tpr")
  
  #### Other output not further considered
  estim.cp <- sapply(fits, function(ee) crossprod(ee$par$W[,ordr], parms$W) %>% abs %>% diag)
  timing <- sapply(fits, function(e) e$meta$time)
  
  #### Calculate outcome training/test error
  outc.pred <- sapply(fits2, 
    function(e) ssq(with(e$par, dat.tst$outc - cbind(1, (Xtst-Xtst%*%Wo%*%t(Wo)) %*% W) %*% coef(lm(dat.tst$outc ~ (Xtst-Xtst%*%Wo%*%t(Wo)) %*% W)) ))/ssq(dat.tst$outc)
  )
  print("outc")
  outc.trainpred <- sapply(fits2, 
    function(e) ssq(with(e$par, dat$outc - cbind(1, (X-X%*%Wo%*%t(Wo)) %*% W) %*% coef(lm(dat$outc ~ (X-X%*%Wo%*%t(Wo)) %*% W)) ))/ssq(dat$outc)
  )
  
  return(list(Pred = preds, trainPred = preds2, Outc = outc.pred, trainOutc = outc.trainpred, TPR = estim.tpr, timing = timing, CP = estim.cp))
  
}

#### Switch between normal simulations, selected design with outcome and rank misspecification
#### Normal is both false, let at most one be true.
doOutcome = FALSE
doRankMisspec = FALSE
doDistr = FALSE
if(doOutcome + doRankMisspec + doDistr > 1) stop("Please choose max 1 scenario")

#### The desired scenarios, only for normal
#### The for loop iterates over all scenarios
normal_scenarios <- merge(
  merge(data.frame(SampleSize = c("low","high")), data.frame(Noise = c("low","high"))), 
  data.frame(Hetero = c("homogen", "heterogen"))
)

for(norm_scen in 1:nrow(normal_scenarios)){ 
  print(paste("Now at", norm_scen, "out of", nrow(normal_scenarios)))
SampleSize = list(high = 1000, low=100)
Dims = list(high = c(2e4, 125), low = c(2e3, 25))
dims = "low" # or high, recommended to leave it on low, or run using mc.cores = 2
Noise = list(high=c(0.95, 0.05, 0.8), low=c(0.4, 0.4, 0.8))
Homogen = list(homogen=TRUE, heterogen=FALSE)

samp = normal_scenarios$SampleSize[norm_scen] 
noise = normal_scenarios$Noise[norm_scen] 
hom = normal_scenarios$Hetero[norm_scen] 

#### For speeding up calculations
#### Original values are nr_steps = 1e3 and nr_testdata = 1e4
nr_steps = 1e2
nr_testdata = 1e3


#### Generate parameters based on scenarios described in Section 3
#### FOR Default
parms <- gen_par(Dims[[dims]][1], Dims[[dims]][2], 5, 5, 5, Noise[[noise]], homogen=Homogen[[hom]])

#### *ONLY FOR OUTCOME*
if(doOutcome) parms <- gen_par(2e4, 1e4, 2, 1, 1, c(0.2,0.4,0.8), homogen=FALSE)

#### True top 25% loadings per component 
true.tpr <- parms$W %>% apply(2, function(eee) order(-eee^2)[1:round(nrow(parms$W)/4)])

#### Ideally, this would be repeated using a computing cluster
#### We provide an alternative, using mclapply for replication on a laptop

#### Run f in normal case
if(!doOutcome && !doDistr) {
  print(paste("Working on", ifelse(doRankMisspec, "Rank misspec", "Default")))
  f_outp <- function(niks) f(parms, N=ifelse(samp == "low", 100, 1000))
  outp <- mclapply(mc.cores=4, 1:16, f_outp)
  save(outp, file=paste0("SIMU",ifelse(doRankMisspec,"rank",""),"_N",samp,"_p",dims,"_n",noise,"_",hom,".RData"))
}

#### Run f for sample size = 25
#### *ONLY FOR OUTCOME*
if(doOutcome) {
  print("Working on Outcome")
  f_outp <- function(niks) f(parms, N=25)
  outp <- mclapply(mc.cores=4, 1:16, f_outp)
  save(outp, file=paste0("SIMUoutcome.RData"))
  break
}

#### Run f for 4 distributions
#### ONLY FOR DISTRIBUTIONS scenario
if(doDistr){
  print("Working on Distribution misspec")
  f_outp <- function(niks) lapply(c(student_t=2, poiss=3, binom=4), 
                function(ii) {print(names(distr_list)[ii]); f(parms, N=ifelse(samp == "low", 100, 1000), distr_list[[ii]])})
  outp <- mclapply(mc.cores=4, 1:16, f_outp)
  save(outp, file=paste0("SIMUdistr_N",samp,"_p",dims,"_n",noise,"_",hom,".RData"))
}

} #######################

proc.time() - tic

gc()
mem_used()

cat("\nENDED : SUCCESS\n")
