
dir <- getwd()
## setwd("~/anclik/")
setwd("~/Documents/Code/R/anclik/")
source("anclik.R")
setwd(dir)
rm(dir)

source("R/epp.R")
source("R/generics.R")

source("R/IMIS.R")


#################
####  Prior  ####
#################

ldinvgamma <- function(x, alpha, beta){
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  return(log.density)
}

logiota.unif.prior <- c(log(1e-14), log(0.0025))
tau2.prior.rate <- 0.5

invGammaParameter <- 0.001   #Inverse gamma parameter for tau^2 prior for spline
ancbias.pr.mean <- 0.15
ancbias.pr.sd <- 1.0
muSS <- 1/11.5               #1/duration for r steady state prior

lprior <- function(theta, fp){

  nk <- fp$numKnots
  tau2 <- exp(theta[nk+3])

  return(sum(dnorm(theta[3:nk], 0, sqrt(tau2), log=TRUE)) +
         dunif(theta[nk+1], logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE) + 
         dnorm(theta[nk+2], ancbias.pr.mean, ancbias.pr.sd, log=TRUE) +
         ldinvgamma(tau2, invGammaParameter, invGammaParameter))
}

################################
####                        ####
####  HH survey likelihood  ####
####                        ####
################################

fnPrepareHHSLikData <- function(hhs, anchor.year = 1970L){
  hhs$W.hhs <- qnorm(hhs$prev)
  hhs$v.hhs <- 2*pi*exp(hhs$W.hhs^2)*hhs$se^2
  hhs$sd.W.hhs <- sqrt(hhs$v.hhs)
  hhs$idx <- hhs$year - (anchor.year - 1)

  hhslik.dat <- subset(hhs, used)
  return(hhslik.dat)
}

fnHHSll <- function(qM, hhslik.dat){
  return(sum(dnorm(hhslik.dat$W.hhs, qM[hhslik.dat$idx], hhslik.dat$sd.W.hhs, log=TRUE)))
}


###########################
####                   ####
####  Full likelihood  ####
####                   ####
###########################

fnCreateLikDat <- function(epp.data){

  likdat <- list(anclik.dat = fnPrepareANCLikelihoodData(epp.data$anc.prev, epp.data$anc.n),
                 hhslik.dat = fnPrepareHHSLikData(epp.data$hhs))
  likdat$lastdata.idx <- max(unlist(likdat$anclik.dat$anc.idx.lst), likdat$hhslik.dat$idx)
  return(likdat)
}

fnCreateParam <- function(theta, fp){

  u <- theta[1:fp$numKnots]
  beta <- numeric(fp$numKnots)
  beta[1] <- u[1]
  beta[2] <- u[1]+u[2]
  for(i in 3:fp$numKnots)
    beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]

  return(list(rvec = as.vector(fp$rvec.spldes %*% beta),
              iota = exp(theta[fp$numKnots+1]),
              ancbias = theta[fp$numKnots+2]))
}

ll <- function(theta, fp, likdat){

  param <- fnCreateParam(theta, fp)
  fp <- update(fp, list=param)

  if(min(fp$rvec)<0 || max(fp$rvec)>20) # Test positivity of rvec
    return(-Inf)
  
  mod <- fnEPP(fp)

  qM.all <- qnorm(prev(mod))
  qM.preg <- if(exists("pregprev", where=fp) && fp$pregprev) qnorm(fnPregPrev(mod, fp)) else qM.all

  if(any(is.na(qM.all)) || any(qM.all[-1] == -Inf))
    return(-Inf)

  ll.anc <- log(fnANClik(qM.preg+fp$ancbias, likdat$anclik.dat))
  ll.hhs <- fnHHSll(qM.all, likdat$hhslik.dat)

  if(exists("equil.rprior", where=fp) && fp$equil.rprior){
    rvec.ann <- fp$rvec[fp$proj.steps %% 1 == 0.5]
    equil.rprior.mean <- muSS/(1-pnorm(qM.all[likdat$lastdata.idx]))
    equil.rprior.sd <- sqrt(mean((muSS/(1-pnorm(qM.all[likdat$lastdata.idx - 10:1])) - rvec.ann[likdat$lastdata.idx - 10:1])^2))  # empirical sd based on 10 previous years
    ll.rprior <- sum(dnorm(rvec.ann[(likdat$lastdata.idx+1L):length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts year after last data
  } else
    ll.rprior <- 0
  
  return(ll.anc+ll.hhs+ll.rprior)
}


##########################
####  IMIS functions  ####
##########################

## Note: requires fp and likdat to be in global environment

sample.prior <- function(n){
  
  mat <- matrix(NA, n, fp$numKnots+3)
  
  ## sample penalty variance
  tau2 <- rexp(n, tau2.prior.rate)                  # variance of second-order spline differences
  
  mat[,1] <- rnorm(n, 1.5, 1)                                                     # u[1]
  mat[,2:fp$numKnots] <- rnorm(n*(fp$numKnots-1), 0, sqrt(tau2))                  # u[2:numKnots]
  mat[,fp$numKnots+1] <-  runif(n, logiota.unif.prior[1], logiota.unif.prior[2])  # iota
  mat[,fp$numKnots+2] <-  rnorm(n, ancbias.pr.mean, ancbias.pr.sd)                # ancbias parameter
  mat[,fp$numKnots+3] <- log(tau2)                                                # tau2
  
  
  return(mat)
}

prior <- function(theta){
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(lprior(theta[i,], fp))))))
}

likelihood <- function(theta){
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(ll(theta[i,], fp, likdat))))))
}
