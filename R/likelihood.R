#################
####  Prior  ####
#################

ldinvgamma <- function(x, alpha, beta){
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  return(log.density)
}

## r-spline prior parameters
logiota.unif.prior <- c(log(1e-14), log(0.0025))
tau2.prior.rate <- 0.5

invGammaParameter <- 0.001   #Inverse gamma parameter for tau^2 prior for spline
muSS <- 1/11.5               #1/duration for r steady state prior

## r-trend prior parameters
t0.unif.prior <- c(1970, 1990)
## t1.unif.prior <- c(10, 30)
## logr0.unif.prior <- c(1/11.5, 10)
t1.pr.mean <- 20.0
t1.pr.sd <- 4.5
logr0.pr.mean <- 0.42
logr0.pr.sd <- 0.23
## rtrend.beta.pr.mean <- 0.0
## rtrend.beta.pr.sd <- 0.2
rtrend.beta.pr.mean <- c(0.46, 0.17, -0.68, -0.038)
rtrend.beta.pr.sd <- c(0.12, 0.07, 0.24, 0.009)


###########################################
####                                   ####
####  Site level ANC data (SS and RT)  ####
####                                   ####
###########################################

ancbias.pr.mean <- 0.15
ancbias.pr.sd <- 1.0
vinfl.prior.rate <- 1/0.015

ancrtsite.beta.pr.mean <- 0
ancrtsite.beta.pr.sd <- 1.0
ancrtsite.vinfl.pr.rate <- 1/0.015


#' Prepare site-level ANC prevalence data for EPP random-effects likelihood
#'
#' @param eppd EPP data object
#' @param anchor.year year in which EPP data inputs start
#' NOTE: requires year to be stored in column names of anc.prev
prepare_ancsite_likdat <- function(eppd, anchor.year=1970L){

  anc.prev <- eppd$anc.prev
  anc.n <- eppd$anc.n
  anc.used <- eppd$anc.used

  ancobs.idx <- mapply(intersect, lapply(as.data.frame(t(!is.na(anc.prev))), which),
                       lapply(as.data.frame(t(!is.na(anc.n))), which), SIMPLIFY=FALSE)
  ## limit to years with both prevalence and N observations (likely input errors in EPP if not)

  nobs <- sapply(ancobs.idx, length)

  anc.years.lst <- lapply(ancobs.idx, function(i) as.integer(colnames(anc.prev)[i]))
  anc.prev.lst <- setNames(lapply(seq_along(ancobs.idx), function(i) as.numeric(anc.prev[i, ancobs.idx[[i]]])), rownames(anc.prev))
  anc.n.lst <- setNames(lapply(seq_along(ancobs.idx), function(i) as.numeric(anc.n[i, ancobs.idx[[i]]])), rownames(anc.n))

  X.lst <- mapply(cbind, Intercept=lapply(nobs, rep, x=1), ancrt=lapply(nobs, rep, x=0), SIMPLIFY=FALSE)

  if(exists("ancrtsite.prev", where=eppd)){
    ancrtsite.prev <- eppd$ancrtsite.prev
    ancrtsite.n <- eppd$ancrtsite.n
    
    ancrtsite.prev <- ancrtsite.prev[anc.used,,drop=FALSE]  # keep only used sites
    ancrtsite.n <- ancrtsite.n[anc.used,,drop=FALSE]        # keep only used sites

    ancrtsiteobs.idx <- mapply(intersect, lapply(as.data.frame(t(!is.na(ancrtsite.prev))), which),
                               lapply(as.data.frame(t(!is.na(ancrtsite.n))), which), SIMPLIFY=FALSE)
    ## limit to years with both prevalence and N observations (likely input errors in EPP if not)

    nobs <- sapply(ancrtsiteobs.idx, length)
    
    ancrtsite.years.lst <- lapply(ancrtsiteobs.idx, function(i) as.integer(colnames(ancrtsite.prev)[i]))
    ancrtsite.prev.lst <- setNames(lapply(seq_along(ancrtsiteobs.idx), function(i) as.numeric(ancrtsite.prev[i, ancrtsiteobs.idx[[i]]])), rownames(ancrtsite.prev))
    ancrtsite.n.lst <- setNames(lapply(seq_along(ancrtsiteobs.idx), function(i) as.numeric(ancrtsite.n[i, ancrtsiteobs.idx[[i]]])), rownames(ancrtsite.n))

    ancrtsite.X.lst <- mapply(cbind, Intercept=lapply(nobs, rep, x=1), ancrt=lapply(nobs, rep, x=1), SIMPLIFY=FALSE)

    ## Combine SS and RT data
    anc.years.lst <- mapply(c, anc.years.lst, ancrtsite.years.lst, SIMPLIFY=FALSE)
    anc.prev.lst <- mapply(c, anc.prev.lst, ancrtsite.prev.lst, SIMPLIFY=FALSE)
    anc.n.lst <- mapply(c, anc.n.lst, ancrtsite.n.lst, SIMPLIFY=FALSE)
    X.lst <- mapply(rbind, X.lst, ancrtsite.X.lst, SIMPLIFY=FALSE)
  }

  ## eliminate records with no observations
  anc.years.lst <- anc.years.lst[sapply(anc.years.lst, length) > 0]
  anc.prev.lst <- anc.prev.lst[sapply(anc.years.lst, length) > 0]
  anc.n.lst <- anc.n.lst[sapply(anc.years.lst, length) > 0]
  X.lst <- X.lst[sapply(anc.years.lst, length) > 0]
    
  x.lst <- mapply(function(p, n) (p*n+0.5)/(n+1), anc.prev.lst, anc.n.lst, SIMPLIFY=FALSE)
  W.lst <- lapply(x.lst, qnorm)
  v.lst <- mapply(function(W, x, n) 2*pi*exp(W^2)*x*(1-x)/n, W.lst, x.lst, anc.n.lst, SIMPLIFY=FALSE)
  anc.idx.lst <- lapply(anc.years.lst, "-", anchor.year-1)  ## index of observations relative to output prevalence vector


  anclik.dat <- list(W.lst = W.lst,
                     v.lst = v.lst,
                     n.lst = anc.n.lst,
                     X.lst = X.lst,
                     anc.idx.lst = anc.idx.lst)

  return(anclik.dat)
}

ll_anc <- function(qM, coef=c(0, 0), vinfl=0, anclik.dat){

  ## linear model offset
  mu <- lapply(lapply(anclik.dat$X.lst, "%*%", coef), c)

  d.lst <- mapply(function(w, mu, idx) w - (qM[idx]+mu), anclik.dat$W.lst, mu, anclik.dat$anc.idx.lst, SIMPLIFY=FALSE)
  v.lst <- lapply(anclik.dat$v.lst, "+", vinfl)

  return(log(anclik::anc_resid_lik(d.lst, v.lst)))
}



#############################################
####                                     ####
####  ANCRT census likelihood functions  ####
####                                     ####
#############################################

## prior parameters for ANCRT census
ancrtcens.bias.pr.mean <- 0
ancrtcens.bias.pr.sd <- 1.0
ancrtcens.vinfl.pr.rate <- 1/0.015

prepare_ancrtcens_likdat <- function(dat, anchor.year){

  x.ancrt <- (dat$prev*dat$n+0.5)/(dat$n+1)
  dat$W.ancrt <- qnorm(x.ancrt)
  dat$v.ancrt <- 2*pi*exp(dat$W.ancrt^2)*x.ancrt*(1-x.ancrt)/dat$n
  dat$idx <- dat$year - anchor.year+1

  return(dat)
}

ll_ancrtcens <- function(qM.preg, ancrtcens.dat, fp){
  sum(dnorm(ancrtcens.dat$W.ancrt,
            qM.preg[ancrtcens.dat$idx] + fp$ancrtcens.bias,
            sqrt(ancrtcens.dat$v.ancrt + fp$ancrtcens.vinfl), log=TRUE))
}



lprior <- function(theta, fp){

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(fp$eppmod %in% c("rspline", "logrspline")){
    epp_nparam <- fp$numKnots+2
    
    nk <- fp$numKnots
    tau2 <- exp(theta[nk+2])
    
    lpr <- sum(dnorm(theta[3:nk], 0, sqrt(tau2), log=TRUE)) +
      dunif(theta[nk+1], logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE) + 
      ldinvgamma(tau2, invGammaParameter, invGammaParameter) + log(tau2)   # + log(tau2): multiply likelihood by jacobian of exponential transformation

  } else { # rtrend
    epp_nparam <- 7

    lpr <- dunif(round(theta[1]), t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
      ## dunif(theta[2], t1.unif.prior[1], t1.unif.prior[2], log=TRUE) +
      dnorm(round(theta[2]), t1.pr.mean, t1.pr.sd, log=TRUE) +
      ## dunif(theta[3], logr0.unif.prior[1], logr0.unif.prior[2], log=TRUE) +
      dnorm(theta[3], logr0.pr.mean, logr0.pr.sd, log=TRUE) +
      sum(dnorm(theta[4:7], rtrend.beta.pr.mean, rtrend.beta.pr.sd, log=TRUE))
  }

  
  lpr <- lpr + dnorm(theta[epp_nparam+1], ancbias.pr.mean, ancbias.pr.sd, log=TRUE)
  if(!exists("v.infl", where=fp)){
    anclik_nparam <- 2
    lpr <- lpr + dexp(exp(theta[epp_nparam+2]), vinfl.prior.rate, TRUE) + theta[epp_nparam+2]         # additional ANC variance
  } else
    anclik_nparam <- 1


  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both")){
    lpr <- lpr + dnorm(theta[paramcurr+1], ancrtcens.bias.pr.mean, ancrtcens.bias.pr.sd, log=TRUE)
    if(!exists("ancrtcens.vinfl", fp)){
      lpr <- lpr + dexp(exp(theta[paramcurr+2]), ancrtcens.vinfl.pr.rate, TRUE) + theta[paramcurr+2]
      paramcurr <- paramcurr+2
    } else 
      paramcurr <- paramcurr+1
  } else if(exists("ancrt", fp) && fp$ancrt %in% c("site", "both")){
    lpr <- lpr +
      dnorm(theta[paramcurr+1], ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd, log=TRUE) ## +
      ## dexp(exp(theta[np]), ancrtsite.vinfl.pr.rate, TRUE) + theta[np]
  }
  
  return(lpr)
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

fnCreateLikDat <- function(epp.data, anchor.year=1970L){

  likdat <- list(anclik.dat = prepare_ancsite_likdat(epp.data, anchor.year=anchor.year),
                 hhslik.dat = fnPrepareHHSLikData(epp.data$hhs, anchor.year=anchor.year))
  if(exists("ancrtcens", where=epp.data))
    likdat$ancrtcens.dat <- prepare_ancrtcens_likdat(epp.data$ancrtcens, anchor.year=anchor.year)

  
  likdat$lastdata.idx <- max(unlist(likdat$anclik.dat$anc.idx.lst),
                             likdat$hhslik.dat$idx,
                             likdat$ancrtcens.dat$idx)
  likdat$firstdata.idx <- min(unlist(likdat$anclik.dat$anc.idx.lst),
                              likdat$hhslik.dat$idx,
                              likdat$ancrtcens.dat$idx)
  return(likdat)
}

fnCreateParam <- function(theta, fp){

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"
  
  if(fp$eppmod %in% c("rspline", "logrspline")){

    epp_nparam <- fp$numKnots+2
    
    u <- theta[1:fp$numKnots]
    beta <- numeric(fp$numKnots)
    beta[1] <- u[1]
    beta[2] <- u[1]+u[2]
    for(i in 3:fp$numKnots)
      beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]
    
    param <- list(rvec = as.vector(fp$rvec.spldes %*% beta),
                  iota = exp(theta[fp$numKnots+1]))

    if(fp$eppmod == "logrspline")
      param$rvec <- exp(param$rvec)

  } else { # rtrend
    epp_nparam <- 7
    param <- list(tsEpidemicStart = fp$proj.steps[which.min(abs(fp$proj.steps - (round(theta[1]-0.5)+0.5)))], # t0
                  rtrend = list(tStabilize = round(theta[1]-0.5)+0.5+round(theta[2]),  # t0 + t1
                                r0 = exp(theta[3]),              # r0
                                beta = theta[4:7]))
  }

  param$ancbias <- theta[epp_nparam+1]
  if(!exists("v.infl", where=fp)){
    anclik_nparam <- 2
    param$v.infl <- exp(theta[epp_nparam+2])
  } else
    anclik_nparam <- 1


  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both")){
    param$ancrtcens.bias <- theta[paramcurr+1]
    if(!exists("ancrtcens.vinfl", fp)){
      param$ancrtcens.vinfl <- exp(theta[paramcurr+2])
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", fp) && fp$ancrt %in% c("site", "both")){
    param$ancrtsite.beta <- theta[paramcurr+1]
    ## param$ancrtsite.vinfl <- exp(theta[length(theta)])
  }

  return(param)
}

ll <- function(theta, fp, likdat){

  param <- fnCreateParam(theta, fp)
  fp <- update(fp, list=param)

  if(!exists("eppmod", where=fp) || fp$eppmod == "rspline")
    if(min(fp$rvec)<0 || max(fp$rvec)>20) # Test positivity of rvec
      return(-Inf)
  
  mod <- simmod(fp)

  qM.all <- suppressWarnings(qnorm(prev(mod)))
  qM.preg <- if(exists("pregprev", where=fp) && !fp$pregprev) qM.all else suppressWarnings(qnorm(fnPregPrev(mod, fp)))

  if(any(is.na(qM.preg)) || any(qM.preg[likdat$firstdata.idx:likdat$lastdata.idx] == -Inf) || any(qM.preg[likdat$firstdata.idx:likdat$lastdata.idx] > 2)) # prevalence not greater than pnorm(2) = 0.977
    return(-Inf)

  ll.anc <- ll_anc(qM.preg, coef=c(fp$ancbias, fp$ancrtsite.beta), vinfl=fp$v.infl, likdat$anclik.dat)
  ll.hhs <- fnHHSll(qM.all, likdat$hhslik.dat)

  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both"))
    ll.ancrt <- ll_ancrtcens(qM.preg, likdat$ancrtcens.dat, fp)
  else
    ll.ancrt <- 0

  if(exists("equil.rprior", where=fp) && fp$equil.rprior){
    rvec.ann <- fp$rvec[fp$proj.steps %% 1 == 0.5]
    equil.rprior.mean <- muSS/(1-pnorm(qM.all[likdat$lastdata.idx]))
    equil.rprior.sd <- sqrt(mean((muSS/(1-pnorm(qM.all[likdat$lastdata.idx - 10:1])) - rvec.ann[likdat$lastdata.idx - 10:1])^2))  # empirical sd based on 10 previous years
    ll.rprior <- sum(dnorm(rvec.ann[(likdat$lastdata.idx+1L):length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts year after last data
  } else
    ll.rprior <- 0
  
  return(ll.anc+ll.hhs+ll.ancrt+ll.rprior)
}


##########################
####  IMIS functions  ####
##########################

sample.prior <- function(n, fp){

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  ## Calculate number of parameters
  if(fp$eppmod %in% c("rspline", "logrspline"))
    epp_nparam <- fp$numKnots+2L
  else
    epp_nparam <- 7

  if(!exists("v.infl", fp))
    anclik_nparam <- 2
  else
    anclik_nparam <- 1

  if(exists("ancrt", fp) && fp$ancrt == "both")
    ancrt_nparam <- 2
  else if(exists("ancrt", fp) && fp$ancrt == "census")
    ancrt_nparam <- 1
  else if(exists("ancrt", fp) && fp$ancrt == "site")
    ancrt_nparam <- 1
  else
    ancrt_nparam <- 0

  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both") && !exists("ancrtcens.vinfl", fp))
    ancrt_nparam <- ancrt_nparam+1


  nparam <- epp_nparam+anclik_nparam+ancrt_nparam

  ## Create matrix for storing samples
  mat <- matrix(NA, n, nparam)

  if(fp$eppmod %in% c("rspline", "logrspline")){
    epp_nparam <- fp$numKnots+2
       
    ## sample penalty variance
    tau2 <- rexp(n, tau2.prior.rate)                  # variance of second-order spline differences

    if(fp$eppmod == "rspline")
      mat[,1] <- rnorm(n, 1.5, 1)                                                   # u[1]
    else # logrspline
      mat[,1] <- rnorm(n, 0.2, 1)                                                   # u[1]
    mat[,2:fp$numKnots] <- rnorm(n*(fp$numKnots-1), 0, sqrt(tau2))                  # u[2:numKnots]
    mat[,fp$numKnots+1] <-  runif(n, logiota.unif.prior[1], logiota.unif.prior[2])  # iota
    mat[,fp$numKnots+2] <- log(tau2)                                                # tau2

  } else { # r-trend

    mat[,1] <- runif(n, t0.unif.prior[1], t0.unif.prior[2])        # t0
    ## mat[,2] <- runif(n, t1.unif.prior[1], t1.unif.prior[2])        # t1
    mat[,2] <- rnorm(n, t1.pr.mean, t1.pr.sd)
    ## mat[,3] <- runif(n, logr0.unif.prior[1], logr0.unif.prior[2])  # r0
    mat[,3] <- rnorm(n, logr0.pr.mean, logr0.pr.sd)  # r0
    mat[,4:7] <- t(matrix(rnorm(4*n, rtrend.beta.pr.mean, rtrend.beta.pr.sd), 4, n))  # beta
  }

  ## sample ANC bias paramters
  mat[,epp_nparam+1] <- rnorm(n, ancbias.pr.mean, ancbias.pr.sd)   # ancbias parameter
  if(!exists("v.infl", where=fp))
    mat[,epp_nparam+2] <- log(rexp(n, vinfl.prior.rate))

  ## sample ANCRT parameters
  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", where=fp) && fp$ancrt %in% c("census", "both")){
    mat[,paramcurr+1] <- rnorm(n, ancrtcens.bias.pr.mean, ancrtcens.bias.pr.sd)
    if(!exists("ancrtcens.vinfl", fp)){
      mat[,paramcurr+2] <- log(rexp(n, ancrtcens.vinfl.pr.rate))
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", where=fp) && fp$ancrt %in% c("site", "both")){
    mat[,paramcurr+1] <- rnorm(n, ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd)
    ## mat[,nparam] <- log(rexp(n, ancrtsite.vinfl.pr.rate))
  }

  return(mat)
}

prior <- function(theta, fp){
  if(is.vector(theta))
    return(exp(lprior(theta, fp)))
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(lprior(theta[i,], fp))))))
}

likelihood <- function(theta, fp, likdat){
  if(is.vector(theta))
    return(exp(ll(theta, fp, likdat)))
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(ll(theta[i,], fp, likdat))))))
}
