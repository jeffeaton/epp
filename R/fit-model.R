fitmod <- function(obj, ..., B0 = 1e5, B = 1e4, B.re = 3000, number_k = 500, D=0, opt_iter=0,
                   sample.prior=epp:::sample.prior,
                   prior=epp:::prior,
                   likelihood=epp:::likelihood){

  ## ... : updates to fixed parameters (fp) object to specify fitting options

  fp <- attr(obj, 'eppfp')
  fp <- update(fp, ...)

  ## Prepare likelihood data
  eppd <- attr(obj, "eppd")
  
  if(exists("ancrt", fp) && fp$ancrt == "none")
    eppd$ancrtcens <- eppd$ancrtsite.prev <- eppd$ancrtsite.n <- NULL
  else if(exists("ancrt", fp) && fp$ancrt == "census")
    eppd$ancrtsite.prev <- eppd$ancrtsite.n <- NULL
  else if(exists("ancrt", fp) && fp$ancrt == "site")
    eppd$ancrtcens <- NULL

  if(is.null(eppd$ancrtcens) && is.null(eppd$ancrtsite.prev)){
    fp$ancrt <- "none"
    fp$ancrtsite.beta <- 0
  } else if(!is.null(eppd$ancrtcens) && is.null(eppd$ancrtsite.prev)){
    fp$ancrt <- "census"
    fp$ancrtsite.beta <- 0
  } else if(is.null(eppd$ancrtcens) && !is.null(eppd$ancrtsite.prev))
    fp$ancrt <- "site"
  else
    fp$ancrt <- "both"

  likdat <- fnCreateLikDat(eppd, floor(fp$proj.steps[1]))
  
  ## If IMIS fails, start again
  fit <- try(stop(""), TRUE)
  while(inherits(fit, "try-error")){
    start.time <- proc.time()
    fit <- try(IMIS(B0, B, B.re, number_k, D, opt_iter, fp=fp, likdat=likdat,
                    sample.prior=sample.prior, prior=prior, likelihood=likelihood))
    fit.time <- proc.time() - start.time
  }
  fit$fp <- fp
  fit$likdat <- likdat
  fit$time <- fit.time

  class(fit) <- "eppfit"

  return(fit)
}


sim_rvec_rwproj <- function(rvec, firstidx, lastidx, dt){
  logrvec <- log(rvec)
  sd <- sqrt(mean(diff(logrvec[firstidx:lastidx])^2))
  projidx <- (lastidx+1):length(rvec)

  ## simulate differences in projection log(rvec)
  ## variance increases with time: sigma^2*(t-t1) [Hogan 2012]

  ## Note: implementation below matches R code from Dan Hogan and Java code
  ## from Tim Brown. Implies that variance grows at rate dt*(t-t1)
  ldiff <- rnorm(length(projidx), 0, sd*sqrt((1+dt*(projidx-lastidx-1))))

  rvec[projidx] <- exp(logrvec[lastidx] + cumsum(ldiff))
  return(rvec)
}


## simulate incidence and prevalence
simfit.eppfit <- function(fit, rwproj=FALSE){
  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))

  if(rwproj){
    if(exists("eppmod", where=fit$fp) && fit$fp$eppmod == "rtrend")
      stop("Random-walk projection is only used with r-spline model")

    fit$rvec.spline <- sapply(fit$param, "[[", "rvec")
    firstidx <- which(fit$fp$proj.steps == fit$fp$tsEpidemicStart)
    lastidx <- (fit$likdat$lastdata.idx-1)/fit$fp$dt+1

    ## replace rvec with random-walk simulated rvec
    fit$param <- lapply(fit$param, function(par){par$rvec <- sim_rvec_rwproj(par$rvec, firstidx, lastidx, fit$fp$dt); par})
  }
  
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list=par))
  mod.list <- lapply(fp.list, simmod)
  
  fit$rvec <- sapply(mod.list, attr, "rvec")
  fit$prev <- sapply(mod.list, prev)
  fit$incid <- mapply(incid, mod = mod.list, fp = fp.list)
  fit$popsize <- sapply(mod.list, rowSums)
  fit$pregprev <- mapply(fnPregPrev.epp, mod.list, fp.list)
  
  return(fit)
}



## Function to do the following:
## (1) Read data, EPP subpopulations, and popualation inputs
## (2) Prepare timestep inputs for each EPP subpopulation

prepare_epp_fit <- function(pjnz, proj.end=2015.5){

  ## epp
  eppd <- read_epp_data(pjnz)
  epp.subp <- read_epp_subpops(pjnz)
  epp.input <- read_epp_input(pjnz)

  epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd)

  ## output
  val <- setNames(vector("list", length(eppd)), names(eppd))

  set.list.attr <- function(obj, attrib, value.lst)
    mapply(function(set, value){ attributes(set)[[attrib]] <- value; set}, obj, value.lst)

  val <- set.list.attr(val, "eppd", eppd)
  val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar, proj.end = proj.end))
  val <- set.list.attr(val, "country", attr(eppd, "country"))
  val <- set.list.attr(val, "region", names(eppd))

  return(val)
}
