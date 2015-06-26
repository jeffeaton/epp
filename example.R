##############################################################################
###                                                                        ###
###  Example of fitting EPP r-spline and r-trend model to data from        ###
###  to data from Botswana.                                                ###
###                                                                        ###
###  Created on 19 June 2015 by Jeff Eaton (jeffrey.eaton@imperial.ac.uk)  ###
###                                                                        ###
##############################################################################

setwd("~/Documents/Code/R/epp/")

source("~/Documents/Code/R/read-epp-spectrum/read-epp-files.R")  # https://github.com/jeffeaton/read-epp-spectrum
source("R/epp.R")
source("R/likelihood.R")


## Function to do the following:
## (1) Read data, EPP subpopulations, and popualation inputs
## (2) Prepare timestep inputs for each EPP subpopulation

prepare.epp.fit <- function(filepath, proj.end=2015.5){

  ## epp
  eppd <- read.epp.data(paste(filepath, ".xml", sep=""))
  epp.subp <- read.epp.subpops(paste(filepath, ".xml", sep=""))
  epp.input <- read.epp.input(filepath)

  epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd)

  ## output
  val <- setNames(vector("list", length(eppd)), names(eppd))

  set.list.attr <- function(obj, attrib, value.lst)
    mapply(function(set, value){ attributes(set)[[attrib]] <- value; set}, obj, value.lst)

  val <- set.list.attr(val, "eppd", eppd)
  val <- set.list.attr(val, "likdat", lapply(eppd, fnCreateLikDat))
  val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar, proj.end = proj.end))
  val <- set.list.attr(val, "country", attr(eppd, "country"))
  val <- set.list.attr(val, "region", names(eppd))

  return(val)
}


## Read Botswana data and prepare fit (available for download: http://apps.unaids.org/spectrum/)

bw.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Botswana 2014/Botswana 2014_Nat 19_06_14-c"   

bw.out <- prepare.epp.fit(bw.path, proj.end=2015.5)


#########################
####  Run EPP model  ####
#########################

## r-spline model: fixed parameter values
theta.rspline <- c(2.16003605, -0.76713859, 0.21682066, 0.03286402, 0.21494412,
                   0.40138627, -0.08235464, -16.32721684, 0.21625028, -2.97511957)

fp <- attr(bw.out$Urban, "eppfp")
param <- fnCreateParam(theta.rspline, fp)
fp.rspline <- update(fp, list=param)
mod.rspline <- fnEPP(fp.rspline)

round(prev(mod.rspline), 3)           # prevalence
round(incid(mod.rspline, fp.rspline), 4)  # incidence

likdat <- attr(bw.out$Urban, "likdat")
qM <- qnorm(prev(mod.rspline))                              # probit-tranformed prevalence
log(fnANClik(qM + fp.rspline$ancbias, likdat$anclik.dat))   # ANC likelihood
fnHHSll(qM, likdat$hhslik.dat)                              # survey likelihood
ll(theta.rspline, fp.rspline, likdat)



## r-trend model: fixed parameter values
param.rtrend <- list(beta = c(0.46, 0.17, -0.68, -0.038),
                     t.stabilize = 1978+20,
                     r0 = exp(0.42))
fp.rtrend <- update(fp, eppmod = "rtrend", iota = 0.0025, tsEpidemicStart = 1978.0,
                    rtrend = par.rtrend, ancbias = 0.21625028)
mod.rtrend <- fnEPP(fp.rtrend)

round(prev(mod.rtrend), 3)           # prevalence
round(incid(mod.rtrend, fp.rtrend), 4)  # incidence

qM <- qnorm(prev(mod.rtrend))                             # probit-tranformed prevalence
log(fnANClik(qM + fp.rtrend$ancbias, likdat$anclik.dat))  # ANC likelihood
fnHHSll(qM, likdat$hhslik.dat)                             # survey likelihood
## ll(theta.rtrend, fp.rtrend, likdat)


#########################
####  Fit EPP model  ####
#########################

fit.mod <- function(obj, ..., B0 = 1e5, B = 1e4, B.re = 3000, number_k = 500){
  ## ... : updates to fixed parameters (fp) object to specify fitting options
  
  likdat <<- attr(obj, 'likdat')  # put in global environment for IMIS functions.
  fp <<- attr(obj, 'eppfp')
  fp <<- update(fp, ...)
  
  fit <- IMIS(B0, B, B.re, number_k)
  fit$fp <- fp
  fit$likdat <- likdat

  rm(fp, likdat, pos=.GlobalEnv)

  return(fit)
}
  
bw.out$Urban <- fit.mod(bw.out$Urban, equil.rprior=TRUE, B0=1e4, B=1e3)
bw.out$Rural <- fit.mod(bw.out$Rural, equil.rprior=TRUE, B0=1e4, B=1e3)
## Note: This crashes if there are fewer than two parameter combinations
##       with non-zero likelihood. In this case run again with different
##       seed, or larger B0.


######################################
####  Simulate posterior outputs  ####
######################################

## simulate incidence and prevalence
sim.mod <- function(fit){
  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list=par))
  fit$mod <- lapply(fp.list, fnEPP)
  fit$prev <- sapply(fit$mod, prev)
  fit$incid <- mapply(incid, mod = fit$mod, fp = fp.list)
  return(fit)
}

bw.out$Urban <- sim.mod(bw.out$Urban)
bw.out$Rural <- sim.mod(bw.out$Rural)


## Plot prevalence, incidence, r(t)
cred.region <- function(x, y, ...)
  polygon(c(x, rev(x)), c(y[1,], rev(y[2,])), border=NA, ...)

transp <- function(col, alpha=0.5)
  return(apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha)))

plot.prev <- function(fit, ylim=c(0, 0.22), col="blue"){
  plot(1970:2015, rowMeans(fit$prev), type="n", ylim=ylim, ylab="", yaxt="n", xaxt="n")
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  cred.region(1970:2015, apply(fit$prev, 1, quantile, c(0.025, 0.975)), col=transp(col, 0.3))
  lines(1970:2015, rowMeans(fit$prev), col=col)
  ##
  points(fit$likdat$hhslik.dat$year, fit$likdat$hhslik.dat$prev, pch=20)
  segments(fit$likdat$hhslik.dat$year,
           y0=pnorm(fit$likdat$hhslik.dat$W.hhs - qnorm(0.975)*fit$likdat$hhslik.dat$sd.W.hhs),
           y1=pnorm(fit$likdat$hhslik.dat$W.hhs + qnorm(0.975)*fit$likdat$hhslik.dat$sd.W.hhs))
}

plot.incid <- function(fit, ylim=c(0, 0.05), col="blue"){
  plot(1970:2015, rowMeans(fit$incid), type="n", ylim=ylim, ylab="", yaxt="n", xaxt="n")
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  cred.region(1970:2015, apply(fit$incid, 1, quantile, c(0.025, 0.975)), col=transp(col, 0.3))
  lines(1970:2015, rowMeans(fit$incid), col=col)
}

plot.rvec <- function(fit, ylim=c(0, 3), col="blue"){
  rvec <- sapply(fit$param, "[[", "rvec")
  plot(fit$fp$proj.steps, rowMeans(rvec), type="n", ylim=ylim, ylab="", yaxt="n")
  axis(2, labels=FALSE)
  cred.region(fit$fp$proj.steps, apply(rvec, 1, quantile, c(0.025, 0.975)), col=transp(col, 0.3))
  lines(fit$fp$proj.steps, rowMeans(rvec), col=col)
}



quartz(h=3.6, w=6, pointsize=8)

par(mfrow=c(2,3), tcl=-0.25, mgp=c(2, 0.5, 0), mar=c(2, 3.5, 2, 1), las=1, cex=1.0)
##
plot.prev(bw.out$Urban, col="darkred", ylim=c(0, 0.3))
axis(2, tick="no")
axis(1, tick="no")
mtext("prevalence", 2, 2.5, las=3)
mtext("Botswana Urban", line=0.5, adj=-1, font=2, cex=1.2)
##
plot.incid(bw.out$Urban, col="darkred", ylim=c(0, 0.06))
axis(2, tick="no")
axis(1, tick="no")
mtext("incidence", 2, 2.5, las=3)
##
plot.rvec(bw.out$Urban, col="darkred")
axis(2, tick="no")
axis(1, tick="no")
mtext("r(t)", 2, 2.5, las=3)
####
plot.prev(bw.out$Rural, col="darkolivegreen", ylim=c(0, 0.3))
axis(2, tick="no")
axis(1, tick="no")
mtext("prevalence", 2, 2.5, las=3)
mtext("Botswana Rural", line=0.5, adj=-1, font=2, cex=1.2)
##
plot.incid(bw.out$Rural, col="darkolivegreen", ylim=c(0, 0.06))
axis(2, tick="no")
axis(1, tick="no")
mtext("incidence", 2, 2.5, las=3)
##
plot.rvec(bw.out$Rural, col="darkolivegreen")
axis(2, tick="no")
axis(1, tick="no")
mtext("r(t)", 2, 2.5, las=3)


##########################################################
####  Simulate ANC posterior predictive distribution  ####
##########################################################

add.b.site <- function(fit){
  qM.mat <- sweep(qnorm(fit$prev), 2, sapply(fit$param, "[[", "ancbias"), "+")
  fit$b.site <- apply(qM.mat, 2, sample.b.site, fit$likdat$anclik.dat)
  return(fit)
}

add.pred.site <- function(fit){
  qM.mat <- sweep(qnorm(fit$prev), 2, sapply(fit$param, "[[", "ancbias"), "+")
  fit$pred.site <- lapply(seq(along=fit$param), function(ii) sample.pred.site(qM.mat[,ii], fit$b.site[,ii], fit$likdat$anclik.dat))
  return(fit)
}

pred.coverage <- function(fit){
  pred.quant <- apply(sapply(fit$pred.site, unlist), 1, quantile, c(0.025, 0.975))
  obs <- pnorm(unlist(fit$likdat$anclik.dat$W.lst))
  return(mean(obs > pred.quant[1,] & obs < pred.quant[2,]))
}

pred.quantile <- function(fit){
  pred.mat <- sapply(fit$pred.site, unlist)
  obs <- pnorm(unlist(fit$likdat$anclik.dat$W.lst))
  pred.quant <- sapply(seq_along(obs), function(i) ecdf(pred.mat[i,])(obs[i]))
  fit$pred.quant <- split(pred.quant, rep(names(fit$likdat$anclik.dat$W.lst), sapply(fit$likdat$anclik.dat$W.lst, length)))
  return(fit)
}

## Sample site-level random effects
bw.out$Urban <- add.b.site(bw.out$Urban)
bw.out$Rural <- add.b.site(bw.out$Rural)

## Sample from clinic posterior predictive distribution
bw.out$Urban <- add.pred.site(bw.out$Urban)
bw.out$Rural <- add.pred.site(bw.out$Rural)

## In-sample coverage of 95% prediction interval
pred.coverage(bw.out$Urban)
pred.coverage(bw.out$Rural)

## Q-Q plot of predicted vs. theoretical quantiles for ANC prevalence
bw.out$Urban <- pred.quantile(bw.out$Urban)
bw.out$Rural <- pred.quantile(bw.out$Rural)

quartz(w=6, h=3, pointsize=9)

par(mfrow=c(1,2), tcl=-0.25, mgp=c(2, 0.5, 0), mar=c(3, 3, 2.5, 1), las=1, cex=1.0)
##
plot(seq(0, 1, length.out=length(unlist(bw.out$Urban$pred.quant))),
     sort(unlist(bw.out$Urban$pred.quant)),
     pch=20, cex=0.5,
     main="Botswana Urban",
     xlab="Theoretical quantiles",
     ylab="Observed quantiles")
abline(a=0, b=1)
##
plot(seq(0, 1, length.out=length(unlist(bw.out$Rural$pred.quant))),
     sort(unlist(bw.out$Rural$pred.quant)),
     pch=20, cex=0.5,
     main="Botswana Rural",
     xlab="Theoretical quantiles",
     ylab="Observed quantiles")
abline(a=0, b=1)
