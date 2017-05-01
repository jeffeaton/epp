read_eppout_posterior <- function(dir, model="rtrend"){
  
  eppout <- readLines(list.files(dir, "_full.csv", full.names=TRUE))

  reg.idx <- grep("Summary of results and resamples for", eppout)
  reg_names <- sub("(.+)_([A-Za-z]*).*", "\\2", eppout[reg.idx])
  
  urban.idx <- grep("Summary of results and resamples for.*Urban", eppout)
  rural.idx <- grep("Summary of results and resamples for.*Rural", eppout)
  min(grep("=========", eppout)[which(grep("=========", eppout) > urban.idx+2)])

  if(model=="rtrend"){
    urban.resample <- read.csv(text=eppout[(urban.idx+27):(min(grep("=========", eppout)[which(grep("=========", eppout) > urban.idx+2)])-1)], header=FALSE)
    rural.resample <- read.csv(text=eppout[(rural.idx+27):(min(grep("=========", eppout)[which(grep("=========", eppout) > rural.idx+2)])-1)], header=FALSE)
    names(urban.resample) <- names(rural.resample) <- c("Sample", "Count", "Weight", "Likelihood", "Prior", "t0", "t1", "logr0", "beta0", "beta1", "beta2", "beta3", "vinfl", "ancbias")
    
    urban.resample <- urban.resample[rep(seq_len(nrow(urban.resample)), urban.resample$Count), c("t0", "t1", "logr0", "beta0", "beta1", "beta2", "beta3", "ancbias", "vinfl")]
    rural.resample <- rural.resample[rep(seq_len(nrow(rural.resample)), rural.resample$Count), c("t0", "t1", "logr0", "beta0", "beta1", "beta2", "beta3", "ancbias", "vinfl")]
    urban.resample$logvinfl <- log(urban.resample$vinfl)
    rural.resample$logvinfl <- log(rural.resample$vinfl)
    urban.resample$vinfl <- NULL
    rural.resample$vinfl <- NULL
  }

  if(model=="rspline"){
    urban.resample <- read.csv(text=eppout[(urban.idx+20):(min(grep("=========", eppout)[which(grep("=========", eppout) > urban.idx+2)])-1)], header=FALSE)
    rural.resample <- read.csv(text=eppout[(rural.idx+20):(min(grep("=========", eppout)[which(grep("=========", eppout) > rural.idx+2)])-1)], header=FALSE)
    names(urban.resample) <- names(rural.resample) <- c("Sample", "Count", "Weight", "Likelihood", "Prior", "logiota", "logtau2", paste0("u", 1:7), "vinfl", "ancbias")
    
    urban.resample <- urban.resample[rep(seq_len(nrow(urban.resample)), urban.resample$Count), c(paste0("u", 1:7), "logiota", "logtau2", "ancbias", "vinfl")]
    rural.resample <- rural.resample[rep(seq_len(nrow(rural.resample)), rural.resample$Count), c(paste0("u", 1:7), "logiota", "logtau2", "ancbias", "vinfl")]
    ## urban.resample$logvinfl <- log(urban.resample$vinfl)
    ## rural.resample$logvinfl <- log(rural.resample$vinfl)
    ## urban.resample$vinfl <- NULL
    ## rural.resample$vinfl <- NULL
  }

  return(list(Urban=urban.resample, Rural=rural.resample))
}


read_rvec <- function(dir){
  eppout <- readLines(list.files(dir, "_CD4.out", full.names=TRUE))
  rvec <- setNames(as.numeric(scan(text=eppout[grep("RValues", eppout)], sep=",", what="raw", quiet=TRUE)[-1]),
                   scan(text=eppout[grep("Times", eppout)], sep=",", what="raw", quiet=TRUE)[-1])
  rvec <- rvec[!is.na(rvec)]
  return(rvec)
}

read_annual_outputs <- function(dir){
  eppout <- readLines(list.files(dir, "_CD4.out", full.names=TRUE))
  out <- as.data.frame(t(read.table(text=eppout[grep("Annual output values from integrator without any calibration:", eppout)+2:13], row.names=1)))

  class(out) <- c("eppdebugann", "data.frame")
  return(out)
}


read_ts_outputs <- function(dir){
  eppout <- readLines(list.files(dir, "_CD4.out", full.names=TRUE))

  breaks.idx <- grep("-----------------", eppout)
  pop.idx <- grep("Populations in various compartments by timestep", eppout)
  hiv.idx <- grep("HIV and ART related values", eppout)

  hiv.out <- read.table(text=eppout[do.call(seq.int, as.list(breaks.idx[which(breaks.idx > hiv.idx)[1:2]]+c(1,-1)))])
  names(hiv.out) <- scan(text=gsub("# ", "", gsub("[ ][ #]+", ";", eppout[hiv.idx+2])), what="character", sep=";", quiet=TRUE)[-1]
  pop.out <- read.table(text=eppout[do.call(seq.int, as.list(breaks.idx[which(breaks.idx > pop.idx)[1:2]]+c(1,-1)))])
  names(pop.out) <- scan(text=gsub("[ ][ ]+", ";", eppout[pop.idx+2]), what="character", sep=";", quiet=TRUE)[-1]

  out <- cbind(pop.out, hiv.out[-1])

  rvec <- read_rvec(dir)
  out <- merge(out, data.frame(Year=names(rvec), rvec=rvec), all.x=TRUE)
  
  class(out) <- c("eppdebug", "data.frame")

  return(out)
}




pop15to49.eppdebug <- function(eppdebug){setNames(eppdebug[["Total pop N"]], eppdebug$Year)}
artcov15to49.eppdebug <- function(eppdebug){setNames(eppdebug[["on ART"]] / eppdebug[["current HIV"]], eppdebug$Year)}
prev.eppdebug <- function(eppdebug){setNames(eppdebug[["current HIV"]] / eppdebug[["Total pop N"]], eppdebug$Year)}
incid.eppdebug <- function(eppdebug){NA}

pop15to49.eppdebugann <- function(eppdebugann){setNames(eppdebugann$Pop, eppdebugann$Year)}
artcov15to49.eppdebugann <- function(eppdebugann){NA}
prev.eppdebugann <- function(eppdebugann){setNames(eppdebugann$Prev / 100, eppdebugann$Year)}
incid.eppdebugann <- function(eppdebugann){setNames(eppdebugann$Incid / 100, eppdebugann$Year)}
