
###########################################################################
####  Function to read inputs from Spectrum to EPP (.ep1, .ep3, .ep4)  ####
###########################################################################

read_epp_input <- function(ep.path){

  ## ep1
  ep1 <- scan(paste(ep.path, ".ep1", sep=""), "character", sep="\n")

  country.idx <- which(sapply(ep1, substr, 1, 7) == "COUNTRY")
  firstprojyr.idx <-  which(sapply(ep1, substr, 1, 11) == "FIRSTPROJYR")
  lastprojyr.idx <-  which(sapply(ep1, substr, 1, 10) == "LASTPROJYR")
  popstart.idx <- which(ep1 == "POPSTART")+1
  popend.idx <- which(ep1 == "POPEND")-1

  country <- as.character(read.csv(text=ep1[country.idx], header=FALSE, as.is=TRUE)[2])
  country.code <- as.integer(read.csv(text=ep1[country.idx], header=FALSE)[3])

  start.year <- as.integer(read.csv(text=ep1[firstprojyr.idx], header=FALSE)[2])
  stop.year <- as.integer(read.csv(text=ep1[lastprojyr.idx], header=FALSE)[2])
  epp.pop <- setNames(read.csv(text=ep1[popstart.idx:popend.idx], header=FALSE, as.is=TRUE),
                      c("year", "pop15to49", "pop15", "pop50", "netmigr"))

  ## ep4
  ep4 <- scan(paste(ep.path, ".ep4", sep=""), "character", sep="\n")

  cd4lim.idx <- which(sapply(ep4, substr, 1, 12) == "CD4LOWLIMITS")
  lambda.idx <- which(sapply(ep4, substr, 1, 6) == "LAMBDA")
  cd4init.idx <- which(sapply(ep4, substr, 1, 13) == "NEWINFECTSCD4")
  mu.idx <- which(sapply(ep4, substr, 1, 3) == "MU_")
  alpha1.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA1")
  alpha2.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA2")
  alpha3.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA3")
  infectreduc.idx <- which(sapply(ep4, substr, 1, 11) == "INFECTREDUC")
  artstart.idx <- which(ep4 == "ARTSTART")+1
  artend.idx <- which(ep4 == "ARTEND")-1

  cd4lim <- as.integer(read.csv(text=ep4[cd4lim.idx], header=FALSE)[-1])
  cd4init <- as.matrix(read.csv(text=ep4[cd4init.idx], header=FALSE, row.names=1))
  lambda <- as.matrix(read.csv(text=ep4[lambda.idx], header=FALSE, row.names=1))
  mu <- as.matrix(read.csv(text=ep4[mu.idx], header=FALSE, row.names=1))
  alpha1 <- as.matrix(read.csv(text=ep4[alpha1.idx], header=FALSE, row.names=1))
  alpha2 <- as.matrix(read.csv(text=ep4[alpha2.idx], header=FALSE, row.names=1))
  alpha3 <- as.matrix(read.csv(text=ep4[alpha3.idx], header=FALSE, row.names=1))
  infectreduc <- as.numeric(read.csv(text=ep4[infectreduc.idx], header=FALSE)[2])

  epp.art <- setNames(read.csv(text=ep4[artstart.idx:artend.idx], header=FALSE, as.is=TRUE),
                      c("year", "m.isperc", "m.val", "f.isperc", "f.val", "cd4thresh", "m.perc50plus", "f.perc50plus", "perc50plus", ""))

  eppin <- list(start.year       = start.year,
                stop.year        = stop.year,
                epp.pop          = epp.pop,
                cd4lowlim        = cd4lim,
                cd4initperc      = cd4init,
                cd4stage.dur     = lambda,
                cd4mort          = mu,
                artmort.less6mos = alpha1,
                artmort.6to12mos = alpha2,
                artmort.after1yr = alpha3,
                infectreduc      = infectreduc,
                epp.art          = epp.art)
  class(eppin) <- "eppin"
  attr(eppin, "country") <- country
  attr(eppin, "country.code") <- country.code

  return(eppin)
}


############################################################################
####  Function to read prevalence data used in EPP fitting (from .xml)  ####
############################################################################

read_epp_data <- function(epp.xml){

  if (!require("XML", quietly = TRUE))
    stop("read_epp_data() requires the package 'XML'. Please install it.", call. = FALSE)
      
  obj <- xmlTreeParse(epp.xml)

  r <- xmlRoot(obj)[[1]]
  eppSetChildren.idx <- which(xmlSApply(r, xmlAttrs) == "eppSetChildren")
  country <- xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetCountry")]][[1]])

  epp.data <- list() # declare list to store output
  attr(epp.data, "country") <- country

  for(eppSet.idx in 1:xmlSize(r[[eppSetChildren.idx]])){

    eppSet <- r[[eppSetChildren.idx]][[eppSet.idx]][[1]]
    eppName <- xmlToList(eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "name")]][["string"]])

    ##  ANC data  ##

    siteNames.idx <- which(xmlSApply(eppSet, xmlAttrs) == "siteNames")
    siteSelected.idx <- which(xmlSApply(eppSet, xmlAttrs) == "siteSelected")
    survData.idx <- which(xmlSApply(eppSet, xmlAttrs) == "survData")
    survSampleSizes.idx <- which(xmlSApply(eppSet, xmlAttrs) == "survSampleSizes")

    siteNames <- xmlSApply(eppSet[[siteNames.idx]][[1]], xmlSApply, xmlToList, FALSE)
    siteIdx <- as.numeric(xmlSApply(eppSet[[siteNames.idx]][[1]], xmlAttrs)) ## 0 based
    anc.used <- as.logical(xmlSApply(eppSet[[siteSelected.idx]][[1]], xmlSApply, xmlToList, FALSE))

    nsites <- length(siteNames)
    nANCyears <- max(as.integer(xmlSApply(eppSet[[survData.idx]][["array"]][[1]][[1]], xmlAttrs))) + 1

    ## ANC prevalence
    anc.prev <- matrix(NA, nsites, nANCyears)
    rownames(anc.prev) <- siteNames
    colnames(anc.prev) <- 1985+0:(nANCyears-1)
    for(clinic.idx in 1:nsites){
      clinic <- eppSet[[survData.idx]][["array"]][[clinic.idx]][[1]]
      prev <- as.numeric(xmlSApply(clinic, xmlSApply, xmlToList, FALSE))
      idx <- as.integer(xmlSApply(clinic, xmlAttrs)) + 1
      anc.prev[clinic.idx, idx] <- prev
    }
    anc.prev[is.na(anc.prev)] <- 0.0 ## NOTE: appears that if value is 0.0, the array index is omitted from XML file, might apply elsewhere.
    anc.prev[anc.prev == -1] <- NA
    anc.prev <- anc.prev/100

    ## ANC sample sizes
    anc.n <- matrix(NA, nsites, nANCyears)
    rownames(anc.n) <- siteNames
    colnames(anc.n) <- 1985+0:(nANCyears-1)
    for(clinic.idx in 1:nsites){
      clinic <- eppSet[[survSampleSizes.idx]][["array"]][[clinic.idx]][[1]]
      n <- as.numeric(xmlSApply(clinic, xmlSApply, xmlToList, FALSE))
      idx <- as.integer(xmlSApply(clinic, xmlAttrs)) + 1
      anc.n[clinic.idx, idx] <- n
    }
    anc.n[anc.n == -1] <- NA


    ##  HH surveys  ##

    hhsUsed.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyIsUsed")
    hhsHIV.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyHIV")
    hhsSampleSize.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveySampleSize")
    hhsSE.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyStandardError")
    hhsYear.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyYears")

    nhhs <- max(xmlSize(eppSet[[hhsYear.idx]][[1]]),
                xmlSize(eppSet[[hhsHIV.idx]][[1]]),
                xmlSize(eppSet[[hhsSE.idx]][[1]]),
                xmlSize(eppSet[[hhsSampleSize.idx]][[1]]),
                xmlSize(eppSet[[hhsUsed.idx]][[1]]))

    hhs <- data.frame(year = rep(NA, nhhs), prev = rep(NA, nhhs), se = rep(NA, nhhs), n = rep(NA, nhhs), used = rep(NA, nhhs))

    hhs$year[as.integer(xmlSApply(eppSet[[hhsYear.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[hhsYear.idx]][[1]], xmlSApply, xmlToList, FALSE))
    hhs$prev[as.integer(xmlSApply(eppSet[[hhsHIV.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[hhsHIV.idx]][[1]], xmlSApply, xmlToList, FALSE))/100
    hhs$se[as.integer(xmlSApply(eppSet[[hhsSE.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[hhsSE.idx]][[1]], xmlSApply, xmlToList, FALSE))/100
    hhs$n[as.integer(xmlSApply(eppSet[[hhsSampleSize.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[hhsSampleSize.idx]][[1]], xmlSApply, xmlToList, FALSE))
    hhs$used[as.integer(xmlSApply(eppSet[[hhsUsed.idx]][[1]], xmlAttrs))+1] <- as.logical(xmlSApply(eppSet[[hhsUsed.idx]][[1]], xmlSApply, xmlToList, FALSE))

    hhs <- subset(hhs, !is.na(prev))

    epp.data[[eppName]] <- list(country=country,
                                region=eppName,
                                anc.used=anc.used,
                                anc.prev=anc.prev,
                                anc.n=anc.n,
                                hhs=hhs)
  }

  class(epp.data) <- "eppd"

  return(epp.data)
}


################################################################################
####  Function to read subpopulation sizes used in EPP fitting (from .xml)  ####
################################################################################

read_epp_subpops <- function(epp.xml){

  obj <- xmlTreeParse(epp.xml)
  r <- xmlRoot(obj)[[1]]
  eppSetChildren.idx <- which(xmlSApply(r, xmlAttrs) == "eppSetChildren")
  country <- xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetCountry")]][[1]])

  epp.pops <- list() # declare list to store output
  attr(epp.pops, "country") <- country

  workset.startyear <- as.integer(xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetStartYear")]][[1]]))
  workset.endyear <- as.integer(xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetEndYear")]][[1]]))
  ## workset.popbaseyear <- as.integer(xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetPopBaseYear")]][[1]])) # not sure what this is...

  pop15to49.idx <- which(xmlSApply(r, xmlAttrs) == "pop15to49")
  pop15.idx <- which(xmlSApply(r, xmlAttrs) == "pop15")
  pop50.idx <- which(xmlSApply(r, xmlAttrs) == "pop50")
  netMigration.idx <- which(xmlSApply(r, xmlAttrs) == "netMigration")

  epp.pops$total <-  data.frame(year = workset.startyear:workset.endyear,
                                pop15to49 = 0,
                                pop15 = 0,
                                pop50 = 0,
                                netmigr = 0)
  epp.pops$total$pop15to49[as.integer(xmlSApply(r[[pop15to49.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(r[[pop15to49.idx]][[1]], xmlSApply, xmlToList))
  epp.pops$total$pop15[as.integer(xmlSApply(r[[pop15.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(r[[pop15.idx]][[1]], xmlSApply, xmlToList))
  epp.pops$total$pop50[as.integer(xmlSApply(r[[pop50.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(r[[pop50.idx]][[1]], xmlSApply, xmlToList))
  epp.pops$total$netmigr[as.integer(xmlSApply(r[[netMigration.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(r[[netMigration.idx]][[1]], xmlSApply, xmlToList))

  epp.pops$subpops <- list()

  for(eppSet.idx in 1:xmlSize(r[[eppSetChildren.idx]])){

    eppSet <- r[[eppSetChildren.idx]][[eppSet.idx]][[1]]
    eppName <- xmlToList(eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "name")]][["string"]])

    pop15to49.idx <- which(xmlSApply(eppSet, xmlAttrs) == "pop15to49")
    pop15.idx <- which(xmlSApply(eppSet, xmlAttrs) == "pop15")
    pop50.idx <- which(xmlSApply(eppSet, xmlAttrs) == "pop50")
    netMigration.idx <- which(xmlSApply(eppSet, xmlAttrs) == "netMigration")

    subp <- data.frame(year = workset.startyear:workset.endyear,
                       pop15to49 = 0,
                       pop15 = 0,
                       pop50 = 0,
                       netmigr = 0)
    subp$pop15to49[as.integer(xmlSApply(eppSet[[pop15to49.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[pop15to49.idx]][[1]], xmlSApply, xmlToList))
    subp$pop15[as.integer(xmlSApply(eppSet[[pop15.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[pop15.idx]][[1]], xmlSApply, xmlToList))
    subp$pop50[as.integer(xmlSApply(eppSet[[pop50.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[pop50.idx]][[1]], xmlSApply, xmlToList))
    subp$netmigr[as.integer(xmlSApply(eppSet[[netMigration.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[netMigration.idx]][[1]], xmlSApply, xmlToList))

    epp.pops$subpops[[eppName]] <- subp
  }

  class(epp.pops) <- "eppsubp"

  return(epp.pops)
}


###################
####  Example  ####
###################

## sa.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/South Africa 2014/SouthAfrica_p_2014June10-H-c"
## sa.epp.input <- read.epp.input(sa.path)
## sa.eppd <- read.epp.data(paste(sa.path, ".xml", sep=""))
## sa.eppsubp <- read.epp.subpops(paste(sa.path, ".xml", sep=""))
