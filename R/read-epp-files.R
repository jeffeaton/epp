
###########################################################################
####  Function to read inputs from Spectrum to EPP (.ep1, .ep3, .ep4)  ####
###########################################################################

read_epp_input <- function(pjnz){

  ## ep1
  ep1file <- grep(".ep1", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  con <- unz(pjnz, ep1file)
  ep1 <- scan(con, "character", sep="\n")
  close(con)

  country.idx <- grep("COUNTRY", ep1)
  firstprojyr.idx <-  which(sapply(ep1, substr, 1, 11) == "FIRSTPROJYR")
  lastprojyr.idx <-  which(sapply(ep1, substr, 1, 10) == "LASTPROJYR")
  popstart.idx <- grep("POPSTART", ep1)+1
  popend.idx <- grep("POPEND", ep1)-1

  country <- as.character(read.csv(text=ep1[country.idx], header=FALSE, as.is=TRUE)[2])
  country.code <- as.integer(read.csv(text=ep1[country.idx], header=FALSE)[3])

  start.year <- as.integer(read.csv(text=ep1[firstprojyr.idx], header=FALSE)[2])
  stop.year <- as.integer(read.csv(text=ep1[lastprojyr.idx], header=FALSE)[2])
  epp.pop <- setNames(read.csv(text=ep1[popstart.idx:popend.idx], header=FALSE, as.is=TRUE),
                      c("year", "pop15to49", "pop15", "pop50", "netmigr"))

  ## ep4
  ep4file <- grep(".ep4", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  con <- unz(pjnz, ep4file)
  ep4 <- scan(con, "character", sep="\n")
  close(con)

  cd4lim.idx <- which(sapply(ep4, substr, 1, 12) == "CD4LOWLIMITS")
  lambda.idx <- which(sapply(ep4, substr, 1, 6) == "LAMBDA")
  cd4init.idx <- which(sapply(ep4, substr, 1, 13) == "NEWINFECTSCD4")
  mu.idx <- which(sapply(ep4, substr, 1, 3) == "MU_")
  alpha1.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA1")
  alpha2.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA2")
  alpha3.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA3")
  infectreduc.idx <- which(sapply(ep4, substr, 1, 11) == "INFECTREDUC")
  artstart.idx <- grep("ARTSTART", ep4)+1
  artend.idx <- grep("ARTEND", ep4)-1

  DS <- 7 # disease stages

  cd4lim <- as.integer(read.csv(text=ep4[cd4lim.idx], header=FALSE)[-1][1:DS])
  cd4init <- as.matrix(read.csv(text=ep4[cd4init.idx], header=FALSE, row.names=1)[,1:DS])
  lambda <- as.matrix(read.csv(text=ep4[lambda.idx], header=FALSE, row.names=1)[,1:(DS-1)])
  mu <- as.matrix(read.csv(text=ep4[mu.idx], header=FALSE, row.names=1)[,1:DS])
  alpha1 <- as.matrix(read.csv(text=ep4[alpha1.idx], header=FALSE, row.names=1)[,1:DS])
  alpha2 <- as.matrix(read.csv(text=ep4[alpha2.idx], header=FALSE, row.names=1)[,1:DS])
  alpha3 <- as.matrix(read.csv(text=ep4[alpha3.idx], header=FALSE, row.names=1)[,1:DS])
  infectreduc <- as.numeric(read.csv(text=ep4[infectreduc.idx], header=FALSE)[2])

  epp.art <- setNames(read.csv(text=ep4[artstart.idx:artend.idx], header=FALSE, as.is=TRUE),
                      c("year", "m.isperc", "m.val", "f.isperc", "f.val", "cd4thresh", "m.perc50plus", "f.perc50plus", "perc50plus", "1stto2ndline"))

  specpop.idx <- grep("SPECPOP", ep4)
  if(length(specpop.idx)){
    art.specpop <- setNames(read.csv(text=ep4[specpop.idx], header=FALSE,
                                     colClasses=c("NULL", "character", "numeric", "integer"))[,1:3],
                            c("specpop", "percelig", "yearelig"))
    art.specpop$percelig <- art.specpop$percelig/100
  } else
    art.specpop <- data.frame(specpop=character(), percelig=numeric(), yearelig=integer())

  cd4median.start.idx <- which(ep4 == "CD4MEDIAN_START")+1
  cd4median.end.idx <- which(ep4 == "CD4MEDIAN_END")-1
  if(length(cd4median.start.idx) > 0)
    epp.pop$cd4median <- read.csv(text=ep4[cd4median.start.idx:cd4median.end.idx], header=FALSE, colClasses=c("NULL", "numeric"))[[1]]
  else
    epp.pop$cd4median <- 0

  hivp15yr.start.idx <- which(ep4 == "HIVPOS_15YEAROLDS")+1
  hivp15yr.end.idx <- which(ep4 == "HIVPOS_15YEAROLDS_END")-1
  if(length(hivp15yr.start.idx) > 0)
    epp.pop$hivp15yr <- read.csv(text=ep4[hivp15yr.start.idx:hivp15yr.end.idx], header=FALSE, colClasses=c("NULL", "numeric"))[[1]]
  else
    epp.pop$hivp15yr <- 0

  art15yr.start.idx <- which(ep4 == "HIVPOS_15YEAROLDSART")+1
  art15yr.end.idx <- which(ep4 == "HIVPOS_15YEAROLDSART_END")-1
  if(length(art15yr.start.idx) > 0)
    epp.art$art15yr <- read.csv(text=ep4[art15yr.start.idx:art15yr.end.idx], header=FALSE, colClasses=c("NULL", "numeric"))[[1]]
  else
    epp.art$art15yr <- 0

  artdropout.start.idx <- which(ep4 == "ARTDROPOUTRATE")+1
  artdropout.end.idx <- which(ep4 == "ARTDROPOUTRATE_END")-1
  if(length(artdropout.start.idx) > 0)
    epp.art$artdropout <- read.csv(text=ep4[artdropout.start.idx:artdropout.end.idx], header=FALSE, colClasses=c("NULL", "numeric"))[[1]]
  else
    epp.art$artdropout <- 0

  hivp15yr.cd4dist.idx <- which(ep4 == "HIVPOS15_CD4")+1
  if(length(hivp15yr.cd4dist.idx) > 0)
    hivp15yr.cd4dist <- as.numeric(read.csv(text=ep4[hivp15yr.cd4dist.idx], header=FALSE))
  else
    hivp15yr.cd4dist <- rep(0, length(cd4lim))

  art15yr.cd4dist.idx <- which(ep4 == "HIVPOS15ART_CD4")+1
  if(length(art15yr.cd4dist.idx) > 0)
    art15yr.cd4dist <- as.numeric(read.csv(text=ep4[art15yr.cd4dist.idx], header=FALSE))
  else
    art15yr.cd4dist <- rep(0, length(cd4lim))


  ## Between 2014 and 2015, Spectrum switched from passing CD4 stage duration to passing
  ## CD4 stage progression rate. This change is unmarked in the .ep4 file.  Guess which
  ## is correct based on the first stage duration, assuming first stage should be > 1 year.
  ## THIS IS NOT IDIOT PROOF!!!!

  if(mean(lambda[,1]) > 1)
    lambda <- lambda
  else
    lambda <- 1/lambda

  ## XML (for epidemic start year)

  r <- get_eppxml_workset(pjnz)

  ## Note: tag "epidemicStartYrVarR" doesn't appear to change...
  ## Use epidemic start from first EPP subpopulation fit
  projsets <- xml_find_all(r, ".//object")
  projsets <- projsets[which(xml_attr(projsets, "class") == "epp2011.core.sets.ProjectionSet")]
  eppSet <- xml_children(projsets[1])
  epidemic.start <- as.integer(xml_double(eppSet[which(xml_attr(eppSet, "property") == "priorT0vr")]))

  eppin <- list(start.year       = start.year,
                stop.year        = stop.year,
                epidemic.start   = epidemic.start,
                epp.pop          = epp.pop,
                cd4lowlim        = cd4lim,
                cd4initperc      = cd4init,
                cd4stage.dur     = lambda,
                cd4mort          = mu,
                artmort.less6mos = alpha1,
                artmort.6to12mos = alpha2,
                artmort.after1yr = alpha3,
                infectreduc      = infectreduc,
                epp.art          = epp.art,
                art.specpop      = art.specpop,
                hivp15yr.cd4dist = hivp15yr.cd4dist,
                art15yr.cd4dist  = art15yr.cd4dist)
  class(eppin) <- "eppin"
  attr(eppin, "country") <- country
  attr(eppin, "country.code") <- country.code

  return(eppin)
}

#' Return workset with names from EPP .xml file
#'
#' @param x `xml_node` object with tag 'array'
#' @import xml2
get_eppxml_workset <- function(pjnz){

  xmlfile <- grep(".xml", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  if(!length(xmlfile)){
    warning(paste0("No EPP .xml file found for ", basename(pjnz)))
    return(NULL)
  }

  con <- unz(pjnz, xmlfile)
  epp.xml <- read_xml(con)

  r <- xml_children(xml_child(epp.xml))
  names(r) <- xml_attr(r, "property")

  return(r)
}

#' Parse java array from EPP .xml file
#'
#' @param x `xml_node` object with tag 'array'
#' @import xml2
.parse_array <- function(x){

  a_length <- as.integer(xml_attr(x, "length"))
  a_mode <- switch(xml_attr(x, "class"),
                   int = "integer",
                   double = "numeric",
                   boolean = "logical",
                   java.lang.String = "character")
  arr <- vector(a_mode, a_length)

  ## Note: default declaration will be 0 for a_mode \in {integer, numeric},
  ##       FALSE for a_mode = logical, and "" for a_mode = character. This
  ##       matches the java array defaults in which indices are omitted for
  ##       these values.

  elem <- xml_children(x)
  idx <- as.integer(xml_attr(elem, "index")) + 1L # java is 0-based
  fn <- switch(a_mode,
               integer = xml_integer,
               numeric = xml_double,
               logical = function(xx) as.logical(xml_text(elem)),
               character = xml_text)
  arr[idx] <- fn(elem)

  return(arr)
}


#' Parse java matrix from EPP .xml file
#'
#' @param x `xml_node` object with tag 'array' and type '[D' or '[I'
#' @import xml2
.parse_matrix <- function(x){

  if(!xml_attr(x, "class") %in% c("[D", "[I"))
    stop("Tried to invoke .parse_matrix() on array node not of class '[D' or '[I]'.")

  m_rows <- as.integer(xml_attr(x, "length"))
  rows <- xml_children(x)
  idx <- as.integer(xml_attr(rows, "index")) + 1L
  rows <- lapply(xml_find_first(rows, "array"), .parse_array)

  m_cols <- length(rows[[1]])
  if(!all(sapply(rows, length) == m_cols))
    stop("Not all rows are same length. This might be a ragged array.")

  mat <- matrix(nrow = m_rows, ncol=m_cols)
  for(i in seq_along(rows))
    mat[idx[i],] <- rows[[i]]

  return(mat)
}


#' Read EPP fitting surveillance data
#'
#' Reads the HIV prevalence from sentinel surveillance and household surveillance
#' data from the EPP .xml file within a PJNZ.
#'
#' @param pjnz file path to Spectrum PJNZ file.
#'
#' @details
#' EPP projection sets are identified in the .xml file by searching the XML tree
#' for tag "object", and then selecting objects with "class" attribute equal to
#' "epp2011.core.sets.ProjectionSet".
#'
#' @import xml2
#' @export
read_epp_data <- function(pjnz){

  r <- get_eppxml_workset(pjnz)
  if(is.null(r))
    return(NULL)

  country <- xml_text(r[["worksetCountry"]])
  country_code <- xml_integer(r[["countryCode"]])

  epp.data <- list() # declare list to store output
  attr(epp.data, "country") <- country
  attr(epp.data, "country_code") <- country_code

  ## ANC/HSS input mode appears only defined in second set if updated, so define global version.
  ## Defaults to "HSS mdoe", which is no ANC-RT data.
  input_mode <- "HSS"

  obj <- xml_find_all(r, ".//object")
  projsets <- obj[which(xml_attr(obj, "class") == "epp2011.core.sets.ProjectionSet")]

  for(eppSet in projsets){

    projset_id <- as.integer(gsub("[^0-9]", "", xml_attr(eppSet, "id")))

    eppSet <- xml_children(eppSet)
    names(eppSet) <- xml_attr(eppSet, "property")

    eppName <- xml_text(eppSet[["name"]])

    ##  ANC data  ##

    if(exists("siteNames", eppSet)) {
      siteNames <- .parse_array(xml_find_first(eppSet[["siteNames"]], "array"))
      nsites <- length(siteNames)

      ## ANC site used
      anc.used <- .parse_array(xml_find_first(eppSet[["siteSelected"]], "array"))

      ## ANC prevalence
      anc.prev <- .parse_matrix(xml_find_first(eppSet[["survData"]], "array"))
      dimnames(anc.prev) <- list(site=siteNames, year=1985+0:(ncol(anc.prev)-1))
      anc.prev[anc.prev == -1] <- NA
      anc.prev <- anc.prev/100

      ## ANC sample sizes
      anc.n <- .parse_matrix(xml_find_first(eppSet[["survSampleSizes"]], "array"))
      dimnames(anc.n) <- list(site=siteNames, year=1985+0:(ncol(anc.n)-1))
      anc.n[anc.n == -1] <- NA

      ## ANC-RT site level

      if(length(eppSet[["dataInputMode"]]) &&
         length(xml_find_first(eppSet[["dataInputMode"]], ".//string")))
        input_mode <- xml_text(xml_find_first(eppSet[["dataInputMode"]], ".//string"))

      if(length(eppSet[["PMTCTData"]]) && input_mode == "ANC"){

        ancrtsite.prev <- .parse_matrix(xml_find_first(eppSet[["PMTCTData"]], "array"))
        dimnames(ancrtsite.prev) <- list(site=siteNames, year=1985+0:(ncol(ancrtsite.prev)-1))
        ancrtsite.prev[ancrtsite.prev == -1] <- NA
        ancrtsite.prev <- ancrtsite.prev/100

        ancrtsite.n <- .parse_matrix(xml_find_first(eppSet[["PMTCTSiteSampleSizes"]], "array"))
        dimnames(ancrtsite.n) <- list(site=siteNames, year=1985+0:(ncol(ancrtsite.n)-1))
        ancrtsite.n[ancrtsite.n == -1] <- NA

      } else {
        ancrtsite.prev <- NULL
        ancrtsite.n <- NULL
      }
    } else {
      anc.used <- NULL
      anc.prev <- NULL
      anc.n <- NULL
      ancrtsite.prev <- NULL
      ancrtsite.n <- NULL
    }


    ## ANC-RT census level
    if(length(eppSet[["censusPMTCTSurvData"]]) && input_mode == "ANC"){

      ancrtcens.prev <- .parse_array(xml_find_first(eppSet[["censusPMTCTSurvData"]], "array"))
      names(ancrtcens.prev) <- 1985+0:(length(ancrtcens.prev)-1)
      ancrtcens.prev[ancrtcens.prev == -1] <- NA
      ancrtcens.prev <- ancrtcens.prev/100

      ancrtcens.n <- .parse_array(xml_find_first(eppSet[["censusPMTCTSampleSizes"]], "array"))
      ancrtcens.n[ancrtcens.n == -1] <- NA

      ancrtcens <- data.frame(year=as.integer(names(ancrtcens.prev)),
                              prev=ancrtcens.prev, n=ancrtcens.n)
      ancrtcens <- subset(ancrtcens, !is.na(prev) | !is.na(n))
    } else {
      ancrtcens <- NULL
    }

    ##  HH surveys  ##
    if("surveyData" %in% names(eppSet)) {
      ## warning(paste("File", basename(pjnz), "uses new EPP data structure for HH survey data. Parsers are not yet implemented"))

      svys <- xml_children(xml_children(eppSet[["surveyData"]]))
      
      parse_survey <- function(ns){
        ns <- xml_children(xml_child(ns))
        attrs <- lapply(ns, xml_attrs)

        val <- lapply(lapply(ns[sapply(attrs,  "%in%", x = "getField")], xml_children), xml_text)

        v2 <- lapply(ns[!sapply(attrs,  "%in%", x = "getField")], xml_children)
        v2 <- lapply(v2, lapply, xml_children)
        v2 <- lapply(v2, lapply, xml_text)

        val <- c(val, unlist(v2, FALSE))
        val <- setNames(sapply(val, "[", 2), sapply(val, "[", 1))

        cols <- c("name" = "name",
                  "year" = "year",
                  "used",
                  "n",
                  "surveyHIV" = "prev",
                  "surveyStandardError" = "se",
                  "incidence" = "incid",
                  "standardError" = "incid_se",
                  "prev_incid_corr",
                  "incidence_cohort",
                  "usingIncidenceData" = "incid_used")

        names(val) <- cols[names(val)]
        val[setdiff(cols, names(val))] <- NA
        val["used"] <- TRUE

        ## Default constant DEFAULT_SURVEY_YEAR = 2009 used to initialise variable
        ## in Java EPP code (Robert Puckett; 12 Sep 2020).
        if(is.na(val[["year"]]))
          val[["year"]] <- 2009
        
        val[cols]
      }
      
      hhs <- lapply(svys, parse_survey)
      hhs <- as.data.frame(do.call(rbind, hhs))
      hhs <- as.data.frame(lapply(hhs, type.convert, as.is = TRUE))

      hhs$prev <- hhs$prev / 100
      hhs$se <- hhs$se / 100
      hhs$incid <- hhs$incid / 100
      hhs$incid_se <- hhs$incid_se / 100

    } else if("surveyYears" %in% names(eppSet)) {
      hhs <- data.frame(year = .parse_array(xml_find_first(eppSet[["surveyYears"]], "array")),
                        prev = .parse_array(xml_find_first(eppSet[["surveyHIV"]], "array"))/100,
                        se = .parse_array(xml_find_first(eppSet[["surveyStandardError"]], "array"))/100,
                        n = NA,
                        used = .parse_array(xml_find_first(eppSet[["surveyIsUsed"]], "array")))

      if(!is.null(eppSet[["inputInc"]])){
        hhs$incid <- .parse_array(xml_find_first(eppSet[["inputInc"]], "array")) / 100
        hhs$incid_se <- .parse_array(xml_find_first(eppSet[["inputIncSE"]], "array")) / 100
        hhs$prev_incid_corr <- .parse_array(xml_find_first(eppSet[["inputIncPrevCorr"]], "array"))
        hhs$incid_cohort <- .parse_array(xml_find_first(eppSet[["incIsCohort"]], "array"))
      } else {
        hhs$incid <- -1
        hhs$incid_se <- NA
        hhs$prev_incid_corr <- NA
        hhs$incid_cohort <- NA
      }

      hhs <- subset(hhs, prev > 0 | used | se != 0.01)
      if(nrow(hhs))
        hhs[hhs$incid < 0, c("incid", "incid_se", "prev_incid_corr", "incid_cohort")] <- NA
    } else {
      hhs <- data.frame()
    }
    

    epp.data[[eppName]] <- list(country=country,
                                region=eppName,
                                projset_id = projset_id,
                                anc.used=anc.used,
                                anc.prev=anc.prev,
                                anc.n=anc.n,
                                ancrtsite.prev=ancrtsite.prev,
                                ancrtsite.n=ancrtsite.n,
                                ancrtcens=ancrtcens,
                                hhs=hhs)
  }

  class(epp.data) <- "eppd"

  return(epp.data)
}


#' Read EPP subpopulation configuration
#'
#' Reads the subpopulation configuration and population sizes from the EPP .xml
#' file within a PJNZ.
#'
#' @param pjnz file path to Spectrum PJNZ file.
#'
#' @details
#' EPP projection sets are identified in the .xml file by searching the XML tree
#' for tag "object", and then selecting objects with "class" attribute equal to
#' "epp2011.core.sets.ProjectionSet".
#'
#' @import xml2
#' @export

read_epp_subpops <- function(pjnz){

  r <- get_eppxml_workset(pjnz)

  epp.pops <- list() # declare list to store output
  attr(epp.pops, "country") <- xml_text(r[["worksetCountry"]])
  attr(epp.pops, "country_code") <- xml_integer(r[["countryCode"]])

  epidemicType <- tolower(xml_text(xml_find_first(r[["epidemicType"]], ".//string")))
  attr(epp.pops, "epidemicType") <- epidemicType

  startyear <- xml_integer(r[["worksetStartYear"]])
  endyear <- xml_integer(r[["worksetEndYear"]])
  popbaseyear <- xml_integer(r[["worksetPopBaseYear"]])

  epp.pops$total <-  data.frame(year = startyear:endyear,
                                pop15to49 = .parse_array(xml_find_first(r[["pop15to49"]], "array")),
                                pop15     = .parse_array(xml_find_first(r[["pop15"]], "array")),
                                pop50     = .parse_array(xml_find_first(r[["pop50"]], "array")),
                                netmigr   = .parse_array(xml_find_first(r[["netMigration"]], "array")))

  epp.pops$subpops <- list()

  obj <- xml_find_all(r, ".//object")
  projsets <- obj[which(xml_attr(obj, "class") == "epp2011.core.sets.ProjectionSet")]

  for(eppSet in projsets){

    projset_id <- as.integer(gsub("[^0-9]", "", xml_attr(eppSet, "id")))

    eppSet <- xml_children(eppSet)
    names(eppSet) <- xml_attr(eppSet, "property")

    eppName <- xml_text(eppSet[["name"]])

    subp <- data.frame(year = startyear:endyear,
                       pop15to49 = .parse_array(xml_find_first(eppSet[["pop15to49"]], "array")),
                       pop15     = .parse_array(xml_find_first(eppSet[["pop15"]], "array")),
                       pop50     = .parse_array(xml_find_first(eppSet[["pop50"]], "array")),
                       netmigr   = .parse_array(xml_find_first(eppSet[["netMigration"]], "array")))

    attr(subp, "projset_id") <- projset_id
    attr(subp, "epidemic.start") <- as.integer(xml_double(eppSet[["priorT0vr"]]))

    if(epidemicType == "concentrated"){

      subpop <- .parse_array(xml_find_first(eppSet[["specSubPop"]], "array"))
      names(subpop) <- c("low_risk", "msm", "msw", "fsw", "clients",
                         "idu", "prisoners", "transgender", "anc")

      if(length(eppSet[["percentageMale"]]))
        percent_male <- xml_double(eppSet[["percentageMale"]])/100
      else
        percent_male <- 0.0

      turnover <- as.logical(length(eppSet[["turnedOver"]])) &&
        as.logical(xml_text(eppSet[["turnedOver"]]))

      if(turnover){
        duration <- xml_double(eppSet[["duration"]])
        assign_id <- xml_attr(xml_find_first(eppSet[["groupToAssignTo"]], ".//object"), "id")
        assign_id <- as.integer(gsub("[^0-9]", "", assign_id))
        assignmentType <- switch(xml_text(xml_find_first(eppSet[["assignmentMethod"]], ".//string")),
                                 ASSIGN_REPLACE_PREVALENCE = "replace",
                                 ASSIGN_ADD_PREVALENCE = "add")

      } else {
        duration <- NA
        assign_id <- NA
        assignmentType <- NA
      }

      attr(subp, "subpop") <-  names(which(subpop))
      attr(subp, "percent_male") <- percent_male
      attr(subp, "turnover") <- turnover
      attr(subp, "duration") <- duration
      attr(subp, "assign_id") <- assign_id
      attr(subp, "assignmentType") <- assignmentType
    }

    epp.pops$subpops[[eppName]] <- subp
  }

  projset_ids <- sapply(epp.pops$subpops, attr, "projset_id")
  assign_ids <- sapply(epp.pops$subpops, attr, "assign_id")
  assign_name <- names(projset_ids)[match(assign_ids, projset_ids)]

  epp.pops$subpops <- Map("attr<-", epp.pops$subpops, "assign_name", assign_name)

  class(epp.pops) <- "eppsubp"

  return(epp.pops)
}


#####################################################
####  Read EPP prevalence and incidence outputs  ####
#####################################################

read_spt <- function(pjnz){

  sptfile <- grep("\\.SPT$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  con <- unz(pjnz, sptfile)
  spt <- scan(con, "character", sep="\n")
  close(con)

  break.rows <- which(spt == "==")
  sub.break.rows <- which(spt == "=")
  n.years <- sub.break.rows[2] - break.rows[1] - 3
  regions <- sapply(strsplit(as.character(spt[break.rows+1]), "\\\\"), function(x) strsplit(x[2], ":")[[1]][1])[-length(break.rows)]
  regions[is.na(regions)] <- "National"

  out <- lapply(sub.break.rows[-1],
                function(idx){
                  dat <- spt[idx-n.years:1]
                  mat <- data.frame(t(sapply(strsplit(dat, ","), as.numeric)), row.names=1)
                  mat[,1:2] <- mat[,1:2]/100
                  names(mat) <- c("prev", "incid", "pop")[1:ncol(mat)]
                  return(mat)
                })

  names(out) <- regions
  return(out)
}

read_spu <- function(pjnz){

  spufile <- grep("\\.SPU$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  spu <- read.csv(unz(pjnz, spufile), header=FALSE)

  n.resamp <- as.numeric(strsplit(as.character(spu[1,1]), " ")[[1]][3])
  break.rows <- which(spu[,1] == "==")
  n.years <- break.rows[2] - break.rows[1] - 2
  count <- sapply(strsplit(as.character(spu[break.rows[-1]-(n.years+1),1]), " "), function(x) as.numeric(x[2]))

  years <- as.numeric(as.character(spu[break.rows[1]-n.years:1,1]))

  incid <- sapply(break.rows[-1], function(idx) spu[idx-n.years:1,3])[,rep(1:length(count), count)]/100
  prev <- sapply(break.rows[-1], function(idx) spu[idx-n.years:1,2])[,rep(1:length(count), count)]/100

  rownames(incid) <- years
  rownames(prev) <- years

  return(list("incid"=incid, "prev"=prev))
}




###################
####  Example  ####
###################

## sa.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/South Africa 2014/SouthAfrica_p_2014June10-H-c"
## sa.epp.input <- read.epp.input(sa.path)
## sa.eppd <- read.epp.data(paste(sa.path, ".xml", sep=""))
## sa.eppsubp <- read.epp.subpops(paste(sa.path, ".xml", sep=""))
