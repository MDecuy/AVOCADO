#' @title Avocado algorithm - Single pixel version
#' @name PLUGPhenAnoRFDPLUS
#' @description Calculate the Anomaly and their likelihood values (based on the difference between the pixel values and reference vegetation) for a numeric vector
#' @author  Roberto O. Chavez, Mathieu Decuyper
#' @param x numeric vector. Time series of vegetation index (e.g. NDVI, EVI, NDMI)
#' @param phen Numeric vector with the values of the reference vegetation (all pixels of the reference area are considered)
#' @param dates A date vector. The number of dates must be equal to the number of values of time series.
#' @param h Numeric. Geographic hemisphere to define the starting date of the growing season. h=1 Northern Hemisphere; h=2 Southern Hemisphere.
#' @param anop Numeric vector with the number of values that are in x. For those values the anomalies and likelihoods will be calculated based on the phen. For example a time series has 450 values, so anop=c(1:450).
#' @param rge A vector containing minimum and maximum values of the response variable used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge =c(0,10000)
#' @return vector of length([x]*2) containing all anomalies, followed by their position within the reference frequency distribution (RFD)
#' @seealso \code{\link{PLUGPhenAnoRFDMapPLUS}}
#' @import npphen
#' @importFrom lubridate yday
#' @examples
#' \donttest{
#' #==========================================================================================================
#' # Reading landsat data and getting scene dates 
#' ndmibrick <- rast("~/Dropbox/GITHUB/avocado_testing_site/Pinus_test1_Landsat_NDMI_2000-01-01_2022-12-31.tif")#Load in your raster stack with the dates
#' #Load in the csv file obtained in GEE - beware column with the dates and "date" as header
#' dates.table <- read.csv('~/Dropbox/GITHUB/avocado_testing_site/Pinus_test1_Date_Landsat_NDMI_2000_2022.csv')
#' lan.dates <- as.Date(dates.table$date, format='%d-%m-%Y') #change order to '%Y-%m-%d'if you only get NA. Make sure that in this case there is a column name 'dates' in the csv file

#' #===================================================================================================================
#' #Create the reference curve
#' ref.core.shp <- vect("~/Dropbox/GITHUB/avocado_testing_site/Pinus_RefFor1.shp")
#' plot(ndmibrick[[1]])
#' plot(ref.core.shp,add=T)

#' ref.brick <- crop(ndmibrick,ref.core.shp)
#' ref.brick <- mask(ref.brick,ref.core.shp)
#' plot(ref.brick[[1]])
#' plot(ref.core.shp,add=T)

#' fin <- nrow(ref.brick)*ncol(ref.brick)
#' p2 <- as.numeric(extract(ref.brick,1))
#' d1 <- lan.dates

#' for(i in 2:fin) {
#'   pp <-as.numeric(extract(ref.brick,i))
#'   p2 <-c(p2,pp)
#'   d1 <-c(d1,lan.dates)
#' }
#' 
#' PhenKplot(x=p2,d1,h=1, xlab="DOY",ylab="NDMI",rge=c(0,8000))

#' #=========================================================================
#' # Task 2, calculating ano-prob + dist-regrowth FOR A SINGLE PIXEL (to test out the parameters before running large areas)
#' #Pixel location
#' plot(ndmibrick[[1]])
#' crs.data <- crs(ndmibrick)
#' xy <- rbind(c(4.918970,51.478009))#Drought detection in 2018
#' p <- vect(xy, crs=crs.data)
#' plot(p, add=T)

#' ts.inter <- as.numeric(extract(ndmibrick,p))
#' plot(ts.inter, ylim=c(-2000,8000))
#' #=========================================================================
#' }
#' @export
PLUGPhenAnoRFDPlus <-
  function(x, phen, dates, h, anop, rge) {
    # a.Preparing dataset
    
    if (length(rge) != 2) {
      stop("rge must be a vector of length 2")
    }
    if (rge[1] > rge[2]) {
      stop("rge vector order must be minimum/maximum")
    }
    if (length(dates) != length(x)) {
      stop("N of dates and files do not match")
    }
    
    # ref.min <- min(refp)
    # ref.max <- max(refp)
    ano.min <- min(anop)
    ano.max <- max(anop)
    ano.len <- ano.max - ano.min + 1
    len2 <- 2 * ano.len
    
    if (ano.min >= ano.max) {
      stop("for anop, lower value > upper value")
    }
    
    if (all(is.na(x))) {
      return(rep(NA, len2))
    }
    
    DOY <- lubridate::yday(dates)
    DOY[which(DOY == 366)] <- 365
    D1 <- cbind(DOY, phen)
    D2 <- cbind(DOY[ano.min:ano.max], x[ano.min:ano.max])
    
    if (length(unique(D1[, 2])) < 10 | (nrow(D1) - sum(is.na(D1))) < (0.1 * length(D1))) {
      return(rep(NA, len2))
    }
    
    if (all(is.na(D2[, 2]))) {
      return(rep(NA, len2))
    }
    
    # b. Kernel calculation using the reference period (D1)
    if (h != 1 && h != 2) {
      stop("Invalid h")
    }
    DOGS <- cbind(seq(1, 365), c(seq(185, 365), seq(1, 184)))
    if (h == 2) {
      for (i in 1:nrow(D1)) {
        D1[i, 1] <- DOGS[which(DOGS[, 1] == D1[i, 1], arr.ind = TRUE), 2]
      }
    }
    
    Hmat <- ks::Hpi(na.omit(D1))
    Hmat[1, 2] <- Hmat[2, 1]
    K1 <- ks::kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), 
                  xmax = c(365, rge[2]), gridsize = c(365, 500))
    K1Con <- K1$estimate
    for (j in 1:365) {
      Kdiv <- sum(K1$estimate[j, ])
      ifelse(Kdiv == 0, K1Con[j, ] <- 0, K1Con[j, ] <- K1$estimate[j, ] / sum(K1$estimate[j, ]))
    }
    first.no.NA.DOY <- min(D1[, 1][which(is.na(D1[, 2]) == FALSE)])
    last.no.NA.DOY <- max(D1[, 1][which(is.na(D1[, 2]) == FALSE)])
    MAXY <- apply(K1Con, 1, max)
    for (i in 1:365) {
      n.select <- which(K1Con[i, ] == MAXY[i], arr.ind = TRUE)
      if (length(n.select) > 1) {
        n <- n.select[1]
        MAXY[i] <- NA
      }
      if (length(n.select) == 1) {
        n <- n.select
        MAXY[i] <- median(K1$eval.points[[2]][n])
      }
      if (i < first.no.NA.DOY) {
        MAXY[i] <- NA
      }
      if (i > last.no.NA.DOY) {
        MAXY[i] <- NA
      }
    }
    
    # c. Calculating cumulative bivariate density distribution
    
    h2d <- list()
    h2d$x <- seq(1, 365)
    h2d$y <- seq(rge[1], rge[2], len = 500)
    h2d$density <- K1Con/sum(K1Con)
    uniqueVals <- rev(unique(sort(h2d$density)))
    cumRFDs <- cumsum(uniqueVals)
    names(cumRFDs) <- uniqueVals
    h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
    h2d$cumDensity[] <- cumRFDs[as.character(h2d$density)]
    na.sta <- first.no.NA.DOY - 1
    na.end <- last.no.NA.DOY + 1
    if(na.sta>=1) {h2d$cumDensity[1:na.sta,] <- NA}
    if(na.end<=365) {h2d$cumDensity[na.end:365,] <- NA}
    
    # c. Anomaly calculation (D2)
    
    if (h == 2) {
      for (i in 1:nrow(D2)) {
        D2[i, 1] <- DOGS[which(DOGS[, 1] == D2[i, 1], arr.ind = TRUE), 2]
      }
    }
    
    # d. Calculating anomalies AND their likelihoods based on D2
    
    # d.1 Anomalies
    Anoma <- rep(NA, ano.len)
    for (i in 1:nrow(D2)) {
      Anoma[i] <- as.integer(D2[i, 2] - MAXY[D2[i, 1]])
    }
    Anoma[1:ano.len]
    
    # d.2 Likelihood
    rowAnom <- matrix(NA, nrow = nrow(D2), ncol = 500)
    for (i in 1:nrow(D2)) {
      rowAnom[i, ] <- abs(h2d$y - D2[i, 2])
    }
    rowAnom2 <- unlist(apply(rowAnom, 1, function(x) {
      if (all(is.na(x))) {
        NA
      }
      else {
        which.min(x)
      }
    }))
    AnomRFD <- rep(NA, nrow(D2))
    for (i in 1:nrow(D2)) {
      AnomRFD[i] <- h2d$cumDensity[D2[i, 1], rowAnom2[i]]
    }
    AnomRFDPerc <- round(100 * (AnomRFD))
    plot(AnomRFDPerc[1:ano.len], xlab = "Time", ylab = "RFD (%)", font.lab = 2)
    abline(h = 90, col = "red")
    abline(h = 95, col = "red")
    abline(h = 99, col = "red")
    AnomRFDPerc[1:ano.len]
    # Two results in a single numeric vector
    AnomPlusRFD <- c(Anoma[1:ano.len], AnomRFDPerc[1:ano.len])
    AnomPlusRFD[1:len2]
  }

#' @title Avocado algorithm - Wall-to-wall map version
#' @name PLUGPhenAnoRFDMapPLUS
#' @description Calculate the Anomaly and their likelihood values (based on the difference between the pixel values and reference vegetation) for a RasterStack
#' @author  Roberto O. Chavez, Mathieu Decuyper
#' @param s RasterStack time series of vegetation index (e.g. NDVI, EVI, NDMI)
#' @param anop Numeric vector with the number of values that are in x. For those values the anomalies and likelihoods will be calculated based on the phen. For example a time series has 450 values, so anop=c(1:450).
#' @param phen Numeric vector with the values of the reference vegetation
#' @param h Numeric. Geographic hemisphere to define the starting date of the growing season. h=1 Northern Hemisphere; h=2 Southern Hemisphere.
#' @param nCluster Numeric. Number of CPU's to be used for the job.
#' @param outname Character vector with the output path and filename with extension or only the filename and extension if work directory was set. More information: See writeRaster
#' @param datatype Character vector that determines the interpretation of values written to disk. More information: See \code{\link{writeFormats}}
#' @param rge A vector containing minimum and maximum values of the response variable used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge =c(0,10000)
#' @return RasterStack of nlayers([s]stack*2) containing all anomalies, followed by their likelihoods

#' @import npphen
#' @import terra
#' @importFrom lubridate yday 
#' @import RColorBrewer
#'
#' @seealso \code{\link{PLUGPhenAnoRFDPLUS}}
#' @examples
#' \donttest{
#' source("PLUGPhenAnoRFDMapPLUS_20190902.R") # Load in the mapping function
#' dates <- lan.dates # The dates from your time-series brick (x)
#' PLUGPhenAnoRFDMapPLUS(s = MDD, dates = dates, h = 1, phen = phen, anop = c(1:n), nCluster = 1, outname = "YourDirectory/Filename.tif", format = "GTiff", datatype = "INT2S", rge = c(0, 10000))
#' }
#' @export
PLUGPhenAnoRFDMapPLUS <-
  function(s, phen, dates, h, anop, nCluster, outname, datatype, rge) {
    ff <- function(x) {
      # a.Preparing dataset
      
      if (length(rge) != 2) {
        stop("rge must be a vector of length 2")
      }
      if (rge[1] > rge[2]) {
        stop("rge vector order must be minimum/maximum")
      }
      if (length(dates) != length(x)) {
        stop("N of dates and files do not match")
      }
      if (length(x) < length(anop)) {
        stop("Inconsistent anop. Argument anop can't be grater than length(x)")
      }
      
      # ref.min <- min(refp)
      # ref.max <- max(refp)
      ano.min <- min(anop)
      ano.max <- max(anop)
      ano.len <- ano.max - ano.min + 1
      len2 <- 2 * ano.len
      
      if (ano.min >= ano.max) {
        stop("for anop, lower value > upper value")
      }
      
      if (all(is.na(x))) {
        return(rep(NA, len2))
      }
      
      DOY <- lubridate::yday(dates)
      DOY[which(DOY == 366)] <- 365
      D1 <- cbind(DOY, phen)
      D2 <- cbind(DOY[ano.min:ano.max], x[ano.min:ano.max])
      
      if (length(unique(D1[, 2])) < 10 | (nrow(D1) - sum(is.na(D1))) < (0.1 * length(D1))) {
        return(rep(NA, len2))
      }
      
      if (all(is.na(D2[, 2]))) {
        return(rep(NA, len2))
      }
      
      # b. Kernel calculation using the reference period (D1)
      if (h != 1 && h != 2) {
        stop("Invalid h")
      }
      DOGS <- cbind(seq(1, 365), c(seq(185, 365), seq(1, 184)))
      if (h == 2) {
        for (i in 1:nrow(D1)) {
          D1[i, 1] <- DOGS[which(DOGS[, 1] == D1[i, 1], arr.ind = TRUE), 2]
        }
      }
      
      Hmat <- ks::Hpi(na.omit(D1))
      Hmat[1, 2] <- Hmat[2, 1]
      K1 <- ks::kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), xmax = c(365, rge[2]), gridsize = c(365, 500))
      K1Con <- K1$estimate
      for (j in 1:365) {
        Kdiv <- sum(K1$estimate[j, ])
        ifelse(Kdiv == 0, K1Con[j, ] <- 0, K1Con[j, ] <- K1$estimate[j, ] / sum(K1$estimate[j, ]))
      }
      
      first.no.NA.DOY <- min(D1[,1][which(is.na(D1[,2])==FALSE)])
      last.no.NA.DOY <- max(D1[,1][which(is.na(D1[,2])==FALSE)])
      
      MAXY <- apply(K1Con, 1, max)
      for (i in 1:365) {
        n.select <- which(K1Con[i, ] == MAXY[i], arr.ind = TRUE)
        if (length(n.select) > 1) {
          n <- n.select[1]
          MAXY[i] <- NA
        }
        if (length(n.select) == 1) {
          n <- n.select
          MAXY[i] <- median(K1$eval.points[[2]][n])
        }
        if (i < first.no.NA.DOY) {MAXY[i] <- NA}
        if (i > last.no.NA.DOY) {MAXY[i] <- NA}
      }
      
      # c. Calculating cumulative bivariate density distribution
      
      h2d <- list()
      h2d$x <- seq(1, 365)
      h2d$y <- seq(rge[1], rge[2], len = 500)
      h2d$density <- K1Con / sum(K1Con)
      uniqueVals <- rev(unique(sort(h2d$density)))
      cumRFDs <- cumsum(uniqueVals)
      names(cumRFDs) <- uniqueVals
      h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
      h2d$cumDensity[] <- cumRFDs[as.character(h2d$density)]
      na.sta <- first.no.NA.DOY-1
      na.end <- last.no.NA.DOY+1
      if(na.sta>=1) {h2d$cumDensity[1:na.sta,] <- NA}
      if(na.end<=365) {h2d$cumDensity[na.end:365,] <- NA}
      
      # d. Calculating anomalies AND their RFD based on D2
      
      if (h == 2) {
        for (i in 1:nrow(D2)) {
          D2[i, 1] <- DOGS[which(DOGS[, 1] == D2[i, 1], arr.ind = TRUE), 2]
        }
      }
      
      # d.1 Anomalies
      Anoma <- rep(NA, ano.len)
      for (i in 1:nrow(D2)) {
        Anoma[i] <- as.integer(D2[i, 2] - MAXY[D2[i, 1]])
      }
      Anoma[1:ano.len]
      
      # d.2 RFD
      rowAnom <- matrix(NA, nrow = nrow(D2), ncol = 500)
      for (i in 1:nrow(D2)) {
        rowAnom[i, ] <- abs(h2d$y - D2[i, 2])
      }
      rowAnom2 <- unlist(apply(rowAnom, 1, function(x) {
        if (all(is.na(x))) {
          NA
        } else {
          which.min(x)
        }
      }))
      AnomRFD <- rep(NA, nrow(D2))
      for (i in 1:nrow(D2)) {
        AnomRFD[i] <- h2d$cumDensity[D2[i, 1], rowAnom2[i]]
      }
      AnomRFDPerc <- round(100 * (AnomRFD))
      AnomRFDPerc[1:ano.len]
      # Two results in a single numeric vector
      AnomPLUSRFD <- c(Anoma[1:ano.len], AnomRFDPerc[1:ano.len])
      AnomPLUSRFD[1:len2]
    }
    
    #----------------------------------------------------------------------------------------
    # cluster processing
    app(s, fun = ff, filename = outname, cores = nCluster, overwrite = T, wopt = list(datatype = datatype))
  }
