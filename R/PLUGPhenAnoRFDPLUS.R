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
#' @import raster
#' @import RColorBrewer
#' @import rgdal
#' @import bfastSpatial
#' @import lubridate
#' @import rts
#' @examples
#' \donttest{
#' #' # load MDD datasets
#' MDD <- stack(system.file("extdata", "MDD_NDMI_1990_2020.gri", package = "AVOCADO"))
#' load(system.file("extdata", "MDDref.RData", package = "AVOCADO"))
#' ############ Section 1: Reference vegetatation curve ###################
#' lan.info <- getSceneinfo(names(MDD))
#' lan.dates <- as.Date(lan.info$date)
#' MDDref <- readOGR(dsn = path.expand("YourDirectory"), layer = "MDDref")
#' # MDDref <- spTransform(MDDref, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #in case of different coordinate systems
#' ref.ext <- extent(MDDref)
#' ref.brick <- crop(MDD, ref.ext)
#' fin <- nrow(ref.brick) * ncol(ref.brick)
#' phen <- extract(ref.brick, 1)
#' d1 <- lan.dates
#' }

#' for(i in 2:fin) {
#'  pp <-extract(ref.brick,i)
#'  phen <-c(phen,pp)
#'  d1 <-c(d1,lan.dates)
#' }
#' PhenKplot(phen,d1,h=1,nGS=365, xlab="DOY",ylab="NDMI",rge=c(0,10000))
#'
#' Alternative: to take the median of all the pixel observations per date (Be aware: smoother curve but loss of natural variation)
#' ts.med <- NULL
#' for (i in 1:nlayers(ref.brick)) {
#'   if (all(is.na(getValues(ref.brick[[i]])))) {
#'     median <- NA
#'   } else {
#'     nn <- ncell(ref.brick[[i]])
#'     ss <- sampleRegular(ref.brick[[i]], size=nn)
#'     median <- median(ss, na.rm=T)
#'   }
#'   ts.med <- c(ts.med,median)
#' }
#' phen <- round(ts.med)
#' PhenKplot(phen,d1,h=1,nGS=365, xlab="DOY",ylab="NDMI",rge=c(0,10000))
#'
#' ############ Section 2:Calculate anomaly and their likelihoods values ##################
#' sample.rts <- rts(MDD,lan.dates)
#' pix.num <- cellFromXY(sample.rts,c(X,Ycoordinates))
#' ts.inter <- extract(sample.rts,pix.num)
#' vec <- as.vector(ts.inter)
#' dd <- lan.dates # NON SORTED DATES
#' nn <- seq(1,length(vec),1)
#' df <-data.frame(cbind(nn,dd))
#' sort <-df[order(as.Date(df$dd)),]
#' sort$vec <- vec
#' sort$date <- as.Date(sort$dd)
#' corr <- seq(1,length(vec),1)
#' sort$corr <- corr
#' x <- sort$vec
#' dates <- as.Date(sort$date)
#' ano_rfd <- PLUGPhenAnoRFDPLUS(x=x,phen=phen,dates=dates,h=1,anop=c(1:nlayers),rge=c(1,10000))
#' }
#' @export
PLUGPhenAnoRFDPLUS <-
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

    DOY <- yday(dates)
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

    Hmat <- Hpi(na.omit(D1))
    Hmat[1, 2] <- Hmat[2, 1]
    K1 <- kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), xmax = c(365, rge[2]), gridsize = c(365, 500))
    K1Con <- K1$estimate
    for (j in 1:365) {
      Kdiv <- sum(K1$estimate[j, ])
      ifelse(Kdiv == 0, K1Con[j, ] <- 0, K1Con[j, ] <- K1$estimate[j, ] / sum(K1$estimate[j, ]))
    }
    MAXY <- apply(K1Con, 1, max)
    for (i in 1:365) {
      MAXY[i] <- median(K1$eval.points[[2]][which(K1Con[i, ] == MAXY[i], arr.ind = TRUE)])
    }

    # c. Calculating cumulative bivariate density distribution

    h2d <- list()
    h2d$x <- seq(1, 365)
    h2d$y <- seq(rge[1], rge[2], len = 500)
    h2d$density <- K1Con / sum(K1Con)
    uniqueVals <- rev(unique(sort(h2d$density)))
    cumProbs <- cumsum(uniqueVals)
    names(cumProbs) <- uniqueVals
    h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
    h2d$cumDensity[] <- cumProbs[as.character(h2d$density)]

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
      } else {
        which.min(x)
      }
    }))
    AnomProb <- rep(NA, nrow(D2))
    for (i in 1:nrow(D2)) {
      AnomProb[i] <- h2d$cumDensity[D2[i, 1], rowAnom2[i]]
    }
    AnomProbPerc <- round(100 * (AnomProb))
    plot(AnomProbPerc[1:ano.len], xlab = "Time", ylab = "Likelihood (%)", font.lab = 2)
    abline(h = 90, col = "red")
    abline(h = 95, col = "red")
    abline(h = 99, col = "red")
    AnomProbPerc[1:ano.len]
    # Two results in a single numeric vector
    AnomPLUSProb <- c(Anoma[1:ano.len], AnomProbPerc[1:ano.len])
    AnomPLUSProb[1:len2]
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
#' @import raster
#' @import RColorBrewer
#' @import rgdal
#' @import bfastSpatial
#' @import lubridate
#' @import rts
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

      DOY <- yday(dates)
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

      Hmat <- Hpi(na.omit(D1))
      Hmat[1, 2] <- Hmat[2, 1]
      K1 <- kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), xmax = c(365, rge[2]), gridsize = c(365, 500))
      K1Con <- K1$estimate
      for (j in 1:365) {
        Kdiv <- sum(K1$estimate[j, ])
        ifelse(Kdiv == 0, K1Con[j, ] <- 0, K1Con[j, ] <- K1$estimate[j, ] / sum(K1$estimate[j, ]))
      }

      MAXY <- apply(K1Con, 1, max)
      for (i in 1:365) {
        MAXY[i] <- median(K1$eval.points[[2]][which(K1Con[i, ] == MAXY[i], arr.ind = TRUE)])
      }

      # c. Calculating cumulative bivariate density distribution

      h2d <- list()
      h2d$x <- seq(1, 365)
      h2d$y <- seq(rge[1], rge[2], len = 500)
      h2d$density <- K1Con / sum(K1Con)
      uniqueVals <- rev(unique(sort(h2d$density)))
      cumProbs <- cumsum(uniqueVals)
      names(cumProbs) <- uniqueVals
      h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
      h2d$cumDensity[] <- cumProbs[as.character(h2d$density)]

      # c. Anomaly calculation (D2)

      if (h == 2) {
        for (i in 1:nrow(D2)) {
          D2[i, 1] <- DOGS[which(DOGS[, 1] == D2[i, 1], arr.ind = TRUE), 2]
        }
      }

      # d. Calculating anomalies AND their likelihoods based on D2

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

      # d.2 Likelihood
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
      AnomProb <- rep(NA, nrow(D2))
      for (i in 1:nrow(D2)) {
        AnomProb[i] <- h2d$cumDensity[D2[i, 1], rowAnom2[i]]
      }
      AnomProbPerc <- round(100 * (AnomProb))
      AnomProbPerc[1:ano.len]
      # Two results in a single numeric vector
      AnomPLUSProb <- c(Anoma[1:ano.len], AnomProbPerc[1:ano.len])
      AnomPLUSProb[1:len2]
    }

    #----------------------------------------------------------------------------------------
    # cluster processing
    app(s, fun = ff, filename = outname, cores = nCluster, overwrite = T, wopt = list(datatype = datatype))
  }
