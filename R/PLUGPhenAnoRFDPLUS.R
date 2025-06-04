#' @title Avocado algorithm - calculate likelihood of vegetation-index–time space for reference vegetation 
#' @name PhenRef2d
#' @description Calculate the cumulative density distribution and their likelihood values based on reference vegetation - input for Wall-to-wall map version (i.e. PLUGPhenAnoRFDPLUS function)
#' @author  Roberto O. Chavez, Mathieu Decuyper
#' @param phen Numeric vector with the values of the reference vegetation
#' @param dates A date vector. The number of dates must be equal to the number of values of time series.
#' @param h Numeric. Geographic hemisphere to define the starting date of the growing season. h=1 Northern Hemisphere; h=2 Southern Hemisphere.
#' @param anop Numeric vector with the number of values that are in x. For those values the anomalies and likelihoods will be calculated based on the phen. For example a time series has 450 values, so anop=c(1:450).
#' @param rge A vector containing minimum and maximum values of the response variable used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge =c(0,10000)
#' @return List which contains cumulative bivariate density distribution (cumDensity), and maximum likelihood of the vegetation-index–time space (MAXY) 
#' @import npphen
#' @import stats
#' @importFrom lubridate yday
#' @seealso \code{\link{PLUGPhenAnoRFDMapPLUS}} \code{\link{PLUGPhenAnoRFDPLUS}}
#' @examples
#' \dontrun{
#' # Loading raster data
#' library(terra)
#' library(npphen)
#' 
#' MDD <- rast(system.file("extdata", "MDD_NDMI_1990_2020.tif", package = "AVOCADO"))
#' # load dates vector
#' load(system.file("extdata", "MDD_dates.RData", package = "AVOCADO"))
#' # load  reference forest shapefile
#' MDDref <- vect(system.file("extdata", "MDDref.gpkg", package = "AVOCADO"))
#' # Create the reference curve
#' ref.ext <- ext(MDDref)
#' ref.brick <- crop(MDD, ref.ext)
#' fin <- nrow(ref.brick) * ncol(ref.brick)
#' phen <- as.numeric(terra::extract(ref.brick, 1))
#' d1 <- MDD_dates
#' for (i in 2:fin) {
#'   pp <- as.numeric(terra::extract(ref.brick, i))
#'   phen <- c(phen, pp)
#'   d1 <- c(d1, MDD_dates)
#' }
#' MDD_fref <- PhenRef2d(phen, d1, h = 1, anop = c(1:1063), rge = c(0, 10000))
#' 
#' # plot reference curve + probabilities (same result as using PhenKplot())
#' image(MDD_fref$x, MDD_fref$y, MDD_fref$cumDensity, xlab = "DOY", ylab = "NMDI", 
#' font.lab = 2, breaks = c(0, 0.5, 0.75, 0.9, 0.95), col = grDevices::heat.colors(n = 4, alpha = 0.6))
#' contour(MDD_fref$x, MDD_fref$y, MDD_fref$cumDensity, levels = c(0, 0.5, 0.75, 0.9, 0.95), add = T, col = grDevices::grey(0.25), labcex = 1)
#' lines(seq(1, 365), MDD_fref$MAXY, lwd = 3, col = "dark red")
#'
#' }
#' 
#' @export
PhenRef2d <-
  function(phen, dates, h, anop, rge) {
    # a.Preparing dataset
    
    if (length(rge) != 2) {
      stop("rge must be a vector of length 2")
    }
    if (rge[1] > rge[2]) {
      stop("rge vector order must be minimum/maximum")
    }
    if (length(dates) != length(phen)) {
      stop("N of dates and files do not match")
    }
    
    ano.min <- min(anop)
    ano.max <- max(anop)
    ano.len <- ano.max - ano.min + 1
    len2 <- 2 * ano.len
    
    if (ano.min >= ano.max) {
      stop("for anop, lower value > upper value")
    }
    
    DOY <- lubridate::yday(dates)
    DOY[which(DOY == 366)] <- 365
    D1 <- cbind(DOY, phen)
    
    if (length(unique(D1[, 2])) < 10 | (nrow(D1) - sum(is.na(D1))) < (0.1 * length(D1))) {
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
    h2d$density <- K1Con / sum(K1Con)
    uniqueVals <- rev(unique(sort(h2d$density)))
    cumRFDs <- cumsum(uniqueVals)
    names(cumRFDs) <- uniqueVals
    h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
    h2d$cumDensity[] <- cumRFDs[as.character(h2d$density)]
    na.sta <- first.no.NA.DOY - 1
    na.end <- last.no.NA.DOY + 1
    if (na.sta >= 1) {
      h2d$cumDensity[1:na.sta, ] <- NA
    }
    if (na.end <= 365) {
      h2d$cumDensity[na.end:365, ] <- NA
    }
    
    # d. prepare output: attach MAXY to h2d list
    
    h2d$MAXY <- MAXY
    h2d
    
  }

#' @title Avocado algorithm - Single pixel version
#' @name PLUGPhenAnoRFDPLUS
#' @description Calculate the Anomaly and their likelihood values (based on the difference between the pixel values and reference vegetation) for a numeric vector
#' @author  Roberto O. Chavez, Mathieu Decuyper
#' @param x numeric vector. Time series of vegetation index (e.g. NDVI, EVI, NDMI)
#' @param phenref List which contains cumulative bivariate density distribution and maximum likelihood of the vegetation-index–time space based on reference vegetation obtained from \code{\link{PhenRef2d}}
#' @param dates A date vector. The number of dates must be equal to the number of values of time series.
#' @param h Numeric. Geographic hemisphere to define the starting date of the growing season. h=1 Northern Hemisphere; h=2 Southern Hemisphere.
#' @param anop Numeric vector with the number of values that are in x. For those values the anomalies and likelihoods will be calculated based on the phen. For example a time series has 450 values, so anop=c(1:450).
#' @return vector of length([x]*2) containing all anomalies, followed by their position within the reference frequency distribution (RFD)
#' @seealso \code{\link{PLUGPhenAnoRFDMapPLUS}} \code{\link{PhenRef2d}}
#' @import npphen
#' @import stats
#' @importFrom lubridate yday
#' @importFrom graphics abline
#' @examples
#' \dontrun{
#' # ===================================================================================================================
#' # Loading raster data
#' library(terra)
#' library(npphen)
#' MDD <- rast(system.file("extdata", "MDD_NDMI_1990_2020.tif", package = "AVOCADO"))
#' # load dates vector
#' load(system.file("extdata", "MDD_dates.RData", package = "AVOCADO"))
#' # load  reference forest data (output from PhenRef2d)
#' load(system.file("extdata", "MDD_forestReference.RData", package = "AVOCADO"))
#' ## time series extraction for a single pixel
#' px <- vect(cbind(-69.265, -12.48))
#' plot(MDD[[1]])
#' plot(px, add = T)
#'
#' # extract series
#' px_series <- as.numeric(terra::extract(MDD, px, ID = F))
#' plot(MDD_dates, px_series, type = "b", xlab = "", ylab = "NDMI")
#' ## Anomaly calculation
#' anom_rfd <- PLUGPhenAnoRFDPLUS(x = px_series, phenref = MDD_fref, dates = MDD_dates, h = 2, anop = c(1:1063))
#' }
#'
#' @export
PLUGPhenAnoRFDPLUS <-
  function(x, phenref, dates, h, anop) {
    # a.Preparing dataset
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
    D2 <- cbind(DOY[ano.min:ano.max], x[ano.min:ano.max])


    if (all(is.na(D2[, 2]))) {
      return(rep(NA, len2))
    }

    # b. Calculating anomalies AND their RFD based on D2
    if (h != 1 && h != 2) {
      stop("Invalid h")
    }
    DOGS <- cbind(seq(1, 365), c(seq(185, 365), seq(1, 184)))
    if (h == 2) {
      for (i in 1:nrow(D2)) {
        D2[i, 1] <- DOGS[which(DOGS[, 1] == D2[i, 1], arr.ind = TRUE), 2]
      }
    }

    
    # b.1 Anomalies 
    MAXY <- phenref$MAXY
    
    Anoma <- rep(NA, ano.len)
    for (i in 1:nrow(D2)) {
      Anoma[i] <- as.integer(D2[i, 2] - MAXY[D2[i, 1]])
    }
    Anoma[1:ano.len]
    
    # b.2 RFD
    rowAnom <- matrix(NA, nrow = nrow(D2), ncol = 500)
    for (i in 1:nrow(D2)) {
      rowAnom[i, ] <- abs(phenref$y - D2[i, 2])
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
      AnomRFD[i] <- phenref$cumDensity[D2[i, 1], rowAnom2[i]]
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
#' @param s SpatRaster time series of vegetation index (e.g. NDVI, EVI, NDMI)
#' @param anop Numeric vector with the number of values that are in x. For those values the anomalies and likelihoods will be calculated based on the phen. For example a time series has 450 values, so anop=c(1:450).
#' @param phenref List which contains cumulative bivariate density distribution and maximum likelihood of the vegetation-index–time space based on reference vegetation obtained from \code{\link{PhenRef2d}}
#' @param dates A date vector. The number of dates must be equal to the number of values of time series.
#' @param h Numeric. Geographic hemisphere to define the starting date of the growing season. h=1 Northern Hemisphere; h=2 Southern Hemisphere.
#' @param nCluster Numeric. Number of CPU's to be used for the job.
#' @param outname Character vector with the output path and filename with extension or only the filename and extension if work directory was set. More information: See writeRaster
#' @param datatype Character vector that determines the interpretation of values written to disk. More information: See \code{\link[terra]{writeRaster}}
#' @return SpatRaster of nlayers([s]stack*2) containing all anomalies, followed by their likelihoods
#' @import npphen
#' @import terra
#' @import stats
#' @importFrom lubridate yday
#' @seealso \code{\link{PLUGPhenAnoRFDPLUS}}
#' @examples
#' \dontrun{
#' # Loading raster data
#' library(terra)
#' library(npphen)
#' 
#' MDD <- rast(system.file("extdata", "MDD_NDMI_1990_2020.tif", package = "AVOCADO"))
#' # load dates vector
#' load(system.file("extdata", "MDD_dates.RData", package = "AVOCADO"))
#' # load  reference forest data (output from PhenRef2d)
#' load(system.file("extdata", "MDD_forestReference.RData", package = "AVOCADO"))
#'
#' ## Anomaly calculation
#' # checking availiable cores and leave one free
#' nc1 <- parallel::detectCores() - 1
#' PLUGPhenAnoRFDMapPLUS(
#'   s = MDD, dates = MDD_dates, h = 1, phenref = MDD_fref, anop = c(1:1063),
#'   nCluster = nc1, outname = "YourDirectory/MDD_AnomalyLikelihood.tif",
#'   datatype = "INT2S")
#' )
#' # The output file contains all the anomalies, followed by all their likelihoods
#' # and thus has twice the number of layers as the input raster stack.
#' }
#' @export
PLUGPhenAnoRFDMapPLUS <-
  function(s, phenref, dates, h, anop, nCluster, outname, datatype) {
    ff <- function(x) {
      # a.Preparing dataset

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
      
      D2 <- cbind(DOY[ano.min:ano.max], x[ano.min:ano.max])

      if (all(is.na(D2[, 2]))) {
        return(rep(NA, len2))
      }

      if (h != 1 && h != 2) {
        stop("Invalid h")
      }
      DOGS <- cbind(seq(1, 365), c(seq(185, 365), seq(1, 184)))

      # b. Calculating anomalies AND their RFD based on D2
      if (h == 2) {
        for (i in 1:nrow(D2)) {
          D2[i, 1] <- DOGS[which(DOGS[, 1] == D2[i, 1], arr.ind = TRUE), 2]
        }
      }

      # b.1 Anomalies
      MAXY <- phenref$MAXY
      
      Anoma <- rep(NA, ano.len)
      for (i in 1:nrow(D2)) {
        Anoma[i] <- as.integer(D2[i, 2] - MAXY[D2[i, 1]])
      }
      Anoma[1:ano.len]

      # b.2 RFD
      rowAnom <- matrix(NA, nrow = nrow(D2), ncol = 500)
      for (i in 1:nrow(D2)) {
        rowAnom[i, ] <- abs(phenref$y - D2[i, 2])
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
        AnomRFD[i] <- phenref$cumDensity[D2[i, 1], rowAnom2[i]]
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
