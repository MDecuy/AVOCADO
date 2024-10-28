#' @title Avocado algorithm - Single pixel version
#' @name dist.reg
#' @description Continuous vegetation change detection
#' @author  Roberto O. Chavez, Mathieu Decuyper
#' @param x Vector of anomaly and likelihood time-series created in using \code{\link{PLUGPhenAnoRFDPLUS}}
#' @param dates The julian dates for each scene in the time series. Should be in the following format: "YYYY-MM-DD"
#' @param rfd The reference frequency distribution. Determines if an anomaly falls outside the 95 percent of the reference frequency distribution. For example a value that fall in a RDF >= 0.95, indicates that the detected anomaly belongs to the 5 percent of lowest values and is a potential disturbance/change.
#' @param dstrb_thr sets a threshold (number of days) in which if a regrowth is detected within n days after a potential disturbance, then the candidate disturbance date is neglected.
#' @param rgrow_thr sets a threshold (number of days) in which if a disturbance is detected within n days after a potential regrowth, then the candidate regrowth date is neglected.
#' @param cdates Sets the number of consecutive data points for a change to be detected (minimum 2 and maximum 5 data points)
#' @return The output will consist of a vector with 16 values: 1st value = year of disturbance detection; 2nd value = corresponding anomalies; 3rd=year of regrowth detection; 4th value = corresponding anomalies band 3; 5th value = 2nd disturbance detection for data points that had a previous disturbance and regrowth; 6th value etc.
#' @seealso \code{\link{dist.reg.map}}
#'
#' @import npphen
#' @import RColorBrewer
#'
#' @examples
#' \dontrun{
#' dist.reg(x = ano_rfd, dates = dates, rfd = 0.95, dstrb_thr = 1, rgrow_thr = 730, cdates = 3)
#' }
#'
#' @export
dist.reg <-
  function(x, dates, rfd, dstrb_thr, rgrow_thr, cdates) { # dstrb_thr = 730, #rgrow_thr =  1. i.e. disturbance will be neglected if there is a regrowth in the upcoming 2 years
    if (length(dates) != length(x) / 2) {
      stop("N of dates and files do not match")
    }

    # Checking for a valid cdates value
    if (cdates == 2) {
      cd2 <- 1
    } else {
      cd2 <- 0
    }
    if (cdates == 3) {
      cd3 <- 1
    } else {
      cd3 <- 0
    }
    if (cdates == 4) {
      cd4 <- 1
    } else {
      cd4 <- 0
    }
    if (cdates == 5) {
      cd5 <- 1
    } else {
      cd5 <- 0
    }
    cd.sum <- cd2 + cd3 + cd4 + cd5
    if (cd.sum != 1) {
      stop("cdates argument must be either 2 or 3 or 4 or 5")
    }

    if (all(is.na(x))) {
      return(as.numeric(rep(NA, 16)))
    }

    vv <- as.vector(x)
    ano.ini <- 1
    ano.fin <- length(vv) / 2
    prob.ini <- ano.fin + 1
    prob.fin <- as.numeric(length(vv))
    ano <- vv[ano.ini:ano.fin]
    prob <- vv[prob.ini:prob.fin]
    rfd <- rfd * 100

    # Preparing data -> sorting if the are not sorted already
    df <- data.frame(cbind(ano, prob, dates))
    df$YY <- substr(dates, 1, 4)
    sort <- df[order(df$dates), ]
    sort <- na.omit(sort)
    corr <- seq(1, length(sort$ano), 1)
    sort$corr <- corr
    vvv.ano <- sort$ano
    vvv.prob <- sort$prob
    vvv.ano[vvv.prob < rfd & vvv.ano < 0] <- NA # delete negative anomalies with < rfd
    text <- paste("grey are negative anomalies with rfd < ", rfd / 100, sep = "")
    plot(sort$ano,
      ylim = c(-5000, 5000), ylab = "Anomaly", xlab = "Scene number", font.sub = 3,
      main = "Disturbance (red) and Regrowth (green) detection", cex.main = 0.9,
      sub = text, col = "grey"
    )
    points(seq(1, length(vvv.ano), 1), vvv.ano)
    abline(h = 0, col = "red")
    fin <- length(sort$ano)

    #----------------------------------------
    # Checking for the last period, according to rgrow_thr and dstrb_thr
    last.ddis <- max(df$dates) - dstrb_thr # last date a disturbance can be found
    last.dreg <- max(df$dates) - rgrow_thr # last date a regrowth can be found

    yy.dist1 <- NA
    val.dist1 <- NA
    corr.dist1 <- NA
    dd.dist1 <- NA
    yy.reg1 <- NA
    val.reg1 <- NA
    corr.reg1 <- NA
    dd.reg1 <- NA

    yy.dist2 <- NA
    val.dist2 <- NA
    corr.dist2 <- NA
    dd.dist2 <- NA
    yy.reg2 <- NA
    val.reg2 <- NA
    corr.reg2 <- NA
    dd.reg2 <- NA

    yy.dist3 <- NA
    val.dist3 <- NA
    corr.dist3 <- NA
    dd.dist3 <- NA
    yy.reg3 <- NA
    val.reg3 <- NA
    corr.reg3 <- NA
    dd.reg3 <- NA

    yy.dist4 <- NA
    val.dist4 <- NA
    corr.dist4 <- NA
    dd.dist4 <- NA
    yy.reg4 <- NA
    val.reg4 <- NA
    corr.reg4 <- NA
    dd.reg4 <- NA

    #----------------------------------------------------------
    # First cycle

    # To catch disturbance 1
    val.dist <- sort$ano[which(sort$prob >= rfd & sort$ano < 0)]
    corr.dist <- sort$corr[which(sort$prob >= rfd & sort$ano < 0)]
    yy.dist <- sort$YY[which(sort$prob >= rfd & sort$ano < 0)]
    dd.dist <- sort$dates[which(sort$prob >= rfd & sort$ano < 0)]

    if (length(val.dist) < cdates) {
      return(as.numeric(rep(NA, 16)))
    }
    regrowth <- "off"
    i <- 1

    while (regrowth == "off") {
      if (i == (length(val.dist) - cdates + 1)) {
        regrowth <- "none"
      } else {
        # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
        if (cdates == 2) {
          cond <- "((corr.dist[i]+1)==corr.dist[i+1])"
        }
        if (cdates == 3) {
          cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"
        }
        if (cdates == 4) {
          cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"
        }
        if (cdates == 5) {
          cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"
        }

        if (eval(parse(text = cond))) { # test the selected condition according to cdates

          ssort <- sort[corr.dist[i]:fin, ]
          corr.reg <- ssort$corr[which(ssort$ano >= 0)]
          dd.reg <- ssort$dates[which(ssort$ano >= 0)]

          rdura <- NULL
          for (t in 1:length(dd.reg)) {
            rresta <- dd.reg[t] - dd.dist[i]
            rdura <- c(rdura, rresta)
            w1y.cor.reg <- corr.reg[which(rdura < dstrb_thr)]
          }

          rcorr <- NULL

          if (length(w1y.cor.reg) <= 1) {
            rcorr <- 0
          } else {
            for (k in 1:(length(w1y.cor.reg) - 1)) {
              if ((w1y.cor.reg[k] + 1) == w1y.cor.reg[k + 1]) {
                rifcor <- 1
              } else {
                rifcor <- 0
              }
              rcorr <- c(rcorr, rifcor)
            }
          }

          if (sum(rcorr) == 0) { ############# Here makes the condition of not having consecutive regrowth within a year
            val.dist1 <- val.dist[i]
            corr.dist1 <- corr.dist[i]
            yy.dist1 <- yy.dist[i]
            dd.dist1 <- dd.dist[i]
            regrowth <- "on"
          }
        }
      }
      i <- i + 1
    }
    abline(v = corr.dist1, col = "red", lty = "dashed")

    # To catch regrowth 1
    ini1 <- corr.dist[i]
    fin <- length(sort$ano)
    sort2 <- sort[ini1:fin, ]
    j <- 1

    val.reg <- sort2$ano[which(sort2$ano >= 0)]
    corr.reg <- sort2$corr[which(sort2$ano >= 0)]
    yy.reg <- sort2$YY[which(sort2$ano >= 0)]
    dd.reg <- sort2$dates[which(sort2$ano >= 0)]


    if (length(val.reg) < cdates) {
      regrowth <- "off"
    }
    while (regrowth == "on") {
      if (j == (length(val.reg) - cdates + 1)) {
        regrowth <- "off"
      } else {
        # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
        if (cdates == 2) {
          cond <- "((corr.reg[j]+1)==corr.reg[j+1])"
        }
        if (cdates == 3) {
          cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"
        }
        if (cdates == 4) {
          cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"
        }
        if (cdates == 5) {
          cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"
        }

        if (eval(parse(text = cond))) { # test the selected condition according to cdates

          ssort2 <- sort[corr.reg[j]:fin, ]
          corr.dist <- ssort2$corr[which(ssort2$prob >= rfd & ssort2$ano < 0)]
          dd.dist <- ssort2$dates[which(ssort2$prob >= rfd & ssort2$ano < 0)]

          ddura <- NULL
          for (t in 1:length(dd.dist)) {
            resta <- dd.dist[t] - dd.reg[j]
            ddura <- c(ddura, resta)
            w1y.cor.dist <- corr.dist[which(ddura < rgrow_thr)]
          }

          dcorr <- NULL

          if (length(w1y.cor.dist) <= 1) {
            dcorr <- 0
          } else {
            for (k in 1:(length(w1y.cor.dist) - 1)) {
              if ((w1y.cor.dist[k] + 1) == w1y.cor.dist[k + 1]) {
                ifcor <- 1
              } else {
                ifcor <- 0
              }
              dcorr <- c(dcorr, ifcor)
            }
          }

          if (sum(dcorr) == 0) { ############# Here makes the condition of not having consecutive disturbances within a year
            val.reg1 <- val.reg[j]
            corr.reg1 <- corr.reg[j]
            yy.reg1 <- yy.reg[j]
            dd.reg1 <- dd.reg[j]
            regrowth <- "off"
          }
        }
      }
      j <- j + 1
    }
    abline(v = corr.reg1, col = "green", lty = "dashed")

    #----------------------------------------------------------
    # Second cycle
    if (is.na(val.reg1) != T) { # Here we continue if there was a regrowth
      ini2 <- corr.reg[j]
      sort3 <- sort[ini2:fin, ]

      # To catch disturbance 2
      val.dist <- sort3$ano[which(sort3$prob >= rfd & sort3$ano < 0)]
      corr.dist <- sort3$corr[which(sort3$prob >= rfd & sort3$ano < 0)]
      yy.dist <- sort3$YY[which(sort3$prob >= rfd & sort3$ano < 0)]
      dd.dist <- sort3$dates[which(sort3$prob >= rfd & sort3$ano < 0)]

      if (length(val.dist) >= cdates) {
        regrowth <- "off"
        i <- 1
        while (regrowth == "off") {
          if (i == (length(val.dist) - cdates + 1)) {
            regrowth <- "none"
          } else {
            # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
            if (cdates == 2) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])"
            }
            if (cdates == 3) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"
            }
            if (cdates == 4) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"
            }
            if (cdates == 5) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"
            }

            if (eval(parse(text = cond))) { # test the selected condition according to cdates

              ssort3 <- sort[corr.dist[i]:fin, ]
              corr.reg <- ssort3$corr[which(ssort3$ano >= 0)]
              dd.reg <- ssort3$dates[which(ssort3$ano >= 0)]

              rdura <- NULL
              for (t in 1:length(dd.reg)) {
                rresta <- dd.reg[t] - dd.dist[i]
                rdura <- c(rdura, rresta)
                w1y.cor.reg <- corr.reg[which(rdura < dstrb_thr)]
              }

              rcorr <- NULL

              if (length(w1y.cor.reg) <= 1) {
                rcorr <- 0
              } else {
                for (k in 1:(length(w1y.cor.reg) - 1)) {
                  if ((w1y.cor.reg[k] + 1) == w1y.cor.reg[k + 1]) {
                    rifcor <- 1
                  } else {
                    rifcor <- 0
                  }
                  rcorr <- c(rcorr, rifcor)
                }
              }

              if (sum(rcorr) == 0) { ############# Here makes the condition of not having consecutive regrowth within a year
                val.dist2 <- val.dist[i]
                corr.dist2 <- corr.dist[i]
                yy.dist2 <- yy.dist[i]
                dd.dist2 <- dd.dist[i]
                regrowth <- "on"
              }
            }
          }
          i <- i + 1
        }
        abline(v = corr.dist2, col = "red", lty = "dashed")

        # To catch regrowth 2
        ini3 <- corr.dist[i]
        sort4 <- sort[ini3:fin, ]
        j <- 1
        val.reg <- sort4$ano[which(sort4$ano >= 0)]
        corr.reg <- sort4$corr[which(sort4$ano >= 0)]
        yy.reg <- sort4$YY[which(sort4$ano >= 0)]
        dd.reg <- sort4$dates[which(sort4$ano >= 0)]

        if (length(val.reg) < cdates) {
          regrowth <- "off"
        }
        while (regrowth == "on") {
          if (j == (length(val.reg) - cdates + 1)) {
            regrowth <- "off"
          } else {
            # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
            if (cdates == 2) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])"
            }
            if (cdates == 3) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"
            }
            if (cdates == 4) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"
            }
            if (cdates == 5) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"
            }

            if (eval(parse(text = cond))) { # test the selected condition according to cdates

              ssort4 <- sort[corr.reg[j]:fin, ]
              corr.dist <- ssort4$corr[which(ssort4$prob >= rfd & ssort4$ano < 0)]
              dd.dist <- ssort4$dates[which(ssort4$prob >= rfd & ssort4$ano < 0)]

              ddura <- NULL
              for (t in 1:length(dd.dist)) {
                resta <- dd.dist[t] - dd.reg[j]
                ddura <- c(ddura, resta)
                w1y.cor.dist <- corr.dist[which(ddura < rgrow_thr)]
              }

              dcorr <- NULL

              if (length(w1y.cor.dist) <= 1) {
                dcorr <- 0
              } else {
                for (k in 1:(length(w1y.cor.dist) - 1)) {
                  if ((w1y.cor.dist[k] + 1) == w1y.cor.dist[k + 1]) {
                    ifcor <- 1
                  } else {
                    ifcor <- 0
                  }
                  dcorr <- c(dcorr, ifcor)
                }
              }

              if (sum(dcorr) == 0) { ############# Here makes the condition of not having consecutive disturbances within a year
                val.reg2 <- val.reg[j]
                corr.reg2 <- corr.reg[j]
                yy.reg2 <- yy.reg[j]
                dd.reg2 <- dd.reg[j]
                regrowth <- "off"
              }
            }
          }
          j <- j + 1
        }
        abline(v = corr.reg2, col = "green", lty = "dashed")
      }
    }
    #----------------------------------------------------------
    # Third cycle
    if (is.na(val.reg2) != T) { # Here we continue if there was a regrowth
      ini4 <- corr.reg[j]
      sort5 <- sort[ini4:fin, ]

      # To catch disturbance 3
      val.dist <- sort5$ano[which(sort5$prob >= rfd & sort5$ano < 0)]
      corr.dist <- sort5$corr[which(sort5$prob >= rfd & sort5$ano < 0)]
      yy.dist <- sort5$YY[which(sort5$prob >= rfd & sort5$ano < 0)]
      dd.dist <- sort5$dates[which(sort5$prob >= rfd & sort5$ano < 0)]

      if (length(val.dist) >= cdates) {
        regrowth <- "off"
        i <- 1
        while (regrowth == "off") {
          if (i == (length(val.dist) - cdates + 1)) {
            regrowth <- "none"
          } else {
            # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
            if (cdates == 2) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])"
            }
            if (cdates == 3) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"
            }
            if (cdates == 4) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"
            }
            if (cdates == 5) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"
            }

            if (eval(parse(text = cond))) { # test the selected condition according to cdates

              ssort5 <- sort[corr.dist[i]:fin, ]
              corr.reg <- ssort5$corr[which(ssort5$ano >= 0)]
              dd.reg <- ssort5$dates[which(ssort5$ano >= 0)]

              rdura <- NULL
              for (t in 1:length(dd.reg)) {
                rresta <- dd.reg[t] - dd.dist[i]
                rdura <- c(rdura, rresta)
                w1y.cor.reg <- corr.reg[which(rdura < dstrb_thr)]
              }

              rcorr <- NULL

              if (length(w1y.cor.reg) <= 1) {
                rcorr <- 0
              } else {
                for (k in 1:(length(w1y.cor.reg) - 1)) {
                  if ((w1y.cor.reg[k] + 1) == w1y.cor.reg[k + 1]) {
                    rifcor <- 1
                  } else {
                    rifcor <- 0
                  }
                  rcorr <- c(rcorr, rifcor)
                }
              }

              if (sum(rcorr) == 0) { ############# Here makes the condition of not having consecutive regrowth within a year
                val.dist3 <- val.dist[i]
                corr.dist3 <- corr.dist[i]
                yy.dist3 <- yy.dist[i]
                dd.dist3 <- dd.dist[i]
                regrowth <- "on"
              }
            }
          }
          i <- i + 1
        }
        abline(v = corr.dist3, col = "red", lty = "dashed")

        # To catch regrowth 3
        ini5 <- corr.dist[i]
        sort6 <- sort[ini5:fin, ]
        j <- 1
        val.reg <- sort6$ano[which(sort6$ano >= 0)]
        corr.reg <- sort6$corr[which(sort6$ano >= 0)]
        yy.reg <- sort6$YY[which(sort6$ano >= 0)]
        dd.reg <- sort6$dates[which(sort6$ano >= 0)]

        if (length(val.reg) < cdates) {
          regrowth <- "off"
        }
        while (regrowth == "on") {
          if (j == (length(val.reg) - cdates + 1)) {
            regrowth <- "off"
          } else {
            # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
            if (cdates == 2) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])"
            }
            if (cdates == 3) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"
            }
            if (cdates == 4) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"
            }
            if (cdates == 5) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"
            }

            if (eval(parse(text = cond))) { # test the selected condition according to cdates

              ssort6 <- sort[corr.reg[j]:fin, ]
              corr.dist <- ssort6$corr[which(ssort6$prob >= rfd & ssort6$ano < 0)]
              dd.dist <- ssort6$dates[which(ssort6$prob >= rfd & ssort6$ano < 0)]

              ddura <- NULL
              for (t in 1:length(dd.dist)) {
                resta <- dd.dist[t] - dd.reg[j]
                ddura <- c(ddura, resta)
                w1y.cor.dist <- corr.dist[which(ddura < rgrow_thr)]
              }

              dcorr <- NULL

              if (length(w1y.cor.dist) <= 1) {
                dcorr <- 0
              } else {
                for (k in 1:(length(w1y.cor.dist) - 1)) {
                  if ((w1y.cor.dist[k] + 1) == w1y.cor.dist[k + 1]) {
                    ifcor <- 1
                  } else {
                    ifcor <- 0
                  }
                  dcorr <- c(dcorr, ifcor)
                }
              }

              if (sum(dcorr) == 0) { ############# Here makes the condition of not having consecutive disturbances within a year
                val.reg3 <- val.reg[j]
                corr.reg3 <- corr.reg[j]
                yy.reg3 <- yy.reg[j]
                dd.reg3 <- dd.reg[j]
                regrowth <- "off"
              }
            }
          }
          j <- j + 1
        }
        abline(v = corr.reg3, col = "green", lty = "dashed")
      }
    }
    #----------------------------------------------------------
    # Fourth cycle
    if (is.na(val.reg3) != T) { # Here we continue if there was a regrowth
      ini6 <- corr.reg[j]
      sort7 <- sort[ini6:fin, ]
      # To catch disturbance 4
      val.dist <- sort7$ano[which(sort7$prob >= rfd & sort7$ano < 0)]
      corr.dist <- sort7$corr[which(sort7$prob >= rfd & sort7$ano < 0)]
      yy.dist <- sort7$YY[which(sort7$prob >= rfd & sort7$ano < 0)]
      dd.dist <- sort7$dates[which(sort7$prob >= rfd & sort7$ano < 0)]

      if (length(val.dist) >= cdates) {
        regrowth <- "off"
        i <- 1
        while (regrowth == "off") {
          if (i == (length(val.dist) - cdates + 1)) {
            regrowth <- "none"
          } else {
            # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
            if (cdates == 2) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])"
            }
            if (cdates == 3) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"
            }
            if (cdates == 4) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"
            }
            if (cdates == 5) {
              cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"
            }

            if (eval(parse(text = cond))) { # test the selected condition according to cdates

              ssort7 <- sort[corr.dist[i]:fin, ]
              corr.reg <- ssort7$corr[which(ssort7$ano >= 0)]
              dd.reg <- ssort7$dates[which(ssort7$ano >= 0)]

              rdura <- NULL
              for (t in 1:length(dd.reg)) {
                rresta <- dd.reg[t] - dd.dist[i]
                rdura <- c(rdura, rresta)
                w1y.cor.reg <- corr.reg[which(rdura < dstrb_thr)]
              }

              rcorr <- NULL

              if (length(w1y.cor.reg) <= 1) {
                rcorr <- 0
              } else {
                for (k in 1:(length(w1y.cor.reg) - 1)) {
                  if ((w1y.cor.reg[k] + 1) == w1y.cor.reg[k + 1]) {
                    rifcor <- 1
                  } else {
                    rifcor <- 0
                  }
                  rcorr <- c(rcorr, rifcor)
                }
              }

              if (sum(rcorr) == 0) { ############# Here makes the condition of not having consecutive regrowth within a year
                val.dist4 <- val.dist[i]
                corr.dist4 <- corr.dist[i]
                yy.dist4 <- yy.dist[i]
                dd.dist4 <- dd.dist[i]
                regrowth <- "on"
              }
            }
          }
          i <- i + 1
        }
        abline(v = corr.dist4, col = "red", lty = "dashed")

        # To catch regrowth 4
        ini7 <- corr.dist[i]
        sort8 <- sort[ini7:fin, ]
        j <- 1
        val.reg <- sort8$ano[which(sort8$ano >= 0)]
        corr.reg <- sort8$corr[which(sort8$ano >= 0)]
        yy.reg <- sort8$YY[which(sort8$ano >= 0)]
        dd.reg <- sort8$dates[which(sort8$ano >= 0)]

        if (length(val.reg) < cdates) {
          regrowth <- "off"
        }
        while (regrowth == "on") {
          if (j == (length(val.reg) - cdates + 1)) {
            regrowth <- "off"
          } else {
            # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
            if (cdates == 2) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])"
            }
            if (cdates == 3) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"
            }
            if (cdates == 4) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"
            }
            if (cdates == 5) {
              cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"
            }

            if (eval(parse(text = cond))) { # test the selected condition according to cdates

              ssort8 <- sort[corr.reg[j]:fin, ]
              corr.dist <- ssort8$corr[which(ssort8$prob >= rfd & ssort8$ano < 0)]
              dd.dist <- ssort8$dates[which(ssort8$prob >= rfd & ssort8$ano < 0)]

              ddura <- NULL
              for (t in 1:length(dd.dist)) {
                resta <- dd.dist[t] - dd.reg[j]
                ddura <- c(ddura, resta)
                w1y.cor.dist <- corr.dist[which(ddura < rgrow_thr)]
              }

              dcorr <- NULL

              if (length(w1y.cor.dist) <= 1) {
                dcorr <- 0
              } else {
                for (k in 1:(length(w1y.cor.dist) - 1)) {
                  if ((w1y.cor.dist[k] + 1) == w1y.cor.dist[k + 1]) {
                    ifcor <- 1
                  } else {
                    ifcor <- 0
                  }
                  dcorr <- c(dcorr, ifcor)
                }
              }

              if (sum(dcorr) == 0) { ############# Here makes the condition of not having consecutive disturbances within a year
                val.reg4 <- val.reg[j]
                corr.reg4 <- corr.reg[j]
                yy.reg4 <- yy.reg[j]
                dd.reg4 <- dd.reg[j]
                regrowth <- "off"
              }
            }
          }
          j <- j + 1
        }
        abline(v = corr.reg4, col = "green", lty = "dashed")
      }
    }
    # ---------------------------------------------------
    # Deleting outputs in the last rgrow_thr y dstrb_thr period
    if (dd.dist1 >= last.ddis | is.na(dd.dist1)) {
      yy.dist1 <- NA
      val.dist1 <- NA
      corr.dist1 <- NA
    }
    if (dd.reg1 >= last.dreg | is.na(dd.reg1)) {
      yy.reg1 <- NA
      val.reg1 <- NA
      corr.reg1 <- NA
    }

    if (dd.dist2 >= last.ddis | is.na(dd.dist2)) {
      yy.dist2 <- NA
      val.dist2 <- NA
      corr.dist2 <- NA
    }
    if (dd.reg2 >= last.dreg | is.na(dd.reg2)) {
      yy.reg2 <- NA
      val.reg2 <- NA
      corr.reg2 <- NA
    }

    if (dd.dist3 >= last.ddis | is.na(dd.dist3)) {
      yy.dist3 <- NA
      val.dist3 <- NA
      corr.dist3 <- NA
    }
    if (dd.reg3 >= last.dreg | is.na(dd.reg3)) {
      yy.reg3 <- NA
      val.reg3 <- NA
      corr.reg3 <- NA
    }

    if (dd.dist4 >= last.ddis | is.na(dd.dist4)) {
      yy.dist4 <- NA
      val.dist4 <- NA
      corr.dist4 <- NA
    }
    if (dd.reg4 >= last.dreg | is.na(dd.reg4)) {
      yy.reg4 <- NA
      val.reg4 <- NA
      corr.reg4 <- NA
    }

    # ---------------------------------------------------
    # New plot in case some dist or reg was deleted according to dstrb_thr and rgrow_thr
    plot(sort$ano,
      ylim = c(-5000, 5000), ylab = "Anomaly", xlab = "Scene number", font.sub = 3,
      main = "Disturbance (red) and Regrowth (green) detection", cex.main = 0.9,
      sub = text, col = "grey"
    )
    points(seq(1, length(vvv.ano), 1), vvv.ano)
    abline(h = 0, col = "black")
    abline(v = corr.dist1, col = "red", lty = "dashed")
    abline(v = corr.reg1, col = "green", lty = "dashed")
    abline(v = corr.dist2, col = "red", lty = "dashed")
    abline(v = corr.reg2, col = "green", lty = "dashed")
    abline(v = corr.dist3, col = "red", lty = "dashed")
    abline(v = corr.reg3, col = "green", lty = "dashed")
    abline(v = corr.dist4, col = "red", lty = "dashed")
    abline(v = corr.reg4, col = "green", lty = "dashed")
    # ---------------------------------------------------
    output <- as.numeric(c(
      yy.dist1, val.dist1, yy.reg1, val.reg1,
      yy.dist2, val.dist2, yy.reg2, val.reg2,
      yy.dist3, val.dist3, yy.reg3, val.reg3,
      yy.dist4, val.dist4, yy.reg4, val.reg4
    ))
    output[1:16]
  }


#' @title Avocado algorithm - Wall-to-wall map version
#' @name dist.reg.map
#' @description Continuous vegetation change detection
#' @author  Roberto O. Chávez, Mathieu Decuyper
#' @param s RasterStack of anomalies and their likelihood created with \code{\link{PLUGPhenAnoRFDMapPLUS}}
#' @param dates The julian dates for each scene in the time series brick. Should be in the following format: "YYYY-MM-DD"
#' @param rfd The reference frequency distribution. Determines if an anomaly falls outside the 95 percent of the reference frequency distribution. For example a value that fall in a RDF >= 0.95, indicates that the detected anomaly belongs to the 5 percent of lowest values and is a potential disturbance/change.
#' @param dstrb_thr sets a threshold (number of days) in which if a regrowth is detected within n days after a potential disturbance, then the candidate disturbance date is neglected.
#' @param rgrow_thr sets a threshold (number of days) in which if a disturbance is detected within n days after a potential regrowth, then the candidate regrowth date is neglected.
#' @param cdates Sets the number of consecutive data points for a change to be detected (minimum 2 and maximum 5 data points).
#' @param nCluster Numeric. Number of CPU's to be used for the job.
#' @param outname Character vector with the output path and filename with extension or only the filename and extension if work directory was set. More information: See writeRaster
#' @param datatype Character vector that determines the interpretation of values written to disk. More information: See \code{\link[terra]{writeRaster}}
#' @return The output will consist of a RasterStack with 16 bands: Band1=disturbance detection; Band2=corresponding anomalies band 1; Band3=regrowth detection; Band4=corresponding anomalies band 3; Band5= 2nd disturbance detection for data points that had a previous disturbance and regrowth; Band 6 etc.
#' @seealso \code{\link{dist.reg}}
#' @import npphen
#' @import terra
#' @import RColorBrewer
#'
#' @examples
#' \dontrun{
#' source("dist.reg.map.R") # Load in the mapping function
#' dates <- lan.dates # The dates from your time-series brick (x)
#' ano.rfd.st <- brick("ano.rfd.st.tif") # Load in the anomaly-rfd brick that you created in section 2
#' dist.reg.map(
#'   s = ano.rfd.st, dates = lan.dates, rfd = 0.95, dstrb_thr = 1, rgrow_thr = 730,
#'   nCluster = 1, cdates = 3, outname = "YourDirectory/ChangeMap.tif", 
#'   format = "GTiff", datatype = "INT2S"
#' )
#' }
#' @export
dist.reg.map <-
  function(s, dates, rfd, dstrb_thr, rgrow_thr, cdates, nCluster, outname, datatype) {
    ff <- function(x) {
      if (length(dates) != length(x) / 2) {
        stop("N of dates and files do not match")
      }

      # Checking for a valid cdates value
      if (cdates == 2) {
        cd2 <- 1
      } else {
        cd2 <- 0
      }
      if (cdates == 3) {
        cd3 <- 1
      } else {
        cd3 <- 0
      }
      if (cdates == 4) {
        cd4 <- 1
      } else {
        cd4 <- 0
      }
      if (cdates == 5) {
        cd5 <- 1
      } else {
        cd5 <- 0
      }
      cd.sum <- cd2 + cd3 + cd4 + cd5
      if (cd.sum != 1) {
        stop("cdates argument must be either 2 or 3 or 4 or 5")
      }

      if (all(is.na(x))) {
        return(as.numeric(rep(NA, 16)))
      }

      vv <- as.vector(x)
      ano.ini <- 1
      ano.fin <- length(vv) / 2
      prob.ini <- ano.fin + 1
      prob.fin <- as.numeric(length(vv))
      ano <- vv[ano.ini:ano.fin]
      prob <- vv[prob.ini:prob.fin]
      rfd <- rfd * 100

      # Preparing data -> sorting if the are not sorted already
      df <- data.frame(cbind(ano, prob, dates))
      df$YY <- substr(dates, 1, 4)
      sort <- df[order(df$dates), ]
      sort <- na.omit(sort)
      corr <- seq(1, length(sort$ano), 1)
      sort$corr <- corr
      vvv.ano <- sort$ano
      vvv.prob <- sort$prob
      vvv.ano[vvv.prob < rfd & vvv.ano < 0] <- NA # delete negative anomalies with < rfd
      fin <- length(sort$ano)
      #----------------------------------------
      # Checking for the last period, according to rgrow_thr and dstrb_thr
      last.ddis <- max(df$dates) - dstrb_thr # last date a disturbance can be found
      last.dreg <- max(df$dates) - rgrow_thr # last date a regrowth can be found

      yy.dist1 <- NA
      val.dist1 <- NA
      corr.dist1 <- NA
      dd.dist1 <- NA
      yy.reg1 <- NA
      val.reg1 <- NA
      corr.reg1 <- NA
      dd.reg1 <- NA

      yy.dist2 <- NA
      val.dist2 <- NA
      corr.dist2 <- NA
      dd.dist2 <- NA
      yy.reg2 <- NA
      val.reg2 <- NA
      corr.reg2 <- NA
      dd.reg2 <- NA

      yy.dist3 <- NA
      val.dist3 <- NA
      corr.dist3 <- NA
      dd.dist3 <- NA
      yy.reg3 <- NA
      val.reg3 <- NA
      corr.reg3 <- NA
      dd.reg3 <- NA

      yy.dist4 <- NA
      val.dist4 <- NA
      corr.dist4 <- NA
      dd.dist4 <- NA
      yy.reg4 <- NA
      val.reg4 <- NA
      corr.reg4 <- NA
      dd.reg4 <- NA

      #----------------------------------------------------------
      # First cycle

      # To catch disturbance 1
      val.dist <- sort$ano[which(sort$prob >= rfd & sort$ano < 0)]
      corr.dist <- sort$corr[which(sort$prob >= rfd & sort$ano < 0)]
      yy.dist <- sort$YY[which(sort$prob >= rfd & sort$ano < 0)]
      dd.dist <- sort$dates[which(sort$prob >= rfd & sort$ano < 0)]

      if (length(val.dist) < cdates) {
        return(as.numeric(rep(NA, 16)))
      }
      regrowth <- "off"
      i <- 1

      while (regrowth == "off") {
        if (i == (length(val.dist) - cdates + 1)) {
          regrowth <- "none"
        } else {
          # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
          if (cdates == 2) {
            cond <- "((corr.dist[i]+1)==corr.dist[i+1])"
          }
          if (cdates == 3) {
            cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"
          }
          if (cdates == 4) {
            cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"
          }
          if (cdates == 5) {
            cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"
          }

          if (eval(parse(text = cond))) { # test the selected condition according to cdates

            ssort <- sort[corr.dist[i]:fin, ]
            corr.reg <- ssort$corr[which(ssort$ano >= 0)]
            dd.reg <- ssort$dates[which(ssort$ano >= 0)]

            rdura <- NULL
            for (t in 1:length(dd.reg)) {
              rresta <- dd.reg[t] - dd.dist[i]
              rdura <- c(rdura, rresta)
              w1y.cor.reg <- corr.reg[which(rdura < dstrb_thr)]
            }

            rcorr <- NULL

            if (length(w1y.cor.reg) <= 1) {
              rcorr <- 0
            } else {
              for (k in 1:(length(w1y.cor.reg) - 1)) {
                if ((w1y.cor.reg[k] + 1) == w1y.cor.reg[k + 1]) {
                  rifcor <- 1
                } else {
                  rifcor <- 0
                }
                rcorr <- c(rcorr, rifcor)
              }
            }

            if (sum(rcorr) == 0) { ############# Here makes the condition of not having consecutive regrowth within a year
              val.dist1 <- val.dist[i]
              corr.dist1 <- corr.dist[i]
              yy.dist1 <- yy.dist[i]
              dd.dist1 <- dd.dist[i]
              regrowth <- "on"
            }
          }
        }
        i <- i + 1
      }
      # To catch regrowth 1
      ini1 <- corr.dist[i]
      fin <- length(sort$ano)
      sort2 <- sort[ini1:fin, ]
      j <- 1

      val.reg <- sort2$ano[which(sort2$ano >= 0)]
      corr.reg <- sort2$corr[which(sort2$ano >= 0)]
      yy.reg <- sort2$YY[which(sort2$ano >= 0)]
      dd.reg <- sort2$dates[which(sort2$ano >= 0)]


      if (length(val.reg) < cdates) {
        regrowth <- "off"
      }
      while (regrowth == "on") {
        if (j == (length(val.reg) - cdates + 1)) {
          regrowth <- "off"
        } else {
          # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
          if (cdates == 2) {
            cond <- "((corr.reg[j]+1)==corr.reg[j+1])"
          }
          if (cdates == 3) {
            cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"
          }
          if (cdates == 4) {
            cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"
          }
          if (cdates == 5) {
            cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"
          }

          if (eval(parse(text = cond))) { # test the selected condition according to cdates

            ssort2 <- sort[corr.reg[j]:fin, ]
            corr.dist <- ssort2$corr[which(ssort2$prob >= rfd & ssort2$ano < 0)]
            dd.dist <- ssort2$dates[which(ssort2$prob >= rfd & ssort2$ano < 0)]

            ddura <- NULL
            for (t in 1:length(dd.dist)) {
              resta <- dd.dist[t] - dd.reg[j]
              ddura <- c(ddura, resta)
              w1y.cor.dist <- corr.dist[which(ddura < rgrow_thr)]
            }

            dcorr <- NULL

            if (length(w1y.cor.dist) <= 1) {
              dcorr <- 0
            } else {
              for (k in 1:(length(w1y.cor.dist) - 1)) {
                if ((w1y.cor.dist[k] + 1) == w1y.cor.dist[k + 1]) {
                  ifcor <- 1
                } else {
                  ifcor <- 0
                }
                dcorr <- c(dcorr, ifcor)
              }
            }

            if (sum(dcorr) == 0) { ############# Here makes the condition of not having consecutive disturbances within a year
              val.reg1 <- val.reg[j]
              corr.reg1 <- corr.reg[j]
              yy.reg1 <- yy.reg[j]
              dd.reg1 <- dd.reg[j]
              regrowth <- "off"
            }
          }
        }
        j <- j + 1
      }
      #----------------------------------------------------------
      # Second cycle
      if (is.na(val.reg1) != T) { # Here we continue if there was a regrowth
        ini2 <- corr.reg[j]
        sort3 <- sort[ini2:fin, ]

        # To catch disturbance 2
        val.dist <- sort3$ano[which(sort3$prob >= rfd & sort3$ano < 0)]
        corr.dist <- sort3$corr[which(sort3$prob >= rfd & sort3$ano < 0)]
        yy.dist <- sort3$YY[which(sort3$prob >= rfd & sort3$ano < 0)]
        dd.dist <- sort3$dates[which(sort3$prob >= rfd & sort3$ano < 0)]

        if (length(val.dist) >= cdates) {
          regrowth <- "off"
          i <- 1
          while (regrowth == "off") {
            if (i == (length(val.dist) - cdates + 1)) {
              regrowth <- "none"
            } else {
              # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
              if (cdates == 2) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])"
              }
              if (cdates == 3) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"
              }
              if (cdates == 4) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"
              }
              if (cdates == 5) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"
              }

              if (eval(parse(text = cond))) { # test the selected condition according to cdates

                ssort3 <- sort[corr.dist[i]:fin, ]
                corr.reg <- ssort3$corr[which(ssort3$ano >= 0)]
                dd.reg <- ssort3$dates[which(ssort3$ano >= 0)]

                rdura <- NULL
                for (t in 1:length(dd.reg)) {
                  rresta <- dd.reg[t] - dd.dist[i]
                  rdura <- c(rdura, rresta)
                  w1y.cor.reg <- corr.reg[which(rdura < dstrb_thr)]
                }

                rcorr <- NULL

                if (length(w1y.cor.reg) <= 1) {
                  rcorr <- 0
                } else {
                  for (k in 1:(length(w1y.cor.reg) - 1)) {
                    if ((w1y.cor.reg[k] + 1) == w1y.cor.reg[k + 1]) {
                      rifcor <- 1
                    } else {
                      rifcor <- 0
                    }
                    rcorr <- c(rcorr, rifcor)
                  }
                }

                if (sum(rcorr) == 0) { ############# Here makes the condition of not having consecutive regrowth within a year
                  val.dist2 <- val.dist[i]
                  corr.dist2 <- corr.dist[i]
                  yy.dist2 <- yy.dist[i]
                  dd.dist2 <- dd.dist[i]
                  regrowth <- "on"
                }
              }
            }
            i <- i + 1
          }
          # To catch regrowth 2
          ini3 <- corr.dist[i]
          sort4 <- sort[ini3:fin, ]
          j <- 1
          val.reg <- sort4$ano[which(sort4$ano >= 0)]
          corr.reg <- sort4$corr[which(sort4$ano >= 0)]
          yy.reg <- sort4$YY[which(sort4$ano >= 0)]
          dd.reg <- sort4$dates[which(sort4$ano >= 0)]

          if (length(val.reg) < cdates) {
            regrowth <- "off"
          }
          while (regrowth == "on") {
            if (j == (length(val.reg) - cdates + 1)) {
              regrowth <- "off"
            } else {
              # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
              if (cdates == 2) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])"
              }
              if (cdates == 3) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"
              }
              if (cdates == 4) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"
              }
              if (cdates == 5) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"
              }

              if (eval(parse(text = cond))) { # test the selected condition according to cdates

                ssort4 <- sort[corr.reg[j]:fin, ]
                corr.dist <- ssort4$corr[which(ssort4$prob >= rfd & ssort4$ano < 0)]
                dd.dist <- ssort4$dates[which(ssort4$prob >= rfd & ssort4$ano < 0)]

                ddura <- NULL
                for (t in 1:length(dd.dist)) {
                  resta <- dd.dist[t] - dd.reg[j]
                  ddura <- c(ddura, resta)
                  w1y.cor.dist <- corr.dist[which(ddura < rgrow_thr)]
                }

                dcorr <- NULL

                if (length(w1y.cor.dist) <= 1) {
                  dcorr <- 0
                } else {
                  for (k in 1:(length(w1y.cor.dist) - 1)) {
                    if ((w1y.cor.dist[k] + 1) == w1y.cor.dist[k + 1]) {
                      ifcor <- 1
                    } else {
                      ifcor <- 0
                    }
                    dcorr <- c(dcorr, ifcor)
                  }
                }

                if (sum(dcorr) == 0) { ############# Here makes the condition of not having consecutive disturbances within a year
                  val.reg2 <- val.reg[j]
                  corr.reg2 <- corr.reg[j]
                  yy.reg2 <- yy.reg[j]
                  dd.reg2 <- dd.reg[j]
                  regrowth <- "off"
                }
              }
            }
            j <- j + 1
          }
        }
      }
      #----------------------------------------------------------
      # Third cycle
      if (is.na(val.reg2) != T) { # Here we continue if there was a regrowth
        ini4 <- corr.reg[j]
        sort5 <- sort[ini4:fin, ]

        # To catch disturbance 3
        val.dist <- sort5$ano[which(sort5$prob >= rfd & sort5$ano < 0)]
        corr.dist <- sort5$corr[which(sort5$prob >= rfd & sort5$ano < 0)]
        yy.dist <- sort5$YY[which(sort5$prob >= rfd & sort5$ano < 0)]
        dd.dist <- sort5$dates[which(sort5$prob >= rfd & sort5$ano < 0)]

        if (length(val.dist) >= cdates) {
          regrowth <- "off"
          i <- 1
          while (regrowth == "off") {
            if (i == (length(val.dist) - cdates + 1)) {
              regrowth <- "none"
            } else {
              # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
              if (cdates == 2) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])"
              }
              if (cdates == 3) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"
              }
              if (cdates == 4) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"
              }
              if (cdates == 5) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"
              }

              if (eval(parse(text = cond))) { # test the selected condition according to cdates

                ssort5 <- sort[corr.dist[i]:fin, ]
                corr.reg <- ssort5$corr[which(ssort5$ano >= 0)]
                dd.reg <- ssort5$dates[which(ssort5$ano >= 0)]

                rdura <- NULL
                for (t in 1:length(dd.reg)) {
                  rresta <- dd.reg[t] - dd.dist[i]
                  rdura <- c(rdura, rresta)
                  w1y.cor.reg <- corr.reg[which(rdura < dstrb_thr)]
                }

                rcorr <- NULL

                if (length(w1y.cor.reg) <= 1) {
                  rcorr <- 0
                } else {
                  for (k in 1:(length(w1y.cor.reg) - 1)) {
                    if ((w1y.cor.reg[k] + 1) == w1y.cor.reg[k + 1]) {
                      rifcor <- 1
                    } else {
                      rifcor <- 0
                    }
                    rcorr <- c(rcorr, rifcor)
                  }
                }

                if (sum(rcorr) == 0) { ############# Here makes the condition of not having consecutive regrowth within a year
                  val.dist3 <- val.dist[i]
                  corr.dist3 <- corr.dist[i]
                  yy.dist3 <- yy.dist[i]
                  dd.dist3 <- dd.dist[i]
                  regrowth <- "on"
                }
              }
            }
            i <- i + 1
          }
          # To catch regrowth 3
          ini5 <- corr.dist[i]
          sort6 <- sort[ini5:fin, ]
          j <- 1
          val.reg <- sort6$ano[which(sort6$ano >= 0)]
          corr.reg <- sort6$corr[which(sort6$ano >= 0)]
          yy.reg <- sort6$YY[which(sort6$ano >= 0)]
          dd.reg <- sort6$dates[which(sort6$ano >= 0)]

          if (length(val.reg) < cdates) {
            regrowth <- "off"
          }
          while (regrowth == "on") {
            if (j == (length(val.reg) - cdates + 1)) {
              regrowth <- "off"
            } else {
              # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
              if (cdates == 2) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])"
              }
              if (cdates == 3) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"
              }
              if (cdates == 4) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"
              }
              if (cdates == 5) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"
              }

              if (eval(parse(text = cond))) { # test the selected condition according to cdates

                ssort6 <- sort[corr.reg[j]:fin, ]
                corr.dist <- ssort6$corr[which(ssort6$prob >= rfd & ssort6$ano < 0)]
                dd.dist <- ssort6$dates[which(ssort6$prob >= rfd & ssort6$ano < 0)]

                ddura <- NULL
                for (t in 1:length(dd.dist)) {
                  resta <- dd.dist[t] - dd.reg[j]
                  ddura <- c(ddura, resta)
                  w1y.cor.dist <- corr.dist[which(ddura < rgrow_thr)]
                }

                dcorr <- NULL

                if (length(w1y.cor.dist) <= 1) {
                  dcorr <- 0
                } else {
                  for (k in 1:(length(w1y.cor.dist) - 1)) {
                    if ((w1y.cor.dist[k] + 1) == w1y.cor.dist[k + 1]) {
                      ifcor <- 1
                    } else {
                      ifcor <- 0
                    }
                    dcorr <- c(dcorr, ifcor)
                  }
                }

                if (sum(dcorr) == 0) { ############# Here makes the condition of not having consecutive disturbances within a year
                  val.reg3 <- val.reg[j]
                  corr.reg3 <- corr.reg[j]
                  yy.reg3 <- yy.reg[j]
                  dd.reg3 <- dd.reg[j]
                  regrowth <- "off"
                }
              }
            }
            j <- j + 1
          }
        }
      }
      #----------------------------------------------------------
      # Fourth cycle
      if (is.na(val.reg3) != T) { # Here we continue if there was a regrowth
        ini6 <- corr.reg[j]
        sort7 <- sort[ini6:fin, ]
        # To catch disturbance 4
        val.dist <- sort7$ano[which(sort7$prob >= rfd & sort7$ano < 0)]
        corr.dist <- sort7$corr[which(sort7$prob >= rfd & sort7$ano < 0)]
        yy.dist <- sort7$YY[which(sort7$prob >= rfd & sort7$ano < 0)]
        dd.dist <- sort7$dates[which(sort7$prob >= rfd & sort7$ano < 0)]

        if (length(val.dist) >= cdates) {
          regrowth <- "off"
          i <- 1
          while (regrowth == "off") {
            if (i == (length(val.dist) - cdates + 1)) {
              regrowth <- "none"
            } else {
              # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
              if (cdates == 2) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])"
              }
              if (cdates == 3) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"
              }
              if (cdates == 4) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"
              }
              if (cdates == 5) {
                cond <- "((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"
              }

              if (eval(parse(text = cond))) { # test the selected condition according to cdates

                ssort7 <- sort[corr.dist[i]:fin, ]
                corr.reg <- ssort7$corr[which(ssort7$ano >= 0)]
                dd.reg <- ssort7$dates[which(ssort7$ano >= 0)]

                rdura <- NULL
                for (t in 1:length(dd.reg)) {
                  rresta <- dd.reg[t] - dd.dist[i]
                  rdura <- c(rdura, rresta)
                  w1y.cor.reg <- corr.reg[which(rdura < dstrb_thr)]
                }

                rcorr <- NULL

                if (length(w1y.cor.reg) <= 1) {
                  rcorr <- 0
                } else {
                  for (k in 1:(length(w1y.cor.reg) - 1)) {
                    if ((w1y.cor.reg[k] + 1) == w1y.cor.reg[k + 1]) {
                      rifcor <- 1
                    } else {
                      rifcor <- 0
                    }
                    rcorr <- c(rcorr, rifcor)
                  }
                }

                if (sum(rcorr) == 0) { ############# Here makes the condition of not having consecutive regrowth within a year
                  val.dist4 <- val.dist[i]
                  corr.dist4 <- corr.dist[i]
                  yy.dist4 <- yy.dist[i]
                  dd.dist4 <- dd.dist[i]
                  regrowth <- "on"
                }
              }
            }
            i <- i + 1
          }
          # To catch regrowth 4
          ini7 <- corr.dist[i]
          sort8 <- sort[ini7:fin, ]
          j <- 1
          val.reg <- sort8$ano[which(sort8$ano >= 0)]
          corr.reg <- sort8$corr[which(sort8$ano >= 0)]
          yy.reg <- sort8$YY[which(sort8$ano >= 0)]
          dd.reg <- sort8$dates[which(sort8$ano >= 0)]

          if (length(val.reg) < cdates) {
            regrowth <- "off"
          }
          while (regrowth == "on") {
            if (j == (length(val.reg) - cdates + 1)) {
              regrowth <- "off"
            } else {
              # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
              if (cdates == 2) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])"
              }
              if (cdates == 3) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"
              }
              if (cdates == 4) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"
              }
              if (cdates == 5) {
                cond <- "((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"
              }

              if (eval(parse(text = cond))) { # test the selected condition according to cdates

                ssort8 <- sort[corr.reg[j]:fin, ]
                corr.dist <- ssort8$corr[which(ssort8$prob >= rfd & ssort8$ano < 0)]
                dd.dist <- ssort8$dates[which(ssort8$prob >= rfd & ssort8$ano < 0)]

                ddura <- NULL
                for (t in 1:length(dd.dist)) {
                  resta <- dd.dist[t] - dd.reg[j]
                  ddura <- c(ddura, resta)
                  w1y.cor.dist <- corr.dist[which(ddura < rgrow_thr)]
                }

                dcorr <- NULL

                if (length(w1y.cor.dist) <= 1) {
                  dcorr <- 0
                } else {
                  for (k in 1:(length(w1y.cor.dist) - 1)) {
                    if ((w1y.cor.dist[k] + 1) == w1y.cor.dist[k + 1]) {
                      ifcor <- 1
                    } else {
                      ifcor <- 0
                    }
                    dcorr <- c(dcorr, ifcor)
                  }
                }

                if (sum(dcorr) == 0) { ############# Here makes the condition of not having consecutive disturbances within a year
                  val.reg4 <- val.reg[j]
                  corr.reg4 <- corr.reg[j]
                  yy.reg4 <- yy.reg[j]
                  dd.reg4 <- dd.reg[j]
                  regrowth <- "off"
                }
              }
            }
            j <- j + 1
          }
        }
      }
      # ---------------------------------------------------
      # Deleting outputs in the last rgrow_thr y dstrb_thr period
      if (dd.dist1 >= last.ddis | is.na(dd.dist1)) {
        yy.dist1 <- NA
        val.dist1 <- NA
        corr.dist1 <- NA
      }
      if (dd.reg1 >= last.dreg | is.na(dd.reg1)) {
        yy.reg1 <- NA
        val.reg1 <- NA
        corr.reg1 <- NA
      }

      if (dd.dist2 >= last.ddis | is.na(dd.dist2)) {
        yy.dist2 <- NA
        val.dist2 <- NA
        corr.dist2 <- NA
      }
      if (dd.reg2 >= last.dreg | is.na(dd.reg2)) {
        yy.reg2 <- NA
        val.reg2 <- NA
        corr.reg2 <- NA
      }

      if (dd.dist3 >= last.ddis | is.na(dd.dist3)) {
        yy.dist3 <- NA
        val.dist3 <- NA
        corr.dist3 <- NA
      }
      if (dd.reg3 >= last.dreg | is.na(dd.reg3)) {
        yy.reg3 <- NA
        val.reg3 <- NA
        corr.reg3 <- NA
      }

      if (dd.dist4 >= last.ddis | is.na(dd.dist4)) {
        yy.dist4 <- NA
        val.dist4 <- NA
        corr.dist4 <- NA
      }
      if (dd.reg4 >= last.dreg | is.na(dd.reg4)) {
        yy.reg4 <- NA
        val.reg4 <- NA
        corr.reg4 <- NA
      }


      # ---------------------------------------------------
      output <- as.numeric(c(
        yy.dist1, val.dist1, yy.reg1, val.reg1,
        yy.dist2, val.dist2, yy.reg2, val.reg2,
        yy.dist3, val.dist3, yy.reg3, val.reg3,
        yy.dist4, val.dist4, yy.reg4, val.reg4
      ))
      output[1:16]
    }

    #----------------------------------------------------------------------------------------
    # Function using n clusters
    app(s, fun = ff, filename = outname, cores = nCluster, overwrite = T, wopt = list(datatype = datatype))
  }
