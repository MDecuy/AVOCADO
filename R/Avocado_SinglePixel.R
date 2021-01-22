#' @title Avocado algorithm - Single pixel version
#' 
#' @description Continuous vegetation change detection 
#' 
#' @author  Roberto O. Chavez, Mathieu Decuyper

#' @section 1: Create the reference vegetatation curve (phen).
#' 
#' @param x RasterBrick. Image time series brick
#' @param dates A date vector. The number of dates must be equal to the number of layers of x.
#' @param loc Location. SpatialPolygons overlapping one or multiple pixels.
#' 
#' @section 2: Calculate the Anomaly and probability values (based on the difference between the pixel values and reference vegetation)
#' 
#' @param x RasterBrick. Image time series brick
#' @param anop Numeric vector with the number of layers that are in x. For those layers the anomalies and probabilities will be calculated based on the phen. For example a raster stack has 450 layers, so anop=c(1:450).
#' @param phen Numeric vector with the values of the reference vegetation
#' @param h Numeric. Geographic hemisphere to define the starting date of the growing season. h=1 Northern Hemisphere; h=2 Southern Hemisphere.

#' @section 3: Automatic detection of the disturbance and regrowth events
#' 
#' @param s Anomaly and probability time-series created in section 2
#' @param dates The julian dates for each scene in the RasterBrick (x). Should be in the following format: "YYYY-MM-DD"
#' @param p The probability (%) of the anomalies being a candidate disturbance/change data point
#' @param cdates Sets the number of consecutive data points for a change to be detected (minimum 2 and maximum 5 data points).
#' @param rth sets a threshold (number of days) in which if a regrowth is detected within n days after a potential disturbance, then the candidate disturbance date is neglected.
#' @param dth sets a threshold (number of days) in which if a disturbance is detected within n days after a potential regrowth, then the candidate regrowth date is neglected.

#' @import npphen
#' @import raster
#' @import RColorBrewer
#' @import rgdal
#' @import bfastSpatial
#' @import lubridate
#' @import rts
#' 
#' @examples
#' # load MDD dataset
#' data(MDD)
#' data(RefForest)
#' 
#' ############ Section 1: Reference vegetatation curve ################### 
#' 
#' x <- brick("MDD")
#' lan.info <- getSceneinfo(names(x))
#' lan.dates <-as.Date(lan.info$date)
#' loc <- readOGR(dsn=path.expand("YourDirectory"), layer="RefForest")
#' #loc <- spTransform(loc, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #in case of different coordinate systems
#' ref.ext<-extent(loc)
#' ref.brick <- crop(x,ref.ext)
#' fin <- nrow(ref.brick)*ncol(ref.brick)
#' phen <- extract(ref.brick,1)
#' d1 <- lan.dates

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
#' #' PhenKplot(phen,d1,h=1,nGS=365, xlab="DOY",ylab="NDMI",rge=c(0,10000))

#' ############ Section 2: Calculate anomaly and probability values ###################
#' #Function - PLUGPhenAnoProbPLUS
#' 
PLUGPhenAnoProbPLUS <-
function(x,phen,dates,h,anop,rge) {
  
  # a.Preparing dataset
  
  if(length(rge)!=2){stop("rge must be a vector of length 2")}
  if(rge[1]>rge[2]){stop("rge vector order must be minimum/maximum")}
  if(length(dates)!=length(x)){stop("N of dates and files do not match")}
  
  #ref.min <- min(refp)
  #ref.max <- max(refp)
  ano.min <- min(anop)
  ano.max <- max(anop)
  ano.len <- ano.max-ano.min+1
  len2 <- 2*ano.len
  
  if(ano.min>=ano.max){stop("for anop, lower value > upper value")}
  
  if (all(is.na(x))) {
    return(rep(NA,len2))
  }
  
  DOY <- yday(dates)
  DOY[which(DOY==366)]<-365	
  D1<-cbind(DOY,phen)
  D2<-cbind(DOY[ano.min:ano.max],x[ano.min:ano.max])
  
  if(length(unique(D1[,2]))<10 | (nrow(D1)-sum(is.na(D1)))<(0.1*length(D1))) {
    return(rep(NA,len2))
  }
  
  if (all(is.na(D2[,2]))) {
    return(rep(NA,len2))
  }
  
  # b. Kernel calculation using the reference period (D1)
  if(h!=1 && h!=2){stop("Invalid h")}
  DOGS<-cbind(seq(1,365),c(seq(185,365),seq(1,184)))
  if(h==2){
    for(i in 1:nrow(D1)){
      D1[i,1]<-DOGS[which(DOGS[,1]==D1[i,1],arr.ind=TRUE),2]}}
  
  Hmat<-Hpi(na.omit(D1))
  Hmat[1,2]<-Hmat[2,1]
  K1<-kde(na.omit(D1),H=Hmat,xmin=c(1,rge[1]),xmax=c(365,rge[2]),gridsize=c(365,500))
  K1Con<-K1$estimate
  for(j in 1:365){
    Kdiv<-sum(K1$estimate[j,])
    ifelse(Kdiv==0,K1Con[j,]<-0,K1Con[j,]<-K1$estimate[j,]/sum(K1$estimate[j,]))}
  MAXY<-apply(K1Con,1,max)
  for(i in 1:365){
    MAXY[i]<-median(K1$eval.points[[2]][which(K1Con[i,]==MAXY[i],arr.ind=TRUE)])}
  
  # c. Calculating cumulative bivariate density distribution
  
  h2d <- list()
  h2d$x <- seq(1,365)
  h2d$y <- seq(rge[1],rge[2],len=500)
  h2d$density <- K1Con/sum(K1Con)
  uniqueVals <- rev(unique(sort(h2d$density)))
  cumProbs <- cumsum(uniqueVals)
  names(cumProbs) <- uniqueVals
  h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
  h2d$cumDensity[] <- cumProbs[as.character(h2d$density)]
  
  # c. Anomaly calculation (D2)
  
  if(h==2){
    for(i in 1:nrow(D2)){
      D2[i,1]<-DOGS[which(DOGS[,1]==D2[i,1],arr.ind=TRUE),2]}}
  
  # d. Calculating anomalies AND probabilities based on D2
  
  # d.1 Anomalies
  Anoma<-rep(NA,ano.len)
  for(i in 1:nrow(D2)){
    Anoma[i]<-as.integer(D2[i,2]-MAXY[D2[i,1]])}
  Anoma[1:ano.len]
  
  # d.2 Probabilities
  rowAnom<-matrix(NA,nrow=nrow(D2),ncol=500)
  for(i in 1:nrow(D2)){
    rowAnom[i,]<-abs(h2d$y-D2[i,2])}
  rowAnom2<-unlist(apply(rowAnom,1,function(x) {if(all(is.na(x))) {NA} else {which.min(x)}}))
  AnomProb<-rep(NA,nrow(D2))
  for( i in 1:nrow(D2)){
    AnomProb[i]<-h2d$cumDensity[D2[i,1],rowAnom2[i]]}
  AnomProbPerc<-round(100*(AnomProb))
  plot(AnomProbPerc[1:ano.len], xlab = "Time", ylab = "Probability (%)", font.lab = 2)
  abline(h = 90, col = "red")
  abline(h = 95, col = "red")
  abline(h = 99, col = "red")
  AnomProbPerc[1:ano.len]
  # Two results in a single numeric vector
  AnomPLUSProb <- c(Anoma[1:ano.len],AnomProbPerc[1:ano.len])
  AnomPLUSProb[1:len2]
}

#' @examples
#sample.rts <- rts(ndmibrick,lan.dates)
#pix.num <- cellFromXY(sample.rts,c(X,Ycoordinates))
#ts.inter <- extract(sample.rts,pix.num)
#vec <- as.vector(ts.inter)
#dd <- lan.dates # NON SORTED DATES
#nn <- seq(1,length(vec),1)
#df <-data.frame(cbind(nn,dd))
#sort <-df[order(as.Date(df$dd)),]
#sort$vec <- vec
#sort$date <- as.Date(sort$dd)
#corr <- seq(1,length(vec),1)
#sort$corr <- corr
#x <- sort$vec
#dates <- as.Date(sort$date)
#ano_prob <- PLUGPhenAnoProbPLUS(x=x,phen=phen,dates=dates,h=1,anop=c(1:nlayers),rge=c(1,10000))

#' ############ Section 3: Automatic detection of the disturbance and regrowth ################### 

#' #Function - dist.reg
dist.reg <-
function(x,dates,p,dth,rth,cdates) {
  if (length(dates) != length(x)/2) {
    stop("N of dates and files do not match")
  }
  
  # Checking for a valid cdates value
  if(cdates==2){cd2<-1}else{cd2<-0}
  if(cdates==3){cd3<-1}else{cd3<-0}
  if(cdates==4){cd4<-1}else{cd4<-0}
  if(cdates==5){cd5<-1}else{cd5<-0}
  cd.sum <- cd2+cd3+cd4+cd5
  if (cd.sum!=1) {
    stop("cdates argument must be either 2 or 3 or 4 or 5")
  }
  
  if (all(is.na(x))) {
    return(as.numeric(rep(NA, 16)))
  }
  
  vv <- as.vector(x)
  ano.ini <-1
  ano.fin <-length(vv)/2
  prob.ini <- ano.fin+1
  prob.fin <- as.numeric(length(vv))
  ano  <- vv[ano.ini:ano.fin]
  prob <- vv[prob.ini:prob.fin]
  
  #Preparing data -> sorting if the are not sorted already
  df <-data.frame(cbind(ano,prob,dates))
  df$YY <- substr(dates,1,4)
  sort <-df[order(df$dates),]
  sort <- na.omit(sort)
  corr <- seq(1,length(sort$ano),1)
  sort$corr <- corr
  vvv.ano <- sort$ano
  vvv.prob <- sort$prob
  vvv.ano[vvv.prob<p&vvv.ano<0]<-NA #delete non-significant negative anomalies
  text <- paste("grey are non-significant negative anomalies with p<",p,sep="")
  plot(sort$ano,  ylim=c(-5000,5000), ylab="Anomaly",xlab="Scene number",font.sub=3,
       main="Disturbance (red) and Regrowth (green) detection", cex.main=0.9,
       sub=text, col="grey")
  points(seq(1,length(vvv.ano),1),vvv.ano)
  abline(h=0, col="red")
  fin <- length(sort$ano)
  #----------------------------------------
  # Checking for the last period, according to dth and rth
  last.ddis <- max(df$dates)-rth
  last.dreg <- max(df$dates)-dth
  
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
  val.dist <- sort$ano[which(sort$prob>=p & sort$ano<0)]
  corr.dist <-sort$corr[which(sort$prob>=p & sort$ano<0)]
  yy.dist <-  sort$YY[which(sort$prob>=p & sort$ano<0)]
  dd.dist <- sort$dates[which(sort$prob>=p & sort$ano<0)]
  
  if (length(val.dist)<cdates) {
    return(as.numeric(rep(NA, 16)))
  }
  regrowth<-"off"
  i <-1
  
  while (regrowth=="off") {
    if(i==(length(val.dist)-cdates+1)){regrowth<-"none"} else {
      
      # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
      if(cdates==2){cond<-"((corr.dist[i]+1)==corr.dist[i+1])"}
      if(cdates==3){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"}
      if(cdates==4){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"}
      if(cdates==5){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"}
      
      if(eval(parse(text=cond))) { # test the selected condition according to cdates
        
        ssort <- sort[corr.dist[i]:fin,]
        corr.reg <- ssort$corr[which(ssort$ano>=0)]
        dd.reg <- ssort$dates[which(ssort$ano>=0)]
        
        rdura <- NULL
        for(t in 1:length(dd.reg)) {
          rresta <-dd.reg[t]-dd.dist[i]
          rdura <-c(rdura,rresta)
          w1y.cor.reg <- corr.reg[which(rdura<rth)]
        }
        
        rcorr <-NULL
        
        if(length(w1y.cor.reg)<=1){rcorr<-0} else {
          for(k in 1:(length(w1y.cor.reg)-1)) {
            if((w1y.cor.reg[k]+1)==w1y.cor.reg[k+1]){
              rifcor <- 1
            } else {
              rifcor <- 0
            }
            rcorr <-c(rcorr,rifcor)
          }
        }
        
        if(sum(rcorr)==0) {############# Here makes the condition of not having consecutive regrowth within a year
          val.dist1 <- val.dist[i]
          corr.dist1 <- corr.dist[i]
          yy.dist1 <- yy.dist[i]
          dd.dist1 <- dd.dist[i]
          regrowth<-"on"
        }
      }
    }
    i <- i+1
  }
  abline(v=corr.dist1, col="red",lty="dashed")
  
  # To catch regrowth 1
  ini1 <- corr.dist[i]
  fin <- length(sort$ano)
  sort2 <- sort[ini1:fin,]
  j <-1
  
  val.reg <- sort2$ano[which(sort2$ano>=0)]
  corr.reg <- sort2$corr[which(sort2$ano>=0)]
  yy.reg <- sort2$YY[which(sort2$ano>=0)]
  dd.reg <- sort2$dates[which(sort2$ano>=0)]
  
  
  if (length(val.reg)<cdates) {
    regrowth<-"off"
  }
  while (regrowth=="on") {
    if(j==(length(val.reg)-cdates+1)){regrowth<-"off"} else {
      
      # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
      if(cdates==2){cond<-"((corr.reg[j]+1)==corr.reg[j+1])"}
      if(cdates==3){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"}
      if(cdates==4){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"}
      if(cdates==5){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"}
      
      if(eval(parse(text=cond))) { # test the selected condition according to cdates
        
        ssort2 <- sort[corr.reg[j]:fin,]
        corr.dist <-ssort2$corr[which(ssort2$prob>=p & ssort2$ano<0)]
        dd.dist <- ssort2$dates[which(ssort2$prob>=p & ssort2$ano<0)]
        
        ddura <- NULL
        for(t in 1:length(dd.dist)) {
          resta <-dd.dist[t]-dd.reg[j]
          ddura <-c(ddura,resta)
          w1y.cor.dist <- corr.dist[which(ddura<dth)]
        }
        
        dcorr <-NULL
        
        if(length(w1y.cor.dist)<=1){dcorr<-0} else {
          for(k in 1:(length(w1y.cor.dist)-1)) {
            if((w1y.cor.dist[k]+1)==w1y.cor.dist[k+1]){
              ifcor <- 1
            } else {
              ifcor <- 0
            }
            dcorr <-c(dcorr,ifcor)
          }
        }
        
        if(sum(dcorr)==0) {############# Here makes the condition of not having consecutive disturbances within a year
          val.reg1 <- val.reg[j]
          corr.reg1 <- corr.reg[j]
          yy.reg1 <- yy.reg[j]
          dd.reg1 <- dd.reg[j]
          regrowth<-"off"
        }
      }
    }
    j <- j+1
  }
  abline(v=corr.reg1, col="green",lty="dashed")
  
  #----------------------------------------------------------
  # Second cycle
  if(is.na(val.reg1)!=T) {      # Here we continue if there was a regrowth
    ini2 <- corr.reg[j]
    sort3 <- sort[ini2:fin,]
    
    # To catch disturbance 2
    val.dist <- sort3$ano[which(sort3$prob>=p & sort3$ano<0)]
    corr.dist <-sort3$corr[which(sort3$prob>=p & sort3$ano<0)]
    yy.dist <-  sort3$YY[which(sort3$prob>=p & sort3$ano<0)]
    dd.dist <- sort3$dates[which(sort3$prob>=p & sort3$ano<0)]
    
    if (length(val.dist)>=cdates) {
      regrowth<-"off"
      i <-1
      while (regrowth=="off") {
        if(i==(length(val.dist)-cdates+1)){regrowth<-"none"} else {
          
          # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
          if(cdates==2){cond<-"((corr.dist[i]+1)==corr.dist[i+1])"}
          if(cdates==3){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"}
          if(cdates==4){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"}
          if(cdates==5){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"}
          
          if(eval(parse(text=cond))) { # test the selected condition according to cdates
            
            ssort3 <- sort[corr.dist[i]:fin,]
            corr.reg <- ssort3$corr[which(ssort3$ano>=0)]
            dd.reg <- ssort3$dates[which(ssort3$ano>=0)]
            
            rdura <- NULL
            for(t in 1:length(dd.reg)) {
              rresta <-dd.reg[t]-dd.dist[i]
              rdura <-c(rdura,rresta)
              w1y.cor.reg <- corr.reg[which(rdura<rth)]
            }
            
            rcorr <-NULL
            
            if(length(w1y.cor.reg)<=1){rcorr<-0} else {
              for(k in 1:(length(w1y.cor.reg)-1)) {
                if((w1y.cor.reg[k]+1)==w1y.cor.reg[k+1]){
                  rifcor <- 1
                } else {
                  rifcor <- 0
                }
                rcorr <-c(rcorr,rifcor)
              }
            }
            
            if(sum(rcorr)==0) {############# Here makes the condition of not having consecutive regrowth within a year
              val.dist2 <- val.dist[i]
              corr.dist2 <- corr.dist[i]
              yy.dist2 <- yy.dist[i]
              dd.dist2 <- dd.dist[i]
              regrowth<-"on"
            }
          }
        }
        i <- i+1
      }
      abline(v=corr.dist2, col="red",lty="dashed")
      
      # To catch regrowth 2
      ini3 <- corr.dist[i]
      sort4 <- sort[ini3:fin,]
      j <-1
      val.reg <- sort4$ano[which(sort4$ano>=0)]
      corr.reg <- sort4$corr[which(sort4$ano>=0)]
      yy.reg <- sort4$YY[which(sort4$ano>=0)]
      dd.reg <- sort4$dates[which(sort4$ano>=0)]
      
      if (length(val.reg)<cdates) {
        regrowth<-"off"
      }
      while (regrowth=="on") {
        if(j==(length(val.reg)-cdates+1)){regrowth<-"off"} else {
          
          # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
          if(cdates==2){cond<-"((corr.reg[j]+1)==corr.reg[j+1])"}
          if(cdates==3){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"}
          if(cdates==4){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"}
          if(cdates==5){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"}
          
          if(eval(parse(text=cond))) { # test the selected condition according to cdates
            
            ssort4 <- sort[corr.reg[j]:fin,]
            corr.dist <-ssort4$corr[which(ssort4$prob>=p & ssort4$ano<0)]
            dd.dist <- ssort4$dates[which(ssort4$prob>=p & ssort4$ano<0)]
            
            ddura <- NULL
            for(t in 1:length(dd.dist)) {
              resta <-dd.dist[t]-dd.reg[j]
              ddura <-c(ddura,resta)
              w1y.cor.dist <- corr.dist[which(ddura<dth)]
            }
            
            dcorr <-NULL
            
            if(length(w1y.cor.dist)<=1){dcorr<-0} else {
              for(k in 1:(length(w1y.cor.dist)-1)) {
                if((w1y.cor.dist[k]+1)==w1y.cor.dist[k+1]){
                  ifcor <- 1
                } else {
                  ifcor <- 0
                }
                dcorr <-c(dcorr,ifcor)
              }
            }
            
            if(sum(dcorr)==0) {############# Here makes the condition of not having consecutive disturbances within a year
              val.reg2 <- val.reg[j]
              corr.reg2 <- corr.reg[j]
              yy.reg2 <- yy.reg[j]
              dd.reg2 <- dd.reg[j]
              regrowth<-"off"
            }
          }
        }
        j <- j+1
      }
      abline(v=corr.reg2, col="green",lty="dashed")
    }
  }
  #----------------------------------------------------------
  # Third cycle
  if(is.na(val.reg2)!=T) {      # Here we continue if there was a regrowth
    ini4 <- corr.reg[j]
    sort5 <- sort[ini4:fin,]
    
    # To catch disturbance 3
    val.dist <- sort5$ano[which(sort5$prob>=p & sort5$ano<0)]
    corr.dist <-sort5$corr[which(sort5$prob>=p & sort5$ano<0)]
    yy.dist <-  sort5$YY[which(sort5$prob>=p & sort5$ano<0)]
    dd.dist <- sort5$dates[which(sort5$prob>=p & sort5$ano<0)]
    
    if (length(val.dist)>=cdates) {
      regrowth<-"off"
      i <-1
      while (regrowth=="off") {
        if(i==(length(val.dist)-cdates+1)){regrowth<-"none"} else {
          
          # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
          if(cdates==2){cond<-"((corr.dist[i]+1)==corr.dist[i+1])"}
          if(cdates==3){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"}
          if(cdates==4){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"}
          if(cdates==5){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"}
          
          if(eval(parse(text=cond))) { # test the selected condition according to cdates
            
            ssort5 <- sort[corr.dist[i]:fin,]
            corr.reg <- ssort5$corr[which(ssort5$ano>=0)]
            dd.reg <- ssort5$dates[which(ssort5$ano>=0)]
            
            rdura <- NULL
            for(t in 1:length(dd.reg)) {
              rresta <-dd.reg[t]-dd.dist[i]
              rdura <-c(rdura,rresta)
              w1y.cor.reg <- corr.reg[which(rdura<rth)]
            }
            
            rcorr <-NULL
            
            if(length(w1y.cor.reg)<=1){rcorr<-0} else {
              for(k in 1:(length(w1y.cor.reg)-1)) {
                if((w1y.cor.reg[k]+1)==w1y.cor.reg[k+1]){
                  rifcor <- 1
                } else {
                  rifcor <- 0
                }
                rcorr <-c(rcorr,rifcor)
              }
            }
            
            if(sum(rcorr)==0) {############# Here makes the condition of not having consecutive regrowth within a year
              val.dist3 <- val.dist[i]
              corr.dist3 <- corr.dist[i]
              yy.dist3 <- yy.dist[i]
              dd.dist3 <- dd.dist[i]
              regrowth<-"on"
            }
          }
        }
        i <- i+1
      }
      abline(v=corr.dist3, col="red",lty="dashed")
      
      # To catch regrowth 3
      ini5 <- corr.dist[i]
      sort6 <- sort[ini5:fin,]
      j <-1
      val.reg <- sort6$ano[which(sort6$ano>=0)]
      corr.reg <- sort6$corr[which(sort6$ano>=0)]
      yy.reg <- sort6$YY[which(sort6$ano>=0)]
      dd.reg <- sort6$dates[which(sort6$ano>=0)]
      
      if (length(val.reg)<cdates) {
        regrowth<-"off"
      }
      while (regrowth=="on") {
        if(j==(length(val.reg)-cdates+1)){regrowth<-"off"} else {
          
          # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
          if(cdates==2){cond<-"((corr.reg[j]+1)==corr.reg[j+1])"}
          if(cdates==3){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"}
          if(cdates==4){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"}
          if(cdates==5){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"}
          
          if(eval(parse(text=cond))) { # test the selected condition according to cdates
            
            ssort6 <- sort[corr.reg[j]:fin,]
            corr.dist <-ssort6$corr[which(ssort6$prob>=p & ssort6$ano<0)]
            dd.dist <- ssort6$dates[which(ssort6$prob>=p & ssort6$ano<0)]
            
            ddura <- NULL
            for(t in 1:length(dd.dist)) {
              resta <-dd.dist[t]-dd.reg[j]
              ddura <-c(ddura,resta)
              w1y.cor.dist <- corr.dist[which(ddura<dth)]
            }
            
            dcorr <-NULL
            
            if(length(w1y.cor.dist)<=1){dcorr<-0} else {
              for(k in 1:(length(w1y.cor.dist)-1)) {
                if((w1y.cor.dist[k]+1)==w1y.cor.dist[k+1]){
                  ifcor <- 1
                } else {
                  ifcor <- 0
                }
                dcorr <-c(dcorr,ifcor)
              }
            }
            
            if(sum(dcorr)==0) {############# Here makes the condition of not having consecutive disturbances within a year
              val.reg3 <- val.reg[j]
              corr.reg3 <- corr.reg[j]
              yy.reg3 <- yy.reg[j]
              dd.reg3 <- dd.reg[j]
              regrowth<-"off"
            }
          }
        }
        j <- j+1
      }
      abline(v=corr.reg3, col="green",lty="dashed")
    }
  }
  #----------------------------------------------------------
  # Fourth cycle
  if(is.na(val.reg3)!=T) {      # Here we continue if there was a regrowth
    ini6 <- corr.reg[j]
    sort7 <- sort[ini6:fin,]
    # To catch disturbance 4
    val.dist <- sort7$ano[which(sort7$prob>=p & sort7$ano<0)]
    corr.dist <-sort7$corr[which(sort7$prob>=p & sort7$ano<0)]
    yy.dist <-  sort7$YY[which(sort7$prob>=p & sort7$ano<0)]
    dd.dist <- sort7$dates[which(sort7$prob>=p & sort7$ano<0)]
    
    if (length(val.dist)>=cdates) {
      regrowth<-"off"
      i <-1
      while (regrowth=="off") {
        if(i==(length(val.dist)-cdates+1)){regrowth<-"none"} else {
          
          # Consecutive dates to flag dist, based on the "cdates" argument # possible values are 2,3,4,5
          if(cdates==2){cond<-"((corr.dist[i]+1)==corr.dist[i+1])"}
          if(cdates==3){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])"}
          if(cdates==4){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])"}
          if(cdates==5){cond<-"((corr.dist[i]+1)==corr.dist[i+1])&((corr.dist[i]+2)==corr.dist[i+2])&((corr.dist[i]+3)==corr.dist[i+3])&((corr.dist[i]+4)==corr.dist[i+4])"}
          
          if(eval(parse(text=cond))) { # test the selected condition according to cdates
            
            ssort7 <- sort[corr.dist[i]:fin,]
            corr.reg <- ssort7$corr[which(ssort7$ano>=0)]
            dd.reg <- ssort7$dates[which(ssort7$ano>=0)]
            
            rdura <- NULL
            for(t in 1:length(dd.reg)) {
              rresta <-dd.reg[t]-dd.dist[i]
              rdura <-c(rdura,rresta)
              w1y.cor.reg <- corr.reg[which(rdura<rth)]
            }
            
            rcorr <-NULL
            
            if(length(w1y.cor.reg)<=1){rcorr<-0} else {
              for(k in 1:(length(w1y.cor.reg)-1)) {
                if((w1y.cor.reg[k]+1)==w1y.cor.reg[k+1]){
                  rifcor <- 1
                } else {
                  rifcor <- 0
                }
                rcorr <-c(rcorr,rifcor)
              }
            }
            
            if(sum(rcorr)==0) {############# Here makes the condition of not having consecutive regrowth within a year
              val.dist4 <- val.dist[i]
              corr.dist4 <- corr.dist[i]
              yy.dist4 <- yy.dist[i]
              dd.dist4 <- dd.dist[i]
              regrowth<-"on"
            }
          }
        }
        i <- i+1
      }
      abline(v=corr.dist4, col="red",lty="dashed")
      
      # To catch regrowth 4
      ini7 <- corr.dist[i]
      sort8 <- sort[ini7:fin,]
      j <-1
      val.reg <- sort8$ano[which(sort8$ano>=0)]
      corr.reg <- sort8$corr[which(sort8$ano>=0)]
      yy.reg <- sort8$YY[which(sort8$ano>=0)]
      dd.reg <- sort8$dates[which(sort8$ano>=0)]
      
      if (length(val.reg)<cdates) {
        regrowth<-"off"
      }
      while (regrowth=="on") {
        if(j==(length(val.reg)-cdates+1)){regrowth<-"off"} else {
          
          # Consecutive dates to flag regrowth, based on the "cdates" argument # possible values are 2,3,4,5
          if(cdates==2){cond<-"((corr.reg[j]+1)==corr.reg[j+1])"}
          if(cdates==3){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])"}
          if(cdates==4){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])"}
          if(cdates==5){cond<-"((corr.reg[j]+1)==corr.reg[j+1])&((corr.reg[j]+2)==corr.reg[j+2])&((corr.reg[j]+3)==corr.reg[j+3])&((corr.reg[j]+4)==corr.reg[j+4])"}
          
          if(eval(parse(text=cond))) { # test the selected condition according to cdates
            
            ssort8 <- sort[corr.reg[j]:fin,]
            corr.dist <-ssort8$corr[which(ssort8$prob>=p & ssort8$ano<0)]
            dd.dist <- ssort8$dates[which(ssort8$prob>=p & ssort8$ano<0)]
            
            ddura <- NULL
            for(t in 1:length(dd.dist)) {
              resta <-dd.dist[t]-dd.reg[j]
              ddura <-c(ddura,resta)
              w1y.cor.dist <- corr.dist[which(ddura<dth)]
            }
            
            dcorr <-NULL
            
            if(length(w1y.cor.dist)<=1){dcorr<-0} else {
              for(k in 1:(length(w1y.cor.dist)-1)) {
                if((w1y.cor.dist[k]+1)==w1y.cor.dist[k+1]){
                  ifcor <- 1
                } else {
                  ifcor <- 0
                }
                dcorr <-c(dcorr,ifcor)
              }
            }
            
            if(sum(dcorr)==0) {############# Here makes the condition of not having consecutive disturbances within a year
              val.reg4 <- val.reg[j]
              corr.reg4 <- corr.reg[j]
              yy.reg4 <- yy.reg[j]
              dd.reg4 <- dd.reg[j]
              regrowth<-"off"
            }
          }
        }
        j <- j+1
      }
      abline(v=corr.reg4, col="green",lty="dashed")
    }
  }
  # ---------------------------------------------------
  # Deleting outputs in the last dth y rth period
  if(dd.dist1 >= last.ddis | is.na(dd.dist1)){
    yy.dist1 <- NA
    val.dist1 <- NA
  }
  if(dd.reg1 >= last.dreg | is.na(dd.reg1)){
    yy.reg1 <- NA
    val.reg1 <- NA
  }
  
  if(dd.dist2 >= last.ddis | is.na(dd.dist2)){
    yy.dist2 <- NA
    val.dist2 <- NA
  }
  if(dd.reg2 >= last.dreg | is.na(dd.reg2)){
    yy.reg2 <- NA
    val.reg2 <- NA
  }
  
  if(dd.dist3 >= last.ddis | is.na(dd.dist3)){
    yy.dist3 <- NA
    val.dist3 <- NA
  }
  if(dd.reg3 >= last.dreg | is.na(dd.reg3)){
    yy.reg3 <- NA
    val.reg3 <- NA
  }
  
  if(dd.dist4 >= last.ddis | is.na(dd.dist4)){
    yy.dist4 <- NA
    val.dist4 <- NA
  }
  if(dd.reg4 >= last.dreg | is.na(dd.reg4)){
    yy.reg4 <- NA
    val.reg4 <- NA
  }
  # ---------------------------------------------------
  output <-as.numeric(c(yy.dist1,val.dist1,yy.reg1,val.reg1,
                        yy.dist2,val.dist2,yy.reg2,val.reg2,
                        yy.dist3,val.dist3,yy.reg3,val.reg3,
                        yy.dist4,val.dist4,yy.reg4,val.reg4))
  output[1:16]
}

#' 
#' @example 
#' dist.reg(x=ano_prob,dates=dates, p=99,rth = 1, dth = 730, cdates=3)
