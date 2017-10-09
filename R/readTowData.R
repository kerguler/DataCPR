if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package 'sp' needed for this function to work. Please install it.",
      call. = FALSE)
}

l.spline <- function(dat,x=NULL) {
    if (is.null(x)) {x<-dat$time}
    f.lon<-splinefun(dat$time,dat$lon)
    f.lat<-splinefun(dat$time,dat$lat)
    return(data.frame(
        "time"=x,
        "lon"=f.lon(x),
        "lat"=f.lat(x)
    ))
}
l.partial<-function(dat,x=NULL) {
    if (is.null(x)) {x<-dat$time}
    # A linear model for each interval separately for longitude and latitude
    l.lons<-NULL
    l.lats<-NULL
    for (i in 1:(nrow(dat)-1)) {
        l.lons[[i]]<-lm(lon~time,data=dat[i:(i+1),])
        l.lats[[i]]<-lm(lat~time,data=dat[i:(i+1),])
    }
    # Find which date corresponds to which data interval
    ids<-findInterval(x,dat$time,all.inside=TRUE)
    # Prepare the arrays for output
    lons<-rep(0,length(x))
    lats<-rep(0,length(x))
    # For each linear model, predict new lon,lat with the set of data points corresponding to the data interval
    for (i in 1:(nrow(dat)-1)) {
        x2<-data.frame("time"=x[ids==i])
        lons[ids==i]<-predict(l.lons[[i]],x2)
        lats[ids==i]<-predict(l.lats[[i]],x2)
    }
    return(data.frame(
        "time"=x,
        "lon"=lons,
        "lat"=lats
    ))
}

read.cprfile <- function(filename) {
    dat <- read.csv(filename,header=TRUE,comment.char="#",as.is=TRUE,blank.lines.skip=TRUE)
    dat <- dat[rowSums(dat=="")!=ncol(dat),]
}

read.TowLog <- function(filename) {
    dat<-read.cprfile(filename)
    d<-data.frame(
        "time" = strptime(paste(as.character(dat[["Date.UTC"]]),as.character(dat[["Time.UTC"]])),"%m/%d/%Y %H:%M",tz="UTC"),
        "lon" = as.numeric(sp::char2dms(as.character(dat[["Lon"]]),chd="o")),
        "lat" = as.numeric(sp::char2dms(as.character(dat[["Lat"]]),chd="o")),
        "dist" = as.numeric(dat[["Distance"]])
    )
    return(d)
}

read.TowAIS <- function(filename) {
    dat<-read.cprfile(filename)
    d<-data.frame(
        "time" = strptime(as.character(dat[["Timestamp..UTC."]]),"%m/%d/%Y %H:%M",tz="UTC"),
        "lon" = as.numeric(dat[["Longitude"]]),
        "lat" = as.numeric(dat[["Latitude"]])
    )
    return(d)
}

read.CTD <- function(filename,time.offset) {
    dat<-read.cprfile(filename)
    d<-data.frame(
        "time" = strptime(as.character(dat$V2),"%d.%m.%y %H:%M:%S",tz="UTC")+3600*time.offset,
        "temp" = as.numeric(gsub(',','.',dat$V3)),
        "pres" = as.numeric(gsub(',','.',dat$V4)),
        "salin" = as.numeric(gsub(',','.',dat$V5))
    )
    return(d)
}

# https://www.programiz.com/r-programming/S4-class

towClass <- setClass("towClass", slots=list(tow.log="character",
                                            tow.ais="character",
                                            tow.ctd="character",
                                            tow.pci="character",
                                            tow.id="numeric",
                                            silk.start="numeric",
                                            silk.end="numeric",
                                            time.offset="numeric",
                                            data.log="data.frame",
                                            data.ais="data.frame",
                                            data.ctd="data.frame",
                                            data.pci="data.frame"))

setMethod("initialize",
          "towClass",
          function(.Object, tow.log=NULL, tow.ais=NULL, tow.ctd=NULL, tow.pci=NULL, tow.id=0, silk.start=0, silk.end=0, time.offset=0) {
              .Object@tow.id <- tow.id
              .Object@time.offset <- time.offset
              .Object@silk.start <- silk.start
              .Object@silk.end <- silk.end
              if (length(tow.log)) {
                  .Object@tow.log <- tow.log
                  .Object@data.log <- read.TowLog(.Object@tow.log)
              }
              if (length(tow.ais)) {
                  .Object@tow.ais <- tow.ais
                  .Object@data.ais <- read.TowAIS(.Object@tow.ais)
              }
              if (length(tow.ctd)) {
                  .Object@tow.ctd <- tow.ctd
                  .Object@data.ctd <- read.CTD(.Object@tow.ctd,.Object@time.offset)
              }
              return(.Object)
          })

setGeneric(name="lonlat", def=function(object,...){standardGeneric("lonlat")})
setMethod("lonlat",
          "towClass",
          function(object,x=NULL) {
              if (length(object@tow.log) && !length(object@tow.ais)) {
                  return(l.partial(object@data.log,x))
              } else if (length(object@tow.ais) && !length(object@tow.log)) {
                  return(l.spline(object@data.ais,x))
              }
              return(NULL)
          })

setMethod("show",
          "towClass",
          function(object) {
              cat("Tow ID...........:", object@tow.id, "\n")
              cat("Tow Log..........:", object@tow.log, "\n")
              cat("Tow AIS..........:", object@tow.ais, "\n")
              cat("Tow CTD..........:", object@tow.ctd, "\n")
              cat("Tow PCI..........:", object@tow.pci, "\n")
              cat("Silk (start, end):", sprintf("%g, %g",object@silk.start,object@silk.end), "\n")
              cat("Time offset......:", object@time.offset, "\n")
          })

setMethod("plot", signature(x="towClass", y="missing"),
          function(x,  y, dates=NULL, add=FALSE, col="black", ...) {
              y <- if (is.null(dates)) {
                       if (length(x@tow.ais)) {
                           x@data.ais$time
                       } else {
                           seq(x@data.log$time[1],x@data.log$time[length(x@data.log$time)],"mins")
                       }
                   } else {
                       dates
                   }
              l.dl<-lonlat(x,y)
              fun <- if (!add) {plot} else {lines}
              fun(l.dl$lon,l.dl$lat,type="l",lwd=1.5,col=col,xlab="Longitude",ylab="Latitude",...)
              if (length(x@tow.ais))
                  points(x@data.ais$lon,x@data.ais$lat,pch=1,cex=1.5,col=rgb(0,0,0,0.5))
              points(x@data.log$lon,x@data.log$lat,pch=16,cex=1,col="red")
          })
