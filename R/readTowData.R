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

trim <- function (x) gsub("^\\s+|\\s+$", "", as.character(x))

parse.date <- function(str) {
    str <- trim(str)
    s <- strsplit(str,"/")[[1]]
    if (length(s) == 3) {
        if (nchar(s[3])==4) return(paste(s[2],s[1],s[3],sep=" "))
        return(paste(s[2],s[1],format(strptime(s[3],format="%y"),"%Y"),sep=" "))
    }
    s <- strsplit(str,"[.]")[[1]]
    if (length(s) == 3) {
        if (nchar(s[3])==4) return(paste(s[1],s[2],s[3],sep=" "))
        return(paste(s[1],s[2],format(strptime(s[3],format="%y"),"%Y"),sep=" "))
    }
    return(NA)
}

parse.time <- function(str) {
    str <- trim(str)
    s <- strsplit(str,":")[[1]]
    if (length(s) == 3) {
        return(paste(s,collapse=" "))
    }
    return(paste(c(s,0),collapse=" "))
}

read.timestamp <- function(dat,time.zone) {
    tm <- NULL
    cn <- colnames(dat)
    fmt <- "%d %m %Y %H %M %S"
    if ("Date_UTC"%in%cn && "Time_UTC"%in%cn) {
        dat.date <- sapply(dat[["Date_UTC"]],parse.date)
        dat.time <- sapply(dat[["Time_UTC"]],parse.time)
        tm <- as.POSIXlt(as.POSIXct(paste(dat.date,dat.time,sep=" "),tz="UTC",format=fmt),tz=time.zone,usetz=TRUE)
    } else if ("Timestamp_UTC"%in%cn) {
        tmp <- strsplit(trim(dat[["Timestamp_UTC"]])," +")
        tmp <- unlist(lapply(tmp,function(x) paste(parse.date(x[1]),parse.time(x[2]),sep=" ")))
        tm <- as.POSIXlt(as.POSIXct(tmp,tz="UTC",format=fmt),tz=time.zone,usetz=TRUE)
    } else if ("Timestamp_Local"%in%cn) {
        tmp <- strsplit(trim(dat[["Timestamp_Local"]])," +")
        tmp <- unlist(lapply(tmp,function(x) paste(parse.date(x[1]),parse.time(x[2]),sep=" ")))
        tm <- as.POSIXlt(tmp,tz=time.zone,format=fmt)
    } else {
        warning("I haven't found the date/time column I was looking for.\nPlease provide one of these:\n[-] Date_UTC / Time_UTC\n[-] Timestamp_UTC\n[-] Timestamp_Local")
    }
    return(tm)
}

parse.lldegree <- function(s) {
    return(as.numeric(sp::char2dms(sprintf(
        "%sd%sm%ss%s",
        trim(s[1]),
        trim(s[2]),
        trim(s[3]),
        trim(s[4])),chd="d",chm="m",chs="s")))
}

read.lonlat <- function(dat) {
    cn <- colnames(dat)
    if ("Lat_Degree"%in%cn && "Lat_Minute"%in%cn && "Lat_Second"%in%cn && "Lat_NS"%in%cn && "Lon_Degree"%in%cn && "Lon_Minute"%in%cn && "Lon_Second"%in%cn && "Lon_EW"%in%cn) {
        return(list(
            apply(cbind(dat[["Lon_Degree"]],dat[["Lon_Minute"]],dat[["Lon_Second"]],dat[["Lon_EW"]]),1,function(s) parse.lldegree(s)),
            apply(cbind(dat[["Lat_Degree"]],dat[["Lat_Minute"]],dat[["Lat_Second"]],dat[["Lat_EW"]]),1,function(s) parse.lldegree(s))
        ))
    } else if ("Latitude"%in%cn && "Longitude"%in%cn) {
        return(list(
            as.numeric(dat[["Longitude"]]),
            as.numeric(dat[["Latitude"]])
        ))
    }
    warning("I haven't found the lon/lat column I was looking for.\nPlease provide one of these:\n[-] Lat_Degree / Lat_Minutes / Lat_Seconds / Lat_NS / Lon_Degree / Lon_Minutes / Lon_Seconds / Lon_EW\n[-] Latitude / Longitude")
}

read.cprfile <- function(filename) {
    dat <- read.csv(filename,header=TRUE,comment.char="#",as.is=TRUE,blank.lines.skip=TRUE)
    dat <- dat[rowSums(dat=="",na.rm=TRUE)!=ncol(dat),]
    return(dat)
}

read.TowLog <- function(filename,time.zone) {
    dat<-read.cprfile(filename)
    ll<-read.lonlat(dat)
    cn<-colnames(dat)
    d<-data.frame(
        "time" = read.timestamp(dat,time.zone),
        "lon" = ll[[1]],
        "lat" = ll[[2]],
        "dist" = if ("Distance"%in%cn) {as.numeric(dat[["Distance"]])} else {NA}
    )
    return(d)
}

read.TowAIS <- function(filename,time.zone) {
    dat<-read.cprfile(filename)
    ll<-read.lonlat(dat)
    d<-data.frame(
        "time" = read.timestamp(dat,time.zone),
        "lon" = ll[[1]],
        "lat" = ll[[2]]
    )
    return(d)
}

read.CTD <- function(filename,time.zone) {
    dat<-read.cprfile(filename)
    cn <- colnames(dat)
    d<-data.frame(
        "time" = read.timestamp(dat,time.zone),
        "temp" = if ("Temperature_C"%in%cn) {as.numeric(gsub(',','.',dat[["Temperature_C"]]))} else {NA},
        "pres" = if ("Pressure_Bar"%in%cn) {as.numeric(gsub(',','.',dat[["Pressure_Bar"]]))} else {NA},
        "salin" = if ("Salinity_psu"%in%cn) {as.numeric(gsub(',','.',dat[["Salinity_psu"]]))} else {NA},
        "depth" = if ("Depth_m"%in%cn) {as.numeric(gsub(',','.',dat[["Depth_m"]]))} else {NA}
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
                                            time.zone="character",
                                            data.log="data.frame",
                                            data.ais="data.frame",
                                            data.ctd="data.frame",
                                            data.pci="data.frame"))

setMethod("initialize",
          "towClass",
          function(.Object, tow.log=NULL, tow.ais=NULL, tow.ctd=NULL, tow.pci=NULL, tow.id=0, silk.start=0, silk.end=0, time.zone="Etc/GMT-2") {
              .Object@tow.id <- tow.id
              .Object@time.zone <- time.zone
              if (!(.Object@time.zone%in%OlsonNames())) {
                  .Object@time.zone <- "Etc/GMT-2"
                  warning("This is not a valid time zone: ",time.zone,"\nPlease see the documentation for Sys.timezone for more information.\nUsing ",.Object@time.zone," as default.")
              }
              .Object@silk.start <- silk.start
              .Object@silk.end <- silk.end
              if (length(tow.log)) {
                  .Object@tow.log <- tow.log
                  .Object@data.log <- read.TowLog(.Object@tow.log,.Object@time.zone)
              }
              if (length(tow.ais)) {
                  .Object@tow.ais <- tow.ais
                  .Object@data.ais <- read.TowAIS(.Object@tow.ais,.Object@time.zone)
              }
              if (length(tow.ctd)) {
                  .Object@tow.ctd <- tow.ctd
                  .Object@data.ctd <- read.CTD(.Object@tow.ctd,.Object@time.zone)
              }
              return(.Object)
          })

setGeneric(name="lonlat", def=function(object,...){standardGeneric("lonlat")})
setMethod("lonlat",
          "towClass",
          function(object,x=NULL) {
              if (length(object@tow.log) && !length(object@tow.ais)) {
                  return(l.partial(object@data.log,x))
              } else if (length(object@tow.ais)) {
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
              cat("Time zone........:", object@time.zone, "\n")
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
