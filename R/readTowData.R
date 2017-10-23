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

plot.coord <- function(x,  y, dates=NULL, add=FALSE, col="black", ...) {
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
}

rescale<-function(a,b=c(0,1),r=NULL) {
    if (is.null(r)) r<-range(a,na.rm=TRUE)
    return(b[1]+(b[2]-b[1])*(a-r[1])/(r[2]-r[1]))
}
rotmat<-function(theta) matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow=2,ncol=2)
data.rot<-function(x,y,r,dir=1)
{
    rot90<-rotmat(-0.5*pi)
    rotm90<-rotmat(0.5*pi)
    ret<-matrix(NA,ncol=2,nrow=length(x))
    mat<-rbind(cbind(x,y),c(x[length(x)-1],y[length(y)-1]))
    for (i in 1:length(x)) {
        m<-t(mat[i+1,])-mat[i,]
        m<-m/sqrt(sum(m^2))
        rot<-if (dir>0) {rot90} else {rotm90}
        if (i==length(x)) {
            rot<-if (dir>0) {rotm90} else {rot90}
        }
        ret[i,]<-(m%*%rot)*r[i] + mat[i,]
    }
    return(data.frame(x=ret[,1],y=ret[,2]))
}
poly2d<-function(ll1,ll2) {
    return(data.frame(x=c(ll1$x,rev(ll2$x)),
                      y=c(ll1$y,rev(ll2$y))))
}
dash2d<-function(ll1,ll2) {
    return(data.frame(x=unlist(lapply(1:nrow(ll1),function(i) c(ll1$x[i],ll2$x[i],NA))),
                      y=unlist(lapply(1:nrow(ll1),function(i) c(ll1$y[i],ll2$y[i],NA)))))
}
plot.ctd <- function(x, xlim=NA, ylim=NA, ...) {
    d<-x@data.ctd
    d<-d[(d$time>=x@start.date) & (d$time<=x@stop.date),]
    ll<-lonlat(x,d$time)
    xy<-data.frame(x=ll$lon,y=ll$lat)
    # ---
    if (any(is.na(xlim))) xlim<-range(ll$lon)+c(-1.0,1.0)
    if (any(is.na(ylim))) ylim<-range(ll$lat)+c(-0.5,0.5)
    m<-map("world",xlim=xlim,ylim=ylim,resolution=0)
    # ---
    lgd<-list(c(
        "black",
        "Time",
        c(strftime(tow@start.date,"%d.%m.%y\n%H:%M:%S"),strftime(tow@stop.date,"%d.%m.%y\n%H:%M:%S")),
        "rev"
    ))
    if (!all(is.na(d$temp))) {
        ll2<-data.rot(xy$x,xy$y,rescale(d$temp,c(0.01,0.25)),dir=1)
        polygon(poly2d(xy,ll2),border="red")
        lgd<-append(lgd,list(c(
            "red",
            expression("Temperature ("~degree~"C)"),
            range(d$temp,na.rm=TRUE),
            "fwd"
            )))
    }
    if (!all(is.na(d$pres))) {
        ll2<-data.rot(xy$x,xy$y,rescale(d$pres,c(0.01,0.25)),dir=1)
        polygon(poly2d(xy,ll2),border="violet")
        lgd<-append(lgd,list(c(
            "violet",
            "Pressure (Bar)",
            range(d$pres,na.rm=TRUE),
            "fwd"
            )))
    }
    if (!all(is.na(d$salin))) {
        ll2<-data.rot(xy$x,xy$y,rescale(d$salin,c(0.01,0.25)),dir=-1)
        polygon(poly2d(xy,ll2),border="green")
        lgd<-append(lgd,list(c(
            "green",
            "Salinity (psu)",
            range(d$salin,na.rm=TRUE),
            "fwd"
            )))
    }
    if (!all(is.na(d$depth))) {
        ll2<-data.rot(xy$x,xy$y,rescale(d$depth,c(0.01,0.25)),dir=-1)
        polygon(poly2d(xy,ll2),border="blue")
        lgd<-append(lgd,list(c(
            "blue",
            "Depth (m)",
            range(d$depth,na.rm=TRUE),
            "fwd"
            )))
    }
    # ---
    lines(ll$lon,ll$lat,
          lw=2,
          col="black")
    ord<-c(which(d$time==min(d$time)),
           which(d$time==max(d$time)))
    l <- if (ord[2]>=ord[1]) {ll} else {ll[nrow(ll):1,]}
    arrows(l$lon[2],l$lat[2],l$lon[1],l$lat[1],
           lw=2,
           col="black",
           length=0.1,
           angle=120)
    # ---
    # x1, y1, x2, y2
    bbox<-c(xlim[1],0.15*(ylim[2]-ylim[1])+ylim[1],0.25*(xlim[2]-xlim[1])+xlim[1],ylim[2])
    ow<-0.2*(bbox[3]-bbox[1])
    oh<-0.075*(bbox[4]-bbox[2])
    bbox[1]<-bbox[1]+ow
    bbox[3]<-bbox[3]+ow
    for (n in length(lgd):1) {
        if (eval(lgd[[n]][5]) == "rev") {
            arrows(
                bbox[3],
                bbox[2]-(n-1)*oh,
                bbox[1],
                bbox[2]-(n-1)*oh,
                lw=2,
                length=0.05,
                angle=120,
                col=eval(lgd[[n]][1])
            )
        } else {
            arrows(
                bbox[1],
                bbox[2]-(n-1)*oh,
                bbox[3],
                bbox[2]-(n-1)*oh,
                lw=2,
                length=0.05,
                angle=30,
                col=eval(lgd[[n]][1])
            )
        }
        text(bbox[1],bbox[2]-(n-1)*oh,lgd[[n]][3],pos=3,cex=0.5)
        text(bbox[3],bbox[2]-(n-1)*oh,lgd[[n]][4],pos=3,cex=0.5)
        text(bbox[3],bbox[2]-(n-1)*oh,lgd[[n]][2],pos=4,cex=0.5)
    }
}

# https://www.programiz.com/r-programming/S4-class

towClass <- setClass("towClass", slots=list(tow.log="character",
                                            tow.ais="character",
                                            tow.ctd="character",
                                            tow.pci="character",
                                            tow.id="numeric",
                                            silk.start="numeric",
                                            silk.stop="numeric",
                                            start.date="POSIXct",
                                            stop.date="POSIXct",
                                            time.zone="character",
                                            data.log="data.frame",
                                            data.ais="data.frame",
                                            data.ctd="data.frame",
                                            data.pci="data.frame"))

setMethod("initialize",
          "towClass",
          function(.Object, tow.log=NULL, tow.ais=NULL, tow.ctd=NULL, tow.pci=NULL, tow.id=0, silk.start=0, silk.stop=0, time.zone="Etc/GMT-2") {
              .Object@tow.id <- tow.id
              .Object@time.zone <- time.zone
              if (!(.Object@time.zone%in%OlsonNames())) {
                  .Object@time.zone <- "Etc/GMT-2"
                  warning("This is not a valid time zone: ",time.zone,"\nPlease see the documentation for Sys.timezone for more information.\nUsing ",.Object@time.zone," as default.")
              }
              .Object@silk.start <- silk.start
              .Object@silk.stop <- silk.stop
              if (length(tow.log)) {
                  .Object@tow.log <- tow.log
                  .Object@data.log <- read.TowLog(.Object@tow.log,.Object@time.zone)
                  .Object@start.date <- min(.Object@data.log$time,na.rm=TRUE)
                  .Object@stop.date <- max(.Object@data.log$time,na.rm=TRUE)
              }
              if (length(tow.ais)) {
                  .Object@tow.ais <- tow.ais
                  .Object@data.ais <- read.TowAIS(.Object@tow.ais,.Object@time.zone)
                  .Object@start.date <- min(.Object@data.ais$time,na.rm=TRUE)
                  .Object@stop.date <- max(.Object@data.ais$time,na.rm=TRUE)
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
              cat("Tow ID............:", object@tow.id, "\n")
              cat("Tow Log...........:", object@tow.log, "\n")
              cat("Tow AIS...........:", object@tow.ais, "\n")
              cat("Tow CTD...........:", object@tow.ctd, "\n")
              cat("Tow PCI...........:", object@tow.pci, "\n")
              cat("Silk (start, stop):", sprintf("%g, %g",object@silk.start,object@silk.stop), "\n")
              cat("Time (start, stop):", sprintf("%g, %g",object@start.date,object@stop.date), "\n")
              cat("Time zone........:", object@time.zone, "\n")
          })

setMethod("plot", signature(x="towClass", y="missing"),
          function(x,  y, dates=NULL, add=FALSE, col="black", ...) {
              plot.coord(x, y, dates=dates, add=add, col=col, ...)
              readline(prompt = "Press <Enter> to continue...")
              plot.ctd(x, ...)
          })

