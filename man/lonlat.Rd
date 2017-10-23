\name{lonlat}
\alias{lonlat}
\title{ Continuous Plankton Recorder (CPR) Data Analyser }
\description{
This package is designed for the purpose of data analysis and visualisation for the Mediterranean Continuous Plankton Recorder Survey (MedCPR) project.
}

\usage{
  lonlat(tow,x=NULL)
}
\arguments{
  \item{tow }{The towClass to use for interpolation }
  \item{x }{ISOdate times to calculate location. If not provided, data.ais or data.log is used}
}
    
\value{
  \item{time }{Timestamp of the location }
  \item{lon }{Longitude }
  \item{lat }{Latidute }
}

\examples{
\dontrun{
tow <- towClass(tow.ais="2MCHL/2_TowAIS_R.csv",
                tow.ctd="2MCHL/2_CTD_R.csv")

lonlat(tow)
lonlat(tow,tow@data.ctd$time)
lonlat(tow,ISOdate(2014,11,05,23,20,30,"Etc/GMT-2"))
}
}
