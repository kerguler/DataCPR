\name{DataCPR}
\alias{DataCPR}
\docType{package}

\title{
Continuous Plankton Recorder (CPR) Data Analyser
}

\description{
This package is designed for the purpose of data analysis and visualisation for the Mediterranean Continuous Plankton Recorder Survey (MedCPR) project.
}

\author{
\enumerate{
\item Kamil Erguler, Associate Research Scientist, EEWRC, The Cyprus Institute \email{k.erguler@cyi.ac.cy}
}
}
\note{ The author is indebted to Rana Abu Alhaija for her invaluable contributions to the development and improvement of this package. }

\usage{
  towClass(...)
}
\arguments{
  \item{tow.id}{The serial number of the specific tow}
  \item{tow.ais}{Data files for the tow}
  \item{tow.log}{Data files for the tow}
  \item{tow.ctd}{Data files for the tow}
  \item{tow.pci}{Data files for the tow}
  \item{silk.start}{Start of silk partitioning}
  \item{silk.end}{End of silk partitioning}
  \item{time.offset}{Time elapsed (in sec) until first data recording}
}
    
\value{
  \item{towClass }{All data and methods are wrapped within this S4 class}
}
\details{
    Cool stuff to wait for...
}
\examples{
\dontrun{
tow <- towClass(tow.log="data/9_TowLog_R.csv",
                tow.ais="data/9_TowAIS_R.csv",
                tow.ctd="data/9_CTD_R.csv",
                tow.pci="data/9_PCI_R.csv",
                tow.id=9,
                silk.start=0,
                silk.end=26.9,
                time.offset=-2)

## Then, we can view the summary of the class contents with the followings:

print(tow)

print(tow@data.log)
print(tow@data.ais)
print(tow@data.ctd)
print(tow@data.pci)

plot(tow)
}
}

\keyword{cpr}
\keyword{plankton}
\keyword{data}


