\name{DateTimeFloX}
\alias{DateTimeFloX}
\title{DateTimeFloX}
\description{Convert FloX date format in as.POSIXct class}
\usage{DateTimeFloX(time, date)}
\arguments{
  \item{time}{numeric vector of time saved from FloX system}
  \item{date}{numeric vector of date saved from FloX system}
}

\value{numeric vector containing the coverted date time in as.POSIXct class.}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

data("FloX_data")
data("wl_FloX")
#define vector of time
time<-FloX_data$time
#define vector of date
date<-FloX_data$date
#convert to as.POSIXTct class
datetime<-DateTimeFloX(time,date)
}
