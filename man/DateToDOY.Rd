\name{DateToDOY}
\alias{DateToDOY}
\title{DateToDOY}
\description{###############################################################################################################################################################
Convert date to julian day of the year and day fraction}
\usage{DateToDOY(datetime)}
\arguments{
  \item{datetime}{numeric vector of POSIXTct class: vector of the date and time}
}

\value{numeric vector containing the coverted date in Day Of the Year (DOY).}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

data("FloX_data")
data("wl_FloX")

#define a vector of POSIXTct class
datetime<-DateTimeFloX(FloX_data$time,FloX_data$date)
#Convert to Julian Day Of the Year and day fraction
doy.dayfract<-DateToDOY(datetime)
}
