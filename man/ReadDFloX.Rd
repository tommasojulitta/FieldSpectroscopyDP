\name{ReadDFloX}
\alias{ReadDFloX}
\title{ReadDFloX}
\description{Read Data from csv files saved by D-FloX}
\usage{ReadDFloX(filename, sep = ";", na.strings = "#N/D", header = FALSE, 
    Ename = "QE_WR", dcEname = "QE_DC_WR", Lname = "QE_VEG", 
    dcLname = "QE_DC_VEG", E2name = "QE_WR2", ITEpos = 6, ITLpos = 8, 
    ITL2pos = 10, datepos = 2, timepos = 3, cycletimepos = 12, 
    cyclenrpos = 1, temp1pos = 14, temp2pos = 16, temp3pos = 18, 
    temp4pos = 20, floxfluoApos = 26, floxfluoBpos = 28, insideH = 22, 
    outsideH = 24, GPStimepos = 31, GPSdatepos = 33, GPSLATpos = 35, 
    GPSLONpos = 37)}
\arguments{
  \item{filename}{character value or vector: names of the file(s) to be opened}
  \item{sep}{the field separator character}
  \item{na.strings}{
}
  \item{header}{logical value indicating whether the file contains the names of the variables as its first line}
  \item{Ename}{character value: string of the name in the ASCII file of the solar radiance vector, if any}
  \item{dcEname}{character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of solar irradiance, if any}
  \item{Lname}{character value: string of the name in the ASCII file of the reflected radiance vector, if any}
  \item{dcLname}{character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of reflected radiance vector, if any}
  \item{E2name}{character value: string of the name in the ASCII file of the reflected radiance vector collected at double IT, if any}
  \item{ITEpos}{numeric value: value of the position (number of column) in the ASCII file of the incoming radiance}
  \item{ITLpos}{numeric value: value of the position (number of column) in the ASCII file of  the reflected radiance}
  \item{ITL2pos}{numeric value: value of the position (number of column) in the ASCII file of  the reflected radiance with double IT}
  \item{datepos}{numeric value: value of the position (number of column) in the ASCII file of the date}
  \item{timepos}{numeric value: value of the position (number of column) in the ASCII file of the time}
  \item{cycletimepos}{numeric value: value of the position (number of column) in the ASCII file of the cycle length}
  \item{cyclenrpos}{numeric value: value of the position (number of column) in the ASCII file  of the cycle number}
  \item{temp1pos}{numeric value: value of the position (number of column) in the ASCII file  of the first temperature sensor}
  \item{temp2pos}{numeric value: value of the position (number of column) in the ASCII file  of the second temperature sensor}
  \item{temp3pos}{numeric value: value of the position (number of column) in the ASCII file  of the third temperature sensor}
  \item{temp4pos}{numeric value: value of the position (number of column) in the ASCII file of the fourth temperature sensor}
  \item{floxfluoApos}{numeric value: value of the position (number of column) in the ASCII file of the FloX O2A fluo estimate}
  \item{floxfluoBpos}{numeric value: value of the position (number of column) in the ASCII file of the FloX O2B fluo estimate}
  \item{insideH}{numeric value: value of the position (number of column) in the ASCII file of the first humidity sensor}
  \item{outsideH}{numeric value: value of the position (number of column) in the ASCII file of the second humidity sensor}
  \item{GPStimepos}{numeric value: value of the position (number of column) in the ASCII file of the GPS time}
  \item{GPSdatepos}{numeric value: value of the position (number of column) in the ASCII file of the GPS date}
  \item{GPSLATpos}{numeric value: value of the position (number of column) in the ASCII file of the GPS LAT}
  \item{GPSLONpos}{numeric value: value of the position (number of column) in the ASCII file of the GPS LON}
}

\value{list containing the FloX data from QE spectrometer and ancillary data collected by the system}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

#Define path and filename
filename<-"..."
#load data
dat<-ReadDFloX(filename = filename)
 
}
