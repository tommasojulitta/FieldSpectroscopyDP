\name{ReadDFloXHR}
\alias{ReadDFloXHR}
\title{ReadDFloXHR}
\description{Read Data from csv files saved by D-FloX}
\usage{ReadDFloXHR(filename, sep = ";", na.strings = "#N/D", header = FALSE, 
    Ename = "HR_WR", dcEname = "HR_DC_WR", Lname = "HR_VEG", 
    dcLname = "HR_DC_VEG", ITEpos = 6, ITLpos = 8, datepos = 2, 
    timepos = 3, cyclenrpos = 1, temp1pos = 12, temp2pos = 14, 
    temp3pos = 16, insideH = 16, outsideH = 18)}
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
  \item{ITEpos}{numeric value: value of the position (number of column) in the ASCII file of integration time of the incoming radiance}
  \item{ITLpos}{numeric value: value of the position (number of column) in the ASCII file of integration time of the reflected radiance}
  \item{datepos}{numeric value: value of the position (number of column) in the ASCII file of integration time of the date}
  \item{timepos}{numeric value: value of the position (number of column) in the ASCII file of integration time of the time}
  \item{cyclenrpos}{numeric value: value of the position (number of column) in the ASCII file of integration time of the cycle number}
  \item{temp1pos}{numeric value: value of the position (number of column) in the ASCII file of integration time of the first temperature sensor}
  \item{temp2pos}{numeric value: value of the position (number of column) in the ASCII file of integration time of the second temperature sensor}
  \item{temp3pos}{numeric value: value of the position (number of column) in the ASCII file of integration time of the third temperature sensor}
  \item{insideH}{numeric value: value of the position (number of column) in the ASCII file of integration time of the first humidity sensor}
  \item{outsideH}{numeric value: value of the position (number of column) in the ASCII file of integration time of the second humidity sensor}
}

\value{list containing the FloX data from HR spectrometer and ancillary data collected by the system}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

#Define path and filename
filename<-"..."
#load data
dat<-ReadDFloXHR(filename = filename)
  
}
