\name{SigTimeInterpolation}
\alias{SigTimeInterpolation}
\title{SigTimeInterpolation}
\description{Compute the linear interpolation of each wavelength of two measuremnts at a specific time occurred in between the two}
\usage{SigTimeInterpolation(FirstTime, SecndTime, TimeBetween, FirstMeas, 
    SecndMeas)}
\arguments{
  \item{FirstTime}{numeric value: time of acquiistion of the first vector}
  \item{SecndTime}{numeric value: time of acquiistion of the second vector}
  \item{TimeBetween}{numeric value: time at which the interpolation want to be done}
  \item{FirstMeas}{numeric vector: spectrum acquire at FirstTime}
  \item{SecndMeas}{numeric vector: spectrum acquire at SecndTime}
}

\value{numeric vector containing spectrum linearly interpolated.}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

data("FloX_data")
data("wl_FloX")

#Define first spectrum
spectrum1<-FloX_data$E[,2]
#Define second spectrum
spectrum2<-FloX_data$E[,3]
#Compute interpolated spectrum
spectrum_int<-SigTimeInterpolation(FirstTime = 10,SecndTime = 25,TimeBetween = 17,FirstMeas = spectrum1,SecndMeas = spectrum2)
#plot
x11()
par(mar=c(5,5,2,2))
plot(wl_FloX,spectrum1,type="l",xlab="WL [nm]",ylab="Digital Counts [-]",ylim = c(0,200000))
lines(wl_FloX,spectrum2,col="green")
lines(wl_FloX,spectrum_int,col="red")
legend("topleft",col=c("black","green","red"),lty=1,cex=1.2,legend=c("spectrum1 ","spectrum2","spectrum_int"),box.col="white")
box()
}
