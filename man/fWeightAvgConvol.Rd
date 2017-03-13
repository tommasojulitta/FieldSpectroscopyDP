\name{fWeightAvgConvol}
\alias{fWeightAvgConvol}
\title{fWeightAvgConvol}
\description{Function to compute the convolution on a specific spectral range given the wavelengths and weights}
\usage{fWeightAvgConvol(Wl = Wl, Ref = Ref, W = W)}
\arguments{
  \item{Wl}{numeric array: wavelenghts spectrometer to be deconvolved}
  \item{Ref}{numeric array: Reflectance}
  \item{W}{numeric array: contains the weights of the user-defined convolution function (same length of Ref)}
}

\value{data frame with the results of the convolution }

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{

data("reflOO")
# define the center wavelength
CWl<-530
# define the Full Widht at Half Maximum
FWHM<-50
# define the Weight for the convolution
zzz <- dnorm(reflOO$Wl, mean = CWl, sd = FWHM/(2*sqrt(2*log(2))))
x.norm = (zzz - min(zzz,na.rm=TRUE))/(max(zzz,na.rm=TRUE) - min(zzz,na.rm=TRUE))
# Exctract the convolved spectrum
res<-fWeightAvgConvol(Wl=reflOO$Wl, Ref=reflOO$Ref, W=x.norm)
#plot results
plot(reflOO$Wl, reflOO$Ref, ylim=c(0,2), xlim=c(400,800), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
lines(reflOO$Wl, x.norm, col="blue")
points(CWl, res$RefConvol, col="red")
legend("topleft",col=c("black","red","blue"),pch=c(NA,1,NA),lty=c(1,NA,1),cex=1.2,legend=c("Reflectance","Convolved value","Convolution function"),box.col="white")
box()
}
