\name{fGaussianConvol}
\alias{fGaussianConvol}
\title{fGaussianConvol}
\description{Function to compute the gaussian convolution on a specific spectral range}
\usage{fGaussianConvol(Wl = Wl, Ref = Ref, CWl = CWl, FWHM = FWHM)}
\arguments{
  \item{Wl}{numeric array: wavelenghts spectrometer to be deconvolved}
  \item{Ref}{numeric array: Reflectance}
  \item{CWl}{scalar: Central Wavelength}
  \item{FWHM}{scalar: Full-Width at Half Maximum of the central Wavelength}
}

\value{numeric value of the gaussian convolution}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{

data("reflOO")
# Select starting wavelenght
CWl<-600
# Select last wavelenght
FWHM<-8
#apply the convolution on the spectra at the defined spectal range
res<-fGaussianConvol(Wl=reflOO$Wl, Ref=reflOO$Ref, CWl=CWl, FWHM=FWHM)
#Plot results
plot(reflOO$Wl, reflOO$Ref, ylim=c(0,0.7),xlim=c(400,900), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
points(CWl, res, col="red")
legend("topleft",col=c("red"),pch=1,cex=1.2,legend=c("convolved value"),box.col="white")
box()
}
