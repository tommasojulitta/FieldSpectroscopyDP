\name{fGaussianConvolMultiWl}
\alias{fGaussianConvolMultiWl}
\title{fGaussianConvolMultiWl}
\description{Function to compute the gaussian convolution on multiple spectral range}
\usage{fGaussianConvolMultiWl(Wl = Wl, Ref = Ref, CWl.v = CWl.v, FWHM.v = FWHM.v, 
    convol = "Gaussian", weights = NULL)}
\arguments{
  \item{Wl}{numeric array: wavelenghts spectrometer to be deconvolved}
  \item{Ref}{numeric array: Reflectance}
  \item{CWl.v}{numeric array: array with the Central Wavelengths}
  \item{FWHM.v}{numeric array: array with the Full-Width at Half Maximum of the central Wavelength}
  \item{convol}{Character: Type of convolution. Default is "Gaussian", alternative "User"}
  \item{weights}{data frame: each column contains the weights of the user-defined convolution function}
}

\value{data frame with the central Wavelength and the result of the gaussian convolution (mean and standard deviation) }

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{

data("reflOO")
#Define the center wavlength of the region to be convolved
CWl.v<-c(400,500,600,700)
#Define the corresponding Full Widht at Half Maximum
FWHM.v<-c(8,8,8,8)
#Extract the convoluted values
res<-fGaussianConvolMultiWl(Wl=reflOO$Wl, Ref=reflOO$Ref, CWl.v=CWl.v, FWHM.v=FWHM.v)
#Plot results
plot(reflOO$Wl, reflOO$Ref, ylim=c(0,0.7),xlim=c(400,900), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
points(res$CWl.v, res$RefConvol, col="red")
legend("topleft",col=c("red"),pch=1,cex=1.2,legend=c("Convolved values"),box.col="white")
box()
}
