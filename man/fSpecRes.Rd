\name{fSpecRes}
\alias{fSpecRes}
\title{fSpecRes}
\description{Function to resample based on gaussian convolution according to the spectral response of different satellite missions}
\usage{fSpecRes(Wl = Wl, Ref = Ref, convol = "Gaussian", Satellite = "MODIS")}
\arguments{
  \item{Wl}{
}
  \item{Ref}{numeric array: Reflectance}
  \item{convol}{Character: Type of convolution. Default is "Gaussian", alternative "Uniform": all the bands in the spectral range are equally weighted,  User"}
  \item{Satellite}{Character: Satellite e.g. "MODIS", "LANDSAT", Default MODIS}
}

\value{data frame with the central Wavelength and the result of the gaussian convolution (mean and standard deviation)}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{

data("reflOO")

#Extract the convolved values on reflectance according to MODIS bands
res<-fSpecRes(Wl=reflOO$Wl, Ref=reflOO$Ref, convol='Gaussian', Satellite='MODIS')
plot(reflOO$Wl, reflOO$Ref, ylim=c(0,0.7),xlim=c(400,900), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
points(res$CWl.v, res$RefConvol, col="red")

#Extract the convolved values on reflectance according to LANDSAT bands
res<-fSpecRes(Wl=reflOO$Wl, Ref=reflOO$Ref, convol='Gaussian', Satellite='LANDSAT')
points(res$CWl.v, res$RefConvol, col="green", pch=19)
legend("topleft",col=c("red","green"),pch=c(1,19),cex=1.2,legend=c("MODIS","LANDSAT"),box.col="white")
box()
}
