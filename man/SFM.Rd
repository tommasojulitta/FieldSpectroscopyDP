\name{SFM}
\alias{SFM}
\title{SFM}
\description{Compute spectral fitting methods on measured data. It returns fluorescence and true reflectance}
\usage{SFM(wl, L, E, fluoFG, O2band, output)}
\arguments{
  \item{wl}{numeric vector: wavelength vector}
  \item{L}{numeric vector: measued reflected radiance}
  \item{E}{numeric vector: measured solar radiance}
  \item{fluoFG}{numeric value: fluorescence estimate derived from iFLD method}
  \item{O2band}{character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation}
  \item{output}{character value: FULL or VALUE referring to output expected. If FULL a data.fame of the spectrum in the considered range of fluorescence and true reflectance is returned. If VALUE the fluorescence and the true reflectance at the selected oxygen band is returned.}
}

\value{List containing: }

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

data("FloX_data")
data("up_coeff")
data("dw_coeff")
data("wl_FloX")

#Get Target Radiance 
L<-GetRadiance(DNSignal=FloX_data$L-FloX_data$dcL,IntegrationTime=FloX_data$IT_L/1000,RadCalCoeff=dw_coeff)
#Get Solar Radiance 
E<-GetRadiance(DNSignal=FloX_data$E-FloX_data$dcE,IntegrationTime=FloX_data$IT_E/1000,RadCalCoeff=up_coeff)
#Estimate fluorescence using iFLD method, used as first guess for Spectral Fitting Methods
iFLD_O2A<-iFLD(wl=wl_FloX,E,L,fwhm =0.3,O2band="A")
iFLD_O2B<-iFLD(wl=wl_FloX,E,L,fwhm =0.4,O2band="B")
#Estimate fluorescence at O2A band using SFM method, used as first guess for Spectral Fitting Methods
res<-SFM(wl = wl_FloX, L = L[,1],E= E[,1],fluoFG = iFLD_O2A$Fluo[1], O2band ="A", output = "FULL")
sfm_FLD_O2A<-res$Results

#Plot results
x11()
range<-which(wl_FloX>=750& wl_FloX<=780)
WL<-wl_FloX[range]
par(mfcol=c(1,2))
plot(wl_FloX,(L[,1]/E[,1]),pch=20,ylim=c(0.5,1.2),xlim=c(755,770),type="l",lwd=3,col="green",xlab="Wavelength [nm]",ylab="Reflectance [-]",cex.lab=1.5,cex.axis=1.5,main = "Reflectance")
lines(WL,sfm_FLD_O2A[,2],pch=20,col="red",lwd=2)
grid()
legend("topleft",col=c("green","red"),lty=1,legend=c("Apparent","Retrieved"),box.col="white",lwd=2,cex=1.5)
box()
plot(WL,sfm_FLD_O2A[,1],ylim=c(0,5),type="l",lwd=3,xlim=c(755,770),col="red",xlab="Wavelength [nm]",ylab="Fluorescence [mW m-2 sr-1 nm-1]",cex.lab=1.5,cex.axis=1.5,main = "Fluorescence")
abline(h=iFLD_O2A$Fluo[1]*1000,lty=2,col="green",lwd=3)
grid()
legend("topleft",col=c("green","red"),lty=c(3,1),legend=c("iFLD","SFM"),box.col="white",lwd=2,cex=1.5)
box()

#Estimate fluorescence at O2B band using SFM method, used as first guess for Spectral Fitting Methods
res<-SFM(wl = wl_FloX, L = L[,1],E= E[,1],fluoFG = iFLD_O2A$Fluo[1], O2band ="B", output = "FULL")
sfm_FLD_O2B<-res$Results
#Plot results
x11()
range<-which(wl_FloX>=684& wl_FloX<=700)
WL<-wl_FloX[range]
par(mfcol=c(1,2))
plot(wl_FloX,(L[,1]/E[,1]),pch=20,ylim=c(0,0.2),xlim=c(686,698),type="l",lwd=3,col="green",xlab="WL [nm]",ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.lab=1.5,cex.axis=1.5,main = "Reflectance")
lines(WL,sfm_FLD_O2B[,2],pch=20,col="red",lwd=2)
grid()
legend("topleft",col=c("green","red"),lty=1,legend=c("Apparent","Retrieved"),box.col="white",lwd=2,cex=1.5)
box()
plot(WL,sfm_FLD_O2B[,1],ylim=c(0,5),type="l",lwd=3,xlim=c(686,698),col="red",xlab="WL [nm]",ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.lab=1.5,cex.axis=1.5,main = "Fluorescence")
abline(h=iFLD_O2B$Fluo[1]*1000,lty=2,col="green",lwd=3)
grid()
legend("topleft",col=c("green","red"),lty=c(3,1),legend=c("iFLD","SFM"),box.col="white",lwd=2,cex=1.5)
box()
}
