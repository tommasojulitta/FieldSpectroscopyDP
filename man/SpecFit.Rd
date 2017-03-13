\name{SpecFit}
\alias{SpecFit}
\title{SpecFit}
\description{Compute the spectral sum of reflected radiance and fluorescence. Inverse mode returns the residual sum of squares. Forward mode returns the modeled radiance obtained by the sum of fluorescence and radiance.}
\usage{SpecFit(x, E, L, fm, wl, run, O2band)}
\arguments{
  \item{x}{numeric vector: parameter for gaussian and spline function}
  \item{E}{numeric vector: measured solar radiance}
  \item{L}{numeric vector: measued reflected radiance}
  \item{fm}{object of class "lm", deriving from first guess functoin}
  \item{wl}{numeric vector: wavelength vector}
  \item{run}{character value. If "inverse" residual sum f squares is returned. If "forward" modelled reflected radiance is returned.}
  \item{O2band}{character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation}
}

\value{if "Inverse" mode the residual sum of squares is returned. If "Forward" mode a numeric vector containing the modeled radiance obtained is returned.}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

data("FloX_data")
data("up_coeff")
data("dw_coeff")
data("wl_FloX")

L<-GetRadiance(DNSignal=FloX_data$L-FloX_data$dcL,IntegrationTime=FloX_data$IT_L/1000,RadCalCoeff=dw_coeff);L<-L[,1]
#Get Solar Radiance 
E<-GetRadiance(DNSignal=FloX_data$E-FloX_data$dcE,IntegrationTime=FloX_data$IT_E/1000,RadCalCoeff=up_coeff);E<-E[,1]
#Estimate fluorescence using iFLD method, used as first guess for Spectral Fitting Methods
iFLD_O2B<-iFLD(wl=wl_FloX,E,L,fwhm =0.4,O2band="B")
#Define range used for =2B band
range<-which(wl_FloX>684& wl_FloX<700)
#Subset the wavelength vector
WL<-wl_FloX[range]
E_sfm<-as.numeric(E[range])
L_sfm<-as.numeric(L[range])
fluoFG<-iFLD_O2B$Fluo
#Compute the first guess parameter both for reflectance and fluorescence
fg<-FirstGuess(wl = WL,L = L_sfm,E = E_sfm, fluo = iFLD_O2B$Fluo,O2band= "B")
#Model the relfected radiance using the first guess as input parameters
modelled_radiance<-SpecFit(fg$first_guess$FG,wl=WL,E=E_sfm,L=L_sfm,fm=fg$fm,run="forward",O2band ="B")
#plot results
x11()
plot(WL,L_sfm,type="l",ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"))
lines(WL,modelled_radiance,col="red",xlab="WL [nm]")
legend("topleft",col=c("black","red"),lty=1,legend=c("Measured reflected radiance","Modeled reflected radiance"),box.col="white",lwd=2,cex=1.5)
box()

  
}
