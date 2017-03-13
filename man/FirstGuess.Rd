\name{FirstGuess}
\alias{FirstGuess}
\title{FirstGuess}
\description{Calculate the first guess of fluorescence and true reflectance based on gaussian function and spline fitting on apparent reflectance}
\usage{FirstGuess(wl, L, E, fluoFG, O2band)}
\arguments{
  \item{wl}{numeric vector: wavelength vector}
  \item{L}{numeric vector: measued reflected radiance}
  \item{E}{numeric vector: measured solar radiance}
  \item{fluoFG}{numeric value: fluorescence estimate derived from iFLD method}
  \item{O2band}{character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation}
}

\value{list of the first guess (including lower and upper boundaries) of the parameters needed for the optimization}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

data("FloX_data")
data("up_coeff")
data("dw_coeff")
data("wl_FloX")

#Get Target Radiance 
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
}
