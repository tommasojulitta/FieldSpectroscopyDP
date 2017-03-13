\name{sFLD}
\alias{sFLD}
\title{sFLD}
\description{Compute the single Fraunhofer Line Discriminator}
\usage{sFLD(wl, E, L, fwhm = 0.5, O2band)}
\arguments{
  \item{wl}{numeric vector: wavelength vector}
  \item{E}{numeric vector or dataframe: vector or dataframe of spectral solar radiance}
  \item{L}{numeric vector or dataframe: vector or dataframe of spectral reflected radiance}
  \item{fwhm}{numeric value: spectral resolution in terms of full width at half maximum of the spectrometer used}
  \item{O2band}{character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation}
}

\value{numeric list containing the computed the estimated true reflectance and fluorescence value from Single Fraunhofer Line Discriminator.}

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
#Estimate Sun induced chlorophyll Fluorescence at O2A band
SFLD_O2A<-sFLD(wl=wl_FloX,E,L,fwhm =0.3,O2band="A")
#Estimate Sun induced chlorophyll Fluorescence at O2B band
SFLD_O2B<-sFLD(wl=wl_FloX,E,L,fwhm =0.3,O2band="B")
#plot results
x11()    
par(mar=c(5,5,2,2))
plot(SFLD_O2A$Fluo*1000,pch=20,xlab="N. measurements",ylab=expression("Radiance [mW m"^-2* "sr"^-1* "nm"^-1*"]"),ylim=c(0,5))    
points(SFLD_O2B$Fluo*1000,pch=20,col="red")
legend("topleft",col=c("black","red"),pch=20,cex=1.2,legend=c("SIF O2A","SIF O2B"),box.col="white")
box()
}
