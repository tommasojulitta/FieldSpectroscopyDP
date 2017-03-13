\name{GetReflectance}
\alias{GetReflectance}
\title{GetReflectance}
\description{Compute the Reflectance Factor from upwelling and downwelling radiance or counts}
\usage{GetReflectance(DownwellingRadiance, UpwellingRadiance)}
\arguments{
  \item{DownwellingRadiance}{numeric vector or data.frame: spectrum/a of of downwelling radiance}
  \item{UpwellingRadiance}{numeric vector or data.frame: spectrum/a of of upwelling radiance}
}

\value{numeric vector or data.frame containing the computed reflectance factor.}

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
#Get Reflectance
Ref<-GetReflectance(UpwellingRadiance = L,DownwellingRadiance=E)
#plot
x11()    
par(mar=c(5,5,2,2))
plot(wl_FloX,Ref[,4],type="l",ylab="Reflectance [-]", xlab="WL [nm]")    

}
