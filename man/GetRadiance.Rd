\name{GetRadiance}
\alias{GetRadiance}
\title{GetRadiance}
\description{Convert Digital number to Radiance}
\usage{GetRadiance(DNSignal, IntegrationTime, RadCalCoeff)}
\arguments{
  \item{DNSignal}{numeric vector or data.frame: spectrum dark current sutracted of downwelling or upwelling channel}
  \item{IntegrationTime}{numeric vector or value: integration time [same unit as the unit used for radcalcoef determination]}
  \item{RadCalCoeff}{numeric vector: wavelength dependent vector of coefficient for calibration}
}

\value{numeric vector or data.frame containing the radiance in physical unit.}

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
#plot
x11()
par(mar=c(5,5,2,2))
plot(wl_FloX,E[,4],type="l",xlab="WL [nm]",ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),ylim = c(0,0.5))
lines(wl_FloX,L[,4],col="green")
legend("topleft",col=c("black","green"),lty=1,cex=1.2,legend=c("Solar Radiance ","Target Radiance"),box.col="white")
box()
}
