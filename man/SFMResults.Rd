\name{SFMResults}
\alias{SFMResults}
\title{SFMResults}
\description{Calculate the optimized fluorescence and true reflectance}
\usage{SFMResults(res, wl, output, fm, O2band)}
\arguments{
  \item{res}{list object: as output from optim.}
  \item{wl}{numeric vector: wavelength vector}
  \item{output}{character value: FULL or VALUE referring to output expected. If FULL a data.fame of the spectrum in the considered range of fluorescence and true reflectance is returned. If VALUE the fluorescence and the true reflectance at the selected oxygen band is returned.}
  \item{fm}{object of class "lm", deriving from first guess functoin}
  \item{O2band}{character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation}
}

\value{numeric vector or data.frame. If output = "FULL" a data.frame containing the spectra of the estimated true reflectance and fluorescence is returned. If output = "VALUE" a vector containing the fluorescence and the reflectance is returned.}

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
#Rescale parameters
pascales<-log10_ceiling(abs(fg$first_guess$FG)/10);  pascales[which(pascales==0)]<-1
#Optimise the modeled reflected radiance
res<-optim(fg$first_guess$FG,fn=SpecFit,wl=WL,E=E_sfm,L=L_sfm,fm=fg$fm,run="inverse", method="L-BFGS-B",lower=fg$first_guess$lb, 
           upper = fg$first_guess$ub,O2band ="B",control = list(parscale=pascales,fnscale=1e-14,factr=1e-1))
### Extract results
Results<-SFMResults(res = res,wl = WL,output = "VALUE",fm=fg$fm,O2band = "B") 
}
