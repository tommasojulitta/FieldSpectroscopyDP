\name{BandMath}
\alias{BandMath}
\title{BandMath}
\description{Compute the mathematic combination of n bands according to a specified formula and to a specified band convlution (average or gaussian)}
\usage{BandMath(wl, spectrum, expr, bands, fun = "mean")}
\arguments{
  \item{wl}{numeric vector: wavelength vector}
  \item{spectrum}{numeric vector or dataframe: spectrum to analyze}
  \item{expr}{expression containing the arithmetic operations between bands. Incremental alphabetich letter are used.}
  \item{bands}{numeric data.frame: center wavelength and fwhm of each bands to be used. Each row of the data.frame is called incremettally considering alphabetic letters. First row a, ..., 26th row z.}
  \item{fun}{function for the band value calculation (e.g. mean, gaussian)}
}

\value{numeric vector or value containing the spectral index selected.}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

data("FloX_data")
data("up_coeffa")
data("dw_coeff")
data("wl_FloX")

#Get Target Radiance 
L<-GetRadiance(DNSignal=FloX_data$L-FloX_data$dcL,IntegrationTime=FloX_data$IT_L/1000,RadCalCoeff=dw_coeff)
#Get Solar Radiance 
E<-GetRadiance(DNSignal=FloX_data$E-FloX_data$dcE,IntegrationTime=FloX_data$IT_E/1000,RadCalCoeff=up_coeff)
#Get Reflectance
Ref<-GetReflectance(UpwellingRadiance = L,DownwellingRadiance=E)

#Define Index expression
expressionND<-"(a-b)/(a+b)"
#Define Index band
ND_bands<-data.frame(cwl= c(780,680),fwhm=c(5,5,5,5))
#Compute index
ND<-BandMath(wl=wl_FloX,spectrum=Ref,expr=parse(text=expressionND),ND_bands,fun="mean")
ND_gaus<-BandMath(wl=wl_FloX,spectrum=Ref,expr=parse(text=expressionND),ND_bands,fun="gaussian")

#plot results    
x11()
par(mar=c(5,5,2,2))
plot(ND,pch=20,xlab="N. measurements",ylab="NDVI - 780nm & 680nm",ylim = c(-1,1),cex=3,col="dark green")
points(ND_gaus,pch=20,col="red")

}
