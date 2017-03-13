\name{fConvol}
\alias{fConvol}
\title{fConvol}
\description{Function to compute the statistics (mean, median, and standard deviation) of a spectrum given a certain spectral range}
\usage{fConvol(Wl = wl, Wl_start = wl[1], Wl_end = wl[length(wl)], Ref = Ref, 
    fun = "mean")}
\arguments{
  \item{Wl}{numeric array: wavelenghts spectrometer}
  \item{Wl_start}{First wavelenght of the spectral range. By default wl[1]}
  \item{Wl_end}{Last wavelenght of the spectral range. By default last element of the array wl}
  \item{Ref}{numeric array: Reflectance}
  \item{fun}{character. statistics to compute (mean, sd, median). By default mean}
}

\value{data frame with mean and standard deviation in the spectral range Wl_start - Wl_end}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{

data("reflOO")
# Select starting wavelenght
Wl_start<-600
# Select last wavelenght
Wl_end<-650
#extract the mean and sd value
res<-fConvol(reflOO$Wl, Wl_start=Wl_start, Wl_end=Wl_end, reflOO$Ref)
#Plot results
plot(reflOO$Wl, reflOO$Ref, ylim=c(0,0.7),xlim=c(400,900), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
points(mean(c(Wl_start,Wl_end)), res$mean, col="red")
legend("topleft",col=c("red"),pch=1,cex=1.2,legend=c("Mean value"),box.col="white")
box()
}
