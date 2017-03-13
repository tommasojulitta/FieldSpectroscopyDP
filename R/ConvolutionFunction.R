#test file: 
#Require:
#library(SDMTools)
#####################################################################################################################################################################

fConvol <- function(
  ### Function to compute the statistics (mean, median, and standard deviation) of a spectrum given a certain spectral range
  Wl=wl,  ##<< numeric array: wavelenghts spectrometer
  Wl_start=wl[1],  ##<< First wavelenght of the spectral range. By default wl[1]
  Wl_end=wl[length(wl)], ##<< Last wavelenght of the spectral range. By default last element of the array wl
  Ref=Ref, ##<< numeric array: Reflectance 
  fun='mean' ##<< character. statistics to compute (mean, sd, median). By default mean
)
{
  spectral_subset<-which(Wl>Wl_start&Wl<Wl_end)
  average_spectrum<-mean(Ref[spectral_subset],na.rm=TRUE)
  sd_spectrum<-sd(Ref[spectral_subset])
  ##value<< data frame with mean and standard deviation in the spectral range Wl_start - Wl_end
  return(data.frame(mean=average_spectrum,sd=sd_spectrum))
}
attr(fConvol,"ex") <- function(){

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


#####################################################################################################################################################################

fGaussianConvol<-function(
  ### Function to compute the gaussian convolution on a specific spectral range
  Wl=Wl,    ##<< numeric array: wavelenghts spectrometer to be deconvolved
  Ref=Ref,   ##<< numeric array: Reflectance
  CWl=CWl,   ##<< scalar: Central Wavelength
  FWHM=FWHM   ##<< scalar: Full-Width at Half Maximum of the central Wavelength
)
{  
  zzz<-dnorm(Wl, mean = CWl, sd = FWHM/(2*sqrt(2*log(2)))) #Conversion FWHM to sigma (2*sqrt(2*log(2)))*sigma
  x.norm = (zzz - min(zzz,na.rm=TRUE))/(max(zzz,na.rm=TRUE) - min(zzz,na.rm=TRUE))
  ##value<< numeric value of the gaussian convolution
  return(weighted.mean(Ref,w=x.norm, na.rm=TRUE))
}
attr(fGaussianConvol,"ex") <- function(){

  data("reflOO")
  # Select starting wavelenght
  CWl<-600
  # Select last wavelenght
  FWHM<-8
  #apply the convolution on the spectra at the defined spectal range
  res<-fGaussianConvol(Wl=reflOO$Wl, Ref=reflOO$Ref, CWl=CWl, FWHM=FWHM)
  #Plot results
  plot(reflOO$Wl, reflOO$Ref, ylim=c(0,0.7),xlim=c(400,900), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
  points(CWl, res, col="red")
  legend("topleft",col=c("red"),pch=1,cex=1.2,legend=c("convolved value"),box.col="white")
  box()
}


#####################################################################################################################################################################

fWeightAvgConvol<-function(
  ### Function to compute the convolution on a specific spectral range given the wavelengths and weights
  Wl=Wl,    ##<< numeric array: wavelenghts spectrometer to be deconvolved
  Ref=Ref,   ##<< numeric array: Reflectance
  W=W   ##<< numeric array: contains the weights of the user-defined convolution function (same length of Ref)
)
{  
  ##value<< data frame with the results of the convolution 
  return(data.frame(RefConvol=weighted.mean(Ref,w=W, na.rm=TRUE),RefConvolSD=wt.sd(Ref,w=W)))
}
attr(fWeightAvgConvol,"ex") <- function(){

  data("reflOO")
  # define the center wavelength
  CWl<-530
  # define the Full Widht at Half Maximum
  FWHM<-50
  # define the Weight for the convolution
  zzz <- dnorm(reflOO$Wl, mean = CWl, sd = FWHM/(2*sqrt(2*log(2))))
  x.norm = (zzz - min(zzz,na.rm=TRUE))/(max(zzz,na.rm=TRUE) - min(zzz,na.rm=TRUE))
  # Exctract the convolved spectrum
  res<-fWeightAvgConvol(Wl=reflOO$Wl, Ref=reflOO$Ref, W=x.norm)
  #plot results
  plot(reflOO$Wl, reflOO$Ref, ylim=c(0,2), xlim=c(400,800), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
  lines(reflOO$Wl, x.norm, col="blue")
  points(CWl, res$RefConvol, col="red")
  legend("topleft",col=c("black","red","blue"),pch=c(NA,1,NA),lty=c(1,NA,1),cex=1.2,legend=c("Reflectance","Convolved value","Convolution function"),box.col="white")
  box()
}
#####################################################################################################################################################################

fGaussianConvolSD<-function(
  ### Function to compute the standard deviation of gaussian convolution on a specific spectral range
  Wl=Wl,    ##<< numeric array: wavelenghts spectrometer to be deconvolved
  Ref=Ref,   ##<< numeric array: Reflectance
  CWl=CWl,   ##<< scalar: Central Wavelength
  FWHM=FWHM   ##<< scalar: Full-Width at Half Maximum of the central Wavelength
)
{  
  zzz<-dnorm(Wl, mean = CWl, sd = FWHM/(2*sqrt(2*log(2)))) #Conversion FWHM to sigma (2*sqrt(2*log(2)))*sigma
  x.norm = (zzz - min(zzz,na.rm=TRUE))/(max(zzz,na.rm=TRUE) - min(zzz,na.rm=TRUE))
  ##value<< numeric value: standard deviation of gaussian convolution on a specific spectral range 
  return(wt.sd(Ref,wt=x.norm))
}
attr(fGaussianConvolSD,"ex") <- function(){

  data("reflOO")
  #Define center wavelength 
  CWl<-600
  #Define Full Width at Half Maximum 
  FWHM<-8
  #Extract SD of the convolved spectrum
  res<-fGaussianConvolSD(Wl=reflOO$Wl, Ref=reflOO$Ref, CWl=CWl, FWHM=FWHM)
  #plot results
  plot(reflOO$Wl, reflOO$Ref, ylim=c(0,0.7),xlim=c(400,900), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
  points(CWl, res, col="red")
  legend("topleft",col=c("red"),pch=1,cex=1.2,legend=c("SD"),box.col="white")
  box()
}
#####################################################################################################################################################################

fGaussianConvolMultiWl<-function(
  ### Function to compute the gaussian convolution on multiple spectral range
  Wl=Wl,       ##<< numeric array: wavelenghts spectrometer to be deconvolved
  Ref=Ref,      ##<< numeric array: Reflectance
  CWl.v=CWl.v,  ##<< numeric array: array with the Central Wavelengths
  FWHM.v=FWHM.v, ##<< numeric array: array with the Full-Width at Half Maximum of the central Wavelength
  convol='Gaussian', ##<< Character: Type of convolution. Default is "Gaussian", alternative "User"
  weights=NULL  ##<< data frame: each column contains the weights of the user-defined convolution function
)
{ 
  RefConvol<-array(NA,length(CWl.v))
  RefSD<-array(NA,length(CWl.v))
  for (i in c(1:length(CWl.v)))
  {
    if ((CWl.v[i] > (min(Wl, na.rm=TRUE)+2*FWHM.v[i])) & (CWl.v[i] < (max(Wl, na.rm=TRUE)-2*FWHM.v[i])))
    {
      RefConvol[i]<-fGaussianConvol(Wl=Wl,Ref=Ref,CWl=CWl.v[i],FWHM=FWHM.v[i])
      RefSD[i]<-fGaussianConvolSD(Wl=Wl,Ref=Ref,CWl=CWl.v[i],FWHM=FWHM.v[i])
    }    
  }
  ##value<< data frame with the central Wavelength and the result of the gaussian convolution (mean and standard deviation) 
  return(data.frame(CWl.v=CWl.v,RefConvol=RefConvol,RefConvolSD=RefSD))
}
attr(fGaussianConvolMultiWl,"ex") <- function(){

  data("reflOO")
  #Define the center wavlength of the region to be convolved
  CWl.v<-c(400,500,600,700)
  #Define the corresponding Full Widht at Half Maximum
  FWHM.v<-c(8,8,8,8)
  #Extract the convoluted values
  res<-fGaussianConvolMultiWl(Wl=reflOO$Wl, Ref=reflOO$Ref, CWl.v=CWl.v, FWHM.v=FWHM.v)
  #Plot results
  plot(reflOO$Wl, reflOO$Ref, ylim=c(0,0.7),xlim=c(400,900), type="l",xlab = "Wl [nm]", ylab= "Reflectance [-]")
  points(res$CWl.v, res$RefConvol, col="red")
  legend("topleft",col=c("red"),pch=1,cex=1.2,legend=c("Convolved values"),box.col="white")
  box()
}
#####################################################################################################################################################################

fSpecRes<-function(
  ### Function to resample based on gaussian convolution according to the spectral response of different satellite missions
  Wl=Wl, 
  Ref=Ref,      ##<< numeric array: Reflectance
  convol='Gaussian', ##<< Character: Type of convolution. Default is "Gaussian", alternative "Uniform": all the bands in the spectral range are equally weighted,  User"
  Satellite='MODIS' ##<< Character: Satellite e.g. "MODIS", "LANDSAT", Default MODIS
)
{ 
  # Load configuration of satellite spectrometers (data frame with name "conf")
  if(Satellite == 'MODIS') conf<-get(data("MODIS"))
  if(Satellite == 'LANDSAT') conf<-get(data("LANDSAT"))
  if(Satellite == 'AVIRIS') conf<-get(data("AVIRIS"))
  
  CWl.v<-conf$CWl
  FWHM.v<-conf$X2FWHM/2     # If the spectral characteristics are reported as 2 x FWHM 
  RefConvol<-array(NA,length(CWl.v))
  RefSD<-array(NA,length(CWl.v))
  
  for (i in c(1:length(CWl.v))){
    if((CWl.v[i]>Wl[1]) & (CWl.v[i]<Wl[length(Wl)])){
      if(convol == 'Gaussian'){    
        RefConvol[i]<-fGaussianConvol(Wl=Wl,Ref=Ref,CWl=CWl.v[i],FWHM=FWHM.v[i])
        RefSD[i]<-fGaussianConvolSD(Wl=Wl,Ref=Ref,CWl=CWl.v[i],FWHM=FWHM.v[i])
      }
      if(convol == 'Uniform'){    
        RefConvol[i]<-fConvol(Wl=Wl,Ref=Ref,Wl_start=CWl.v[i]-FWHM.v[i],Wl_end=CWl.v[i]+FWHM.v[i])$mean
        RefSD[i]<-fConvol(Wl=Wl,Ref=Ref,Wl_start=CWl.v[i]-FWHM.v[i],Wl_end=CWl.v[i]+FWHM.v[i])$sd
      }
    }
  }
  ##value<< data frame with the central Wavelength and the result of the gaussian convolution (mean and standard deviation)
  return(data.frame(CWl.v=CWl.v,RefConvol=RefConvol,RefConvolSD=RefSD))
}
attr(fSpecRes,"ex") <- function(){

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
#####################################################################################################################################################################

#fCreateConfSatellite <- function(
#  ### Function to write the Rdata object with the configuration of the satellites
#  path=path,  ##<< character. Path with the table containing the sensor characteristics: Ch (Channels), CWl (Wavelength at the centre of the channel), SRange (spectral range)
#  confName=confName ##<< character. String with the name of the satellite configuration (The R data will be saved with this name): e.g. MODIS 
#)
#{
#  conf<-read.csv(path, sep=";")
#  save(conf, file=paste("data/",confName,".Rda", sep=""))
#}