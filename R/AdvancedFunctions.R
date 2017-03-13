GetReflectance<-function(
  ### Compute the Reflectance Factor from upwelling and downwelling radiance or counts
  DownwellingRadiance		##<< numeric vector or data.frame: spectrum/a of of downwelling radiance
  ,UpwellingRadiance		##<< numeric vector or data.frame: spectrum/a of of upwelling radiance
)
{
  if (class(DownwellingRadiance)!=class(UpwellingRadiance)){print("error")}
  Reflectance<-UpwellingRadiance/DownwellingRadiance
  ##value<< numeric vector or data.frame containing the computed reflectance factor.
  return(Reflectance)
}

attr(GetReflectance,"ex") <- function(){
  
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

##################################################################################################################################################################

GetRadiance<-function(
  ### Convert Digital number to Radiance
  DNSignal		##<< numeric vector or data.frame: spectrum dark current sutracted of downwelling or upwelling channel
  ,IntegrationTime		##<< numeric vector or value: integration time [same unit as the unit used for radcalcoef determination]
  ,RadCalCoeff  ##<< numeric vector: wavelength dependent vector of coefficient for calibration
)
{
  DNSignal<-data.frame(DNSignal)
  DNIT<-sweep(DNSignal, 2,IntegrationTime , `/`)  
  Rad<-sweep(DNIT, 1,RadCalCoeff , `*`)  
  ##value<< numeric vector or data.frame containing the radiance in physical unit.
  return(Rad)
}

attr(GetRadiance,"ex") <- function(){
  

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

##################################################################################################################################################################

SigTimeInterpolation<-function(
  ### Compute the linear interpolation of each wavelength of two measuremnts at a specific time occurred in between the two
  FirstTime  ##<< numeric value: time of acquiistion of the first vector
  ,SecndTime  ##<< numeric value: time of acquiistion of the second vector
  ,TimeBetween  ##<< numeric value: time at which the interpolation want to be done
  ,FirstMeas  ##<< numeric vector: spectrum acquire at FirstTime
  ,SecndMeas  ##<< numeric vector: spectrum acquire at SecndTime
)
{
  if(TimeBetween<FirstTime|TimeBetween>SecndTime) print("error")
  Int<-numeric(length(FirstMeas))  
  for (n in 1:length(FirstMeas))
  {
    if(FirstMeas[n]== SecndMeas[n]) {Int[n]=FirstMeas[n]}else{
    x<-c(FirstTime,SecndTime)
    y<-c(FirstMeas[n],SecndMeas[n])
    x_pred<-TimeBetween
    new <- data.frame(x = x_pred)
    out<-predict(lm(y ~ x), new,se.fit = TRUE);out<-out$fit
    Int[n]<-out
    }
  }
  ##value<< numeric vector containing spectrum linearly interpolated.
  return(Int)
}

attr(SigTimeInterpolation,"ex") <- function(){
  

    data("FloX_data")
    data("wl_FloX")
    
    #Define first spectrum
    spectrum1<-FloX_data$E[,2]
    #Define second spectrum
    spectrum2<-FloX_data$E[,3]
    #Compute interpolated spectrum
    spectrum_int<-SigTimeInterpolation(FirstTime = 10,SecndTime = 25,TimeBetween = 17,FirstMeas = spectrum1,SecndMeas = spectrum2)
    #plot
    x11()
    par(mar=c(5,5,2,2))
    plot(wl_FloX,spectrum1,type="l",xlab="WL [nm]",ylab="Digital Counts [-]",ylim = c(0,200000))
    lines(wl_FloX,spectrum2,col="green")
    lines(wl_FloX,spectrum_int,col="red")
    legend("topleft",col=c("black","green","red"),lty=1,cex=1.2,legend=c("spectrum1 ","spectrum2","spectrum_int"),box.col="white")
    box()
}




##################################################################################################################################################################

BandMath<-function(
  ### Compute the mathematic combination of n bands according to a specified formula and to a specified band convlution (average or gaussian)
  wl  ##<< numeric vector: wavelength vector
  ,spectrum  ##<< numeric vector or dataframe: spectrum to analyze
  ,expr  ##<< expression containing the arithmetic operations between bands. Incremental alphabetich letter are used. 
  ,bands  ##<< numeric data.frame: center wavelength and fwhm of each bands to be used. Each row of the data.frame is called incremettally considering alphabetic letters. First row a, ..., 26th row z. 
  ,fun = "mean" ##<< function for the band value calculation (e.g. mean, gaussian)
)
{
  spectrum<-data.frame(spectrum)
  names(bands)<-c("CWL","FWHM")
  n_bands<-length(bands$CWL)
  if(fun =="mean")
  {
    for(n in 1:n_bands)
    {    
      band_sel<-which(wl>=bands$CWL[n]-((bands$FWHM[n])/2) & wl<=bands$CWL[n]+((bands$FWHM[n])/2))
      subset<-data.frame(spectrum[band_sel,])
      if(n==1) {band_value<-data.frame(band_value<-apply(subset,2,mean,na.rm=TRUE))}else{
        band_value<-cbind(band_value,apply(subset,2,mean,na.rm=TRUE))}
    }
  }
  
  if(fun =="gaussian")
  {  
    gauss_convol <- apply(spectrum, 2, function(y, x) {fGaussianConvolMultiWl(x,y,CWl.v=as.array(bands$CWL),FWHM.v=bands$FWHM)},x=wl)
    reorderd<-data.frame(matrix(unlist(gauss_convol),ncol=n_bands*3,byrow = TRUE))
    band_value<-reorderd[,(n_bands+1):(n_bands+1+n_bands-1)]
  }
  
  names(band_value)<-make.unique(rep(letters, length.out = length(bands$CWL)), sep='')
  index<-data.frame(apply(band_value,1,function(x) eval(expr, band_value)));index<-as.numeric(index[,1])
  ##value<< numeric vector or value containing the spectral index selected.
  return(index)
}

attr(BandMath,"ex") <- function(){
  

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