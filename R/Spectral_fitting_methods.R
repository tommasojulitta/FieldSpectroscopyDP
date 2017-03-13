
#################################################################################################################################

SpecFit<-function(
  ### Compute the spectral sum of reflected radiance and fluorescence. Inverse mode returns the residual sum of squares. Forward mode returns the modeled radiance obtained by the sum of fluorescence and radiance.
  x,  ##<< numeric vector: parameter for gaussian and spline function
  E,  ##<< numeric vector: measured solar radiance
  L,  ##<< numeric vector: measued reflected radiance
  fm, ##<< object of class "lm", deriving from first guess functoin
  wl,  ##<< numeric vector: wavelength vector
  run,  ##<< character value. If "inverse" residual sum f squares is returned. If "forward" modelled reflected radiance is returned.
  O2band  ##<< character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation
  )
{ 
  #wavelength definition according to oxygen absorption band selected
  if(O2band =="B"){sd<-8;ctr<-684-wl[1]}
  if(O2band =="A"){sd<-24;ctr<-740-wl[1]}
  # Fs 
  dwvl = wl-wl[1]
  Fi <- x[1]*exp(-((dwvl-ctr)^2/(2*sd^2))) 
  
  # Re
  fm$coefficients<-x[2:length(x)]
  R<-predict(fm, data.frame(x = wl))
  
  #modeled reflected radiance 
  L_mod = Fi + (R*E)
  
  RSS<-sum((L-L_mod)^2)

  if (run == 'inverse') return(RSS)
  if (run == 'forward') return(L_mod)
  ##value<< if "Inverse" mode the residual sum of squares is returned. If "Forward" mode a numeric vector containing the modeled radiance obtained is returned.
}



attr(SpecFit,"ex") <- function(){
  

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



#################################################################################################################################

FirstGuess<-function(
  ### Calculate the first guess of fluorescence and true reflectance based on gaussian function and spline fitting on apparent reflectance
  wl,  ##<< numeric vector: wavelength vector
  L,  ##<< numeric vector: measued reflected radiance
  E,  ##<< numeric vector: measured solar radiance
  fluoFG,  ##<< numeric value: fluorescence estimate derived from iFLD method
  O2band  ##<< character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation
)
{
  # First guess on reflectance derived from spline fitting on the apparent reflectance
  Rapp<-L/E
  spectral_region_to_exclude<-data.frame(matrix(c(686.5,690.0,759.3,768.0),ncol=2,byrow=TRUE))
  Rapp<-ExcludeSpectralRegions(wl=wl,Rapp,spectral_region_to_exclude)
  dat_sfm<-data.frame(wl,Rapp);names(dat_sfm)<-c("x","y")
  fm <- lm(y ~ bs(x,df=6), data = dat_sfm)
  sp_coef<-as.numeric(fm$coefficients)
  
  # First guess on fluorecsnce according to oxygen absorption bands
  if(O2band =="A")
  {
    sd = 24  # [lb FG ub]
    ctr = 740; ctr=ctr-wl[1]
    fsPeak <- fluoFG*(((ctr-(740-wl[1]))/sd)^2+1)
  }
  
  if(O2band =="B")
  {
    sd = 8
    ctr = 684; ctr=ctr-wl[1]
    fsPeak = fluoFG*(((ctr-(684-wl[1]))/sd)^2+1);
  }
  # First guess vector and bounds limits
  FG =c(fsPeak,sp_coef)
  lb = c(0,c(sp_coef-0.9*sp_coef))
  ub = c(15,c(sp_coef+0.9*sp_coef))
  first_guess<-data.frame(FG,lb,ub)
  first_guess<-list(fm,first_guess);names(first_guess)<-c("fm","first_guess")
  return(first_guess)
  ##value<< list of the first guess (including lower and upper boundaries) of the parameters needed for the optimization
  
}


attr(FirstGuess,"ex") <- function(){
  

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


#################################################################################################################################

SFMResults<-function(
  ### Calculate the optimized fluorescence and true reflectance
  res, ##<< list object: as output from optim.  
  wl,  ##<< numeric vector: wavelength vector
  output,  ##<< character value: FULL or VALUE referring to output expected. If FULL a data.fame of the spectrum in the considered range of fluorescence and true reflectance is returned. If VALUE the fluorescence and the true reflectance at the selected oxygen band is returned.
  fm,  ##<< object of class "lm", deriving from first guess functoin
  O2band  ##<< character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation
  )
{
  x_out<-res$par 
  dwl = wl-wl[1]
  if(O2band =="B"){sd<-8;ctr<-684-wl[1]}
  if(O2band =="A"){sd<-24;ctr<-740-wl[1]}
  G1 = x_out[1]*exp(-((dwl-ctr)^2/(2*sd^2)));
  #Fluorescence
  f_wvl = G1
  #Reflectance
  fm$coefficients<-x_out[2:length(x_out)]
  r_wvl<-predict(fm, data.frame(x = wl))
  out<-data.frame(f_wvl,r_wvl);names(out)<-c("Fluo","TrueRef")
  if(O2band =="B"){ Fluo<-c(mean(f_wvl[which(round(wl)==687)]),mean(r_wvl[which(round(wl)==687)]));names(Fluo)<-c("Fluo","TrueRef")}
  if(O2band =="A"){ Fluo<-c(mean(f_wvl[which(round(wl)==760)]),mean(r_wvl[which(round(wl)==760)]));names(Fluo)<-c("Fluo","TrueRef")}
  if(output=="FULL") return(out)
  if(output=="VALUE") return(Fluo)
  
  ##value<< numeric vector or data.frame. If output = "FULL" a data.frame containing the spectra of the estimated true reflectance and fluorescence is returned. If output = "VALUE" a vector containing the fluorescence and the reflectance is returned.
  
}




attr(SFMResults,"ex") <- function(){
  

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


#################################################################################################################################

SFM<-function(
  ### Compute spectral fitting methods on measured data. It returns fluorescence and true reflectance
  wl,  ##<< numeric vector: wavelength vector
  L,  ##<< numeric vector: measued reflected radiance
  E,  ##<< numeric vector: measured solar radiance
  fluoFG,  ##<< numeric value: fluorescence estimate derived from iFLD method
  O2band,  ##<< character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation
  output  ##<< character value: FULL or VALUE referring to output expected. If FULL a data.fame of the spectrum in the considered range of fluorescence and true reflectance is returned. If VALUE the fluorescence and the true reflectance at the selected oxygen band is returned.
)
{
  if(O2band =="B")
  { 
    range<-which(wl>684& wl<700)
    fnsc=1e-14
    fac=1e-1
  }
  if(O2band =="A")
  { 
    range<-which(wl>750& wl<780)
    fnsc=1e-12
    fac=1e-15
  }
  
  WL<-wl[range]
  E_sfm<-as.numeric(E[range])*1000
  L_sfm<-as.numeric(L[range])*1000
  fluoFG<-fluoFG*1000
  fg<-FirstGuess(wl = WL,L = L_sfm,E = E_sfm, fluo = fluoFG,O2band= O2band)
  pascales<-log10_ceiling(abs(fg$first_guess$FG)/10);  pascales[which(pascales==0)]<-1
  res<-optim(fg$first_guess$FG,fn=SpecFit,wl=WL,E=E_sfm,L=L_sfm,fm=fg$fm,run="inverse", method="L-BFGS-B",lower=fg$first_guess$lb, 
             upper = fg$first_guess$ub,O2band =O2band,control = list(parscale=pascales,fnscale=fnsc,factr=fac))
  controllb<-as.numeric(res$par) - as.numeric(fg$first_guess$lb);controllb<-which(controllb==0)
  controlub<-as.numeric(res$par) - as.numeric(fg$first_guess$ub);controlub<-which(controlub==0)

  warning<-0
  if(length(controllb)>0){warning<-1}
  if(length(controlub)>0){warning<-2}
  if(length(controllb)>0 & length(controlub)>0){warning<-3}

  convergence<-res$convergence
  costF<-res$value
  
  ###parametri conv
  ### value da optim per cost function
  Results<-SFMResults(res = res,wl = WL,output = output,fm=fg$fm,O2band = O2band) 
  #list including warnings
  out<-list(Results,convergence,costF,warning);names(out)<-c("Results","convergence","costF","warning")
  ##value<< List containing: 
  # - Results of the spectral fitting methods, (see \code{\link{SFMResults}})
  # - Convergence of the optimization, (see \code{\link{optim}})
  # - Cost function, (see \code{\link{SpecFit}})
  # - warning. 0, no warning. 1 results recalc the lower bounday first guess. 2 results recalc the upper bounday first gues. 3 results recalc the lower and upper bounday first gues. 
  return(out)
}


attr(SFM,"ex") <- function(){
  

    data("FloX_data")
    data("up_coeff")
    data("dw_coeff")
    data("wl_FloX")
    
    #Get Target Radiance 
    L<-GetRadiance(DNSignal=FloX_data$L-FloX_data$dcL,IntegrationTime=FloX_data$IT_L/1000,RadCalCoeff=dw_coeff)
    #Get Solar Radiance 
    E<-GetRadiance(DNSignal=FloX_data$E-FloX_data$dcE,IntegrationTime=FloX_data$IT_E/1000,RadCalCoeff=up_coeff)
    #Estimate fluorescence using iFLD method, used as first guess for Spectral Fitting Methods
    iFLD_O2A<-iFLD(wl=wl_FloX,E,L,fwhm =0.3,O2band="A")
    iFLD_O2B<-iFLD(wl=wl_FloX,E,L,fwhm =0.4,O2band="B")
    #Estimate fluorescence at O2A band using SFM method, used as first guess for Spectral Fitting Methods
    res<-SFM(wl = wl_FloX, L = L[,1],E= E[,1],fluoFG = iFLD_O2A$Fluo[1], O2band ="A", output = "FULL")
    sfm_FLD_O2A<-res$Results
    
    #Plot results
    x11()
    range<-which(wl_FloX>=750& wl_FloX<=780)
    WL<-wl_FloX[range]
    par(mfcol=c(1,2))
    plot(wl_FloX,(L[,1]/E[,1]),pch=20,ylim=c(0.5,1.2),xlim=c(755,770),type="l",lwd=3,col="green",xlab="Wavelength [nm]",ylab="Reflectance [-]",cex.lab=1.5,cex.axis=1.5,main = "Reflectance")
    lines(WL,sfm_FLD_O2A[,2],pch=20,col="red",lwd=2)
    grid()
    legend("topleft",col=c("green","red"),lty=1,legend=c("Apparent","Retrieved"),box.col="white",lwd=2,cex=1.5)
    box()
    plot(WL,sfm_FLD_O2A[,1],ylim=c(0,5),type="l",lwd=3,xlim=c(755,770),col="red",xlab="Wavelength [nm]",ylab="Fluorescence [mW m-2 sr-1 nm-1]",cex.lab=1.5,cex.axis=1.5,main = "Fluorescence")
    abline(h=iFLD_O2A$Fluo[1]*1000,lty=2,col="green",lwd=3)
    grid()
    legend("topleft",col=c("green","red"),lty=c(3,1),legend=c("iFLD","SFM"),box.col="white",lwd=2,cex=1.5)
    box()
    
    #Estimate fluorescence at O2B band using SFM method, used as first guess for Spectral Fitting Methods
    res<-SFM(wl = wl_FloX, L = L[,1],E= E[,1],fluoFG = iFLD_O2A$Fluo[1], O2band ="B", output = "FULL")
    sfm_FLD_O2B<-res$Results
    #Plot results
    x11()
    range<-which(wl_FloX>=684& wl_FloX<=700)
    WL<-wl_FloX[range]
    par(mfcol=c(1,2))
    plot(wl_FloX,(L[,1]/E[,1]),pch=20,ylim=c(0,0.2),xlim=c(686,698),type="l",lwd=3,col="green",xlab="WL [nm]",ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.lab=1.5,cex.axis=1.5,main = "Reflectance")
    lines(WL,sfm_FLD_O2B[,2],pch=20,col="red",lwd=2)
    grid()
    legend("topleft",col=c("green","red"),lty=1,legend=c("Apparent","Retrieved"),box.col="white",lwd=2,cex=1.5)
    box()
    plot(WL,sfm_FLD_O2B[,1],ylim=c(0,5),type="l",lwd=3,xlim=c(686,698),col="red",xlab="WL [nm]",ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.lab=1.5,cex.axis=1.5,main = "Fluorescence")
    abline(h=iFLD_O2B$Fluo[1]*1000,lty=2,col="green",lwd=3)
    grid()
    legend("topleft",col=c("green","red"),lty=c(3,1),legend=c("iFLD","SFM"),box.col="white",lwd=2,cex=1.5)
    box()
}

#################################################################################################################################

log10_ceiling <- function(
  ### Round to the next order fo magnitude
  x  ##<< numeric value or vector: number to be rounded
)
{
  10^(ceiling(log10(x)))
  ##value<< numeric value or vector. rounded to the next order of magnitude 
}