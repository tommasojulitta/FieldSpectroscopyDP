sFLD<-function(
  ### Compute the single Fraunhofer Line Discriminator
  wl  ##<< numeric vector: wavelength vector
  ,E  ##<< numeric vector or dataframe: vector or dataframe of spectral solar radiance
  ,L  ##<< numeric vector or dataframe: vector or dataframe of spectral reflected radiance
  ,fwhm = 0.5  ##<< numeric value: spectral resolution in terms of full width at half maximum of the spectrometer used 
  ,O2band  ##<< character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation
)
{
  E<-data.frame(E)
  L<-data.frame(L)
  ### Convert input to dataframe
  if(O2band=="A")
  {
    ### In case of O2band equal to A 
    wl_in<-760
    bufferIn<-5
    bufferOut<-1
    out_in<-0.7535*fwhm+2.8937
    ### parameter to identify the spectral region where to apply the the fluorescence estimation
    Ein<-StatsOnSpectra(wl,wlStart = (wl_in-bufferIn),wlEnd = (wl_in+bufferIn),E,fun="min",margin=2)
    ### Calculation of the solar radiance in the oxygen absorption feature
    wl_out<-numeric(length(Ein))
    for (n in 1: length(Ein)){wl_out[n]<-wl[which(E[,n]==Ein[n])]-out_in}
    ### In case of multiple spectra in input loop to identify the spectral region where to calculate the solar and reflected radiance outside the absorption band
    Eout<-StatsOnSpectra(wl,wlStart=(mean(wl_out)-bufferOut),wlEnd=(mean(wl_out)),E,fun="mean",margin=2)
    Lin<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),L,fun="min",margin=2)
    Lout<-StatsOnSpectra(wl,wlStart=(mean(wl_out)-bufferOut),wlEnd=(mean(wl_out)),L,fun="mean",margin=2)
    ### Calculation of the solar or reflected radiance inside and outside the oxygen absorption feature
    Fluo<-(Eout*Lin - Lout*Ein)/(Eout-Ein)
    TrueRef=(Lout-Lin)/(Eout-Ein)
    ### Calculation of the fluorescence value and the true reflectance
    out<-list(Fluo,TrueRef);names(out)<-c("Fluo","TrueRef")
    ### Lits the output
  }
  if(O2band=="B")
  {
    ### In case of O2band equal to B 
    wl_in<-687
    bufferIn<-5
    bufferOut<-1
    out_in<-0.697*fwhm+1.245
    ### parameter to identify the spectral region where to apply the the fluorescence estimation
    Ein<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),E,fun="min",margin=2)
    ### Calculation of the solar radiance in the oxygen absorption feature
    wl_out<-numeric(length(Ein))
    for (n in 1: length(Ein)){wl_out[n]<-wl[which(E[,n]==Ein[n])]-out_in}
    ### In case of multiple spectra in input loop to identify the spectral region where to calculate the solar and reflected radiance outside the absorption band
    Eout<-StatsOnSpectra(wl,wlStart=(wl_out-bufferOut),wlEnd=wl_out,E,fun="mean",margin=2)
    Lin<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),L,fun="min",margin=2)
    Lout<-StatsOnSpectra(wl,wlStart=(wl_out-bufferOut),wlEnd=wl_out,L,fun="mean",margin=2)
    ### Calculation of the solar or reflected radiance inside and outside the oxygen absorption feature
    Fluo<-(Eout*Lin - Lout*Ein)/(Eout-Ein)
    TrueRef=(Lout-Lin)/(Eout-Ein)
    ### Calculation of the fluorescence value and the true reflectance
    out<-list(Fluo,TrueRef);names(out)<-c("Fluo","TrueRef")
    ### Lits the output
  }
  ##value<< numeric list containing the computed the estimated true reflectance and fluorescence value from Single Fraunhofer Line Discriminator.
  return(out)
}


attr(sFLD,"ex") <- function(){
  

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



#####################################################################################################################################################################

FLD3<-function(
  ### Compute the 3 Fraunhofer Line Discriminator
  wl  ##<< numeric vector: wavelength vector
  ,E  ##<< numeric vector or dataframe: vector or dataframe of spectral solar radiance
  ,L  ##<< numeric vector or dataframe: vector or dataframe of spectral reflected radiance
  ,fwhm = 0.5  ##<< numeric value: spectral resolution in terms of full width at half maximum of the spectrometer used 
  ,O2band  ##<< character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation
)
{
  E<-data.frame(E)
  L<-data.frame(L)
  ### Convert input to dataframe
  if(O2band=="A")
  {
    ### In case of O2band equal to A 
    wl_in<-760
    bufferIn<-5
    bufferOut<-1
    out_in_first<-0.7535*fwhm+2.8937
    out_in_secnd<-10
    ### parameter to identify the spectral region where to apply the the fluorescence estimation
    Ein<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),E,fun="min",margin=2)
    ### Calculation of the solar radiance in the oxygen absorption feature
    wl_out_first<-numeric(length(Ein))
    wl_out_secnd<-numeric(length(Ein))
    for (n in 1: length(Ein)){wl_out_first[n]<-wl[which(E[,n]==Ein[n])]-out_in_first;wl_out_secnd[n]<-wl[which(E[,n]==Ein[n])]+out_in_secnd}
    ### In case of multiple spectra in input loop to identify the spectral region where to calculate the solar and reflected radiance outside the absorption band
    Eout_first<-StatsOnSpectra(wl,wlStart=(wl_out_first-bufferOut),wlEnd=wl_out_first,E,fun="mean",margin=2)
    Eout_secnd<-StatsOnSpectra(wl,wlStart=wl_out_secnd,wlEnd=(wl_out_secnd+bufferOut),E,fun="mean",margin=2)
    Eout<-rowMeans(cbind(Eout_first, Eout_secnd), na.rm=TRUE)
    Lin<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),L,fun="min",margin=2)
    Lout_first<-StatsOnSpectra(wl,wlStart=(mean(wl_out_first)-bufferOut),wlEnd=(mean(wl_out_first)),L,fun="mean",margin=2)
    Lout_secnd<-StatsOnSpectra(wl,wlStart=(mean(wl_out_secnd)),wlEnd=(mean(wl_out_secnd)+bufferOut),L,fun="mean",margin=2)
    Lout<-rowMeans(cbind(Lout_first, Lout_secnd), na.rm=TRUE)
    ### Calculation of the solar or reflected radiance inside and outside the oxygen absorption feature
    Fluo<-(Eout*Lin - Lout*Ein)/(Eout-Ein)
    TrueRef=(Lout-Lin)/(Eout-Ein)
    ### Calculation of the fluorescence value and the true reflectance
    out<-list(Fluo,TrueRef);names(out)<-c("Fluo","TrueRef")
    ### Lits the output
  }
  if(O2band=="B")
  {
    ### In case of O2band equal to B 
    wl_in<-687
    bufferIn<-5
    bufferOut<-1
    out_in_first<-0.697*fwhm+1.245
    out_in_secnd<-8
    ### parameter to identify the spectral region where to apply the the fluorescence estimation
    Ein<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),E,fun="min",margin=2)
    ### Calculation of the solar radiance in the oxygen absorption feature
    wl_out_first<-numeric(length(Ein))
    wl_out_secnd<-numeric(length(Ein))
    for (n in 1: length(Ein)){wl_out_first[n]<-wl[which(E[,n]==Ein[n])]-out_in_first; wl_out_secnd[n]<-wl[which(E[,n]==Ein[n])]+out_in_secnd}
    ### In case of multiple spectra in input loop to identify the spectral region where to calculate the solar and reflected radiance outside the absorption band
    Eout_first<-StatsOnSpectra(wl,wlStart=(mean(wl_out_first)-bufferOut),wlEnd=(mean(wl_out_first)),E,fun="mean",margin=2)
    Eout_secnd<-StatsOnSpectra(wl,wlStart=(mean(wl_out_secnd)),wlEnd=(mean(wl_out_secnd)+bufferOut),E,fun="mean",margin=2)
    Eout<-rowMeans(cbind(Eout_first, Eout_secnd), na.rm=TRUE)
    Lin<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),L,fun="min",margin=2)
    Lout_first<-StatsOnSpectra(wl,wlStart=(mean(wl_out_first)-bufferOut),wlEnd=(mean(wl_out_first)),L,fun="mean",margin=2)
    Lout_secnd<-StatsOnSpectra(wl,wlStart=(mean(wl_out_secnd)),wlEnd=(mean(wl_out_secnd)+bufferOut),L,fun="mean",margin=2)
    Lout<-rowMeans(cbind(Lout_first, Lout_secnd), na.rm=TRUE)
    ### Calculation of the solar or reflected radiance inside and outside the oxygen absorption feature
    Fluo<-(Eout*Lin - Lout*Ein)/(Eout-Ein)
    TrueRef=(Lout-Lin)/(Eout-Ein)
    ### Calculation of the fluorescence value and the true reflectance
    out<-list(Fluo,TrueRef);names(out)<-c("Fluo","TrueRef")
    ### Lits the output
  }
  ##value<< numeric list containing the computed the estimated true reflectance and fluorescence value from 3 Fraunhofer Line Discriminator.
  return(out)
  ###  Return the output
}



attr(FLD3,"ex") <- function(){
  

    data("FloX_data")
    data("up_coeff")
    data("dw_coeff")
    data("wl_FloX")
    
    #Get Target Radiance 
    L<-GetRadiance(DNSignal=FloX_data$L-FloX_data$dcL,IntegrationTime=FloX_data$IT_L/1000,RadCalCoeff=dw_coeff)
    #Get Solar Radiance 
    E<-GetRadiance(DNSignal=FloX_data$E-FloX_data$dcE,IntegrationTime=FloX_data$IT_E/1000,RadCalCoeff=up_coeff)
    #Estimate Sun induced chlorophyll Fluorescence at O2A band
    FLD3_O2A<-FLD3(wl=wl_FloX,E,L,fwhm =0.3,O2band="A")
    #Estimate Sun induced chlorophyll Fluorescence at O2B band
    FLD3_O2B<-FLD3(wl=wl_FloX,E,L,fwhm =0.3,O2band="B")
    #plot results
    x11()    
    par(mar=c(5,5,2,2))
    plot(FLD3_O2A$Fluo*1000,pch=20,xlab="N. measurements",ylab=expression("Radiance [mW m"^-2* "sr"^-1* "nm"^-1*"]"),ylim=c(-10,5))    
    points(FLD3_O2B$Fluo*1000,pch=20,col="red")
    legend("topleft",col=c("black","red"),pch=20,cex=1.2,legend=c("SIF O2A","SIF O2B"),box.col="white")
    box()
  
}


#####################################################################################################################################################################

iFLD<-function(
  ### Compute the improved Fraunhofer Line Discriminator
  wl  ##<< numeric vector: wavelength vector
  ,E  ##<< numeric vector or dataframe: vector or dataframe of spectral solar radiance
  ,L  ##<< numeric vector or dataframe: vector or dataframe of spectral reflected radiance
  ,fwhm = 0.5  ##<< numeric value: spectral resolution in terms of full width at half maximum of the spectrometer used 
  ,O2band  ##<< character value: A or B referring to the oxygen absorption band where to compute the fluorescence estimation
)
{
  L=data.frame(L)
  E=data.frame(E)
  R=data.frame(L/E)
  
  if(is.null(dim(L)) ==TRUE){L[which(is.finite(L)==FALSE)]<-NA;L<-data.frame(L);Ls<-L}else{
    L<-as.vector(unlist(L));L[which(is.finite(L)==FALSE)]<-NA;L<-matrix(L,nrow=dim(E)[1],ncol=dim(E)[2]);L<-data.frame(L);Ls<-L}
  
  if(is.null(dim(E)) ==TRUE){E[which(is.finite(E)==FALSE)]<-NA;E<-data.frame(E);Es<-E}else{
    E<-as.vector(unlist(E));E[which(is.finite(E)==FALSE)]<-NA;E<-matrix(E,nrow=dim(L)[1],ncol=dim(L)[2]);E<-data.frame(E);Es<-E}
  
  if(is.null(dim(R)) ==TRUE){R[which(is.finite(R)==FALSE)]<-NA;R<-data.frame(R);Rs<-R}else{
    R<-as.vector(unlist(R));R[which(is.finite(R)==FALSE)]<-NA;R<-matrix(R,nrow=dim(E)[1],ncol=dim(E)[2]);R<-data.frame(R);Rs<-R}
  ### Convert input to dataframe and create dataframe same as input to remove the oxygen absorption regions
  wl_2_s<-686;wl_2_e<-688;Rs[wl> wl_2_s & wl< wl_2_e,]<-NA
  #--------------------------------------------------------------------------------
  wl_3_s<-757;wl_3_e<-768;Rs[wl> wl_3_s & wl< wl_3_e,]<-NA
  #--------------------------------------------------------------------------------
  ### Remove the region affecting the reflectance 
  R_smoothed <- apply(Rs, 2, function(y, x) { SplineSmoothGapfilling(x, y,df=80)}, x = wl)
  ### Calculate the smoothed reflectance using smoot.spline
  wl_2_s<-680;wl_2_e<-711;Es[wl> wl_2_s & wl< wl_2_e,]<-NA
  #--------------------------------------------------------------------------------
  wl_4_s<-753;wl_4_e<-784;Es[wl> wl_4_s & wl< wl_4_e,]<-NA
  #--------------------------------------------------------------------------------
  ### Remove the region affecting the solar radiance 
  E_smoothed <- apply(Es, 2, function(y, x) { SplineSmoothGapfilling(x, y,df=80)}, x = wl)

  if(O2band=="A")
  {  
    ### In case of O2band equal to A
    wl_in<-760
    bufferIn<-5
    bufferOut<-1
    out_in<-0.7535*fwhm+2.8937
    ### parameter to identify the spectral region where to apply the the fluorescence estimation
    Ein<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),E,fun="min",margin=2)
    Lin<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),L,fun="min",margin=2)
    #wl[which(E[,12]==Ein[12])]
    R_in_smoothed<-numeric(length(Ein))
    for (n in 1: length(Ein)){R_in_smoothed[n]<-R_smoothed[which(E[,n]==Ein[n]),n]}
    E_in_smoothed=StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),E_smoothed,fun="mean",margin=2)
    ### Calculation of the solar radiance, the reflected radiance and the smothed reflectance and the smoothed radiance in the oxygen absorption feature
    wl_out<-numeric(length(Ein))
    out<-numeric(length(Ein))
    for (n in 1: length(Ein)){wl_out[n]<-wl[which(E[,n]==Ein[n])]-out_in}
    for (n in 1: length(Ein)){out[n]<-which(abs(wl-wl_out[n])==min(abs(wl-wl_out[n])))}
    ### In case of multiple spectra in input loop to identify the spectral region where to calculate the solar and reflected radiance outside the absorption band
    Rout<-numeric(length(Ein))
    for (n in 1: length(Ein)){Rout[n]<-R[out[n],n]}
    Eout<-StatsOnSpectra(wl,wlStart=(mean(wl_out)-bufferOut),wlEnd=(mean(wl_out)),E,fun="mean",margin=2)
    Lout<-StatsOnSpectra(wl,wlStart=(mean(wl_out)-bufferOut),wlEnd=(mean(wl_out)),L,fun="mean",margin=2)
    ### Calculation of the solar or reflected radiance outside and outside the oxygen absorption feature
    alpha_r=Rout/R_in_smoothed
    alpha_f=Eout/E_in_smoothed*alpha_r
    #### Calculate correction factors
    Fluo=(alpha_r*Eout*Lin - Ein*Lout)/(alpha_r*Eout-alpha_f*Ein)
    TrueRef=(Lin-Fluo)/Ein
    ### Calculation of the fluorescence value and the true reflectance
    out<-list(Fluo,TrueRef);names(out)<-c("Fluo","TrueRef")
    ### Lits the output
  }
  
  if(O2band=="B")
  { 
   ### In case of O2band equal to B
    wl_in<-687
    bufferIn<-5
    bufferOut<-2
    out_in<-0.697*fwhm+1.245
    ### parameter to identify the spectral region where to apply the the fluorescence estimation
    Ein<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),E,fun="min",margin=2)
    Lin<-StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),L,fun="min",margin=2)
    #wl[which(E[,12]==Ein[12])]
    R_in_smoothed<-numeric(length(Ein))
    for (n in 1: length(Ein)){R_in_smoothed[n]<-R_smoothed[which(E[,n]==Ein[n]),n]}
    E_in_smoothed=StatsOnSpectra(wl,wlStart=(wl_in-bufferIn),wlEnd=(wl_in+bufferIn),E_smoothed,fun="mean",margin=2)
    ### Calculation of the solar radiance, the reflected radiance and the smothed reflectance and the smoothed radiance in the oxygen absorption feature
    wl_out<-numeric(length(Ein))
    out<-numeric(length(Ein))
    for (n in 1: length(Ein)){wl_out[n]<-wl[which(E[,n]==Ein[n])]-out_in}
    for (n in 1: length(Ein)){out[n]<-which(abs(wl-wl_out[n])==min(abs(wl-wl_out[n])))}
    ### In case of multiple spectra in input loop to identify the spectral region where to calculate the solar and reflected radiance outside the absorption band
    Rout<-numeric(length(Ein))
    for (n in 1: length(Ein)){Rout[n]<-R[out[n],n]}
    Eout<-StatsOnSpectra(wl,wlStart=(mean(wl_out)-bufferOut),wlEnd=(mean(wl_out)),E,fun="mean",margin=2)
    Lout<-StatsOnSpectra(wl,wlStart=(mean(wl_out)-bufferOut),wlEnd=(mean(wl_out)),L,fun="mean",margin=2)
      ### Calculation of the solar or reflected radiance outside and outside the oxygen absorption feature
    alpha_r=Rout/R_in_smoothed
    alpha_f=Eout/E_in_smoothed*alpha_r
    #### Calculate correction factors
    Fluo=(alpha_r*Eout*Lin - Ein*Lout)/(alpha_r*Eout-alpha_f*Ein)
    TrueRef=(Lin-Fluo)/Ein
    ### Calculation of the fluorescence value and the true reflectance
    out<-list(Fluo,TrueRef);names(out)<-c("Fluo","TrueRef")
    ### Lits the output
  }
  ##value<< numeric list containing the computed the estimated true reflectance and fluorescence value from improved Fraunhofer Line Discriminator.
  return(out)
  ###  Return the output
}



attr(iFLD,"ex") <- function(){
  

    data("FloX_data")
    data("up_coeff")
    data("dw_coeff")
    data("wl_FloX")
    
    #Get Target Radiance 
    L<-GetRadiance(DNSignal=FloX_data$L-FloX_data$dcL,IntegrationTime=FloX_data$IT_L/1000,RadCalCoeff=dw_coeff)
    #Get Solar Radiance 
    E<-GetRadiance(DNSignal=FloX_data$E-FloX_data$dcE,IntegrationTime=FloX_data$IT_E/1000,RadCalCoeff=up_coeff)
    #Estimate Sun induced chlorophyll Fluorescence at O2A band
    iFLD_O2A<-iFLD(wl=wl_FloX,E,L,fwhm =0.3,O2band="A")
    #Estimate Sun induced chlorophyll Fluorescence at O2B band
    iFLD_O2B<-iFLD(wl=wl_FloX,E,L,fwhm =0.3,O2band="B")
    #plot results
    x11()    
    par(mar=c(5,5,2,2))
    plot(iFLD_O2A$Fluo*1000,pch=20,xlab="N. measurements",ylab=expression("Radiance [mW m"^-2* "sr"^-1* "nm"^-1*"]"),ylim=c(0,5))    
    points(iFLD_O2B$Fluo*1000,pch=20,col="red")
    legend("topleft",col=c("black","red"),pch=20,cex=1.2,legend=c("SIF O2A","SIF O2B"),box.col="white")
    box()
}
