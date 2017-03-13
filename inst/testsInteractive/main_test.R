#packageDir <- system.file(package='FieldSpectroscopyCC')	# only works if package has been installed
#dataDir <- file.path(packageDir,"data") 

.tmp.f <- function(){
  source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyCC\\R\\AdvancedFunctions.R")
  source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyCC\\R\\BaseFunctions.R")
  source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyCC\\R\\IOFunctions.R")
  source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyDP\\R\\AdvancedFunctions.R")
  source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyDP\\R\\IOFunctions.R")
  source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyDP\\R\\ConversionFunctions.R")
  #source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyDP\\R\\SIFSYSTEMFunctions.R")
  source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyDP\\R\\FluorescenceFunctions.R")
  source("C:/Users/Tommi/Documents/R_workspace/FieldSpectroscopyDP/R/ConvolutionFunction.R")
  setwd("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyDP\\data\\")
  dataDir <- "C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyDP\\data\\"
  return(dataDir)
}
dataDir<-.tmp.f()
############################################################################################################################################################

require(RAtmosphere)

############################################################################################################################################################
#-------------------------------------------------TEST TIME INTERPOLATION-----------------------------------------------------------------------------------
############################################################################################################################################################
#select and open file 1 for test
file_1<-paste(dataDir,"file1.txt",sep="")
dc1<-DataLoadVector(file_1,sep="\t",rowSkip = 16 ,nrows =1044, dec = "." ,ncol=2)
#select and open file 2 for test
file_2<-paste(dataDir,"file2.txt",sep="")
dc2<-DataLoadVector(file_2,sep="\t",rowSkip = 16 ,nrows =1044, dec = "." ,ncol=2)
#inteprolate spectrum at a specific time
dc_int<-SigTimeInterpolation(FirstTime = 10,SecndTime = 25,TimeBetween = 11,FirstMeas = dc1,SecndMeas = dc2)
#plot results
x11()
plot(dc1,ylim=c(1600,1650))
points(dc2,col="red")
points(dc_int,pch=20)
############################################################################################################################################################
#-------------------------------------------------TEST TIME BANDMATH-----------------------------------------------------------------------------------
############################################################################################################################################################
#read data
#open file and name variable
filename<-paste(dataDir,"file3.csv",sep="")
dat<-read.csv(filename,sep=";")
wl<-dat$Wl
ref<-dat$Ref
#define wl bands and fwhm
ND_bands<-data.frame(cwl= c(800,680),fwhm=c(10,10))
#define expression to use
expression<-"(a-b)/(a+b)"
#calculate index
NDVI<-BandMath(wl,ref,expr=parse(text=expression),ND_bands,fun="mean")
############################################################################################################################################################
#-------------------------------------------------TEST CALIBRATION FUNCTION----------------------------------------------------------------------------
############################################################################################################################################################
#read calibration files
filename<-paste(dataDir,"radcalch1.csv",sep="")
caldw<-read.csv(filename,sep=";",na.strings = "Infinity",header= FALSE);caldw<-as.numeric(unlist(caldw))
filename<-paste(dataDir,"radcalch2.csv",sep="")
calup<-read.csv(filename,sep=";",na.strings = "Infinity",header= FALSE);calup<-as.numeric(unlist(calup))
filename<-paste(dataDir,"wl_flox.csv",sep="")
wl<-read.csv(filename,sep=";",na.strings = "Infinity",header= TRUE);wl<-as.numeric(unlist(wl))
#read data from SIF system
filename<-paste(dataDir,"FloXEX.csv",sep="")
dat<-ReadFloX(filename = filename)
datetime<-DateTimeFloX(dat$time,dat$date)
#Convert date time to doy.dayfract
doy.dayfract<-DateToDOY(datetime)
#Calculate SZA using function from RAtmospher packages
SZA<-SZA(timein = datetime, Lat =39.940189, Lon = 5.763964)
#Subract dark current from spectra
Downwelling<-DCSubtraction(signal=dat$E,DarkSignal = dat$dcE,type=1)
Upwelling<-DCSubtraction(signal=dat$L,DarkSignal =dat$dcL,type=1)
#Extract Integration time
IntegrationTime<-dat$IT_E/1000
#Apply calibration coefficients to obtain radiance spectra
DWrad<-as.data.frame(GetRadiance(DNSignal = Downwelling,IntegrationTime = IntegrationTime,RadCalCoeff = caldw))
IntegrationTime<-dat$IT_L/1000
Uprad<-as.data.frame(GetRadiance(DNSignal = Upwelling,IntegrationTime = IntegrationTime,RadCalCoeff = calup))
#Calculate reflectance
refrad<-GetReflectance(DWrad,Uprad);refrad[mapply(is.infinite, refrad)] <- NA
gc()
#plot ex reflectance
x11();plot(wl,refrad[,1],type="l",ylim=c(0,10))
############################################################################################################################################################
#-------------------------------------------------TEST READ SPECTRA FUNCTION--------------------------------------------------------------------------------
############################################################################################################################################################
filenames<-list.files(path = paste(dataDir,"lamp\\",sep=""), pattern = ".txt",full.names=TRUE)
IT<-ReadSpectra(filenames,Epos =  2,skip=16,nrow=1044,sep="\t")
############################################################################################################################################################
#-------------------------------------------------TEST BANDMATH FUNCTION------------------------------------------------------------------------------------
############################################################################################################################################################
#Calculate NDVI indices using bandmathover the whole matrix
ND_bands<-data.frame(cwl= c(780,680),fwhm=c(2,2))
expression<-"(a-b)/(a+b)"
#using mean convolution of the bands
NDVI<-BandMath(wl=wl,spectrum=refrad,expr=parse(text=expression),ND_bands,fun="mean");NDVI[which(NDVI<=-1)]<-NA;NDVI[which(NDVI>=1)]<-NA
#using gaussian convolution of the bands
NDVI_g<-BandMath(wl=wl,spectrum=refrad,expr=parse(text=expression),ND_bands,fun="gaussian");NDVI_g[which(NDVI_g<=-1)]<-NA;NDVI_g[which(NDVI_g>=1)]<-NA
#plot results
x11()
plot(NDVI,pch=20,ylim=c(0.5,1))
points(NDVI_g,pch=20,ylim=c(1.5,3),col="red")
#Calculate MTCI indices using bandmathover the whole matrix
MTCI_bands<-data.frame(cwl= c(753,708,681),fwhm=c(5,5,5))
expression<-"(a-b)/(b-c)"
#using mean convolution of the bands
MTCI<-BandMath(wl=wl,spectrum=refrad,expr=parse(text=expression),MTCI_bands,fun="mean");MTCI[which(MTCI<=-10)]<-NA;MTCI[which(MTCI>=10)]<-NA
#using gaussian convolution of the bands
MTCI_g<-BandMath(wl=wl,spectrum=refrad,expr=parse(text=expression),MTCI_bands,fun="gaussian");MTCI_g[which(MTCI_g<=-10)]<-NA;MTCI_g[which(MTCI_g>=10)]<-NA
#plot results
x11()
plot(MTCI,pch=20,ylim=c(1.5,3))
points(MTCI_g,pch=20,ylim=c(1.5,3),col="red")
############################################################################################################################################################
#-------------------------------------------------TEST FLUORESCENCE FUNCTION----------------------------------------------------------------------------
############################################################################################################################################################
#Single FLD
SFLD_O2A<-sFLD(wl=wl,E=DWrad,L=Uprad,fwhm =0.3,O2band="A")
SFLD_O2B<-sFLD(wl=wl,E=DWrad,L=Uprad,fwhm =0.3,O2band="B")
#3 FLD
FLD3_O2A<-FLD3(wl=wl,E=DWrad,L=Uprad,fwhm =0.3,O2band="A")
FLD3_O2B<-FLD3(wl=wl,E=DWrad,L=Uprad,fwhm =0.3,O2band="B")
#iFLD
iFLD_O2A<-iFLD(wl=wl,E=DWrad,L=Uprad,fwhm =0.3,O2band="A")
iFLD_O2B<-iFLD(wl=wl,E=DWrad,L=Uprad,fwhm =0.3,O2band="B")
#plot results
x11()
par(mfcol=c(2,3),mar=c(5,6,2,2))

plot(doy.dayfract,SFLD_O2A$Fluo*1000,pch=20,ylim=c(-100,100),ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.axis=1.4,cex.lab=1.4,main="SFLD O2A")
plot(doy.dayfract,SFLD_O2B$Fluo*1000,pch=20,ylim=c(-100,100),ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.axis=1.4,cex.lab=1.4,main="SFLD O2B")
plot(doy.dayfract,FLD3_O2A$Fluo*1000,pch=20,ylim=c(-100,100),ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.axis=1.4,cex.lab=1.4,main="3FLD O2A")
plot(doy.dayfract,FLD3_O2B$Fluo*1000,pch=20,ylim=c(-100,100),ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.axis=1.4,cex.lab=1.4,main="3FLD O2B")
plot(doy.dayfract,iFLD_O2A$Fluo*1000,pch=20,ylim=c(-100,100),ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.axis=1.4,cex.lab=1.4,main="iFLD O2A")
plot(doy.dayfract,iFLD_O2B$Fluo*1000,pch=20,ylim=c(-100,100),ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),cex.axis=1.4,cex.lab=1.4,main="iFLD O2B")