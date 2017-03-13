ReadDFloX<-function(
  ### Read Data from csv files saved by D-FloX
  filename  ##<< character value or vector: names of the file(s) to be opened
  ,sep=";"  ##<< the field separator character
  ,na.strings = "#N/D"
  ,header=FALSE ##<< logical value indicating whether the file contains the names of the variables as its first line
  ,Ename = "QE_WR" ##<< character value: string of the name in the ASCII file of the solar radiance vector, if any
  ,dcEname = "QE_DC_WR" ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of solar irradiance, if any
  ,Lname ="QE_VEG" ##<< character value: string of the name in the ASCII file of the reflected radiance vector, if any
  ,dcLname = "QE_DC_VEG"  ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of reflected radiance vector, if any
  ,E2name ="QE_WR2"  ##<< character value: string of the name in the ASCII file of the reflected radiance vector collected at double IT, if any
  #,dcL2name = "QE_DC_VEG2"  ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of L2, if any
  ,ITEpos = 6  ##<< numeric value: value of the position (number of column) in the ASCII file of the incoming radiance
  ,ITLpos = 8  ##<< numeric value: value of the position (number of column) in the ASCII file of  the reflected radiance
  ,ITL2pos = 10  ##<< numeric value: value of the position (number of column) in the ASCII file of  the reflected radiance with double IT
  ,datepos = 2  ##<< numeric value: value of the position (number of column) in the ASCII file of the date
  ,timepos = 3  ##<< numeric value: value of the position (number of column) in the ASCII file of the time
  ,cycletimepos = 12  ##<< numeric value: value of the position (number of column) in the ASCII file of the cycle length
  ,cyclenrpos = 1  ##<< numeric value: value of the position (number of column) in the ASCII file  of the cycle number
  ,temp1pos = 14  ##<< numeric value: value of the position (number of column) in the ASCII file  of the first temperature sensor
  ,temp2pos = 16  ##<< numeric value: value of the position (number of column) in the ASCII file  of the second temperature sensor
  ,temp3pos = 18  ##<< numeric value: value of the position (number of column) in the ASCII file  of the third temperature sensor
  ,temp4pos = 20  ##<< numeric value: value of the position (number of column) in the ASCII file of the fourth temperature sensor
  ,floxfluoApos=26  ##<< numeric value: value of the position (number of column) in the ASCII file of the FloX O2A fluo estimate
  ,floxfluoBpos=28  ##<< numeric value: value of the position (number of column) in the ASCII file of the FloX O2B fluo estimate
  ,insideH=22  ##<< numeric value: value of the position (number of column) in the ASCII file of the first humidity sensor
  ,outsideH=24  ##<< numeric value: value of the position (number of column) in the ASCII file of the second humidity sensor
  ,GPStimepos = 31##<< numeric value: value of the position (number of column) in the ASCII file of the GPS time
  ,GPSdatepos = 33##<< numeric value: value of the position (number of column) in the ASCII file of the GPS date
  ,GPSLATpos = 35##<< numeric value: value of the position (number of column) in the ASCII file of the GPS LAT
  ,GPSLONpos = 37##<< numeric value: value of the position (number of column) in the ASCII file of the GPS LON
)
{
  system_data<-read.csv(filename,sep=";",na.strings = "#N/D",header=FALSE,stringsAsFactors=FALSE)
  
  where_E<-which(system_data[,1]==Ename)
  where_E2<-which(system_data[,1]==E2name)
  where_dcE<-which(system_data[,1]==dcEname)
  IT_E<-system_data[(where_E-1),ITEpos]
  
  where_L<-which(system_data[,1]==Lname)
  where_dcL<-which(system_data[,1]==dcLname)
  IT_L<-system_data[(where_E-1),ITLpos]
  
  E<-system_data[where_E,2:dim(system_data)[2]];
  classes<-lapply(E, class);classes<- as.vector(unlist(classes))
  whernot<-which(classes !="integer")
  E[,whernot]<-as.numeric(as.character(unlist(E[,whernot])))
  
  E<- data.frame(matrix(unlist(E), nrow=length(where_E), byrow=F));E<-t(E);E<-na.omit(E)
  dcE<-system_data[where_dcE,2:dim(system_data)[2]];dcE[,whernot]<-as.numeric(as.character(unlist(dcE[,whernot])));dcE<- data.frame(matrix(unlist(dcE), nrow=length(where_dcE), byrow=F));dcE<-t(dcE);dcE<-na.omit(dcE)
  L<-system_data[where_L,2:dim(system_data)[2]];L[,whernot]<-as.numeric(as.character(unlist(L[,whernot])));L<- data.frame(matrix(unlist(L), nrow=length(where_L), byrow=F));L<-t(L);L<-na.omit(L)
  dcL<-system_data[where_dcL,2:dim(system_data)[2]];dcL[,whernot]<-as.numeric(as.character(unlist(dcL[,whernot])));dcL<- data.frame(matrix(unlist(dcL), nrow=length(where_dcL), byrow=F));dcL<-t(dcL);dcL<-na.omit(dcL)
  
  E2<-system_data[where_E2,2:dim(system_data)[2]];
  classes<-lapply(E2, class);classes<- as.vector(unlist(classes))
  whernot<-which(classes !="integer")
  E2[,whernot]<-as.numeric(as.character(unlist(E2[,whernot])))
  E2<- data.frame(matrix(unlist(E2), nrow=length(where_E2), byrow=F));E2<-t(E2);E2<-na.omit(E2)
  
  date<-system_data[(where_E-1),datepos]
  time<-system_data[(where_E-1),timepos]
  GPStime<-system_data[(where_E-1),GPStimepos]
  GPSdate<-system_data[(where_E-1),GPSdatepos]
  GPSLAT<-system_data[(where_E-1),GPSLATpos]
  GPSLON<-system_data[(where_E-1),GPSLONpos]
  cycleduration<-system_data[(where_E-1),cycletimepos]
  
  cyclenr<-system_data[(where_E-1),cyclenrpos]
  temp1<-system_data[(where_E-1),temp1pos]
  temp2<-system_data[(where_E-1),temp2pos]
  temp3<-system_data[(where_E-1),temp3pos]
  temp4<-system_data[(where_E-1),temp4pos]
  Hin<-system_data[(where_E-1),insideH]
  Hout<-system_data[(where_E-1),outsideH]
  
  floxfluoA<-system_data[(where_E-1),floxfluoApos]
  floxfluoB<-system_data[(where_E-1),floxfluoBpos]
  
  rawData<-list(E,E2,dcE,L,dcL,IT_E,IT_L,date,time,cycleduration,cyclenr,temp1,temp2,temp3,temp4,Hin,Hout,GPStime,GPSdate,GPSLAT,GPSLON);names(rawData)<-c("E","E2","dcE","L","dcL","IT_E","IT_L","date","time","cycleduration","cyclenr","temp1","temp2","temp3","temp4","Hin","Hout","GPStime","GPSdate","GPSLAT","GPSLON")
  ##value<< list containing the FloX data from QE spectrometer and ancillary data collected by the system
  return(rawData)
} 



attr(ReadDFloX,"ex") <- function(){
  

    #Define path and filename
    filename<-"..."
    #load data
    dat<-ReadDFloX(filename = filename)
 
}

#####################################################################################################################################################################


ReadDFloXHR<-function(
  ### Read Data from csv files saved by D-FloX
  filename  ##<< character value or vector: names of the file(s) to be opened
  ,sep=";"  ##<< the field separator character
  ,na.strings = "#N/D"
  ,header=FALSE ##<< logical value indicating whether the file contains the names of the variables as its first line
  ,Ename = "HR_WR" ##<< character value: string of the name in the ASCII file of the solar radiance vector, if any
  ,dcEname = "HR_DC_WR" ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of solar irradiance, if any
  ,Lname ="HR_VEG" ##<< character value: string of the name in the ASCII file of the reflected radiance vector, if any
  ,dcLname = "HR_DC_VEG"  ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of reflected radiance vector, if any
  ,ITEpos = 6  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the incoming radiance
  ,ITLpos = 8  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the reflected radiance
  ,datepos = 2  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the date
  ,timepos = 3  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the time
  ,cyclenrpos = 1  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the cycle number
  ,temp1pos = 12  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the first temperature sensor
  ,temp2pos = 14  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the second temperature sensor
  ,temp3pos = 16  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the third temperature sensor
  ,insideH=16  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the first humidity sensor
  ,outsideH=18  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the second humidity sensor
)
{
  system_data<-read.csv(filename,sep=";",na.strings = "#N/D",header=FALSE,stringsAsFactors=FALSE)
  
  where_E<-which(system_data[,1]==Ename)
  where_dcE<-which(system_data[,1]==dcEname)
  IT_E<-system_data[(where_E-1),ITEpos]
  
  where_L<-which(system_data[,1]==Lname)
  where_dcL<-which(system_data[,1]==dcLname)
  IT_L<-system_data[(where_E-1),ITLpos]
  
  E<-system_data[where_E,2:dim(system_data)[2]];
  classes<-lapply(E, class);classes<- as.vector(unlist(classes))
  whernot<-which(classes !="integer")
  E[,whernot]<-as.numeric(as.character(unlist(E[,whernot])))
  
  E<- data.frame(matrix(unlist(E), nrow=length(where_E), byrow=F));E<-t(E);E<-na.omit(E)
  dcE<-system_data[where_dcE,2:dim(system_data)[2]];dcE[,whernot]<-as.numeric(as.character(unlist(dcE[,whernot])));dcE<- data.frame(matrix(unlist(dcE), nrow=length(where_dcE), byrow=F));dcE<-t(dcE);dcE<-na.omit(dcE)
  L<-system_data[where_L,2:dim(system_data)[2]];L[,whernot]<-as.numeric(as.character(unlist(L[,whernot])));L<- data.frame(matrix(unlist(L), nrow=length(where_L), byrow=F));L<-t(L);L<-na.omit(L)
  dcL<-system_data[where_dcL,2:dim(system_data)[2]];dcL[,whernot]<-as.numeric(as.character(unlist(dcL[,whernot])));dcL<- data.frame(matrix(unlist(dcL), nrow=length(where_dcL), byrow=F));dcL<-t(dcL);dcL<-na.omit(dcL)
  

  date<-system_data[(where_E-1),datepos]
  time<-system_data[(where_E-1),timepos]

  cyclenr<-system_data[(where_E-1),cyclenrpos]
  temp1<-system_data[(where_E-1),temp1pos]
  temp2<-system_data[(where_E-1),temp2pos]
  temp3<-system_data[(where_E-1),temp3pos]
  Hin<-system_data[(where_E-1),insideH]
  Hout<-system_data[(where_E-1),outsideH]
  

  rawData<-list(E,dcE,L,dcL,IT_E,IT_L,date,time,cyclenr,temp1,temp2,temp3,Hin,Hout);names(rawData)<-c("E","dcE","L","dcL","IT_E","IT_L","date","time","cyclenr","temp1","temp2","temp3","Hin","Hout")
  ##value<< list containing the FloX data from HR spectrometer and ancillary data collected by the system
  return(rawData)
} 


attr(ReadDFloXHR,"ex") <- function(){
  

    #Define path and filename
    filename<-"..."
    #load data
    dat<-ReadDFloXHR(filename = filename)
  
}


#####################################################################################################################################################################


ReadRoX<-function(
  ### Read Data from csv files saved by D-FloX
  filename  ##<< character value or vector: names of the file(s) to be opened
  ,sep=";"  ##<< the field separator character
  ,na.strings = "#N/D"
  ,header=FALSE ##<< logical value indicating whether the file contains the names of the variables as its first line
  ,Ename = "WR" ##<< character value: string of the name in the ASCII file of the solar radiance vector, if any
  ,dcEname = "DC_WR" ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of solar irradiance, if any
  ,Lname ="VEG" ##<< character value: string of the name in the ASCII file of the reflected radiance vector, if any
  ,dcLname = "DC_VEG"  ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of reflected radiance vector, if any
  ,ITEpos = 6  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the incoming radiance
  ,ITLpos = 8  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the reflected radiance
  ,datepos = 2  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the date
  ,timepos = 3  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the time
  ,cyclenrpos = 1  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the cycle number
  ,temp1pos = 12  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the first temperature sensor
  ,temp2pos = 14  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the second temperature sensor
  ,temp3pos = 16  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the third temperature sensor
  ,insideH=16  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the first humidity sensor
  ,outsideH=18  ##<< numeric value: value of the position (number of column) in the ASCII file of integration time of the second humidity sensor
)
{
  system_data<-read.csv(filename,sep=";",na.strings = "#N/D",header=FALSE,stringsAsFactors=FALSE)
  
  where_E<-which(system_data[,1]==Ename)
  where_dcE<-which(system_data[,1]==dcEname)
  IT_E<-system_data[(where_E-1),ITEpos]
  
  where_L<-which(system_data[,1]==Lname)
  where_dcL<-which(system_data[,1]==dcLname)
  IT_L<-system_data[(where_E-1),ITLpos]
  
  E<-system_data[where_E,2:dim(system_data)[2]];
  classes<-lapply(E, class);classes<- as.vector(unlist(classes))
  whernot<-which(classes !="integer")
  E[,whernot]<-as.numeric(as.character(unlist(E[,whernot])))
  
  E<- data.frame(matrix(unlist(E), nrow=length(where_E), byrow=F));E<-t(E);E<-na.omit(E)
  dcE<-system_data[where_dcE,2:dim(system_data)[2]];dcE[,whernot]<-as.numeric(as.character(unlist(dcE[,whernot])));dcE<- data.frame(matrix(unlist(dcE), nrow=length(where_dcE), byrow=F));dcE<-t(dcE);dcE<-na.omit(dcE)
  L<-system_data[where_L,2:dim(system_data)[2]];L[,whernot]<-as.numeric(as.character(unlist(L[,whernot])));L<- data.frame(matrix(unlist(L), nrow=length(where_L), byrow=F));L<-t(L);L<-na.omit(L)
  dcL<-system_data[where_dcL,2:dim(system_data)[2]];dcL[,whernot]<-as.numeric(as.character(unlist(dcL[,whernot])));dcL<- data.frame(matrix(unlist(dcL), nrow=length(where_dcL), byrow=F));dcL<-t(dcL);dcL<-na.omit(dcL)
  
  
  date<-system_data[(where_E-1),datepos]
  time<-system_data[(where_E-1),timepos]
  
  cyclenr<-system_data[(where_E-1),cyclenrpos]
  temp1<-system_data[(where_E-1),temp1pos]
  temp2<-system_data[(where_E-1),temp2pos]
  temp3<-system_data[(where_E-1),temp3pos]
  Hin<-system_data[(where_E-1),insideH]
  Hout<-system_data[(where_E-1),outsideH]
  
  
  rawData<-list(E,dcE,L,dcL,IT_E,IT_L,date,time,cyclenr,temp1,temp2,temp3,Hin,Hout);names(rawData)<-c("E","dcE","L","dcL","IT_E","IT_L","date","time","cyclenr","temp1","temp2","temp3","Hin","Hout")
  ##value<< list containing the RoX data and ancillary data collected by the system
  return(rawData)
} 


attr(ReadRoX,"ex") <- function(){
  

    #Define path and filename
    filename<-"..."
    #load data
    dat<-ReadRoX(filename = filename)
  
}

