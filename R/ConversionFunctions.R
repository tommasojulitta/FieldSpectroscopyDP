##################################################################################################################################################################
DateToDOY<-function(
  ### Convert date to julian day of the year and day fraction
  datetime  ##<< numeric vector of POSIXTct class: vector of the date and time 
)
{
  datetime<-as.character(datetime)
  doy<-sapply(datetime,function(x) {
    year<-strsplit(x,"")
    year<-paste(unlist(year)[1],unlist(year)[2],unlist(year)[3],unlist(year)[4],sep="")
    doy<-difftime(x,as.POSIXct(as.Date(paste(year,"-01-01 00:00",sep=""), tzone="GMT")),units='days')
    doy<-as.numeric(doy)+1
    return(doy)}
  )
  doy<-as.numeric(doy)
  ##value<< numeric vector containing the coverted date in Day Of the Year (DOY).
  return(doy)
}

attr(DateToDOY,"ex") <- function(){
  

    data("FloX_data")
    data("wl_FloX")
    
    #define a vector of POSIXTct class
    datetime<-DateTimeFloX(FloX_data$time,FloX_data$date)
    #Convert to Julian Day Of the Year and day fraction
    doy.dayfract<-DateToDOY(datetime)
}




##################################################################################################################################################################

DateTimeFloX<-function(
  ### Convert FloX date format in as.POSIXct class
  time ##<< numeric vector of time saved from FloX system
  ,date ##<< numeric vector of date saved from FloX system
)
{
  split<-sapply(time, function(x) strsplit(as.character(x), ""))
  times<-lapply(split,function(x) {
    if(length(unlist(x))==6)
    {
      hour<-paste(unlist(x)[1],unlist(x)[2],sep="")
      min<-paste(unlist(x)[3],unlist(x)[4],sep="")
      sec<-paste(unlist(x)[5],unlist(x)[6],sep="")
    }    
    if(length(unlist(x))==5)
    {
      hour<-paste(0,unlist(x)[1],sep="")
      min<-paste(unlist(x)[2],unlist(x)[3],sep="")
      sec<-paste(unlist(x)[4],unlist(x)[5],sep="")
    }    
    if(length(unlist(x))==4)
    {
      hour<-"00"
      min<-paste(unlist(x)[1],unlist(x)[2],sep="")
      sec<-paste(unlist(x)[3],unlist(x)[4],sep="")  
    }
    if(length(unlist(x))==3)
    {
      hour<-"00"
      min<-paste("0",unlist(x)[1],sep="")
      sec<-paste(unlist(x)[2],unlist(x)[3],sep="")  
    }
    if(length(unlist(x))==2)
    {
      hour<-"00"
      min<-"00"
      sec<-paste(unlist(x)[1],unlist(x)[2],sep="")  
    }
    if(length(unlist(x))==1)
    {
      hour<-"00"
      min<-"00"
      sec<-paste("0",unlist(x)[1],sep="")  
    }
    dates<-paste(hour, ":",min, ":",sec,sep="")
    return(dates)
  })
  
  split<-sapply(date, function(x) strsplit(as.character(x), ""))
  dates<-lapply(split,function(x) { year<-paste(unlist(x)[1],unlist(x)[2],sep="")
  month<-paste(unlist(x)[3],unlist(x)[4],sep="")
  day<-paste(unlist(x)[5],unlist(x)[6],sep="")
  dates<-paste("20",year, "-", month, "-", day,sep="")
  return(dates)
  })
  merged<-mapply(c, dates, times, SIMPLIFY=FALSE)
  datetime<-sapply(merged, function(x) paste(x[1],"  ",x[2],sep=""));datetime<-as.POSIXct(datetime)
  ##value<< numeric vector containing the coverted date time in as.POSIXct class.
  return(datetime)
}



attr(DateTimeFloX,"ex") <- function(){
  

    data("FloX_data")
    data("wl_FloX")
    #define vector of time
    time<-FloX_data$time
    #define vector of date
    date<-FloX_data$date
    #convert to as.POSIXTct class
    datetime<-DateTimeFloX(time,date)
}