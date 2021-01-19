
#this scrip reads the Lizards GPS data from 2015 and calc HR size, overlap etc
#rm(list=ls())

#important params
HomePc=0; #for home PC 0 for lab
removeInterppoint=1;#to remove interpolated points?
SubsetData=0; #To set as 0 for not subsetting#the whole data  points is 985771 points for 2015. 


#"file:///C:/Users/ors/Dropbox/R codes/Liz_GPS_Data_R_Format.rdata"

rm(LizXYdata)
###### Helper functions ######
matlab2POS = function(x, timez = "UTC") {
  # Convert between MATLAB datenum values and R POSIXt time values.
  # 
  # Author: Luke Miller   Feb 20, 2011
  #Convert a numeric  MATLAB datenum (days since 0000-1-1 00:00) to seconds in 
  #the Unix epoch (seconds since 1970-1-1 00:00). Specify a time zone if the 
  #input datenum is anything other than the GMT/UTC time zone. 
  days = x - 719529   # 719529 = days from 1-1-0000 to 1-1-1970
  secs = days * 86400 # 86400 seconds in a day
  # This next string of functions is a complete disaster, but it works.
  # It tries to outsmart R by converting the secs value to a POSIXct value
  # in the UTC time zone, then converts that to a time/date string that 
  # should lose the time zone, and then it performs a second as.POSIXct()
  # conversion on the time/date string to get a POSIXct value in the user's 
  # specified timezone. Time zones are a goddamned nightmare.
  return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
                                        tz = 'UTC'), format = '%Y-%m-%d %H:%M', 
                             tz = 'UTC', usetz = FALSE), tz = timez))
}

#plot.SpatialPolygons + downloaded helper functions- https://github.com/edzer/sp/blob/master/R/SpatialPolygons-displayMethods.R
{plot.SpatialPolygons <- function(x, col, border = par("fg"), add=FALSE, 
                                 xlim=NULL, ylim=NULL, xpd = NULL, density = NULL, angle = 45, 
                                 pbg=NULL, axes = FALSE, lty = par("lty"), ..., setParUsrBB=FALSE,
                                 usePolypath=NULL, rule=NULL, bgMap = NULL) {
  
  if (is.null(pbg))
    pbg = par("bg") # transparent!
  if (!is(x, "SpatialPolygons")) 
    stop("Not a SpatialPolygons object")
  if (is.null(usePolypath)) usePolypath <- get_Polypath()
  if (is.null(rule)) rule <- get_PolypathRule()
  
  if (! add) 
    plot(as(x, "Spatial"), xlim=xlim, ylim=ylim, axes = axes, 
         ..., setParUsrBB=setParUsrBB, bgMap = bgMap)
  
  n <- length(slot(x, "polygons"))
  if (length(border) != n)
    border <- rep(border, n, n)
  polys <- slot(x, "polygons")
  pO <- slot(x, "plotOrder")
  if (!is.null(density)) {
    if (missing(col)) col <- par("fg")
    if (length(col) != n) col <- rep(col, n, n)
    if (length(density) != n)
      density <- rep(density, n, n)
    if (length(angle) != n)
      angle <- rep(angle, n, n)
    for (j in pO) 
      .polygonRingHoles(polys[[j]], border = border[j], 
                        xpd = xpd, density = density[j], angle = angle[j], 
                        col = col[j], pbg = pbg, lty=lty, ...) 
  } else {
    if (missing(col)) col <- NA
    if (length(col) != n) col <- rep(col, n, n)
    for (j in pO) 
      .polygonRingHoles(polys[[j]], col=col[j], 
                        border=border[j], xpd = xpd, pbg = pbg, lty=lty, ...,
                        usePolypath=usePolypath, rule=rule)
  }
}

setMethod("plot", signature(x = "SpatialPolygons", y = "missing"),
          function(x, y, ...) plot.SpatialPolygons(x, ...))

.polygonRingHoles <- function(Sr, col=NA, border=NULL, xpd=NULL, density=NULL,
                              angle=45, pbg, lty = par("lty"), ..., usePolypath=NULL,
                              rule=NULL) {
  if (!is(Sr, "Polygons")) 
    stop("Not an Polygons object")
  if (is.null(usePolypath)) usePolypath <- get_Polypath()
  if (is.null(rule)) rule <- get_PolypathRule()
  if (!is.null(density)) hatch <- TRUE
  else hatch <- FALSE
  pO <- slot(Sr, "plotOrder")
  polys <- slot(Sr, "Polygons")
  
  if (hatch) {
    for (i in pO) {
      if (!slot(polys[[i]], "hole"))
        .polygon(slot(polys[[i]], "coords"), 
                 border = border, xpd = xpd, 
                 density = density, angle = angle,
                 col=col, hatch=TRUE, lty=lty, ...)
      else .polygon(slot(polys[[i]], "coords"), 
                    border = border, xpd = xpd, col=pbg, 
                    density = NULL, lty=lty, ...)
    } 
  } else if (exists("polypath") && usePolypath) {
    Srl <- as(Sr, "Lines")
    crds <- coordinates(Srl)
    if (length(crds) == 1) mcrds <- crds[[1]]
    else {
      NAr <- as.double(c(NA, NA))
      crds1 <- lapply(crds, function(x) rbind(x, NAr))
      mcrds <- do.call(rbind, crds1)
      mcrds <- mcrds[-nrow(mcrds),]
      rownames(mcrds) <- NULL
    }
    polypath(x=mcrds[,1], y=mcrds[,2], border=border, col=col,
             lty=lty, rule=rule, xpd=xpd, ...)
  } else {
    for (i in pO) {
      if (!slot(polys[[i]], "hole"))
        .polygon(slot(polys[[i]], "coords"), 
                 border = border, xpd = xpd, 
                 col=col, lty=lty, ...)
      else .polygon(slot(polys[[i]], "coords"), 
                    border = border, xpd = xpd, col=pbg, lty=lty,
                    ...)
    }
  }
}


.polygon = function(x, y = NULL, density = NULL, angle = 45,
                    border = NULL, col = NA, lty = NULL, xpd = NULL, hatch=NA, ...) {
  if (is.na(hatch)) polygon(x = x, y = y, border = border, 
                            col = col, lty = lty, xpd = xpd, ...)
  else polygon(x = x, y = y, density = density, angle = angle, 
               border = border, lty = lty, xpd = xpd, col=col, ...)
}
}



#Required packages #######
require(R.matlab);require(CircStats); #require(boot) ; require(MASS) ;   require(fields); 
require(adehabitatHR);require(plyr);#require(lme4);require(AICcmodavg);require(ggplot2);
require(rgeos);require(sp);#require(spatstat);   require(scales);


##### Loading Lizards Data from a mat file and some filtering ####
if (HomePc==1) {
  m.input.path <- "D:\\Users\\Orr\\OneDrive\\MATLAB\\Lizards\\MatFiles\\" 
  setwd("D:\\Dropbox\\R codes\\HRandMovement codes")
  #r.out.path <- "D:\\Dropbox\\Matlab\\Lizards\\matFiles\\fromR\\" 
}else if (HomePc==0) {
  m.input.path <- "C:\\Users\\Ors\\OneDrive\\MATLAB\\Lizards\\MatFiles\\" 
  setwd("C:\\Users\\ors\\Dropbox\\R codes\\HRandMovement codes")
  #r.out.path <- "D:\\Dropbox\\Matlab\\Lizards\\matFiles\\fromR\\" 
}



#
#files = list.files(m.input.path, full.names=TRUE)
#if(! (length(files)>0)){print("no files in folder or wrong path")}  
File.name="LizardDataMatrix2015forR"; #for the 2015 data, not 2010 2009
Data.FromMat <-readMat(paste(m.input.path,File.name,".mat", sep = ""))
#str(Data.FromMat)
#columns in 1=tag 2=date&time 3=date only 4=time only 5=Lat 6=Long 7=elevation 8=speed 
#        9=Numof Sattlt 10=DOP param 11=TGSV_val 12=UTM_Easting 13=UTM_Northing 
#        14 distance from center 15 burst number, 16 PointTimeDiff 17  PointDistDiff
#        18 CalcSpeed %19 interpolated


##### processing mat file and saving it ######
#converting into a data.frame and giving names, changing date format
MetaDataOnLiz=data.frame(Data.FromMat$MetaDataOnLiz)
names(MetaDataOnLiz)=c('PropOfFilteredFix','DaysWithData','TotalBursts','AvgBurstsperday','propfilteredpoints1','filteredpoints2','firstDay','LastDay')
MetaDataOnLiz$firstDay=matlab2POS(MetaDataOnLiz$firstDay);#converting dates to r format
MetaDataOnLiz$LastDay=matlab2POS(MetaDataOnLiz$LastDay);#converting dates to r format

LizMatData=data.frame(Data.FromMat$MergedLizMatrixByBurst)
LizMatData[LizMatData=="NaN"]<-NA#replacing NaN with NA, takes time.
names(LizMatData)=c('names','dateWtime','date_only','time_only','Lat','Long','elevation','speed','NumofSattlt',
                    'DOP_param','TGSV_val','UTM_Easting','UTM_Northing','DistFromCenter','burst_number','PointTimeDiff',
                    'PointDistDiff','CalcSpeed','interpolated')

#saving days as burst with a different name for each indiv
LizMatData$DaysAsBurst=as.factor(paste(1+LizMatData$date_only-min(LizMatData$date_only),LizMatData$names,2015,sep="_"));

LizMatData$POSIXctTime=matlab2POS(LizMatData$dateWtime);#converting dates to r format
LizMatData$date_only=matlab2POS(LizMatData$date_only);#converting dates to r format
LizMatData$time_only=matlab2POS(1+LizMatData$time_only);#converting dates to r format must be more than 1 for min date

LizMatData$NameFactor= as.factor(paste(LizMatData$names,2015,sep="_"))#name_year as factor
LizMatData$names= as.factor(LizMatData$names)#name as factor

#throwing all interpolated points  if removeInterppoint= 1 deleting interpolated points #
if (removeInterppoint==1){LizMatData=LizMatData[ which(LizMatData$interpolated==0), ]}#for 2015 no interpolation was done anyway

head(LizMatData)
rm(Data.FromMat, m.input.path)

## saving the whole dataset
Name='Liz_2015GPS_Data_R_FormatAll.rdata'
save(list=c('LizMatData','MetaDataOnLiz'),file=Name)

  
#with(Data.FromMat,data.frame(names,years,Lines,GPSlinesAfterFilter,MissingDaysPerLizard, DaysRange));
#LizMatData=with(Data.FromMat,data.frame(UUTM.Easting ,UUTM.Northing,StepUnHorizAcc,StepUHorizDilPrec,UCalcSpeed,UGPSDateAsNum,UGPSTimeAsNum,UGPSyear,UGPSname,UPointDistDiff,UPointTimeDiff,UisInterpolated,Lat,Lon,UDayFirstPoint,UDayLastPoint,DistFromRoad,DistFromCenter,MindistFromDams,UObsthreatCont,UCurisotyCont,UConspesAgress,RecentSteps));


#which(tt$StepUnHorizAcc==NaN)
#str(LizMatData); View(LizMatData)


# ##### working on a substet for the HR analsis set  SubsetData to 0 to skip #####  
#LizMatData2=LizMatData[seq(from=1, to=100000,by=100),]
if (SubsetData>0){#working on data subset for the HR analysis
  LizMatData2=LizMatData[seq(from=1,to=dim(LizMatData)[1],by=SubsetData),c(1,12,13,20,21,22)]#the subset has also only some of the columns!!!
  head(LizMatData2)
  LizMatData2 <- droplevels(LizMatData2); 
#   #name_year as factor
#   LizMatData2$NameFactor= as.factor(paste(LizMatData2$UGPSname,LizMatData2$UGPSyear,sep="_"))
#   #saving days as burst with a different name for each indiv
#   LizMatData2$DaysAsBurst=as.factor(paste(1+LizMatData2$UGPSDateAsNum-min(LizMatData2$UGPSDateAsNum),LizMatData2$UGPSname,LizMatData2$UGPSyear,sep="_"));
#   
}else{LizMatData2=LizMatData}#working on a substet


##### as.ltraj class #####
print('now converting to track ltraj ')
#conversion to class of adehabitatLT
LizXYdataSubset=as.ltraj(xy = LizMatData2[,c("UTM_Easting","UTM_Northing")], 
                     date = LizMatData2$POSIXctTime, 
                     id = LizMatData2$NameFactor,
                     burst = LizMatData2$DaysAsBurst,
                     typeII=TRUE,
                     infolocs =LizMatData2);

LizXYdataFull=as.ltraj(xy = LizMatData[,c("UTM_Easting","UTM_Northing")], 
                         date = LizMatData$POSIXctTime, 
                         id = LizMatData$NameFactor,
                         burst = LizMatData$DaysAsBurst,
                         typeII=TRUE,
                         infolocs =LizMatData);
#in addiotn the infolocs include the distance and R2N and dx dy and rel.angle turnign angle
head(LizXYdataSubset[[2]])
#just to work on a small file with bursts (days) from 10 lizards
LizXYdataSample=LizXYdataSubset[id=id(LizXYdataSubset[10])]
#LizXYdataSample=LizXYdataFull  [id=id(LizXYdataFull)[1:10]]

#this will plot the tracks of all lizards. heavy
#plot(LizXYdata)
#plot(LizXYdataSample)

## Saving processed file first time, before HR analysis ####
Name='Liz_2015GPS_Data_R_FormatW_xy.rdata'
save(list=ls(),file=Name)
print('saved  Lizard GPS data  ')


## Home Range analysis (output in hactare) #######
  #head(LizMatData2);str(LizMatData2);View(LizMatData2)
#XYind_log2$Indiv2=as.factor(paste('Ind',XYind_log2$ProbSwitch1to2,XYind_log2$ItemsPerPatch,XYind_log2$iteration,XYind_log2$indiv,sep="_"))
LizDataForHRsbst2015=with(LizMatData2,data.frame(NameFactor,UTM_Easting,UTM_Northing))
coordinates(LizDataForHRsbst2015) = c("UTM_Easting", "UTM_Northing") # specify column names
HRlogSbst2015=data.frame(unique(LizMatData2$names),2015,unique(LizDataForHRsbst2015$NameFactor)) ; 
names(HRlogSbst2015)=c('Name','Year','Name_year')

# #creating a subset of just a single year
# LizDataForHRsbst_2009=with(subset(LizMatData2, UGPSDateAsNum < 734200),data.frame(NameFactor,UUTM.Easting,UUTM.Northing))
# LizDataForHRsbst_2010=with(subset(LizMatData2, UGPSDateAsNum > 734300),data.frame(NameFactor,UUTM.Easting,UUTM.Northing))
# LizDataForHRsbst_2009 <- droplevels(LizDataForHRsbst_2009); 
# LizDataForHRsbst_2010 <- droplevels(LizDataForHRsbst_2010); #removing non needed levels
# coordinates(LizDataForHRsbst_2009) = c("UUTM.Easting", "UUTM.Northing") # specify column names
# coordinates(LizDataForHRsbst_2010) = c("UUTM.Easting", "UUTM.Northing") # specify column names

#MCP
print('working on MCP100')
MCP=mcp(LizDataForHRsbst2015, percent = c(100) , unin = "m", unout = "ha")
HRlogSbst2015$MCP100=MCP$area;#rm(MCP)#mean(HRlog$MCP100)/(10^5)

#UD KERNELS 50% and 95% and their ratio in 3 ways

print('working on UD2 href not same4all')
#UD2=kernelUD(LizDataForHRsbst, h = "href", grid = 100,same4all = FALSE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 1,boundary = NULL)
#UD2=kernelUD(LizDataForHRsbst, h = "href", grid = 100,same4all = FALSE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 3,boundary = NULL)
#image(UD2)

#calculating kernels for a single year
print('working on UD2 href  same4all=TRUE')
#UD2_2009=kernelUD(LizDataForHRsbst_2009, h = "href", grid = 100,same4all = TRUE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 0.1,boundary = NULL)
proj4string(LizDataForHRsbst2015) <- CRS("+proj=utm +zone=56") # this gives the spatial reference point (UTM zone) and projection typ
#this worked:  UD2_2015=kernelUD(LizDataForHRsbst2015, h = "href", grid = 1000,same4all = TRUE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 3,boundary = NULL)
UD2_2015=kernelUD(LizDataForHRsbst2015, h = "href", grid = 1000,same4all = TRUE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 3,boundary = NULL)
#UD2_2015=kernelUD(LizDataForHRsbst2015, h = "href", grid = 100,same4all =FALSE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 3,boundary = NULL)
#UD2=kernelUD(LizDataForHRsbst, h = "href", grid = 100,same4all = FALSE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 1,boundary = NULL)

image(UD2_2015)

# #another solution to try to set the correct grid value for 
 x1 <- seq( 341000, 346000,by=200)  # where resolution is the pixel size you desire 
 y1 <- seq(6248000,6250000,by=200) 
 x1y1 <- expand.grid(x=x1,y=y1) 
 coordinates(x1y1) <- ~x+y 
 gridded(x1y1) <- TRUE 
 class(x1y1) 
kUD2_2015=kernelUD(LizDataForHRsbst2015, h = "href", grid = x1y1,same4all = FALSE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 5,boundary = NULL)
getverticeshr(kUD2_2015,percent = 95)
# 
# kud=kernelUD(LizDataForHRsbst2015,h="href", grid=x1y1, kern=c("bivnorm")) 
 

#plotting the core HR/ HR95 and exporting them to GIS ####
#Cores2009<- getverticeshr(UD2_2009,percent = 50)
Cores2015<- getverticeshr(UD2_2015,percent = 50)
HR95_2015<- getverticeshr(UD2_2015,percent = 95)
#class(Cores2009)
#tt=(Cores2009$id)
# Sex2009=c(2,	1,	2,	1,	1,	1,	1,	2,	1,	2,	2,	1,	1,	1,	2,	2,	2,	1,	1,	2,	2,	2,	2,	1,	1,	1,	1,	2,	2,	
#           2,	1,	1,	1,	2,	1,	2,	1,	1,	1,	1,	2,	1,	2,	2,	1,	1,	1,	1,	2,	2);
# Sex2010=c(2,	1,	2,	1,	1,	1,	2,	1,	2,	1,	2,	2,	2,	1,	2,	2,	2,	2,	1,	1,	1,	
#           1,	2,	2,	2,	1,	1,	2,	2,	2,	1,	1,	1,	2,	2,	2,	2,	1,	2,	2,	2,	1,	1,	1,	2,	1,	2,	2,	1,
#           1,	1,	1,	1,	1,	2,	2,	1,	1,	1,	2);
#loading data from EricP dataframe:1m0f in my file 2 is fem
load("C:/Users/ors/Dropbox/R codes/Lizard_Ticks/BT_ticksAVG.Rdata")
Cores2015$Sex2015=c(2,	2,	2,	1,	2,	2,	1,	2,	1,	2,	2,	1,	1,	2,	2,	1,	2,	1,	2,	1,	2,	1,	2,	1,	1,
          1,	1,	2,	2,	1,	2,	1,	1,	2,	2,	1,	1,	1,	2,	1,	2,	1,	1,	2,	2,	1,	2,	2,	1,	2,
          2,	2,	1,	2,	1,	2,	2,	2,	2,	1,	2,	2,	2,	1,	1,	2,	2,	2,  2,	1,	1,	2,	2,	1,	1,	2);#Eric Had two NA 41040,41227, the rest are now the same

#plot(Cores2009, col=5-Sex2009) #blue male green female
plot(Cores2015, col=5-Cores2015$Sex2015) #blue male green female
plot(HR95_2015, col=5-Cores2015$Sex2015) #blue male green female
plot.SpatialPolygons(HR95_2015, col=5-Cores2015$Sex2015)
spplot(Cores2015, zcol="Sex2015",scales=list(draw=T),main="cores2015",do.log = TRUE)
spplot(Cores2015, zcol="id",scales=list(draw=T),main="cores2015")
spplot(HR95_2015, zcol="id",scales=list(draw=T),main="HR95_2015")


## plotting only infested lizards #####
rm(Cores2015Infested)
Cores2015Infested=Cores2015[Cores2015$id %in% c('12434_2015','41211_2015','40019_2015','41207_2015','12847_2015','11885_2015',
                              '2408_2015','40174_2015','10029_2015','41221_2015','41222_2015','40722_2015',	
                              '9367_2015','40772_2015'),]
Cores2015Infested$id <- factor(Cores2015Infested$id) #removing non needed levels
summary(Cores2015Infested)
spplot(Cores2015Infested, zcol="id",scales=list(draw=T),main="Cores2015Infested")

HR95_2015Infested=HR95_2015[HR95_2015$id %in% c('12434_2015','41211_2015','40019_2015','41207_2015','12847_2015','11885_2015',
                                                '2408_2015','40174_2015','10029_2015','41221_2015','41222_2015','40722_2015',	
                                                '9367_2015','40772_2015'),]
HR95_2015Infested$id <- factor(HR95_2015Infested$id) #removing non needed levels
summary(HR95_2015Infested)
spplot(HR95_2015Infested, zcol="id",scales=list(draw=T),main="HR95_2015Infested")

## Exporting to Shapfile ######
library(rgdal)
dir.create("exportsToGIS")
writeOGR(obj=Cores2015Infested, dsn="exportsToGIS", layer="core2015Infested",  driver="ESRI Shapefile")
writeOGR(obj=HR95_2015Infested, dsn="exportsToGIS", layer="HR95_2015Infested", driver="ESRI Shapefile")
writeOGR(obj=Cores2015,         dsn="exportsToGIS", layer="core2015all",       driver="ESRI Shapefile")
writeOGR(obj=HR95_2015,         dsn="exportsToGIS", layer="HR95_2015all",      driver="ESRI Shapefile")


## finding HR center ######
summary(HR95_2015)
UD2_2015
LizYearID <- unique(HR95_2015@data$id)
HRcenterUTMCoor <- sapply(1:length(LizYearID), function(i){gCentroid(HR95_2015[HR95_2015$id==LizYearID[i],])})
head(HRcenterUTMCoor) # a list of SpatialPoints objects, each with the utm coord of the center of the HR
HRcenterUTMCoor[[5]]#they apear twice but its the same point (xmin==xmax, ymin== ymax)

#### #####

UD2b=kernel.area(UD2_2015, percent = c(50,95), standardize = FALSE,unin='m',unout='ha')
#same as above: tt=getverticeshr(UD2, percent = 95,unin='m',unout='km2')
HRlogSbst2015$kde50=(as.numeric(UD2b[1,]))
HRlogSbst2015$kde95=(as.numeric(UD2b[2,]))
HRlogSbst2015$ratio50to95=HRlogSbst2015$kde50/HRlogSbst2015$kde95

# #same4all = TRUE- gives the same resutls so not needed. there is a perfect correaltion with the above method
# print('working on UD2 href same4all=T')
# UD2=kernelUD(LizDataForHRsbst, h = "href", grid = 500,same4all = TRUE, hlim = c(0.1, 1.5), kern = "bivnorm", extent = 3,boundary = NULL)
# UD2b=kernel.area(UD2, percent = c(50,95), standardize = FALSE,unin='m',unout='ha')
# HRlogSbst$kde50s4a=(as.numeric(UD2b[1,]))
# HRlogSbst$kde95s4a=(as.numeric(UD2b[2,]))
# HRlogSbst$ratio50to95s4a=HRlogSbst$kde50s4a/HRlogSbst$kde95s4a
# hist(HRlogSbst$kde95s4a,100)
# hist(HRlogSbst$kde95,100)
# plot(HRlogSbst$kde95s4a,HRlogSbst$kde95)
# plot(HRlogSbst$kde50s4a,HRlogSbst$kde50)

#LSCV
print('working on UD2 LSCV')
UD3=kernelUD(LizDataForHRsbst, h = "LSCV",hlim = c(0.1, 2),grid = 500,kern = "bivnorm", extent = 1.5,boundary = NULL)
plotLSCV(UD3)# Diagnostic of the cross-validation
UD3b=kernel.area(UD3, percent = c(50,95), standardize = FALSE,unin='m',unout='ha')
HRlogSbst$kde95lscv=(as.numeric(UD3b[2,]))
HRlogSbst$kde50lscv=(as.numeric(UD3b[1,]))
HRlogSbst$ratio50to95lscv=HRlogSbst$kde50lscv/HRlogSbst$kde95lscv
hist(HRlogSbst$kde95lscv,100)
 
save(list=ls(),file=Name)
print('saved again with MCP, kde ')


##### plot Hist of HR sizes ##########
par(mfrow=c(2,2)) # 1 row, 3 columns
  hist(HRlogSbst2015$MCP100,breaks =10)# HR in Hactare
  hist(HRlogSbst2015$kde95,breaks =10)# HR in Hactare
  hist(HRlogSbst2015$MCP100,breaks =10)# HR in Hactare
  hist(HRlogSbst2015$kde50,breaks =10)# HR in Hactare
par(mfrow=c(1,1)) # 1 row, 3 columns


##HR sex and years effects male have bigger HR #####
## adding sex column to HRlogSbst
load('LizSexLogger.Rda')
for (i in 1:length(HRlogSbst$Name)){
  HRlogSbst$Sex[i]=LizSexLogger$gender[LizSexLogger$Lizard_ID==HRlogSbst$Name[i]]
  print(HRlogSbst$Name[i])
}

##does sex affect HR size1?
Ttest2009KDE95bysex=t.test(x=HRlogSbst$kde95[HRlogSbst$Year==2009 & HRlogSbst$Sex==1],
                               y=HRlogSbst$kde95[HRlogSbst$Year==2009 & HRlogSbst$Sex==2] )
Ttest2010KDE95bysex=t.test(x=HRlogSbst$kde95[HRlogSbst$Year==2010 & HRlogSbst$Sex==1],
                           y=HRlogSbst$kde95[HRlogSbst$Year==2010 & HRlogSbst$Sex==2] )
TtestKDE95MalebyYear=t.test(x=HRlogSbst$kde95[HRlogSbst$Year==2009 & HRlogSbst$Sex==1],
                            y=HRlogSbst$kde95[HRlogSbst$Year==2010 & HRlogSbst$Sex==1] )
TtestKDE95FembyYear=t.test(x=HRlogSbst$kde95[HRlogSbst$Year==2009 & HRlogSbst$Sex==2],
                            y=HRlogSbst$kde95[HRlogSbst$Year==2010 & HRlogSbst$Sex==2] )
TtestKDE95byYear=t.test(x=HRlogSbst$kde95[HRlogSbst$Year==2009 ],
                        y=HRlogSbst$kde95[HRlogSbst$Year==2010 ] )
TtestBothYearsKDE95bysex=t.test(x=HRlogSbst$kde95[ HRlogSbst$Sex==1],
                                y=HRlogSbst$kde95[ HRlogSbst$Sex==2] )


Ttest2009KDE50bysex=t.test(x=HRlogSbst$kde50[HRlogSbst$Year==2009 & HRlogSbst$Sex==1],
                           y=HRlogSbst$kde50[HRlogSbst$Year==2009 & HRlogSbst$Sex==2] )
Ttest2010KDE50bysex=t.test(x=HRlogSbst$kde50[HRlogSbst$Year==2010 & HRlogSbst$Sex==1],
                           y=HRlogSbst$kde50[HRlogSbst$Year==2010 & HRlogSbst$Sex==2] )
TtestKDE50MalebyYear=t.test(x=HRlogSbst$kde50[HRlogSbst$Year==2009 & HRlogSbst$Sex==1],
                            y=HRlogSbst$kde50[HRlogSbst$Year==2010 & HRlogSbst$Sex==1] )
TtestKDE50FembyYear=t.test(x=HRlogSbst$kde50[HRlogSbst$Year==2009 & HRlogSbst$Sex==2],
                           y=HRlogSbst$kde50[HRlogSbst$Year==2010 & HRlogSbst$Sex==2] )
TtestKDE50byYear=t.test(x=HRlogSbst$kde50[HRlogSbst$Year==2009 ],
                        y=HRlogSbst$kde50[HRlogSbst$Year==2010 ] )
TtestBothYearsKDE50bysex=t.test(x=HRlogSbst$kde50[ HRlogSbst$Sex==1],
                                y=HRlogSbst$kde50[ HRlogSbst$Sex==2] )

Ttest2009MCPbysex=t.test(x=HRlogSbst$MCP[HRlogSbst$Year==2009 & HRlogSbst$Sex==1],
                           y=HRlogSbst$MCP[HRlogSbst$Year==2009 & HRlogSbst$Sex==2] )
Ttest2010MCPbysex=t.test(x=HRlogSbst$MCP[HRlogSbst$Year==2010 & HRlogSbst$Sex==1],
                           y=HRlogSbst$MCP[HRlogSbst$Year==2010 & HRlogSbst$Sex==2] )

HRlogSbst$Name=as.factor(HRlogSbst$Name);
HRlogSbst$Year=as.factor(HRlogSbst$Year);
ModelKde95bySex0 =lmer(kde95 ~ 1 +        (1|Name) + (1|Year) ,data=HRlogSbst, REML = FALSE)
ModelKde95bySex1 =lmer(kde95 ~ 1 +  Sex + (1|Name) + (1|Year) ,data=HRlogSbst, REML = FALSE)

mean(HRlogSbst$kde95[ HRlogSbst$Sex==1 & HRlogSbst$Year==2009]);sd(HRlogSbst$kde95[ HRlogSbst$Sex==1 & HRlogSbst$Year==2009]);range(HRlogSbst$kde95[ HRlogSbst$Sex==1 & HRlogSbst$Year==2009])
mean(HRlogSbst$kde95[ HRlogSbst$Sex==2 & HRlogSbst$Year==2009]);sd(HRlogSbst$kde95[ HRlogSbst$Sex==2 & HRlogSbst$Year==2009]);range(HRlogSbst$kde95[ HRlogSbst$Sex==2 & HRlogSbst$Year==2009])

mean(HRlogSbst$kde50[ HRlogSbst$Sex==1 & HRlogSbst$Year==2009]);sd(HRlogSbst$kde50[ HRlogSbst$Sex==1 & HRlogSbst$Year==2009]);range(HRlogSbst$kde50[ HRlogSbst$Sex==1 & HRlogSbst$Year==2009])
mean(HRlogSbst$kde50[ HRlogSbst$Sex==2 & HRlogSbst$Year==2009]);sd(HRlogSbst$kde50[ HRlogSbst$Sex==2 & HRlogSbst$Year==2009]);range(HRlogSbst$kde50[ HRlogSbst$Sex==2 & HRlogSbst$Year==2009])

mean(HRlogSbst$kde95[ HRlogSbst$Sex==1 & HRlogSbst$Year==2010]);sd(HRlogSbst$kde95[ HRlogSbst$Sex==1 & HRlogSbst$Year==2010]);range(HRlogSbst$kde95[ HRlogSbst$Sex==1 & HRlogSbst$Year==2010])
mean(HRlogSbst$kde95[ HRlogSbst$Sex==2 & HRlogSbst$Year==2010]);sd(HRlogSbst$kde95[ HRlogSbst$Sex==2 & HRlogSbst$Year==2010]);range(HRlogSbst$kde95[ HRlogSbst$Sex==2 & HRlogSbst$Year==2010])

mean(HRlogSbst$kde50[ HRlogSbst$Sex==1 & HRlogSbst$Year==2010]);sd(HRlogSbst$kde50[ HRlogSbst$Sex==1 & HRlogSbst$Year==2010]);range(HRlogSbst$kde50[ HRlogSbst$Sex==1 & HRlogSbst$Year==2010])
mean(HRlogSbst$kde50[ HRlogSbst$Sex==2 & HRlogSbst$Year==2010]);sd(HRlogSbst$kde50[ HRlogSbst$Sex==2 & HRlogSbst$Year==2010]);range(HRlogSbst$kde50[ HRlogSbst$Sex==2 & HRlogSbst$Year==2010])

#p.adjust (TTestVolUD_mm_ObsVsShfl$p.value,method="bonferroni",n=N_multiTest_coorection);


##### HR overlap among individuals 2015 ######
#see https://www.rdocumentation.org/packages/adehabitatHR/versions/0.4.14/topics/kerneloverlap
HRoverlap95_2015= kerneloverlaphr(UD2_2015,            method = 'HR',percent = 95, conditional = TRUE) #c("HR", "PHR", "VI", "BA", "UDOI", "HD"),
HRoverlap50_2015= kerneloverlaphr(UD2_2015,            method = 'HR',percent = 50, conditional = TRUE) #c("HR", "PHR", "VI", "BA", "UDOI", "HD"),
#another alternative. works good as well, just make sure grid is ~500 otherwise it creates a column of 1s to two of the lizards: 
HRoverlap95_2015b=kerneloverlap(LizDataForHRsbst2015,  method = "HR" ,percent = 95, grid=500)#conditional = TRUE not relevant for method=HR
#just make sure grid is ~500 otherwise  lower grid values result in some lizards having less than 5 pixles or it creates a column of 1s to two of the lizards:
#method="HR" computes the proportion of the home range of one animal covered by the home range of another 
# to remove diag (all 1), and zero values:
HRoverlap95_2015[HRoverlap95_2015  < .0000001]   <- NA; diag(HRoverlap95_2015)  <- NA
HRoverlap95_2015b[HRoverlap95_2015b  < .0000001] <- NA; diag(HRoverlap95_2015b) <- NA
HRoverlap50_2015[HRoverlap50_2015  < .0000001]   <- NA; diag(HRoverlap50_2015)  <- NA

View(HRoverlap95_2015);View(HRoverlap95_2015b);View(HRoverlap50_2015);

#checking if the two methods agree. almost perfect 
plot(x=as.vector(HRoverlap95_2015b),y=as.vector(HRoverlap95_2015))
plot(x=as.vector(HRoverlap95_2015),y=as.vector(HRoverlap50_2015))

#converting to a dataframe with names:
HRoverlap95_2015b_df=data.frame(rep(colnames(HRoverlap95_2015b),each=length(colnames(HRoverlap95_2015b))),
                                rep(colnames(HRoverlap95_2015b),times=length(colnames(HRoverlap95_2015b))),
                                as.vector(HRoverlap95_2015b),
                                rep(colSums(!is.na(HRoverlap95_2015b))   ,each=length(colnames(HRoverlap95_2015b))),
                                rep(colMeans(HRoverlap95_2015b,na.rm = T),each=length(colnames(HRoverlap95_2015b)))
                                );names(HRoverlap95_2015b_df) =c('ID1col','ID2rowl','HRoverlap','NofOverlappingLiz','MeanOverLap')
#just non nans
HRoverlap95_2015b_df=HRoverlap95_2015b_df[complete.cases(HRoverlap95_2015b_df),]#removing those NAs
#another dataframe of summary stats for each lizard
HRoverlap_2015_df_condns=data.frame(colnames(HRoverlap95_2015b));names(HRoverlap95_2015_df_condns) =c('ID')
HRoverlap_2015_df_condns$N_OvrlppLiz95=colSums(!is.na(HRoverlap95_2015b)) ;
HRoverlap_2015_df_condns$MeanOverLap95=colMeans(HRoverlap95_2015b,na.rm = T);
HRoverlap_2015_df_condns$N_OvrlppLiz50=colSums(!is.na(HRoverlap50_2015)) ;
HRoverlap_2015_df_condns$MeanOverLap50=colMeans(HRoverlap50_2015,na.rm = T);

View(HRoverlap95_2015b_df)
View(HRoverlap95_2015_df_condns)
save(list=ls(),file=Name)
print('saved again with HR overlap ')

#histograms of overlaps
par(mfrow=c(2,2)) # 1 row, 3 columns
  hist(HRoverlap_2015_df_condns$N_OvrlppLiz95 , breaks= seq(from=0, to=51, by=3))
  hist(HRoverlap_2015_df_condns$MeanOverLap95, breaks= seq(from=0, to=1, by=0.05))
  hist(HRoverlap_2015_df_condns$N_OvrlppLiz50, breaks= seq(from=0, to=21, by=3))
  hist(HRoverlap_2015_df_condns$MeanOverLap50, breaks= seq(from=0, to=1, by=0.05))
  
par(mfrow=c(1,1)) 

# to remove diag (all 1), and zero values:
HRoverlap50_2015[HRoverlap50_2015  < .0000001] <- NA 
diag(HRoverlap50_2015) <- NA
View(HRoverlap50_2015)
#converting to a dataframe with names:
HRoverlap50_2015=data.frame(as.vector(HRoverlap50_2015),
                            rep(colnames(HRoverlap50_2015),times=length(colnames(HRoverlap50_2015))),
                            rep(colnames(HRoverlap50_2015),each=length(colnames(HRoverlap50_2015)))
);
names(HRoverlap50_2009) =c('HRoverlap','ID1row','ID2col')
HRoverlap50_2009=HRoverlap50_2009[complete.cases(HRoverlap50_2009),]#removing those NAs
hist(HRoverlap50_2009$HRoverlap)
##### HR overlap among individuals 2009 -to do ######
HRoverlap50_2009=kerneloverlaphr(UD2_2009, method = 'HR',percent = 50) #c("HR", "PHR", "VI", "BA", "UDOI", "HD"),
# to remove diag (all 1), and zero values:
HRoverlap50_2009[HRoverlap50_2009  < .0000001] <- NA 
diag(HRoverlap50_2009) <- NA
View(HRoverlap50_2009)
#converting to a dataframe with names:
HRoverlap50_2009=data.frame(as.vector(HRoverlap50_2009),
                            rep(colnames(HRoverlap50_2009),times=length(colnames(HRoverlap50_2009))),
                            rep(colnames(HRoverlap50_2009),each=length(colnames(HRoverlap50_2009)))
                            );
names(HRoverlap50_2009) =c('HRoverlap','ID1row','ID2col')
HRoverlap50_2009=HRoverlap50_2009[complete.cases(HRoverlap50_2009),]#removing those NAs
hist(HRoverlap50_2009$HRoverlap)


# Plot the home ranges... And the relocations for 
# INndPlot=5 #will plot the first XX individuals
# MCPplot=MCP[1:INndPlot,];XYind_log2plot=XYind_log2[1:(INndPlot*N_tmStp),]
# XYind_log2plot@data$CumIndv=as.factor(as.numeric(XYind_log2plot@data$CumIndv))
# levels(XYind_log2plot@data$CumIndv)
# plot(MCPplot);plot(XYind_log2plot, col=as.data.frame(XYind_log2plot)[,5], add=TRUE)
# MCP2=as.data.frame(MCP)
# ud3=kernelUD(XYind_log2plot[,5], h = "LSCV",same4all=T)
# plotLSCV(ud3)
# rm(MCP2,MCPplot,XYind_log2plot,INndPlot)


##### looking on movement data ######
require(adehabitatLT)
is.regular(LizXYdataSubset) #is the data regular? i know it is not
#just to have a small file for training
windows();plotltr(LizXYdataSample, "dt/3600/24")
head(infolocs(LizXYdataSample))
LizXYdataSample[[3]]
summary.ltraj(LizXYdataSample)
LizXYdataSample$dist



data(capreotf)

## computes the average speed of the roe deer in a moving window of width
## equal to 60 minutes
toto <- sliwinltr(capreotf, function(x) mean(x$dist/x$dt, na.rm = TRUE),
                  step = 30, type = "time", units = "min")

## zoom before the peak
head(toto[[1]])
plot(toto[[1]][1:538,], ty="l")
windows();
plotltr(LizXYdataSample, "dist")


##### mapping interactions to their cumulative UD 
#lower UD values are more peripherial kde100% is the most distant one





#### Codes from other sections ###################################################

##### working on the trajectory of each individual divided to days (bursts)
# ov <- over(XYind_log2, geometry(map)) #map has to be an object from class SpatialPixelsDataFrame" (package sp)
# return(all(!is.na(ov)))

da=as.POSIXct(strptime(paste(as.character(XYind_log2$Day),"1",as.character(XYind_log2$StepInDay),sep=" "),"%j %H %M" ))
#validation: length(unique(Burst))-(DaysToSimulate*N_Iter*N_indv)


XYind_log2b=as.ltraj(xy = XYind_log2[,c("x","y")], #conversion to class of adehabitatLT
                     date = da, burst=XYind_log2$burst , 
                     id = XYind_log2$CumIndv,                     
                     infolocs =XYind_log2[,1:10]);rm(da);

# head(XYind_log2b[[1]])
# plot(XYind_log2b)
# tt=ld(XYind_log2b)
# ttt=XYind_log2b[id='2']
# sum(ttt[[1]]$dist)
# plotltr(ttt, "dist")

#simple trajectory analysis daily distance and max displacement (  daily or by animal)
#cummulative travel distance  daily or by animal 
CumDyDst=unlist(lapply(1:(DaysToSimulate*N_Iter*N_indv), function(brst) 
{sum(XYind_log2b[burst=brst][[1]]$dist,na.rm=T)}))

CumDyDstMtrx= as.matrix(sapply( seq(from=1,to=(DaysToSimulate*N_Iter*N_indv), by=DaysToSimulate), function(i){
  (CumDyDst[i:(i+DaysToSimulate-1)]) }  ) )#each colum is an individual and the rows are the days (bursts)

CumDstByIndv=unlist(lapply(1:(N_Iter*N_indv), function(indv) 
{sum(bindltraj(XYind_log2b)[id=indv][[1]]$dist,na.rm=T)}))

#colMeans(CumDyDstMtrx)

#Max displacement from origin  
MaxDyDsplc=unlist(lapply(1:(DaysToSimulate*N_Iter*N_indv), function(brst) 
{max(sqrt(XYind_log2b[burst=brst][[1]]$R2n))}))

MaxDyDsplcMtrx= sapply( seq(from=1,to=(DaysToSimulate*N_Iter*N_indv), by=DaysToSimulate), function(i){
  (MaxDyDsplc[i:(i+DaysToSimulate-1)]) }  ) #each colum is an individual and the rows are the days (bursts)

MaxDsplcDstByIndv=unlist(lapply(1:(N_Iter*N_indv), function(indv) 
{max(sqrt(bindltraj(XYind_log2b)[id=indv][[1]]$R2n))}))

#net displacement from start to last locaiton
NetDyDsplc=unlist(lapply(1:(DaysToSimulate*N_Iter*N_indv), function(brst) 
{sqrt(XYind_log2b[burst=brst][[1]]$R2n)[DayLength]}))

NetDyDsplcMtrx= sapply( seq(from=1,to=(DaysToSimulate*N_Iter*N_indv), by=DaysToSimulate), function(i){
  (NetDyDsplc[i:(i+DaysToSimulate-1)]) }  ) #each colum is an individual and the rows are the days (bursts)

NetDsplcDstByIndv=unlist(lapply(1:(N_Iter*N_indv), function(indv) 
{sqrt(bindltraj(XYind_log2b)[id=indv][[1]]$R2n)[(DayLength*DaysToSimulate)]}))

par(mfrow=c(3,1)) # 1 row, 3 columns
hist(MaxDyDsplc)
hist(NetDyDsplc)
hist(CumDyDst)

par(mfrow=c(1,1)) # 1 row, 3 columns
save(list=ls(),file=Name)




# Plot the home ranges... And the relocations for 
# INndPlot=5 #will plot the first XX individuals
# MCPplot=MCP[1:INndPlot,];XYind_log2plot=XYind_log2[1:(INndPlot*N_tmStp),]
# XYind_log2plot@data$CumIndv=as.factor(as.numeric(XYind_log2plot@data$CumIndv))
# levels(XYind_log2plot@data$CumIndv)
# plot(MCPplot);plot(XYind_log2plot, col=as.data.frame(XYind_log2plot)[,5], add=TRUE)
# MCP2=as.data.frame(MCP)
# ud3=kernelUD(XYind_log2plot[,5], h = "LSCV",same4all=T)
# plotLSCV(ud3)
# rm(MCP2,MCPplot,XYind_log2plot,INndPlot)





# debugging
#tt=XYind_log2[which(with(XYind_log2@ data,ProbSwitch1to2==0.8 & ItemsPerPatch==9)),]#getting all steps in this value combinaiton
#unique(tt@data$CumIndv)#10 individuals
#tt@data$CumIndv=as.factor(as.numeric(tt@data$CumIndv))
#LSCV
#UD2tt=kernelUD(tt[,5], h = "LSCV")
#UD2btt=kernel.area(UD2tt, percent = c(50,95), standardize = FALSE)
# mean(as.numeric(UD2btt[1,]))
#mean(HRlog$ kde50lscv[with(HRlog,ProbSwitch1to2==0.8 & ItemsPerPatch==9)])
#meanHRlog$ kde50lscv[with(meanHRlog,ProbSwitch1to2==0.8 & ItemsPerPatch==9)]


save(list=ls(),file=Name)



#### what is the first day? #####
UniqueLiz=unique(LizMatData$name)
NetDyDsplcMtrx= sapply( seq(from=1,to=length(UniqueLiz) ), function(i){
   }  ) #each colum is an individual and the rows are the days (bursts)
StartingDateIndx= sapply(seq(from=1,to=length(UniqueLiz) ), function(i){which(as.numeric(LizMatData$names)==as.numeric(UniqueLiz[i]))[1]})
StartingDates=LizMatData$date_only[StartingDateIndx]
