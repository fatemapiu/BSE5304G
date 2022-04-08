objects()
rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
               data.table,foreign,maptools,dataRetrieval,gdistance)

# Need to figure out which data to download. 
##############################################
#0205551460 LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
##############################################
make_usgs_gage_list=function(siteNo = "0205551460",
                             parameterCd = c("00060","00065"),
                             start.date = "2017-05-01",  # Not frozen to not frozen
                             end.date = "2017-11-01"){    # to still not frozen
  
  USGSlist=list()   # Organize the data in a nice list as in previous labs
  USGSlist[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
  head(USGSlist$flowdata)  # Note that we have 00060 and 00065...
  # And of course we want to work in SI units so:
  USGSlist$flowdata$depth_m=USGSlist$flowdata$X_00065_00000*0.3048
  # m/ft depth
  USGSlist$flowdata$cms=USGSlist$flowdata$X_00060_00000*.02832
  # m3/ft3 flow
  USGSlist[["site"]]=readNWISsite(siteNo)
  head(USGSlist$site)
  class(USGSlist$site$dec_lat_va)
  #
  # Set the Manning Coefficient in the USGS Gage's Site Table
  #
  USGSlist$site$man_n=.035/1.49
  #
  # Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
  coordinates(USGSlist$site)=~dec_long_va+dec_lat_va
  #
  return(USGSlist)
}
USGS02056000=make_usgs_gage_list(siteNo = "02056000")
USGS0205551460=make_usgs_gage_list(siteNo ="0205551460" )
USGS02055100=make_usgs_gage_list(siteNo ="02055100" )
USGS02055000=make_usgs_gage_list(siteNo ="02055000" )
USGS07010038=make_usgs_gage_list(siteNo ="07010038" )

ab_ll=rbind(USGS02056000$site,
            USGS0205551460$site,
            USGS02055100$site,
            USGS02055000$site,
            USGS07010038$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                   trunc((180+coordinates(USGS02055000$site)[1])/6+1), 
                   " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
                      z = 12, prj = proj4_utm,expand=1)

# Lets plot the DEM and the gage locations so we can guess 
# what gages connect with what gages
#
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
# From Lab02, I know I can get an overview of streams with the 
# USGS H
Sys.unsetenv("https_proxy")
url='https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/Shape/NHD_H_03010101_HU8_Shape.zip'
curl_download(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
streams_utm=spTransform(streams,crs_utm)
plot(streams_utm,col="blue",add=T)
#zoom(mydem)

USGS02056000$flowdata=USGS02056000$flowdata[,c(1,2,3,4,5,8,10)]
head(USGS02056000$flowdata,2)  # Note that we have 00060 but missing 00065...

USGS02056000[["rating"]]=readNWISrating(USGS02056000$site$site_no)
plot(USGS02056000$rating$DEP,USGS02056000$rating$INDEP,xlab="DEP",ylab="INDEP")
USGS02056000$flowdata$X_00065_00000=approx(USGS02056000$rating$DEP,
                                           USGS02056000$rating$INDEP, xout = USGS02056000$flowdata$X_00060_00000, ties = min)$y
points(USGS02056000$flowdata$X_00060_00000,USGS02056000$flowdata$X_00065_00000,
       col="red")
#
USGS02056000$flowdata$depth_m=USGS02056000$flowdata$X_00065_00000*0.3048 #ft to m
# m/ft depth
A=SpatialPoints(USGS0205551460$site)# Up gradient site Lick Run
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
#plot(streams_utm,col="blue",add=T)
#plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS0205551460$site$L=SpatialLinesLengths(AtoB) # km to m
USGS0205551460$site$L # reach length in m
USGS0205551460$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/USGS0205551460$site$L
USGS0205551460$site$slope
USGS0205551460$flowdata$B=(USGS0205551460$site$man_n*
                             USGS0205551460$flowdata$cms)/(USGS0205551460$flowdata$depth_m^(5/3)*
                                                             sqrt(USGS0205551460$site$slope))
a=((USGS0205551460$site$man_n*
      USGS0205551460$flowdata$cms)/(1.49*USGS0205551460$flowdata$B*
                                      sqrt(USGS0205551460$site$slope)))^3/5
head(USGS0205551460$flowdata)
plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$B, main="LICKRUN TO ROANOKE RIVER AT NIAGARA, VA") 
plot(USGS0205551460$flowdata$cms,USGS0205551460$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
par(mar=c(4,4,4,4))
plot(USGS0205551460$flowdata$cms,a, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA",col='red',xlim=c(-1, 16),ylim=c(0, 1.8),xlab="Discharge [m3/s]", ylab="Flow depth[m]")

# ck
USGS0205551460$flowdata$ck =5/3*sqrt(USGS0205551460$site$slope)/USGS0205551460$site$man_n*
  (USGS0205551460$flowdata$depth_m^(2/3))

USGS0205551460$flowdata$dt =USGS0205551460$site$L/USGS0205551460$flowdata$ck

plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$dt)
USGS0205551460$flowdata$outTime=USGS0205551460$flowdata$dateTime+USGS0205551460$flowdata$dt
USGS0205551460$flowdata$newwave=
  USGS0205551460$flowdata$cms *1.1 <
  data.table::shift(USGS0205551460$flowdata$cms)
summary(USGS0205551460$flowdata$newwave)
# Add plot of the point found
len=length(USGS0205551460$flowdata$newwave)
USGS0205551460$flowdata$newwave[is.na(USGS0205551460$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS0205551460$flowdata$newwave[i]==T &
     USGS0205551460$flowdata$newwave[i-1]==T){
    USGS0205551460$flowdata$newwave[i]=F
  }
}
plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms,type="l")
points(USGS0205551460$flowdata$dateTime[USGS0205551460$flowdata$newwave],
       USGS0205551460$flowdata$cms[USGS0205551460$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS0205551460$flowdata$newwave == TRUE)
plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms,
     type="l",xlim=c(USGS0205551460$flowdata$dateTime[1109],
                     USGS0205551460$flowdata$dateTime[1109+200]))
lines(USGS0205551460$flowdata$outTime,USGS0205551460$flowdata$cms,col=2)

#........................................................................................
#USGS 02055100 TINKER CREEK NEAR DALEVILLE, VA
A=SpatialPoints(USGS02055100$site)# Up gradient site Lick Run
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
#plot(streams_utm,col="blue",add=T)
#plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS02055100$site$L=SpatialLinesLengths(AtoB) # km to m
USGS02055100$site$L # reach length in m
#.rs.unloadPackage("tidyr")
USGS02055100$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/USGS02055100$site$L
USGS02055100$site$slope
USGS02055100$flowdata$B=(USGS02055100$site$man_n*
                           USGS02055100$flowdata$cms)/(USGS02055100$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS02055100$site$slope))
a=((USGS02055100$site$man_n*
      USGS02055100$flowdata$cms)/(1.49*USGS02055100$flowdata$B*
                                    sqrt(USGS02055100$site$slope)))^3/5
head(USGS02055100$flowdata)
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$B, main="LICKRUN TO ROANOKE RIVER AT NIAGARA, VA") 
plot(USGS02055100$flowdata$cms,USGS02055100$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
par(mar=c(4,4,4,4))
plot(USGS02055100$flowdata$cms,a, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA",col='red',xlim=c(-1, 16),ylim=c(0, 1.8),xlab="Discharge [m3/s]", ylab="Flow depth[m]")

# ck
USGS02055100$flowdata$ck =5/3*sqrt(USGS02055100$site$slope)/USGS02055100$site$man_n*
  (USGS02055100$flowdata$depth_m^(2/3))

USGS02055100$flowdata$dt =USGS02055100$site$L/USGS02055100$flowdata$ck

plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$dt)
USGS02055100$flowdata$outTime=USGS02055100$flowdata$dateTime+USGS02055100$flowdata$dt
USGS02055100$flowdata$newwave=
  USGS02055100$flowdata$cms *1.1 <
  data.table::shift(USGS02055100$flowdata$cms)
summary(USGS02055100$flowdata$newwave)
# Add plot of the point found
len=length(USGS02055100$flowdata$newwave)
USGS02055100$flowdata$newwave[is.na(USGS02055100$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS02055100$flowdata$newwave[i]==T &
     USGS02055100$flowdata$newwave[i-1]==T){
    USGS02055100$flowdata$newwave[i]=F
  }
}
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,type="l")
points(USGS02055100$flowdata$dateTime[USGS02055100$flowdata$newwave],
       USGS02055100$flowdata$cms[USGS02055100$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS02055100$flowdata$newwave == TRUE)

plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,
     type="l",main="TINKER CREEK to ROANOKE RIVER AT NIAGARA,VA ",xlab="Time", ylab="Discharge[m3/s]",xlim=c(USGS02055100$flowdata$dateTime[1109],
                     USGS02055100$flowdata$dateTime[1109+200]))
lines(USGS02055100$flowdata$outTime,USGS02055100$flowdata$cms,col=2)

#...............................................................................................
#USGS 02054530 ROANOKE RIVER AT GLENVAR, VA
A=SpatialPoints(USGS02054530$site)# Up gradient site Lick Run
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
#plot(streams_utm,col="blue",add=T)
#plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS02054530$site$L=SpatialLinesLengths(AtoB) # km to m
USGS02054530$site$L # reach length in m
#.rs.unloadPackage("tidyr")
USGS02054530$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/USGS02054530$site$L
USGS02054530$site$slope
USGS02054530$flowdata$B=(USGS02054530$site$man_n*
                           USGS02054530$flowdata$cms)/(USGS02054530$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS02054530$site$slope))
a=((USGS02054530$site$man_n*
      USGS02054530$flowdata$cms)/(1.49*USGS02054530$flowdata$B*
                                    sqrt(USGS02054530$site$slope)))^3/5
head(USGS02054530$flowdata)
plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$B, main="LICKRUN TO ROANOKE RIVER AT NIAGARA, VA") 
plot(USGS02054530$flowdata$cms,USGS02054530$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
par(mar=c(4,4,4,4))
plot(USGS02054530$flowdata$cms,a, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA",col='red',xlim=c(-1, 16),ylim=c(0, 1.8),xlab="Discharge [m3/s]", ylab="Flow depth[m]")

# ck
USGS02054530$flowdata$ck =5/3*sqrt(USGS02054530$site$slope)/USGS02054530$site$man_n*
  (USGS02054530$flowdata$depth_m^(2/3))

USGS02054530$flowdata$dt =USGS02054530$site$L/USGS02054530$flowdata$ck

plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$dt)
USGS02054530$flowdata$outTime=USGS02054530$flowdata$dateTime+USGS02054530$flowdata$dt
USGS02054530$flowdata$newwave=
  USGS02054530$flowdata$cms *1.1 <
  data.table::shift(USGS02054530$flowdata$cms)
summary(USGS02054530$flowdata$newwave)
# Add plot of the point found
len=length(USGS02054530$flowdata$newwave)
USGS02054530$flowdata$newwave[is.na(USGS02054530$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS02054530$flowdata$newwave[i]==T &
     USGS02054530$flowdata$newwave[i-1]==T){
    USGS02054530$flowdata$newwave[i]=F
  }
}
plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,type="l")
points(USGS02054530$flowdata$dateTime[USGS02054530$flowdata$newwave],
       USGS02054530$flowdata$cms[USGS02054530$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS02054530$flowdata$newwave == TRUE)

plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,
     type="l",main="ROANOKE RIVER AT GLENVAR to ROANOKE RIVER AT NIAGARA,VA ",xlab="Time", ylab="Discharge[m3/s]",xlim=c(USGS02054530$flowdata$dateTime[1109],
                                                                                                             USGS02054530$flowdata$dateTime[1109+200]))
lines(USGS02054530$flowdata$outTime,USGS02054530$flowdata$cms,col=2)

#....................................................................................
#USGS 02055000 ROANOKE RIVER AT ROANOKE, VA
A=SpatialPoints(USGS02055000$site)# Up gradient site Lick Run
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
#plot(streams_utm,col="blue",add=T)
#plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS02055000$site$L=SpatialLinesLengths(AtoB) # km to m
USGS02055000$site$L # reach length in m
#.rs.unloadPackage("tidyr")
USGS02055000$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/USGS02055000$site$L
USGS02055000$site$slope
USGS02055000$flowdata$B=(USGS02055000$site$man_n*
                           USGS02055000$flowdata$cms)/(USGS02055000$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS02055000$site$slope))
a=((USGS02055000$site$man_n*
      USGS02055000$flowdata$cms)/(1.49*USGS02055000$flowdata$B*
                                    sqrt(USGS02055000$site$slope)))^3/5
head(USGS02055000$flowdata)
plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$B, main="LICKRUN TO ROANOKE RIVER AT NIAGARA, VA") 
plot(USGS02055000$flowdata$cms,USGS02055000$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
par(mar=c(4,4,4,4))
plot(USGS02055000$flowdata$cms,a, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA",col='red',xlim=c(-1, 16),ylim=c(0, 1.8),xlab="Discharge [m3/s]", ylab="Flow depth[m]")

# ck
USGS02055000$flowdata$ck =5/3*sqrt(USGS02055000$site$slope)/USGS02055000$site$man_n*
  (USGS02055000$flowdata$depth_m^(2/3))

USGS02055000$flowdata$dt =USGS02055000$site$L/USGS02055000$flowdata$ck

plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$dt)
USGS02055000$flowdata$outTime=USGS02055000$flowdata$dateTime+USGS02055000$flowdata$dt
USGS02055000$flowdata$newwave=
  USGS02055000$flowdata$cms *1.1 <
  data.table::shift(USGS02055000$flowdata$cms)
summary(USGS02055000$flowdata$newwave)
# Add plot of the point found
len=length(USGS02055000$flowdata$newwave)
USGS02055000$flowdata$newwave[is.na(USGS02055000$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS02055000$flowdata$newwave[i]==T &
     USGS02055000$flowdata$newwave[i-1]==T){
    USGS02055000$flowdata$newwave[i]=F
  }
}
plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,type="l")
points(USGS02055000$flowdata$dateTime[USGS02055000$flowdata$newwave],
       USGS02055000$flowdata$cms[USGS02055000$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS02055000$flowdata$newwave == TRUE)

plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,
     type="l",main="ROANOKE RIVER AT ROANOKE to ROANOKE RIVER AT NIAGARA,VA ",xlab="Time", ylab="Discharge[m3/s]",xlim=c(USGS02055000$flowdata$dateTime[1109],
                                                                                                                         USGS02055000$flowdata$dateTime[1109+200]))
lines(USGS02055000$flowdata$outTime,USGS02055000$flowdata$cms,col=2)

#.............................................grad...................................
objects()
rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
               data.table,foreign,maptools,dataRetrieval,gdistance)

# Need to figure out which data to download. 
##############################################
#0205551460 LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
##############################################
make_usgs_gage_list=function(siteNo = "0205551460",
                             parameterCd = c("00060","00065"),
                             start.date = "2017-05-01",  # Not frozen to not frozen
                             end.date = "2017-11-01"){    # to still not frozen
  
  USGSlist=list()   # Organize the data in a nice list as in previous labs
  USGSlist[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
  head(USGSlist$flowdata)  # Note that we have 00060 and 00065...
  # And of course we want to work in SI units so:
  USGSlist$flowdata$depth_m=USGSlist$flowdata$X_00065_00000*0.3048
  # m/ft depth
  USGSlist$flowdata$cms=USGSlist$flowdata$X_00060_00000*.02832
  # m3/ft3 flow
  USGSlist[["site"]]=readNWISsite(siteNo)
  head(USGSlist$site)
  class(USGSlist$site$dec_lat_va)
  #
  # Set the Manning Coefficient in the USGS Gage's Site Table
  #
  USGSlist$site$man_n=.035/1.49
  #
  # Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
  coordinates(USGSlist$site)=~dec_long_va+dec_lat_va
  #
  return(USGSlist)
}

#..................................................................................
#grad hw
#.................................................................................

USGS07010000=make_usgs_gage_list(siteNo = "07010000")
USGS07010035=make_usgs_gage_list(siteNo ="07010035" )
USGS07010030=make_usgs_gage_list(siteNo ="07010030" )
USGS07010022=make_usgs_gage_list(siteNo ="07010022" )

ab_ll=rbind(USGS07010000$site,
            USGS07010035$site,
            USGS07010030$site,
            USGS07010022$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                   trunc((180+coordinates(USGS07010000$site)[1])/6+1), 
                   " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
                      z = 12, prj = proj4_utm,expand=1)

# Lets plot the DEM and the gage locations so we can guess 
# what gages connect with what gages
#
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
# From Lab02, I know I can get an overview of streams with the 
# USGS H
Sys.unsetenv("https_proxy")
url='https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/Shape/NHD_H_07140101_HU8_Shape.zip'
curl_download(url,"NHD_H_07140101_HU8_Shape.zip")
unzip("NHD_H_07140101_HU8_Shape.zip",exdir="07140101")
streams=readOGR("07140101/Shape/NHDFlowline.dbf")
streams_utm=spTransform(streams,crs_utm)
plot(streams_utm,col="blue",add=T)
#zoom(mydem)

USGS07010000$flowdata=USGS07010000$flowdata[,c(1,2,3,4,5,8,10)]
head(USGS07010000$flowdata,2)  # Note that we have 00060 but missing 00065...

USGS07010000[["rating"]]=readNWISrating(USGS07010000$site$site_no)
plot(USGS07010000$rating$DEP,USGS07010000$rating$INDEP,xlab="DEP",ylab="INDEP")
USGS07010000$flowdata$X_00065_00000=approx(USGS07010000$rating$DEP,
                                           USGS07010000$rating$INDEP, xout = USGS07010000$flowdata$X_00060_00000, ties = min)$y
points(USGS07010000$flowdata$X_00060_00000,USGS07010000$flowdata$X_00065_00000,
       col="red")
#
USGS07010000$flowdata$depth_m=USGS07010000$flowdata$X_00065_00000*0.3048 #ft to m

#...................................................................
#07010035
A=SpatialPoints(USGS07010035$site)# Up gradient site Lick Run
B=SpatialPoints(USGS07010000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
#plot(streams_utm,col="blue",add=T)
#plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS07010035$site$L=SpatialLinesLengths(AtoB) # km to m
USGS07010035$site$L # reach length in m
#.rs.unloadPackage("tidyr")
USGS07010035$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/USGS07010035$site$L
USGS07010035$site$slope
USGS07010035$flowdata$B=(USGS07010035$site$man_n*
                           USGS07010035$flowdata$cms)/(USGS07010035$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS07010035$site$slope))
a=((USGS07010035$site$man_n*
      USGS07010035$flowdata$cms)/(1.49*USGS07010035$flowdata$B*
                                    sqrt(USGS07010035$site$slope)))^3/5
head(USGS07010035$flowdata)
plot(USGS07010035$flowdata$dateTime,USGS07010035$flowdata$B, main="LICKRUN TO ROANOKE RIVER AT NIAGARA, VA") 
plot(USGS07010035$flowdata$cms,USGS07010035$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
par(mar=c(4,4,4,4))
plot(USGS07010035$flowdata$cms,a, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA",col='red',xlim=c(-1, 16),ylim=c(0, 1.8),xlab="Discharge [m3/s]", ylab="Flow depth[m]")

# ck
USGS07010035$flowdata$ck =5/3*sqrt(USGS07010035$site$slope)/USGS07010035$site$man_n*
  (USGS07010035$flowdata$depth_m^(2/3))

USGS07010035$flowdata$dt =USGS07010035$site$L/USGS07010035$flowdata$ck

plot(USGS07010035$flowdata$dateTime,USGS07010035$flowdata$dt)
USGS07010035$flowdata$outTime=USGS07010035$flowdata$dateTime+USGS07010035$flowdata$dt
USGS07010035$flowdata$newwave=
  USGS07010035$flowdata$cms *1.1 <
  data.table::shift(USGS07010035$flowdata$cms)
summary(USGS07010035$flowdata$newwave)
# Add plot of the point found
len=length(USGS07010035$flowdata$newwave)
USGS07010035$flowdata$newwave[is.na(USGS07010035$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS07010035$flowdata$newwave[i]==T &
     USGS07010035$flowdata$newwave[i-1]==T){
    USGS07010035$flowdata$newwave[i]=F
  }
}
plot(USGS07010035$flowdata$dateTime,USGS07010035$flowdata$cms,type="l")
points(USGS07010035$flowdata$dateTime[USGS07010035$flowdata$newwave],
       USGS07010035$flowdata$cms[USGS07010035$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS07010035$flowdata$newwave == TRUE)

plot(USGS07010035$flowdata$dateTime,USGS07010035$flowdata$cms,
     type="l",main=" Engelholm Creek to Mississippi River at St. Louis, MO,Mo ",xlab="Time", ylab="Discharge[m3/s]",xlim=c(USGS07010035$flowdata$dateTime[1109],
                                                                                                                         USGS07010035$flowdata$dateTime[1109+200]))
lines(USGS07010035$flowdata$outTime,USGS07010035$flowdata$cms,col=2)

#...................................
#07010030
A=SpatialPoints(USGS07010030$site)# Up gradient site Lick Run
B=SpatialPoints(USGS07010000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
#plot(streams_utm,col="blue",add=T)
#plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS07010030$site$L=SpatialLinesLengths(AtoB) # km to m
USGS07010030$site$L # reach length in m
#.rs.unloadPackage("tidyr")
USGS07010030$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/USGS07010030$site$L
USGS07010030$site$slope
USGS07010030$flowdata$B=(USGS07010030$site$man_n*
                           USGS07010030$flowdata$cms)/(USGS07010030$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS07010030$site$slope))
a=((USGS07010030$site$man_n*
      USGS07010030$flowdata$cms)/(1.49*USGS07010030$flowdata$B*
                                    sqrt(USGS07010030$site$slope)))^3/5
head(USGS07010030$flowdata)
plot(USGS07010030$flowdata$dateTime,USGS07010030$flowdata$B, main="LICKRUN TO ROANOKE RIVER AT NIAGARA, VA") 
plot(USGS07010030$flowdata$cms,USGS07010030$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
par(mar=c(4,4,4,4))
plot(USGS07010030$flowdata$cms,a, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA",col='red',xlim=c(-1, 16),ylim=c(0, 1.8),xlab="Discharge [m3/s]", ylab="Flow depth[m]")

# ck
USGS07010030$flowdata$ck =5/3*sqrt(USGS07010030$site$slope)/USGS07010030$site$man_n*
  (USGS07010030$flowdata$depth_m^(2/3))

USGS07010030$flowdata$dt =USGS07010030$site$L/USGS07010030$flowdata$ck

plot(USGS07010030$flowdata$dateTime,USGS07010030$flowdata$dt)
USGS07010030$flowdata$outTime=USGS07010030$flowdata$dateTime+USGS07010030$flowdata$dt
USGS07010030$flowdata$newwave=
  USGS07010030$flowdata$cms *1.1 <
  data.table::shift(USGS07010030$flowdata$cms)
summary(USGS07010030$flowdata$newwave)
# Add plot of the point found
len=length(USGS07010030$flowdata$newwave)
USGS07010030$flowdata$newwave[is.na(USGS07010030$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS07010030$flowdata$newwave[i]==T &
     USGS07010030$flowdata$newwave[i-1]==T){
    USGS07010030$flowdata$newwave[i]=F
  }
}
plot(USGS07010030$flowdata$dateTime,USGS07010030$flowdata$cms,type="l")
points(USGS07010030$flowdata$dateTime[USGS07010030$flowdata$newwave],
       USGS07010030$flowdata$cms[USGS07010030$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS07010030$flowdata$newwave == TRUE)

plot(USGS07010030$flowdata$dateTime,USGS07010030$flowdata$cms,
     type="l",main=" River Des Peres Tributary to Mississippi River at St. Louis, MO,Mo ",xlab="Time", ylab="Discharge[m3/s]",xlim=c(USGS07010030$flowdata$dateTime[1109],
                                                                                                                           USGS07010030$flowdata$dateTime[1109+200]))
lines(USGS07010030$flowdata$outTime,USGS07010030$flowdata$cms,col=2)

#.........................................................................
#07010022
A=SpatialPoints(USGS07010022$site)# Up gradient site Lick Run
B=SpatialPoints(USGS07010000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
#plot(streams_utm,col="blue",add=T)
#plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS07010022$site$L=SpatialLinesLengths(AtoB) # km to m
USGS07010022$site$L # reach length in m
#.rs.unloadPackage("tidyr")
USGS07010022$site$slope=(extract(mydem,A_utm)-extract(mydem,B_utm))/USGS07010022$site$L
USGS07010022$site$slope
USGS07010022$flowdata$B=(USGS07010022$site$man_n*
                           USGS07010022$flowdata$cms)/(USGS07010022$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS07010022$site$slope))
a=((USGS07010022$site$man_n*
      USGS07010022$flowdata$cms)/(1.49*USGS07010022$flowdata$B*
                                    sqrt(USGS07010022$site$slope)))^3/5
head(USGS07010022$flowdata)
plot(USGS07010022$flowdata$dateTime,USGS07010022$flowdata$B, main="LICKRUN TO ROANOKE RIVER AT NIAGARA, VA") 
plot(USGS07010022$flowdata$cms,USGS07010022$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
par(mar=c(4,4,4,4))
plot(USGS07010022$flowdata$cms,a, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA",col='red',xlim=c(-1, 16),ylim=c(0, 1.8),xlab="Discharge [m3/s]", ylab="Flow depth[m]")

# ck
USGS07010022$flowdata$ck =5/3*sqrt(USGS07010022$site$slope)/USGS07010022$site$man_n*
  (USGS07010022$flowdata$depth_m^(2/3))

USGS07010022$flowdata$dt =USGS07010022$site$L/USGS07010022$flowdata$ck

plot(USGS07010022$flowdata$dateTime,USGS07010022$flowdata$dt)
USGS07010022$flowdata$outTime=USGS07010022$flowdata$dateTime+USGS07010022$flowdata$dt
USGS07010022$flowdata$newwave=
  USGS07010022$flowdata$cms *1.1 <
  data.table::shift(USGS07010022$flowdata$cms)
summary(USGS07010022$flowdata$newwave)
# Add plot of the point found
len=length(USGS07010022$flowdata$newwave)
USGS07010022$flowdata$newwave[is.na(USGS07010022$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS07010022$flowdata$newwave[i]==T &
     USGS07010022$flowdata$newwave[i-1]==T){
    USGS07010022$flowdata$newwave[i]=F
  }
}
plot(USGS07010022$flowdata$dateTime,USGS07010022$flowdata$cms,type="l")
points(USGS07010022$flowdata$dateTime[USGS07010022$flowdata$newwave],
       USGS07010022$flowdata$cms[USGS07010022$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS07010022$flowdata$newwave == TRUE)

plot(USGS07010022$flowdata$dateTime,USGS07010022$flowdata$cms,
     type="l",main=" River Des Peres near University to Mississippi River at St. Louis, MO,Mo ",xlab="Time", ylab="Discharge[m3/s]",xlim=c(USGS07010022$flowdata$dateTime[1109],
                                                                                                                                     USGS07010022$flowdata$dateTime[1109+200]))
lines(USGS07010022$flowdata$outTime,USGS07010022$flowdata$cms,col=2)

