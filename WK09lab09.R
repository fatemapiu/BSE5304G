# Cleaning up
objects()
rm(list=objects())

Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")
options(repos ="http://cran.r-project.org")  # required to get latest libs
if (!require("pacman")) install.packages("pacman")
pacman::p_load(aqp,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,soilDB,circlize,topmodel,DEoptim)
install.packages("EcoHydRology", repos="http://R-Forge.R-project.org")
pacman::p_load(EcoHydRology)
# Lets go back to Lick Run as it has a nice Urban to Forest mix for our P Loss Model
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2022-03-24")

# We want Q in mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

#
# But, we are going to build on last weeks lab separating out 
declat=myflowgage$declat
declon=myflowgage$declon
WXData=FillMissWX(declat, declon,30,
                       date_min="2010-01-01",
                       date_max="2022-03-16")
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")

trunc((180+declon)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)

# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"

# Now we will build our proj4strings which define our “Coordinate 
# Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)
#
# Double check against
area=myflowgage$area   # area in km2
# If the watershed was square, which it is not, the size would be the 
# square root of the area. Also, the gage/pour point is not in the center
# so we will search around the gage.
# Build sp point for USGS gage location, in both _ll and _utm
latlon <- cbind(declon,declat)
gagepoint_ll <- SpatialPoints(latlon)
proj4string(gagepoint_ll)=proj4_ll
gagepoint_utm=spTransform(gagepoint_ll,crs_utm)
# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/place/",declat,",",declon,"/@",
             declat,",",declon,",18z")
browseURL(url)
pourpoint=SpatialPoints(gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=gagepoint_utm@coords
# We are familiar with this basin, and know it is mostly to the NW 
searchlength=sqrt(area*4)*1000 
bboxpts
bboxpts=rbind(bboxpts,bboxpts)
bboxpts[2,1]=bboxpts[2,1]-searchlength
bboxpts[2,2]=bboxpts[2,2]+searchlength
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                        z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydem)
plot(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
Sys.getenv("PATH")
if(!grepl("TauDEM",Sys.getenv("PATH"))){
  old_path <- Sys.getenv("PATH")
  old_path
  Sys.setenv(PATH = paste(old_path,
            paste0(Sys.getenv("HOME"),"/TauDEM/bin"), 
            sep = ":"))
}
# Test to make sure TauDEM runs
system("mpirun aread8")

writeRaster(mydem,filename = "mydem.tif",overwrite=T)
# remember our intro to terminal
# ls; cd ~; pwd;  #Linux/Mac
# dir; cd ; 

z=raster("mydem.tif")
plot(z)

# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel)

plot(z-fel)
# D8 flow directions
system("mpiexec -n 8 d8flowdir -p mydemp.tif -sd8 mydemsd8.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
p=raster("mydemp.tif")
plot(p)
sd8=raster("mydemsd8.tif")
plot(sd8)

# Contributing area
system("mpiexec -n 8 aread8 -p mydemp.tif -ad8 mydemad8.tif")
ad8=raster("mydemad8.tif")
plot(log(ad8))

# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)

# DInf flow directions
system("mpiexec -n 8 dinfflowdir -ang mydemang.tif -slp mydemslp.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
ang=raster("mydemang.tif")
plot(ang)
slp=raster("mydemslp.tif")
plot(slp)

# Dinf contributing area
system("mpiexec -n 8 areadinf -ang mydemang.tif -sca mydemsca.tif")
sca=raster("mydemsca.tif")
plot(log(sca))
#
# This area could be useful in 
threshold=area*10^6/(res(mydem)^2)[1]/10
threshold
# Threshold
system(paste0("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh ",threshold))
src=raster("mydemsrc.tif")
plot(src,xlim=c(pourpoint@coords[1]-5*res(src)[1],pourpoint@coords[1]+5*res(src)[1]),
         ylim=c(pourpoint@coords[2]-5*res(src)[2],pourpoint@coords[2]+5*res(src)[2]))
plot(pourpoint,add=T)

pourpointSPDF=SpatialPointsDataFrame(pourpoint,data.frame(Outlet=row(pourpoint@coords)[1]))
writeOGR(pourpointSPDF, dsn=".", "ApproxOutlets", driver="ESRI Shapefile", overwrite_layer = TRUE)

## a quick R function to write a shapefile
#makeshape.r=function(sname="shape",n=1)
#{
#  xy=locator(n=n)
#  points(xy)
#  
#  #Point
#  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
#  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
#  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
#  write.shapefile(ddShapefile, sname, arcgis=T)
#}
# makeshape.r("ApproxOutlets")

# Move Outlets
system("mpiexec -n 8  moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=readOGR("Outlet.shp")
points(outpt,col=2,add=T)

# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 

# Threshold
system(paste0("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh ",threshold))
src1=raster("mydemsrc1.tif")
plot(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
mydemw=raster("mydemw.tif")
plot(mydemw)
# Plot streams using stream order as width
plot(readOGR("mydemnet.dbf"),add=T)
streams=readOGR("mydemnet.dbf")
# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi,add=T)

mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)

pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
#
# A series of plots to show all of the components
#
par(mfrow = c(2, 2))
plot(TI)
plot(TIC)
plot(mybasinsca)
plot(mybasinslp)
dev.off()
plot(TIC,col=rainbow(5))
# Make a poly with command line gdal (fast)
# 
system("gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp")
mydemw_poly=readOGR("mydemw_poly_gdal.shp")
plot(mydemw_poly,add=T,border="black")
mydemw_poly
writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)
# We will use this ESRI shape file, a zipped group of it, to download 
# our soil extent from the WebSoilSurvey Website
#zip("mydemw.zip",list.files(pattern="mydemw[:.:]"))
# Download to your local machine mydemw.zip from the "Files" tab
# Open the WebSoilSurvey site to: 
#browseURL("https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx")
# "Creat AOI from a zipped shapefile"
# Open "Download Soils Data" Tab
# "Create Download Link" in lower right hand corner
# Right-Click on download link and "Copy Link Address" and 
# paste into a url object
#url=""
#download.file(url,"wss_aoi.zip")
#unzip("wss_aoi.zip")
# Or, if you want to make your life easier, 
mysoil_ll <- mapunit_geom_by_ll_bbox(bbox(spTransform(mydemw_poly,crs_ll)))
crs(mysoil_ll)=crs_ll
mysoil_utm = spTransform(mysoil_ll,crs_utm)
writeOGR(mysoil_utm,dsn=".",layer="mysoil",driver="ESRI Shapefile", overwrite_layer=TRUE)
#m=readOGR("mysoil.shp")
#plot(m)
plot(mysoil_utm,add=T)

# Get the NLCD (USA ONLY)
# Returns a raster
# About time for a break!?!?
pacman::p_load(devtools,terra)
Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")
devtools::install_github("ropensci/FedData")
pacman::p_load(FedData)

NLCD <- get_nlcd(template=TIC, label='TIC',force.redo = TRUE)

# Plot with raster::plot
plot(NLCD)
NLCDutm=projectRaster(NLCD,TIC)
NLCDmask=mask(NLCDutm,mybasinmask)
plot(NLCDmask)


rmysoil_utm=rasterize(mysoil_utm,TIC,field=as.numeric(mysoil_utm$mukey))
unique(rmysoil_utm)
pacman::p_load(circlize)
plot(rmysoil_utm,col=rand_color(length(unique(values(rmysoil_utm)))))
rmysoil_utm=mask(rmysoil_utm,mybasinmask)
plot(rmysoil_utm,col=rand_color(length(unique(values(rmysoil_utm)))))
unique(rmysoil_utm)
#
ratify( TIC*10^9)
unique(ratify(round(mybasinslp*10+1)))
unique(ratify((rmysoil_utm*10^3)))


hru=ratify(TIC*10^9 + (rmysoil_utm*10^3) + round(mybasinslp*10+1))
unique(values(hru))
length(unique(values(hru)))
plot(hru,col=rand_color(length(unique(values(hru)))))


hru_table = levels(hru)[[1]]
origID = hru_table$ID # preserve data order for later reordering
# metadata parameters from a string... this will make more sense
# after the next "head()" command
hru_table$TIclass = as.numeric(substr(sprintf("%10.0f", hru_table$ID), 1,1))
hru_table$mukey = as.numeric(substr(sprintf("%10.0f", hru_table$ID), 2,7))
hru_table$slp = (as.numeric(substr(sprintf("%10.0f", 
                                             hru_table$ID), 8,10))-1)/10
#
# Calculate the area for each unique soil (mukey) X TIClass combination
# using res(raster) for x and y resolution in units of m
# Note that table() function returns the count of the occurrences of
# unique values in the hru raster cells.
hru_table$areaSQKM = as.vector(round(res(hru)[1]*res(hru)[2]*
                                         table(values(hru))/10^6, 3))
#
# To better understand what happened, look at the new hru_table
head(hru_table)
summary(hru_table)

rm("mu2co")
for (mymukey in unique(hru_table$mukey)){
  print(mymukey)
  mukey_statement = format_SQL_in_statement(mymukey)
  q_mu2co = paste("SELECT mukey,cokey FROM component 
           WHERE mukey IN ", mukey_statement, sep="")
  if(!exists("mu2co")){
    mu2co=SDA_query(q_mu2co)} 
  else{
    mu2co=rbind(mu2co,SDA_query(q_mu2co))
  } 
}
View(mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
# cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
# q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r,frag3to10_r  
# FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
# co2ch = SDA_query(q_co2ch)
rm("co2ch")
for (mycokey in unique(mu2co$cokey)){
  print(mycokey)
  cokey_statement = format_SQL_in_statement(mycokey)
  q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r,frag3to10_r FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
  print(q_co2ch)
  if(!exists("co2ch")){
    co2ch=SDA_query(q_co2ch)
  } else{
    try((co2ch=rbind(co2ch,SDA_query(q_co2ch))))
  } 
}
View(co2ch)
rm("co2co")
for (mycokey in unique(mu2co$cokey)){
  print(mycokey)
  cokey_statement = format_SQL_in_statement(mycokey)
  q_co2co = paste("SELECT cokey,slopelenusle_r FROM component WHERE cokey IN ", cokey_statement, sep="")
  print(q_co2co)
  if(!exists("co2co")){
    co2co=SDA_query(q_co2co)} 
  else{
    try((co2co=rbind(co2co,SDA_query(q_co2co))))} 
}
View(co2co)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
mu2ch=merge(mu2ch,co2co)
View(mu2ch)
#
#
#source CNmodel function
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")
pacman::p_load(data.table)
# Use the estimated S for our watershed (Lab06)
Sest = 157
# We will split into 5 VSA areas represented by 5 TI Classes
nTIclass=5
VSAsol=data.table(TIClass=seq(from=nTIclass,to=1),
                    As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
#
VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1
# Calculate TI Class localized sigma and Curve Number
VSAsol[,sigma:=Sest*sSratio]
VSAsol[,CN:=25400/(sigma+254)]
VSAsol

TIC01=modeldata
TIC02=modeldata
TIC03=modeldata
TIC04=modeldata
TIC05=modeldata
# For TIC01 CNavg=VSAParams$CN[1] but confirm
TIC01 = CNmodel(CNmodeldf = TIC01, CNavg=VSAsol$CN[1], 
                declat=declat,declon=declon)
TIC02$P=TIC01$Excess+TIC02$P
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAsol$CN[2], 
                declat=declat,declon=declon)
TIC03$P=TIC02$Excess+TIC03$P
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAsol$CN[3], 
                declat=declat,declon=declon)
TIC04$P=TIC03$Excess+TIC04$P
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAsol$CN[4], 
                declat=declat,declon=declon)
TIC05$P=TIC04$Excess+TIC05$P
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAsol$CN[5], 
                declat=declat,declon=declon)

# Build a dataframe for your DP Load model in TI Class 05,(you will need to do 
#this 5 times for TIC 1-4 where you will need the date, Qpred, and Average 
# Temperature
DPTI05=data.frame(date=TIC05$date,
                  Rt=TIC05$Qpred/1000, 
                  Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
DPTI04=data.frame(date=TIC04$date,
                  Rt=TIC04$Qpred/1000, 
                  Tavg=(TIC04$MaxTemp+TIC04$MinTemp)/2)
DPTI03=data.frame(date=TIC03$date,
                  Rt=TIC03$Qpred/1000, 
                  Tavg=(TIC03$MaxTemp+TIC03$MinTemp)/2)
DPTI02=data.frame(date=TIC02$date,
                  Rt=TIC02$Qpred/1000, 
                  Tavg=(TIC02$MaxTemp+TIC02$MinTemp)/2)
DPTI01=data.frame(date=TIC01$date,
                  Rt=TIC01$Qpred/1000, 
                  Tavg=(TIC01$MaxTemp+TIC01$MinTemp)/2)

# Take it to volumetric water content rather # than volumetric soil content, 
# should be in m of water
tau=9.3  # days
dt=1     # days time step
kF=.015  # Table 2
# Initialize MF and DF
DPTI05$MF=0
DPTI05$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI05$MF[(format(DPTI05$date,"%j") %in% c(121,213,274))]=5.4*10^-4
# Remember what we say about attaching! 
attach(DPTI05)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI05$MF=MF
DPTI05$DF=DF
detach(DPTI05)
rm(list=c("MF","DF")) # Clean up the environment 
dev.off() #reset graphics device
plot(DPTI05$date,DPTI05$MF)
plot(DPTI05$date,DPTI05$DF)


#for TIC4

DPTI04$MF=0
DPTI04$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI04$MF[(format(DPTI04$date,"%j") %in% c(121,213,274))]=5.4*10^-4
# Remember what we say about attaching! 
attach(DPTI04)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI04$MF=MF
DPTI04$DF=DF
detach(DPTI04)
rm(list=c("MF","DF")) # Clean up the environment 
plot(DPTI04$date,DPTI04$MF)
plot(DPTI04$date,DPTI04$DF)

#TIC 3
DPTI03$MF=0
DPTI03$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI03$MF[(format(DPTI03$date,"%j") %in% c(121,213,274))]=5.4*10^-4
# Remember what we say about attaching! 
attach(DPTI03)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI03$MF=MF
DPTI03$DF=DF
detach(DPTI03)
rm(list=c("MF","DF")) # Clean up the environment 
plot(DPTI03$date,DPTI03$MF)
plot(DPTI03$date,DPTI03$DF)

#TIC2
DPTI02$MF=0
DPTI02$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI02$MF[(format(DPTI02$date,"%j") %in% c(121,213,274))]=5.4*10^-4
# Remember what we say about attaching! 
attach(DPTI02)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI02$MF=MF
DPTI02$DF=DF
detach(DPTI02)
rm(list=c("MF","DF")) # Clean up the environment 

plot(DPTI02$date,DPTI02$MF)
plot(DPTI02$date,DPTI02$DF)

#TIC1

DPTI01$MF=0
DPTI01$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI01$MF[(format(DPTI01$date,"%j") %in% c(121,213,274))]=5.4*10^-4
# Remember what we say about attaching! 
attach(DPTI01)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI01$MF=MF
DPTI01$DF=DF
detach(DPTI01)
rm(list=c("MF","DF")) # Clean up the environment 

plot(DPTI01$date,DPTI01$MF)
plot(DPTI01$date,DPTI01$DF)

#from soil

# Calculate your Export Coef from Easton et al 2007 Figure 2 using TI Class 5
# and note that the bold gives the TI Class. This figure gives a range of 
# 0 - 520 micrograms/litre 
# For TIC=5
muTS_TI05=(((520-0)/5)*VSAsol$TIClass[1]+0) # use VSAsol$TIClass table to 
# assign TIC to calc (remember TIC05 is in location 1 in the VSAsol$TIClass 
# table
# Setting range of Soil P values (MS) using the range given on page 7 of 
# Easton et al. 2007 (3.7-18.5mg/kg), then 
# Moore 1993… assume soil density is 2000kg/m^3 of soil
MS_TI05=(((18.5-3.7)/5)*VSAsol$TIClass[1]+3.7)*2000  # range from Easton et 
# al. 2007, pg 7. Moore suggests linear relationship
# We will take care of all of TIClass 05 now as will so 
# it makes sense when you repeat for TI Classes 1-4
# You will use muTS_TI01 thru muTS_TI04 to finish this lab
#
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTI05$muS= muTS_TI05*QS^((DPTI05$Tavg-TR)/10)  # Eq. 5
DPTI05$DS=(DPTI05$muS*MS_TI05*DPTI05$Rt)/10^6          # Eq. 4
plot(DPTI05$date,DPTI05$DS)

#TIC 4
muTS_TI04=(((520-0)/5)*VSAsol$TIClass[2]+0) # use VSAsol$TIClass table to 
MS_TI04=(((18.5-3.7)/5)*VSAsol$TIClass[2]+3.7)*2000  # range from Easton et 
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTI04$muS= muTS_TI04*QS^((DPTI04$Tavg-TR)/10)  # Eq. 5
DPTI04$DS=(DPTI04$muS*MS_TI04*DPTI04$Rt)/10^6          # Eq. 4
plot(DPTI04$date,DPTI04$DS)

#TIC 3
muTS_TI03=(((520-0)/5)*VSAsol$TIClass[3]+0) # use VSAsol$TIClass table to 
MS_TI03=(((18.5-3.7)/5)*VSAsol$TIClass[3]+3.7)*2000  # range from Easton et 
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTI03$muS= muTS_TI03*QS^((DPTI03$Tavg-TR)/10)  # Eq. 5
DPTI03$DS=(DPTI03$muS*MS_TI03*DPTI03$Rt)/10^6          # Eq. 4
plot(DPTI03$date,DPTI03$DS)

#TIC 2
muTS_TI02=(((520-0)/5)*VSAsol$TIClass[4]+0) # use VSAsol$TIClass table to 
MS_TI02=(((18.5-3.7)/5)*VSAsol$TIClass[4]+3.7)*2000  # range from Easton et 
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTI02$muS= muTS_TI02*QS^((DPTI02$Tavg-TR)/10)  # Eq. 5
DPTI02$DS=(DPTI02$muS*MS_TI02*DPTI02$Rt)/10^6          # Eq. 4
plot(DPTI02$date,DPTI02$DS)

#TIC 1
muTS_TI01=(((520-0)/5)*VSAsol$TIClass[5]+0) # use VSAsol$TIClass table to 
MS_TI01=(((18.5-3.7)/5)*VSAsol$TIClass[5]+3.7)*2000  # range from Easton et 
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTI01$muS= muTS_TI01*QS^((DPTI01$Tavg-TR)/10)  # Eq. 5
DPTI01$DS=(DPTI01$muS*MS_TI01*DPTI01$Rt)/10^6          # Eq. 4
plot(DPTI01$date,DPTI01$DS)

# We will assume a simple base flow model, where the stream baseflow 
# B equals the minimum measured flow of the basin. We will build the 
# Base Flow solution directly into the Stream Load dataframe (DPLT).
# Note that we only consider TIC05 here because TIC05 is the class that 
# directly contributes to the stream.
DPLT=data.frame(date=TIC05$date,
                Rt=TIC05$Qpred/1000.0,
                Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
plot(DPLT$date,DPLT$LB)
summary(DPLT)

TIC05$area=myflowgage$area/5
TIC04$area=myflowgage$area/5
TIC03$area=myflowgage$area/5
TIC02$area=myflowgage$area/5
TIC01$area=myflowgage$area/5

DPLT$LT=DPLT$LB +
  TIC05$area*(DPTI05$DF + DPTI05$DS)+
  TIC04$area*(DPTI04$DF + DPTI04$DS)+
  TIC03$area*(DPTI03$DF + DPTI03$DS)+
  TIC02$area*(DPTI02$DF + DPTI02$DS)+
  TIC01$area*(DPTI01$DF + DPTI01$DS)
plot(DPLT$date,DPLT$LT, type="l",ylim=c(1,10))


#.....................................................................
#....................................question1.......................
a=mean(c(18.5,3.7))/10^6*2000*.1*.2
DPLT$total_DS=DPTI01$DS+DPTI02$DS+DPTI03$DS+DPTI04$DS+DPTI05$DS
#d=mean(DPTI01$DS)+mean(DPTI02$DS)+mean(DPTI03$DS)+mean(DPTI04$DS)+mean(DPTI05$DS)
t1=a/(mean(DPLT$total_DS))
print(t1)

#........................................question2.......................

#increasing temperature

#P from soil will be affected

# For TIC=5
DPTI05$muS_5= muTS_TI05*QS^((DPTI05$Tavg+5-TR)/10)  # Eq. 5
DPTI05$DS_5=(DPTI05$muS_5*MS_TI05*DPTI05$Rt)/10^6          # Eq. 4
plot(DPTI05$date,DPTI05$DS_5)

#TIC 4
DPTI04$muS_5= muTS_TI04*QS^((DPTI04$Tavg+5-TR)/10)  # Eq. 5
DPTI04$DS_5=(DPTI04$muS_5*MS_TI04*DPTI04$Rt)/10^6          # Eq. 4
plot(DPTI04$date,DPTI04$DS_5)

#TIC 3
DPTI03$muS_5= muTS_TI03*QS^((DPTI03$Tavg+5-TR)/10)  # Eq. 5
DPTI03$DS_5=(DPTI03$muS_5*MS_TI03*DPTI03$Rt)/10^6          # Eq. 4
plot(DPTI03$date,DPTI03$DS_5)

#TIC 2
DPTI02$muS_5= muTS_TI02*QS^((DPTI02$Tavg+5-TR)/10)  # Eq. 5
DPTI02$DS_5=(DPTI02$muS_5*MS_TI02*DPTI02$Rt)/10^6          # Eq. 4
plot(DPTI02$date,DPTI02$DS_5)

#TIC 1
DPTI01$muS_5= muTS_TI01*QS^((DPTI01$Tavg+5-TR)/10)  # Eq. 5
DPTI01$DS_5=(DPTI01$muS_5*MS_TI01*DPTI01$Rt)/10^6          # Eq. 4
plot(DPTI01$date,DPTI01$DS_5)

DPLT$muB_5=muTB*QB^((DPLT$Tavg+5-TB)/10)  # Easton eq. 10
DPLT$LB_5=DPLT$muB_5*DPLT$B 

DPLT$LT_5=DPLT$LB_5 +
  TIC05$area*(DPTI05$DF + DPTI05$DS_5)+
  TIC04$area*(DPTI04$DF + DPTI04$DS_5)+
  TIC03$area*(DPTI03$DF + DPTI03$DS_5)+
  TIC02$area*(DPTI02$DF + DPTI02$DS_5)+
  TIC01$area*(DPTI01$DF + DPTI01$DS_5)
plot(DPLT$date,DPLT$LT_5, type="l",ylim=c(1,10),add=T,col='red')

summary(DPLT)

pacman::p_load(lattice,ggplot2,dplyr,patchwork,hrbrthemes,tidyverse)

colors <- c("Before Temperature change" = "darkgreen", "After Temperature Change" = "red")
p1<-
  ggplot(data=DPLT, aes(x=date) )+
  geom_line( aes(y=LT_5,color='After Temperature Change'), size=.5)+
  geom_line(aes(y=LT, color="Before Temperature change"), size=.5)+
  
  xlab('Year')+
  ylab('LT[kg/day]')+
  theme(
    axis.title.y = element_text(color = 'black', size=12),
    axis.title.y.right = element_text(color = 'black', size=12),
    plot.title = element_text(hjust = 0.5,size=14),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.position ='top',
    legend.key = element_rect(fill = NA),
    legend.title = element_text(size=10),
    legend.text = element_text(size=10)
  ) +
  ggtitle("P export at Lick Run, VA")+
  scale_color_manual('Parameters:',values = colors)
p1

#..............................question 3......................
#....................................................................
tau=2  # days
dt=1     # days time step
kF=.015  # Table 2F
DPTI05$MF_2=0
DPTI05$DF_2=0
DPTI05$MF_2[(format(DPTI05$date,"%j") %in% c(121,213,274))]=5.4*10^-4
attach(DPTI05)
#
for (i in 2:length(date)){
  if(MF_2[i]<=MF_2[i-1]){
    MF_2[i]=MF_2[i-1]*exp(-dt/tau)-DF_2[i-1]
  }
  DF_2[i]=MF_2[i]*(kF*MF_2[i]*Rt[i]/(1+kF*MF_2[i]*Rt[i]))
}
DPTI05$MF_2=MF_2
DPTI05$DF_2=DF_2
detach(DPTI05)
rm(list=c("MF_2","DF_2")) # Clean up the environment 

#for TIC4
DPTI04$MF_2=0
DPTI04$DF_2=0
DPTI04$MF_2[(format(DPTI04$date,"%j") %in% c(121,213,274))]=5.4*10^-4
# Remember what we say about attaching! 
attach(DPTI04)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF_2[i]<=MF_2[i-1]){
    MF_2[i]=MF_2[i-1]*exp(-dt/tau)-DF_2[i-1]
  }
  DF_2[i]=MF_2[i]*(kF*MF_2[i]*Rt[i]/(1+kF*MF_2[i]*Rt[i]))
}
DPTI04$MF_2=MF_2
DPTI04$DF_2=DF_2
detach(DPTI04)
rm(list=c("MF_2","DF_2")) # Clean up the environment 


#TIC 3
DPTI03$MF_2=0
DPTI03$DF_2=0

DPTI03$MF_2[(format(DPTI03$date,"%j") %in% c(121,213,274))]=5.4*10^-4

attach(DPTI03)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF_2[i]<=MF_2[i-1]){
    MF_2[i]=MF_2[i-1]*exp(-dt/tau)-DF_2[i-1]
  }
  DF_2[i]=MF_2[i]*(kF*MF_2[i]*Rt[i]/(1+kF*MF_2[i]*Rt[i]))
}
DPTI03$MF_2=MF_2
DPTI03$DF_2=DF_2
detach(DPTI03)



#TIC2
DPTI02$MF_2=0
DPTI02$DF_2=0
DPTI02$MF_2[(format(DPTI02$date,"%j") %in% c(121,213,274))]=5.4*10^-4
attach(DPTI02)
#
for (i in 2:length(date)){
  if(MF_2[i]<=MF_2[i-1]){
    MF_2[i]=MF_2[i-1]*exp(-dt/tau)-DF_2[i-1]
  }
  DF_2[i]=MF_2[i]*(kF*MF_2[i]*Rt[i]/(1+kF*MF_2[i]*Rt[i]))
}
DPTI02$MF_2=MF_2
DPTI02$DF_2=DF_2
detach(DPTI02)



#TIC1

DPTI01$MF_2=0
DPTI01$DF_2=0
DPTI01$MF_2[(format(DPTI01$date,"%j") %in% c(121,213,274))]=5.4*10^-4
attach(DPTI01)
for (i in 2:length(date)){
  if(MF_2[i]<=MF_2[i-1]){
    MF_2[i]=MF_2[i-1]*exp(-dt/tau)-DF_2[i-1]
  }
  DF_2[i]=MF_2[i]*(kF*MF_2[i]*Rt[i]/(1+kF*MF_2[i]*Rt[i]))
}
DPTI01$MF_2=MF_2
DPTI01$DF_2=DF_2
detach(DPTI01)
rm(list=c("MF_2","DF_2")) # Clean up the environment 

DPLT$LT_t2=DPLT$LB +
  TIC05$area*(DPTI05$DF_2 + DPTI05$DS)+
  TIC04$area*(DPTI04$DF_2 + DPTI04$DS)+
  TIC03$area*(DPTI03$DF_2 + DPTI03$DS)+
  TIC02$area*(DPTI02$DF_2 + DPTI02$DS)+
  TIC01$area*(DPTI01$DF_2 + DPTI01$DS)
summary(DPLT)

colors <- c("Tau=2 days" = "darkgreen", "Tau=9.3days" = "red")
p2<-
  ggplot(data=DPLT, aes(x=date) )+
  geom_line( aes(y=LT_t2,color='Tau=2 days'), size=.5)+
  geom_line(aes(y=LT, color="Tau=9.3days"), size=.5)+
  
  xlab('Year')+
  ylab('LT[kg/day]')+
  theme(
    axis.title.y = element_text(color = 'black', size=12),
    axis.title.y.right = element_text(color = 'black', size=12),
    plot.title = element_text(hjust = 0.5,size=14),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.position ='top',
    legend.key = element_rect(fill = NA),
    legend.title = element_text(size=10),
    legend.text = element_text(size=10)
  ) +
  ggtitle("P export at Lick Run, VA")+
  scale_color_manual('Parameters:',values = colors)
p2

colors <- c("Tau=2 days" = "darkgreen", "Tau=9.3days" = "red")

#.......................................................
#.......................grad.............................
p_loss<-function(RunoffModel=TIC05,tau=9.3, dt=1, kF=.015, TI=1){
  
  DPTI05=data.frame(date=TIC05$date,
                    Rt=TIC05$Qpred/1000, 
                    Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
  

  DPTI05$MF=0
  DPTI05$DF=0
  DPTI05$MF[(format(DPTI05$date,"%j") %in% c(121,213,274))]=5.4*10^-4
  attach(DPTI05)
  for (i in 2:length(date)){
    if(MF[i]<=MF[i-1]){
      MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
    }
    DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
  }
  DPTI05$MF=MF
  DPTI05$DF=DF
  detach(DPTI05)
  rm(list=c("MF","DF")) # Clean up the environment 
  
  muTS_TI05=(((520-0)/5)*VSAsol$TIClass[6-TI]+0) # use VSAsol$TIClass table to 
  
  MS_TI05=(((18.5-3.7)/5)*VSAsol$TIClass[6-TI]+3.7)*2000  # range from Easton et 
 
  QS= 3.0 # A guess using the middle of the range 1-5
  TR=20   # reference Temperature from Table 2.
  DPTI05$muS= muTS_TI05*QS^((DPTI05$Tavg-TR)/10)  # Eq. 5
  DPTI05$DS=(DPTI05$muS*MS_TI05*DPTI05$Rt)/10^6          # Eq. 4
  
  DPTI05$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
  muTB=2.1*10^(-5) # Easton Table 2
  QB=2.2           # Easton Table 2
  TB=17            # Easton Table 2
  DPTI05$muB=muTB*QB^((DPTI05$Tavg-TB)/10)  # Easton eq. 10
  DPTI05$LB=DPTI05$muB*DPTI05$B     # Easton eq. 9
  TIC05$area=myflowgage$area/5
  TIC04$area=myflowgage$area/5
  TIC03$area=myflowgage$area/5
  TIC02$area=myflowgage$area/5
  TIC01$area=myflowgage$area/5
  DPLT$LT=DPLT$LB +
    TIC05$area*(DPTI05$DF + DPTI05$DS)+
    TIC04$area*(DPTI04$DF + DPTI04$DS)+
    TIC03$area*(DPTI03$DF + DPTI03$DS)+
    TIC02$area*(DPTI02$DF + DPTI02$DS)+
    TIC01$area*(DPTI01$DF + DPTI01$DS)
  return(DPTI05)
  
}

TIC01=p_loss(RunoffModel=TIC01,tau=9.3, dt=1, kF=.015, TI=1)
TIC02=p_loss(RunoffModel=TIC02,tau=9.3, dt=1, kF=.015, TI=2)
TIC03=p_loss(RunoffModel=TIC03,tau=9.3, dt=1, kF=.015, TI=3)
TIC04=p_loss(RunoffModel=TIC04,tau=9.3, dt=1, kF=.015, TI=4)
TIC05=p_loss(RunoffModel=TIC05,tau=9.3, dt=1, kF=.015, TI=5)


l1=mean(TIC01$area*(DPTI01$DF + DPTI01$DS))*365
l2=mean(TIC02$area*(DPTI02$DF + DPTI02$DS))*365
l3=mean(TIC03$area*(DPTI03$DF + DPTI03$DS))*365
l4=mean(TIC04$area*(DPTI04$DF + DPTI04$DS))*365
l5=mean(TIC05$area*(DPTI05$DF + DPTI05$DS))*365

L <- data.frame(b1 = c(l1,l2,l3,l4,l5))
library(raster)
layer<- raster(nrow=5, ncol=1, crs="+proj=utm +zone=17 +datum=WGS84 +units=m")
values(layer) <- L$b1
new_TIC <- reclassify(TIC, layer)
par(mar = c(1.5, 1.5, 1.5, 1.5),xaxs = "i")
plot(new_TIC,xlim=c(585000,596000),main='P_export')
