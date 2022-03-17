#....................................................................................
#.......................................HW2..........................................
#.....................................................................................
# Cleaning up
objects()
rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,classInt)
# Getting standard flow data from USGS 0205551460 `
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2019-01-01")
# We want Q in mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

# We can determine our UTM zone and build our proj4 projection string
trunc((180+myflowgage$declon)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)

# Lat/Lon (_ll)
proj4_ll = "+proj=longlat"
# Now we will build our proj4strings which define our “Coordinate
# Reference Systems” or CRS in future geographic manipulations.
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)
myflowgage$area # area in km2

latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)
proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)
# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/@",
             myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)
# We are going to over estimate our area
sqrt(myflowgage$area) # guestimating square watershed
# For our search we are going to multiply the area by 6 and
# to get the distance
sqrt(myflowgage$area*8)
searchlength=sqrt(myflowgage$area*8)*1000
pourpoint=SpatialPoints(myflowgage$gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=myflowgage$gagepoint_utm@coords
bboxpts=rbind(bboxpts,bboxpts+searchlength)
bboxpts=rbind(bboxpts,bboxpts-searchlength)
bboxpts
bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
bboxpts
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# my DEM
mydem=get_aws_terrain(locations=bboxpts@coords,
                        z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydem)
par(mar = c(1.5, 1.5, 1.5, 1.5),xaxs = "i")
plot(mydem,xlim=c(580000,630000))
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
# Write our raster to a geotiff file that can be used with
# OS level hydrological models
writeRaster(mydem,filename = "mydem.tif",overwrite=T)
# remember our intro to terminal
# ls; cd ~; pwd; #Linux/Mac
# dir; cd ;
#zoom(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path,
                        paste0(Sys.getenv("HOME"),"/TauDEM/bin"),
                        sep = ":"))
system("mpirun aread8")


writeRaster(mydem,filename = "mydem.tif",overwrite=T)
z=raster("mydem.tif")
plot(z)

# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel)


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
#zoom(log(ad8))


# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)
#zoom(gord)

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
#zoom(log(sca))

# Threshold
system("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh 2000")
src=raster("mydemsrc.tif")
plot(src)
plot(pourpoint,add=T,col="red")
zoom(src)

# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

makeshape.r("ApproxOutlets")

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")
approxpt=read.shp("ApproxOutlets.shp")

plot(src)
points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)


# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa,xlim=c(580000,600000)) 


# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 4000")
src1=raster("mydemsrc1.tif")
plot(src1)
#zoom(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
plot(raster("mydemord.tif"))
#zoom(raster("mydemord.tif"))
plot(raster("mydemw.tif"))

# Plot streams using stream order as width
snet=read.shapefile("mydemnet")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Peuker Douglas stream definition
system("mpiexec -n 8 peukerdouglas -fel mydemfel.tif -ss mydemss.tif")
ss=raster("mydemss.tif")
plot(ss)
#zoom(ss)

#  Accumulating candidate stream source cells
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif -wg mydemss.tif")
ssa=raster("mydemssa.tif")
plot(ssa)

#  Drop Analysis
system("mpiexec -n 8 dropanalysis -p mydemp.tif -fel mydemfel.tif -ad8 mydemad8.tif -ssa mydemssa.tif -drp mydemdrp.txt -o Outlet.shp -par 5 500 10 0")

# Deduce that the optimal threshold is 300 
# Stream raster by threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc2.tif -thresh 4000")
plot(raster("mydemsrc2.tif"))

# Stream network
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc2.tif -ord mydemord2.tif -tree mydemtree2.dat -coord mydemcoord2.dat -net mydemnet2.shp -w mydemw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)

plot(raster("mydemw2.tif"))
#zoom(raster("mydemw2.tif"))
snet=read.shapefile("mydemnet2")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi)

# Distance Down
system("mpiexec -n 8 dinfdistdown -ang mydemang.tif -fel mydemfel.tif -src mydemsrc2.tif -m ave v -dd mydemdd.tif",show.output.on.console=F,invisible=F)
plot(raster("mydemdd.tif"))


#Trimming, Cropping, and Masking to make life prettier and easier
mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem,add=T)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)
mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)
TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
#zoom(TI)
pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method
#various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)

#Trimming, Cropping, and Masking to make life prettier and easier
mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(z,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem)

mybasin_fel=crop(fel,mybasinmask)
mybasin_fel=mask(mybasin_fel,mybasinmask)

#plotting Filled DEM, Difference between DEM and Filled DEM, TI, and TI Classes
par(mar=c(2.5, 2.5,2.5,2.5))
par(mfrow = c(2, 2))
plot(mybasin_fel,main='Filled DEM')

plot(mybasindem-mybasin_fel,main='DEM - Filled DEM')
plot(TI,main='TI')
plot(TIC,main='TIC')

#.................................................................................
# Cleaning up
objects()
rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,classInt)
# Getting standard flow data from USGS 01483700 `
myflowgage_id="01483700"
myflowgage=get_usgs_gage(myflowgage_id,begin_date="2019-01-01",end_date="2021-02-01")
# We want Q in mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

# We can determine our UTM zone and build our proj4 projection string
trunc((180+myflowgage$declon)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)

# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
# Now we will build our proj4strings which define our “Coordinate
# Reference Systems” or CRS in future geographic manipulations.
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)
myflowgage$area # area in km2

latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)
proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)
# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/@",
           myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)
# We are going to over estimate our area
sqrt(myflowgage$area) # guestimating square watershed
# For our search we are going to multiply the area by 6 and
# to get the distance
sqrt(myflowgage$area*4)
searchlength=sqrt(myflowgage$area*4)*1000
pourpoint=SpatialPoints(myflowgage$gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=myflowgage$gagepoint_utm@coords
bboxpts=rbind(bboxpts,bboxpts+searchlength)
bboxpts=rbind(bboxpts,bboxpts-searchlength)
bboxpts
bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
bboxpts
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# my DEM
mydem=get_aws_terrain(locations=bboxpts@coords,
                      z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydem)
#par(mar = c(1.5, 1.5, 1.5, 1.5))
plot(mydem,xlim=c(580000,630000))
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
# Write our raster to a geotiff file that can be used with
# OS level hydrological models
writeRaster(mydem,filename = "mydem.tif",overwrite=T)
# remember our intro to terminal
# ls; cd ~; pwd; #Linux/Mac
# dir; cd ;
#zoom(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path,
                        paste0(Sys.getenv("HOME"),"/TauDEM/bin"),
                        sep = ":"))
system("mpirun aread8")


writeRaster(mydem,filename = "mydem.tif",overwrite=T)
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
#zoom(log(ad8))


# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)
#zoom(gord)

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
#zoom(log(sca))

# Threshold
system("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh 6000")
src=raster("mydemsrc.tif")
plot(src)
plot(pourpoint,add=T,col="red")
zoom(src)

# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

makeshape.r("ApproxOutlets")

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")
approxpt=read.shp("ApproxOutlets.shp")

plot(src)
points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)

#zoom(src)


# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 


# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 6000")
src1=raster("mydemsrc1.tif")
plot(src1)
zoom(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
plot(raster("mydemord.tif"))
#zoom(raster("mydemord.tif"))
plot(raster("mydemw.tif"))

# Plot streams using stream order as width
snet=read.shapefile("mydemnet")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Peuker Douglas stream definition
system("mpiexec -n 8 peukerdouglas -fel mydemfel.tif -ss mydemss.tif")
ss=raster("mydemss.tif")
plot(ss)
#zoom(ss)

#  Accumulating candidate stream source cells
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif -wg mydemss.tif")
ssa=raster("mydemssa.tif")
plot(ssa)

#  Drop Analysis
system("mpiexec -n 8 dropanalysis -p mydemp.tif -fel mydemfel.tif -ad8 mydemad8.tif -ssa mydemssa.tif -drp mydemdrp.txt -o Outlet.shp -par 5 500 10 0")

# Deduce that the optimal threshold is 300 
# Stream raster by threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc2.tif -thresh 4000")
plot(raster("mydemsrc2.tif"))

# Stream network
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc2.tif -ord mydemord2.tif -tree mydemtree2.dat -coord mydemcoord2.dat -net mydemnet2.shp -w mydemw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)

plot(raster("mydemw2.tif"))
#zoom(raster("mydemw2.tif"))
snet=read.shapefile("mydemnet2")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi)

# Distance Down
system("mpiexec -n 8 dinfdistdown -ang mydemang.tif -fel mydemfel.tif -src mydemsrc2.tif -m ave v -dd mydemdd.tif",show.output.on.console=F,invisible=F)
plot(raster("mydemdd.tif"))


#Trimming, Cropping, and Masking to make life prettier and easier
mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
#plot(mybasindem,add=T)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)
mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)
TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
#zoom(TI)
pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method
#various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)

par(mfrow = c(2, 2))
plot(TI)
plot(TIC)
plot(mybasinsca)
plot(mybasinslp)
dev.off()
plot(TIC,col=rainbow(5))


#Trimming, Cropping, and Masking to make life prettier and easier
mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(z,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem,add=T)

mybasin_fel=crop(fel,mybasinmask)
mybasin_fel=mask(mybasin_fel,mybasinmask)

#plotting Filled DEM, Difference between DEM and Filled DEM, TI, and TI Classes
par(mar=c(2.5, 2.5,2.5,2.5))
par(mfrow = c(2, 2))
plot(mybasin_fel,main='Filled DEM')

plot(mybasindem-mybasin_fel,main='DEM - Filled DEM')
plot(TI,main='TI')
plot(TIC,main='TIC')
#..........................................................................................
#.........................................HW3.............................................
#..........................................................................................
#soil data
# Make a poly with raster library (slow)
# or from thee command line gdal (fast)
# gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp
mydemw_poly=rasterToPolygons(mydemw,dissolve = T,na.rm = T)
plot(mydemw_poly,add=T,border="red")
mydemw_poly
writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile",overwrite_layer=TRUE)
zip("mydemw.zip",list.files(pattern="mydemw[:.:]"))
browseURL("https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx")
url='https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/ot2r32ewxdxl0mhqlhan5rgo/wss_aoi_2022-03-17_13-51-35.zip'
download.file(url,"mysoil.zip")
unzip("mysoil.zip")

mysoil=readOGR("wss_aoi_2022-03-17_13-51-35/spatial/soilmu_a_aoi.shp")    
# Exploring the mysoil dataset
mybbox=c(mysoil@bbox)
# First associating mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  #renaming the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associating cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
# Last, bringing them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)
aws<-merge(mysoil,mu2chmax,by="mukey")
#nlcd
nlcd=raster("NLCD.tiff")
nlcd=projectRaster(nlcd,crs=crs_utm)
mybasin_nl=crop(nlcd,mybasinmask)
mybasin_nl=mask(mybasin_nl,mybasinmask)
nl=resample(nlcd,mybasinmask)
mybasin_nl=mask(nl,mybasinmask)
plot(mybasin_nl)

#plotting Filled DEM, Difference between DEM and Filled DEM, TI, and TI Classes
par(mar=c(1.5, 1.5,1.5,1.5))
plot(mybasindem,main='DEM')
plot(TIC,main='TIC')
plot(mybasin_nl,main='Land Cover')
spplot(aws, zcol="awc_r",main="AWC",
       colorkey=T, col.regions=topo.colors(100), checkEmptyRC=T,add=T)

#.........................................................................................
#.........................................HW4.............................................
#...........................................................................................
# Cleaning up
objects()
rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,classInt)
# Getting standard flow data from USGS 01102345 `
#myflowgage_id="01483700"
#myflowgage=get_usgs_gage(myflowgage_id,begin_date="2007-01-01",end_date="2008-02-01")
# We want Q in mm/day for the basin
#myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
lon=-80.447565
lat=37.20767
# We can determine our UTM zone and build our proj4 projection string
trunc((180+lon)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+lon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)

# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
# Now we will build our proj4strings which define our “Coordinate
# Reference Systems” or CRS in future geographic manipulations.
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)
#myflowgage$area # area in km2

latlon <- cbind(lon,lat)
gagepoint_ll <- SpatialPoints(latlon)
proj4string(gagepoint_ll)=proj4_ll
gagepoint_utm=spTransform(gagepoint_ll,crs_utm)
# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/@",
           lat,",",lon,",18z")
browseURL(url)
# We are going to over estimate our area
#sqrt(myflowgage$area) # guestimating square watershed
# For our search we are going to multiply the area by 6 and
# to get the distance
sqrt(6*6)
searchlength=sqrt(6*6)*1000
pourpoint=SpatialPoints(gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=gagepoint_utm@coords
bboxpts=rbind(bboxpts,bboxpts+searchlength)
bboxpts=rbind(bboxpts,bboxpts-searchlength)
bboxpts
bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
bboxpts
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# my DEM
mydem=get_aws_terrain(locations=bboxpts@coords,
                      z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydem)
#par(mar = c(1.5, 1.5, 1.5, 1.5),xaxs = "i")
plot(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
# Write our raster to a geotiff file that can be used with
# OS level hydrological models
writeRaster(mydem,filename = "mydem.tif",overwrite=T)
# remember our intro to terminal
# ls; cd ~; pwd; #Linux/Mac
# dir; cd ;
#zoom(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path,
                        paste0(Sys.getenv("HOME"),"/TauDEM/bin"),
                        sep = ":"))
system("mpirun aread8")


writeRaster(mydem,filename = "mydem.tif",overwrite=T)
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
#zoom(log(ad8))


# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)
#zoom(gord)

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
#zoom(log(sca))

# Threshold
system("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh 5000")
src=raster("mydemsrc.tif")
plot(src)
plot(pourpoint,add=T,col="red")
zoom(src)

# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

makeshape.r("ApproxOutlets")

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")
approxpt=read.shp("ApproxOutlets.shp")

plot(src)
points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)

#zoom(src)


# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 


# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 5000")
src1=raster("mydemsrc1.tif")
plot(src1)
zoom(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
plot(raster("mydemord.tif"))
#zoom(raster("mydemord.tif"))
plot(raster("mydemw.tif"))

# Plot streams using stream order as width
snet=read.shapefile("mydemnet")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Peuker Douglas stream definition
system("mpiexec -n 8 peukerdouglas -fel mydemfel.tif -ss mydemss.tif")
ss=raster("mydemss.tif")
plot(ss)
#zoom(ss)

#  Accumulating candidate stream source cells
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif -wg mydemss.tif")
ssa=raster("mydemssa.tif")
plot(ssa)

#  Drop Analysis
system("mpiexec -n 8 dropanalysis -p mydemp.tif -fel mydemfel.tif -ad8 mydemad8.tif -ssa mydemssa.tif -drp mydemdrp.txt -o Outlet.shp -par 5 500 10 0")

# Deduce that the optimal threshold is 300 
# Stream raster by threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc2.tif -thresh 4000")
plot(raster("mydemsrc2.tif"))

# Stream network
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc2.tif -ord mydemord2.tif -tree mydemtree2.dat -coord mydemcoord2.dat -net mydemnet2.shp -w mydemw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)

plot(raster("mydemw2.tif"))
#zoom(raster("mydemw2.tif"))
snet=read.shapefile("mydemnet2")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi)

# Distance Down
system("mpiexec -n 8 dinfdistdown -ang mydemang.tif -fel mydemfel.tif -src mydemsrc2.tif -m ave v -dd mydemdd.tif",show.output.on.console=F,invisible=F)
plot(raster("mydemdd.tif"))


#Trimming, Cropping, and Masking to make life prettier and easier
mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem,add=T)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)
mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)
TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
#zoom(TI)
pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method
#various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)

#Trimming, Cropping, and Masking to make life prettier and easier
mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(z,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem,add=T)

mybasin_fel=crop(fel,mybasinmask)
mybasin_fel=mask(mybasin_fel,mybasinmask)


#plotting Filled DEM, Difference between DEM and Filled DEM, TI, and TI Classes
par(mar=c(2.5, 2.5,2.5,2.5))
par(mfrow = c(2, 2))
plot(mybasin_fel,main='Filled DEM')

plot(mybasindem-mybasin_fel,main='DEM - Filled DEM')
plot(TI,main='TI')
plot(TIC,main='TIC')




