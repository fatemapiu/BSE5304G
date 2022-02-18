#cleaning environment
old=c(1,2,3,4,5) # Numeric vector
c=6:10 # Numeric vector
iwanna=function(a,b){print(old*c)} # a simple function
remove=iwanna(old,c)
objects() # list of objects that I have
rm(list=objects()) # Removes ALL the objects
objects() # environment is cleaned now

#lab04 starts
options(repos ="http://cran.us.r-project.org")  # required to get latest library
# Installing the packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster)
pacman::p_load(EcoHydRology,rnoaa,curl,httr)
# 3 Functions to calculate SWE and excess when soil is drying, 
#wetting, and wetting above capacity

browseURL("https://github.com/vtdrfuka/BSE5304_2022/tree/main/functions")
browseURL("https://github.com/vtdrfuka/BSE5304_2022/blob/main/functions/TMWBmodel.R")
browseURL("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/NSE.R")
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")

# Downloading a soils dataset for my basin based on the WebSoilSurvey method 
url='https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/1fzi4r2xg2q1x3sq4zfgk0xs/wss_aoi_2022-02-17_08-38-07.zip'
download.file(url,"mysoil.zip")
unzip("mysoil.zip")

# using the function to get data from USGS 01102345
#SAUGUS RIVER AT SAUGUS IRONWORKS AT SAUGUS, MA
myflowgage_id="01102345"
myflowgage=get_usgs_gage(myflowgage_id,begin_date="2000-1-1",end_date="2022-02-01")
# converting flow  m3/day to mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
# the soil data by using the soilDB package
mysoil=readOGR("wss_aoi_2022-02-17_08-38-07/spatial/soilmu_a_aoi.shp")    
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

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.001)

summary(terrain(mydem, opt='slope',unit = "radians"))

#gathering weather data
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius =40,
  limit = NULL
)
WXData=meteo_pull_monitors(
  monitors=stns[201,1],
  keep_flags = FALSE,
  date_min = "2000-1-1",
  date_max ="2022-02-01",
  var = c("TMAX","TMIN","PRCP")
)
summary(WXData)
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
summary(modeldata)  #
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
# Compare your precipitation to the flow out of your basin
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MaxTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)]=0 #removal of NAs
summary(modeldata)
TMWB=modeldata
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")

TopSlope=TMWB
MidSlope=TMWB
BotSlope=TMWB  

top=TMWBmodel(TMWB=TopSlope,Slope=0.07093194)
MidSlope$P=top$Excess+TMWB$P
mid=TMWBmodel(TMWB=MidSlope,Slope =0.14762225)
BotSlope$P=mid$Excess+TMWB$P
bot=TMWBmodel(TMWB=BotSlope,Slope=0)
top1= top[top$date >"2021-01-15" & top$date < "2021-04-30", ]
mid1= mid[mid$date >"2021-01-15" & mid$date < "2021-04-30", ]
bot1= bot[bot$date >"2021-01-15" & bot$date < "2021-04-30", ]


colors <- c("Top" = "red", "Mid" = "blue","Bottom"="darkgreen")

ggplot()+
  
  geom_line(data=top1, aes(x=date,y=Excess,color='Top'), size=.5) + 
  geom_line(data=mid1, aes(x=date,y=Excess,color='Mid'), size=.5)+
  geom_line(data=bot1, aes(x=date,y=Excess,color='Bottom'), size=.5)+
  
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
  xlab('Year')+
  ylab('Excess[mm]')+
  ggtitle('Excess at SAUGUS RIVER,MO')+
  scale_color_manual('Parameters:',values = colors)



ggplot()+
  
  geom_line(data=top1, aes(x=date,y=AW,color='Top'), size=.5) + 
  geom_line(data=mid1, aes(x=date,y=AW,color='Mid'), size=.5)+
  geom_line(data=bot1, aes(x=date,y=AW,color='Bottom'), size=.5)+
  
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
  xlab('Year')+
  ylab('AW[mm]')+
  ggtitle('AW at SAUGUS RIVER,MO')+
  scale_color_manual('Parameters:',values = colors)


#.........................................................
#graduate homework problem B

# Installing the packages we will play with today
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster,aqp,sp,rgeos,wk,ggplot2,gridExtra,viridisLite,knitr)

p <- fetchKSSL(bbox=c(-71.500, 42, -71, 42.5))
coordinates(p) <- ~ x + y
proj4string(p) <- '+proj=longlat +datum=WGS84'
p.sp <- as(p, 'SpatialPointsDataFrame')

# perform SDA query on collection of points
mu1.data <- SDA_spatialQuery(p.sp, what = 'geom')

# use local spatial intersection to link source + SDA data
p.sp$mukey <- over(p.sp, mu1.data)$mukey

# join results to original SoilProfileCollection using 'pedlabsampnum'
site(p) <- as.data.frame(p.sp)[, c('pedlabsampnum', 'mukey')]
par(mar=c(0,2,5,2))
groupedProfilePlot(p, groups='mukey', group.name.cex=0.65, color='clay', name='hzn_desgn', id.style='side', label='pedon_id', max.depth=200)
# describe IDs
mtext('user pedon ID', side=2, line=-1.5)
mtext('mukey', side=3, line=-1, at = c(0,0), adj = -0.8)
dev.off()
#sand
par(mar=c(0,2,5,2))
groupedProfilePlot(p, groups='mukey', group.name.cex=0.65, color='sand', name='hzn_desgn', id.style='side', label='pedon_id', max.depth=200)
# describe IDs
mtext('user pedon ID', side=2, line=-1.5)
mtext('mukey', side=3, line=-1, at = c(0,0), adj = -0.8)

#silt
par(mar=c(0,2,5,2))
groupedProfilePlot(p, groups='mukey', group.name.cex=0.65, color='silt', name='hzn_desgn', id.style='side', label='pedon_id', max.depth=200)
# describe IDs
mtext('user pedon ID', side=2, line=-1.5)
mtext('mukey', side=3, line=-1, at = c(0,0), adj = -0.8)

