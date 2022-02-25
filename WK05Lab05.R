objects()
rm(list=objects())
options(repos ="http://cran.us.r-project.org")  # required to get latest libs
# Installing the packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster,ggplot2,patchwork,EcoHydRology,rnoaa,curl,httr)

source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")

# Downloading a soils dataset for Lick Run basin based on the WebSoilSurvey method 

url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/bky5q4i2wkdjmvm3t0d1gniq/wss_aoi_2022-02-11_09-31-24.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")

#Data from USGS 0205551460 
#LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2022-03-01")
# Note that flow returned is in m3/day, converting it mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
# the soil data using the soilDB package
mysoil=readOGR("wss_aoi_2022-02-11_09-31-24/spatial/soilmu_a_aoi.shp")    
#Exploring mysoil dataset 
mybbox=c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
#bringing them back together, and aggregating based on max values of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.001)

summary(terrain(mydem, opt='slope',unit = "degrees"))

stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# have to chose stations with elements that have PRCP, TMAX and TMIN and current data (i.e. Year 2021). 
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
summary(WXData)  #

# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
summary(modeldata)  #
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
#handling missing values
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)]=0 # A Quick removal of NAs
TMWB=modeldata
# Last weeks homework example
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
# Calibrating the parameters one at a time
for (fcres in seq(.1,.5,.1)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=fcres)
  print(paste(fcres,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
# fcres=.3 is the highest NSE
for (SFTmp in seq(-5,20)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.3,SFTmp = SFTmp)
  print(paste(SFTmp,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
for(AWCval in seq(50,350,50)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.3,SFTmp = 11,Tlag = .5,AWCval = AWCval)
  print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
# Best result for "LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA" NSE = .34 . 
TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.3,SFTmp = 11,Tlag = .5,AWCval = 100)
print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))

# For a simple spatial model, use TMWB, to initialize 3 
# slope components for Top, Mid, and Bottom of the hillside. 
summary(modeldata)
TopSlope=modeldata
MidSlope=modeldata
BotSlope=modeldata
# Low slope but highest ksat
summary(mu2chmax)
#topslope has the lowest Z, lowest slope as botslope, midslope has the
#highest slope, we assume that awc in % is
#equal for all three models
awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
#for TopSlope model lowest depth
TopslopeZ=min(mu2chmax$hzdepb_r)*10 #in mm
MidslopeZ=mean(mu2chmax$hzdepb_r)*10
BotslopeZ=max(mu2chmax$hzdepb_r)*10
#AWCVAL is mm of awc= depth*awc(%)
AWCvalTop=TopslopeZ*awcpercent
AWCvalMid=MidslopeZ*awcpercent
AWCvalBot=BotslopeZ*awcpercent
#calculation of slope from terrain
# unit is in degree
summary(terrain(mydem, opt='slope',unit = "degrees"))
SlopeTop=1.837328 #degree
SlopeBot=1.837328 #degree
SlopeMid=40.600004 #degree
# Running the model with the three HRU dataframes
# Low slope but highest ksat
TopSlope = TMWBmodel(TMWB = TopSlope,SFTmp = 1,
                     AWCval = AWCvalTop,
                     Tlag = .5,fcres=.3,Slope = atan(SlopeTop/100))
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5
MidSlope = TMWBmodel(TMWB = MidSlope,SFTmp = 1,
                     AWCval = AWCvalMid,
                     Tlag = .5,fcres=0.5,Slope = atan(SlopeMid/100))
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlope = TMWBmodel(TMWB = BotSlope,SFTmp = 1,
                     AWCval = AWCvalBot,
                     Tlag = .5,fcres=0.2,Slope = atan(SlopeBot/100))
NSeff(BotSlope$Qmm,BotSlope$Qpred)

#............CNmodel........
CNmodel<-function(CNmodeldf, CNavg = 75,IaFrac = 0.05,fnc_slope=0,fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  # Energy Balance based Snow Accumulation
  # and Melt model from the EcoHydRology package.
  attach(CNmodeldf)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat,
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1,
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2,
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(CNmodeldf)
  CNmodeldf$SNO=SNO_Energy$SnowWaterEq_mm
  CNmodeldf$SNOmlt=SNO_Energy$SnowMelt_mm
  CNmodeldf$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  CNmodeldf$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(CNmodeldf)
  CNmodeldf$Albedo=.23
  CNmodeldf$Albedo[CNmodeldf$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  CNmodeldf$PET=PET
  detach(CNmodeldf)
  rm(list="PET")
  CNmodeldf$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  CNmodeldf$dP = 0 # Initializing Net Precipitation
  CNmodeldf$ET = 0 # Initializing ET
  CNmodeldf$AW = 0 # Initializing AW
  CNmodeldf$Excess = 0 # Initializing Excess
  CNmodeldf$S =0 # Initializing S
  CNmodeldf$Qpred=0 # Initializing Qpred
  attach(CNmodeldf)
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN
  CNmodeldf$CNavg = CNavg
  CNmodeldf$SSCNavg = SSCNavg
  CNmodeldf$SSCN = SSCN
  detach(CNmodeldf)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  CNmodeldf$Ia = Ia_init
  attach(CNmodeldf)
  # Those processes that are dependant on prior days conditions, we run as a
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] +
      SNO_Energy$SnowMelt_mm[t] # CN Solution
    # Is the soil saturated, and thus can't take more dP?
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction?
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  CNmodeldf$ET=ET
  CNmodeldf$dP=dP
  CNmodeldf$AW=AW
  CNmodeldf$Excess=Excess
  CNmodeldf$S=S
  CNmodeldf$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(CNmodeldf)
  return(CNmodeldf)
}
#
# Like before, initializing the 3 hillslope classes
#
TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata
# Call the new CNmodel() function with Top,Mid,BotSlope HRU objects,
# passing the Qpred into the lower HRUs HillslopeAboveExcess (as area scaledflow)
TopSlopeCN=CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)
NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)
for(IaFrac in seq(0.01,0.4,0.01)){
  TMWBnew=CNmodel(CNmodeldf=BotSlopeCN,CNavg = 75,IaFrac = IaFrac,fnc_slope=0,fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.2)
  print(paste(IaFrac,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
#
colors <- c("Q_Actual" = "red","Q_CN"="darkgreen")
p1=ggplot()+
  
  geom_line(data=BotSlopeCN, aes(x=date,y=Qpred,color='Q_CN'), size=.5) + 
  geom_line(data=BotSlopeCN, aes(x=date,y=Qmm,color='Q_Actual'), size=.5)+

  
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
  ylab('Q[mm]')+
  annotate("text", label = "NSE=0.1520626", x = as.Date("2017-07-01"), y = 40)+
  ggtitle('Comparing Actual and CN Modeled Discharge at Lick Run,VA')+
  scale_color_manual('Parameters:',values = colors)

colors <- c("Q_Actual" = "red", "Q_TMWB" = "blue")

p2=ggplot()+
  
  geom_line(data=BotSlope, aes(x=date,y=Qmm,color='Q_Actual'), size=.5)+
  geom_line(data=BotSlope, aes(x=date,y=Qpred,color='Q_TMWB'), size=.5)+
  
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
  ylab('Q[mm]')+
  annotate("text", label = "NSE=-6.366674", x = as.Date("2017-07-01"), y = 70)+
  ggtitle('Comparing Actual and TMWB Modeled Discharge at Lick Run,VA')+
  scale_color_manual('Parameters:',values = colors)

p1+p2+ plot_layout(ncol = 1, widths = c(1, 1))

#....................................................
#..............HW3....................................
#......................................................
objects()
rm(list=objects())
options(repos ="http://cran.us.r-project.org")  # required to get latest libs
# Installing the packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster,ggplot2,patchwork,EcoHydRology,rnoaa,curl,httr)

source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")

# Downloading a soils dataset for Saugus River basin based on the WebSoilSurvey method 

url='https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/1fzi4r2xg2q1x3sq4zfgk0xs/wss_aoi_2022-02-17_08-38-07.zip'
download.file(url,"mysoil.zip")
unzip("mysoil.zip")

# using the function to get data from USGS 01102345
#SAUGUS RIVER AT SAUGUS IRONWORKS AT SAUGUS, MA
myflowgage_id="01102345"
myflowgage=get_usgs_gage(myflowgage_id,begin_date="2010-1-1",end_date="2022-02-01")
# Note that flow returned is in m3/day, converting it mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
# the soil data using the soilDB package
mysoil=readOGR("wss_aoi_2022-02-17_08-38-07/spatial/soilmu_a_aoi.shp")    
#Exploring mysoil dataset 
mybbox=c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
#bringing them back together, and aggregating based on max values of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.001)


stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 40,
  limit = NULL
)
# have to chose stations with elements that have PRCP, TMAX and TMIN and current data (i.e. Year 2021). 
WXData=meteo_pull_monitors(
  monitors=stns[201,1],
  keep_flags = FALSE,
  date_min = "2010-1-1",
  date_max ="2022-02-01",
  var = c("TMAX","TMIN","PRCP")
)

# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
summary(modeldata)  #
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
#handling missing values
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)]=0 # A Quick removal of NAs
TMWB=modeldata

source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
# Calibrating the parameters one at a time
for (fcres in seq(.1,.5,.1)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=fcres)
  print(paste(fcres,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
# fcres=.2 is the highest NSE
for (SFTmp in seq(-5,20)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.2,SFTmp = SFTmp)
  print(paste(SFTmp,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
#SFTmp=4 is the highest NSE
for(AWCval in seq(50,350,50)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.2,SFTmp = 4,Tlag = .5,AWCval = AWCval)
  print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
#AWCval=100 is the highest NSE
# Best result for "SAUGUS River, MO" NSE = -0.32462
TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.2,SFTmp = 4,Tlag = .5,AWCval = 100)
print(paste(NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
# For a simple spatial model, use TMWB, to initialize 3 
# slope components for Top, Mid, and Bottom of the hillside. 
summary(modeldata)
TopSlope=modeldata
MidSlope=modeldata
BotSlope=modeldata
# Low slope but highest ksat
#############################################
##GETTING ARGUMENT VALUES BASED ON EACH MODEL
#############################################
summary(mu2chmax)
##topslope has the lowest Z, lowest slope as botslope, midslope has the
#highest slope, we assume that awc in % is
#equal for all three models
awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
### for TopSlope model lowest depth
TopslopeZ=min(mu2chmax$hzdepb_r)*10 #in mm
MidslopeZ=mean(mu2chmax$hzdepb_r)*10
BotslopeZ=max(mu2chmax$hzdepb_r)*10
#####################AWCVAL is mm of awc= depth*awc(%)
AWCvalTop=TopslopeZ*awcpercent
AWCvalMid=MidslopeZ*awcpercent
AWCvalBot=BotslopeZ*awcpercent

summary(terrain(mydem, opt='slope',unit = "radians"))

SlopeTop=0 #radians
SlopeBot=0 #radians
SlopeMid=0.147622 #radians
# Running the model with the three HRU dataframes
# Low slope but highest ksat
TopSlope = TMWBmodel(TMWB = TopSlope,SFTmp =20,
                     AWCval = AWCvalTop,
                     Tlag = .5,fcres=0.2,Slope =SlopeTop)
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5
MidSlope = TMWBmodel(TMWB = MidSlope,SFTmp =20,
                     AWCval = AWCvalMid,
                     Tlag = .5,fcres=0.2,Slope = SlopeMid)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlope = TMWBmodel(TMWB = BotSlope,SFTmp =20,
                     AWCval = AWCvalBot,
                     Tlag = .5,fcres=0.2,Slope = SlopeBot)
NSeff(BotSlope$Qmm,BotSlope$Qpred)
#............CNmodel........
CNmodel<-function(CNmodeldf, CNavg = 75,IaFrac = 0.05,fnc_slope=0,fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  # Energy Balance based Snow Accumulation
  # and Melt model from the EcoHydRology package.
  attach(CNmodeldf)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat,
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1,
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2,
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(CNmodeldf)
  CNmodeldf$SNO=SNO_Energy$SnowWaterEq_mm
  CNmodeldf$SNOmlt=SNO_Energy$SnowMelt_mm
  CNmodeldf$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  CNmodeldf$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(CNmodeldf)
  CNmodeldf$Albedo=.23
  CNmodeldf$Albedo[CNmodeldf$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  CNmodeldf$PET=PET
  detach(CNmodeldf)
  rm(list="PET")
  CNmodeldf$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  CNmodeldf$dP = 0 # Initializing Net Precipitation
  CNmodeldf$ET = 0 # Initializing ET
  CNmodeldf$AW = 0 # Initializing AW
  CNmodeldf$Excess = 0 # Initializing Excess
  CNmodeldf$S =0 # Initializing S
  CNmodeldf$Qpred=0 # Initializing Qpred
  attach(CNmodeldf)
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN
  CNmodeldf$CNavg = CNavg
  CNmodeldf$SSCNavg = SSCNavg
  CNmodeldf$SSCN = SSCN
  detach(CNmodeldf)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  CNmodeldf$Ia = Ia_init
  attach(CNmodeldf)
  # Those processes that are dependant on prior days conditions, we run as a
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] +
      SNO_Energy$SnowMelt_mm[t] # CN Solution
    # Is the soil saturated, and thus can't take more dP?
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction?
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  CNmodeldf$ET=ET
  CNmodeldf$dP=dP
  CNmodeldf$AW=AW
  CNmodeldf$Excess=Excess
  CNmodeldf$S=S
  CNmodeldf$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(CNmodeldf)
  return(CNmodeldf)
}

# Like before, initializing the 3 hillslope classes
#
TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata
# Call the new CNmodel() function with Top,Mid,BotSlope HRU objects,
# passing the Qpred into the lower HRUs HillslopeAboveExcess (as area scaledflow)
TopSlopeCN=CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)
NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)

colors <- c("Q_Actual" = "red","Q_CN"="darkgreen")

p1=ggplot()+
  
  geom_line(data=BotSlopeCN, aes(x=date,y=Qpred,color='Q_CN'), size=.5) + 
  geom_line(data=BotSlopeCN, aes(x=date,y=Qmm,color='Q_Actual'), size=.5)+
  
  
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
  ylab('Q[mm]')+
  annotate("text", label = "NSE=-0.04094479", x = as.Date("2012-07-01"), y = 30)+
  ggtitle('Comparing Actual and CN Modeled Discharge at SAUGUS River, MO')+
  scale_color_manual('Parameters:',values = colors)

colors <- c("Q_Actual" = "red", "Q_TMWB" = "blue")

p2=ggplot()+
  
  geom_line(data=BotSlope, aes(x=date,y=Qmm,color='Q_Actual'), size=.5)+
  geom_line(data=BotSlope, aes(x=date,y=Qpred,color='Q_TMWB'), size=.5)+
  
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
  ylab('Q[mm]')+
  annotate("text", label = "NSE=-0.7050565", x = as.Date("2012-07-01"), y = 50)+
  ggtitle('Comparing Actual and TMWB Modeled Discharge at SAUGUS River, MO')+
  scale_color_manual('Parameters:',values = colors)

p1+p2+ plot_layout(ncol = 1, widths = c(1, 1))

#snow plot
colors <- c("Snow_melt" = "red","SWE"="darkgreen")

p3=ggplot()+
  
  geom_line(data=BotSlopeCN, aes(x=date,y=SnowMelt_mm,color='Snow_melt'), size=.5) + 
  geom_line(data=BotSlopeCN, aes(x=date,y=SnowfallWatEq_mm,color='SWE'), size=.5)+
  
  
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
  ylab('[mm]')+
  ggtitle('CN Modeled Snowmelt and SWE at SAUGUS River, MO')+
  scale_color_manual('Parameters:',values = colors)

colors <- c("Snow_melt" = "red","SWE"="darkgreen")

p4=ggplot()+
  
  geom_line(data=BotSlope, aes(x=date,y=SNOmlt,color='Snow_melt'), size=.5)+
  geom_line(data=BotSlope, aes(x=date,y=SNO,color='SWE'), size=.5)+
  
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
  ylab('[mm]')+
  ggtitle('TMWB Modeled Snowmelt and SWE at SAUGUS River, MO')+
  scale_color_manual('Parameters:',values = colors)

p3+p4+plot_layout(ncol = 1, widths = c(1, 1))

#.........................................................
#...................HW4...................................
#.........................................................
objects()
rm(list=objects())
options(repos ="http://cran.us.r-project.org")  # required to get latest libs
# Installing the packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster,ggplot2,patchwork,EcoHydRology,rnoaa,curl,httr)

source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")

# Downloading a soils dataset for TOWN BROOK NY basin based on the WebSoilSurvey method 
url='https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/x2qg1plc5ik5fu44pft4w1pu/wss_aoi_2022-02-24_20-47-58.zip'
download.file(url,"mysoil.zip")
unzip("mysoil.zip")

# using the function to get data from USGS 01421618
#TOWN BROOK NY
myflowgage_id="01421618"
myflowgage=get_usgs_gage(myflowgage_id,begin_date="2010-1-1",end_date="2022-02-01")
# Note that flow returned is in m3/day, converting it mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
# the soil data using the soilDB package
mysoil=readOGR("wss_aoi_2022-02-24_20-47-58/spatial/soilmu_a_aoi.shp")    
#Exploring mysoil dataset 
mybbox=c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
#bringing them back together, and aggregating based on max values of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.001)


stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 40,
  limit = NULL
)
# have to chose stations with elements that have PRCP, TMAX and TMIN and current data (i.e. Year 2021). 
WXData=meteo_pull_monitors(
  monitors=stns[150,1],
  keep_flags = FALSE,
  date_min = "2010-1-1",
  date_max ="2022-02-01",
  var = c("TMAX","TMIN","PRCP")
)
# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
summary(modeldata)  #
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
#handling missing values
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)]=0 # A Quick removal of NAs
TMWB=modeldata
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
# Calibrating the parameters one at a time
for (fcres in seq(.1,.2,.1)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=fcres)
  print(paste(fcres,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
# fcres=.2 is the highest NSE
for (SFTmp in seq(-5,-4)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.2,SFTmp = SFTmp)
  print(paste(SFTmp,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
#SFTmp=20 is the highest NSE
for(AWCval in seq(50,200,50)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.2,SFTmp =-4,Tlag = .5,AWCval = AWCval)
  print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
#AWCval=200 is the highest NSE
# Best result for "SAUGUS River, MO" NSE = -0.32462
TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.2,SFTmp =-4,Tlag = .5,AWCval = 200)
print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))

CNmodel<-function(CNmodeldf, CNavg = 75,IaFrac = 0.05,fnc_slope=0,fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  # Energy Balance based Snow Accumulation
  # and Melt model from the EcoHydRology package.
  attach(CNmodeldf)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat,
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1,
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2,
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(CNmodeldf)
  CNmodeldf$SNO=SNO_Energy$SnowWaterEq_mm
  CNmodeldf$SNOmlt=SNO_Energy$SnowMelt_mm
  CNmodeldf$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  CNmodeldf$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(CNmodeldf)
  CNmodeldf$Albedo=.23
  CNmodeldf$Albedo[CNmodeldf$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  CNmodeldf$PET=PET
  detach(CNmodeldf)
  rm(list="PET")
  CNmodeldf$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  CNmodeldf$dP = 0 # Initializing Net Precipitation
  CNmodeldf$ET = 0 # Initializing ET
  CNmodeldf$AW = 0 # Initializing AW
  CNmodeldf$Excess = 0 # Initializing Excess
  CNmodeldf$S =0 # Initializing S
  CNmodeldf$Qpred=0 # Initializing Qpred
  attach(CNmodeldf)
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN
  CNmodeldf$CNavg = CNavg
  CNmodeldf$SSCNavg = SSCNavg
  CNmodeldf$SSCN = SSCN
  detach(CNmodeldf)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  CNmodeldf$Ia = Ia_init
  attach(CNmodeldf)
  # Those processes that are dependant on prior days conditions, we run as a
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] +
      SNO_Energy$SnowMelt_mm[t] # CN Solution
    # Is the soil saturated, and thus can't take more dP?
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction?
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  CNmodeldf$ET=ET
  CNmodeldf$dP=dP
  CNmodeldf$AW=AW
  CNmodeldf$Excess=Excess
  CNmodeldf$S=S
  CNmodeldf$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(CNmodeldf)
  return(CNmodeldf)
}


for(IaFrac in seq(0.01,0.01,0.01)){
  TMWBnew=CNmodel(CNmodeldf=TMWB,CNavg = 75,IaFrac = IaFrac,fnc_slope=0,fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.2)
  print(paste(IaFrac,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
#Iafrac=0.01
for(CNavg in seq(10,95,5)){
  TMWBnew=CNmodel(CNmodeldf=TMWB,CNavg =CNavg,IaFrac = 0.01,fnc_slope=0,
                  fnc_aspect=0,func_DAWC=.3,
                  func_z=500,fnc_fcres=.2)
  print(paste(CNavg,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
#CN=95

