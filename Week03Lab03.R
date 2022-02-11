# data from https://goo.gl/Cb8zGn

source("https://goo.gl/Cb8zGn")

# get_usgs_gage() is a function returned by this source() call
get_usgs_gage

# using the function to get data from USGS 01102345
#SAUGUS RIVER AT SAUGUS IRONWORKS AT SAUGUS, MA
myflowgage_id="01102345"
myflowgage=get_usgs_gage(myflowgage_id,begin_date="2000-1-1",end_date="2022-02-01")

#weather data from NOAA

install.packages("rnoaa");
library(rnoaa)
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

#merging flowdata to weather data
modeldata=merge(WXData, myflowgage$flowdata, by.x="date", by.y="mdate")
View(modeldata)

#coverting discgarge from m3/day to mm/day
modeldata$flowmm=myflowgage$flowdata$flow/myflowgage$area/10^3
modeldata$Qmm=modeldata$flowmm

#converting weather parameters
modeldata$MaxTemp_C=modeldata$tmax/10 # Converting to C
modeldata$MinTemp_C=modeldata$tmin/10 # Converting to C
modeldata$P_mm=modeldata$prcp/10 # Converting to mm
View(modeldata)

# Comparing precipitation to the flow out of the basin
mean(modeldata$Qmm)
mean(modeldata$P_mm)
#ignoring nan values
modeldata$P_mm[is.na(modeldata$P_mm)]=0
mean(modeldata$P_mm)
mean(modeldata$MaxTemp_C)
#ignoring nan values
modeldata$MaxTemp_C[is.na(modeldata$MaxTemp_C)]=0
mean(modeldata$MaxTemp_C)
mean(modeldata$MinTemp_C)
#ignoring nan values
modeldata$MinTemp_C[is.na(modeldata$MinTemp_C)]=0
mean(modeldata$MinTemp_C)
summary(modeldata)

#creating new dataset for Thornthwaite Mather model
TMWB=modeldata

#Energy Balance based Snow Accumulation and Melt model 
#EcoHydRology package snowmelt will be used

attach(TMWB)
#for flat slope
SNO_Energy_flat=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                      slope = 0,
                      aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                      startingSnowDensity_kg_m3=450)
#for 10% North
SNO_Energy_N=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                         slope =atan(10/100),
                         aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                         SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                         startingSnowDensity_kg_m3=450)
#for 10% South
SNO_Energy_S=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                         slope = atan(10/100),
                         aspect = 180*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                         SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                         startingSnowDensity_kg_m3=450)
#for 45% NW
SNO_Energy_NW=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                         slope =atan(45/100),
                         aspect = 300*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                         SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                         startingSnowDensity_kg_m3=450)
#for45% SW
SNO_Energy_SW=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                         slope =atan(45/100),
                         aspect = 145*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                         SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                         startingSnowDensity_kg_m3=450)

detach(TMWB)

#plotting SWE for 2001 to 2022
if(!require("pacman"))install.packages("pacman")
pacman::p_load(EcoHydRology,lattice,ggplot2,dplyr,patchwork,hrbrthemes,tidyverse)#loadlibrary

colors <- c('Flat' = 'black', '10%North' = 'red','10%South'='blue','45%NW'='darkgreen','45%SW'='darkorange')

ggplot()+
  
  geom_line(data=SNO_Energy_flat, aes(x=Date,y=SnowWaterEq_mm,color='Flat'), size=.5) + 
  geom_line(data=SNO_Energy_N, aes(x=Date,y=SnowWaterEq_mm,color='10%North'), size=.5)+
  geom_line(data=SNO_Energy_S, aes(x=Date,y=SnowWaterEq_mm,color='10%South'), size=.5)+
  geom_line(data=SNO_Energy_NW, aes(x=Date,y=SnowWaterEq_mm,color='45%NW'), size=.5)+
  geom_line(data=SNO_Energy_SW, aes(x=Date,y=SnowWaterEq_mm,color='45%SW'), size=.5)+

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
  ylab('SWE[mm]')+
  ggtitle('SWE for different Slope and Aspect')+
  scale_color_manual('Parameters:',values = colors)

#kind of messy plotting let's plot for one year
#TMWB1=subset(TMWB, format(as.Date(date),"%Y")==2019)
class(TMWB$date) 
TMWB1= TMWB[TMWB$date >"2019-01-15" & TMWB$date < "2019-03-31", ]
attach(TMWB1)
#for flat slope
SNO_Energy_flat1=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                         slope = 0,
                         aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                         SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                         startingSnowDensity_kg_m3=450)
#for 10% North
SNO_Energy_N1=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                      slope =atan(10/100),
                      aspect = 18*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                      startingSnowDensity_kg_m3=450)
#for 10% South
SNO_Energy_S1=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                      slope = atan(10/100),
                      aspect = 190*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                      startingSnowDensity_kg_m3=450)
#for 45% NW
SNO_Energy_NW1=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                       slope =atan(45/100),
                       aspect = 330*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                       SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                       startingSnowDensity_kg_m3=450)
#for45% SW
SNO_Energy_SW1=SnowMelt(date, P_mm, MaxTemp_C,MinTemp_C, myflowgage$declat,
                       slope =atan(45/100),
                       aspect = 210*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                       SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                       startingSnowDensity_kg_m3=450)

detach(TMWB1)
colors <- c('Flat' = 'black', '10%North' = 'red','10%South'='blue','45%NW'='darkgreen','45%SW'='darkorange')

ggplot()+
  
  geom_line(data=SNO_Energy_flat1, aes(x=Date,y=SnowWaterEq_mm,color='Flat'), size=.5) + 
  geom_line(data=SNO_Energy_N1, aes(x=Date,y=SnowWaterEq_mm,color='10%North'), size=.5)+
  geom_line(data=SNO_Energy_S1, aes(x=Date,y=SnowWaterEq_mm,color='10%South'), size=.5)+
  geom_line(data=SNO_Energy_NW1, aes(x=Date,y=SnowWaterEq_mm,color='45%NW'), size=.5)+
  geom_line(data=SNO_Energy_SW1, aes(x=Date,y=SnowWaterEq_mm,color='45%SW'), size=.5)+
  
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
  xlab('Month')+
  ylab('SWE[mm]')+
  ggtitle('SWE for different Slope and Aspect for 2019')+
  scale_color_manual('Parameters:',values = colors)
sum(SNO_Energy_flat$SnowWaterEq_mm)
sum(SNO_Energy_N$SnowWaterEq_mm)
sum(SNO_Energy_NW$SnowWaterEq_mm)
sum(SNO_Energy_S$SnowWaterEq_mm)
sum(SNO_Energy_SW$SnowWaterEq_mm)

#.........................................................................................................
#HW2
TMWB=modeldata
if(!require("pacman"))install.packages("pacman")
pacman::p_load(EcoHydRology,lattice,ggplot2,dplyr,patchwork,hrbrthemes,tidyverse)#loadlibrary

soilwetting<-function(AWprev,dP_func,AWC_func){
  AW_func<-AWprev+dP_func
  excess_func<-0.0
  c(AW_func,excess_func)
}
# soildrying function
soildrying<-function(AWprev,dP_func,AWC_func){
  AW_func=AWprev*exp(dP_func/AWC_func)
  excess_func<-0.0
  c(AW_func,excess_func)
}

soil_wetting_above_capacity<-function(AWprev,dP_func,AWC_func){
  AW_func<-AWC_func
  excess_func<-AWprev+dP_func-AWC_func
  c(AW_func,excess_func)
}

SFTmp = 3 # referred to as SFTMP in SWAT input (Table 1)
bmlt6 = 7 # referred to as SMFMX in SWAT input (Table 1)
bmlt12 =1.4 # referred to as SMFMN in SWAT input adjusted for season
Tmlt = SFTmp # Assumed to be same as SnowFall Temperature
Tlag = 0 # referred to as TIMP in SWAT input (Table 1)
TMWB$AvgTemp=(TMWB$MaxTemp_C+TMWB$MinTemp_C)/2
TMWB$bmlt = (bmlt6 + bmlt12)/2 + (bmlt6 - bmlt12)/2 *
  sin(2*pi/365*(julian(TMWB$date,origin = as.Date("2000-01-01"))-81))
# Initializing SNO, Tsno as well as the first values of each
TMWB$SNO = 0 # Snow Depth (mm)
TMWB$Tsno = 0 # Snow Temp (C)
TMWB$SNOmlt = 0 # Snow Melt (mm)
attach(TMWB)
for (t in 2:length(date)){
  Tsno[t]= Tsno[t-1] * (1.0-Tlag) + AvgTemp[t] * Tlag
  if(AvgTemp[t] < SFTmp){
    SNO[t]= SNO[t-1] + P_mm[t]
  } else {
    SNOmlt[t]= bmlt[t] * SNO[t-1] * ((Tsno[t]+MaxTemp_C[t])/2 - Tmlt)
    SNOmlt[t]= min(SNOmlt[t],SNO[t-1])
    SNO[t]= SNO[t-1] -SNOmlt[t]
  }
  print(t)
}
plot(date,SNO,type="l")
detach(TMWB)
TMWB$Albedo=.23
TMWB$Albedo[TMWB$SNO>0]=.95
attach(TMWB)

PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp_C,Tmin_C =
                   MinTemp_C,albedo=Albedo,lat_radians = myflowgage$declat*pi/180) * 1000
TMWB$PET=PET


detach(TMWB)
rm(list=c("PET"))
mean((TMWB$PET))
TMWB$PET[is.na(TMWB$PET)]=0
TMWB$AWC=150
TMWB$dP = 0 # Initializing Net Precipitation
TMWB$ET = 0 # Initializing ET
TMWB$AW = 0 # Initializing AW
TMWB$Excess = 0 # Initializing Excess

# Loop to calculate AW and Excess
attach(TMWB)
for (t in 2:length(AW)){
  
  ET[t] = min (AW[t-1],PET[t])
  ET[t] = (AW[t-1]/AWC[t-1])*PET[t] # New Model
  if(AvgTemp[t] >= SFTmp){
    dP[t] = P_mm[t] - ET[t] + SNOmlt[t]
  } else {
    dP[t] = ET[t]
  }
  
  
  if (dP[t]<=0) {
    values<-soildrying(AW[t-1],dP[t],AWC[t])
  } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
    values<-soilwetting(AW[t-1],dP[t],AWC[t])
  } else {
    values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
  }
  AW[t]<-values[1]
  Excess[t]<-values[2]
  print(t)
}
TMWB$AW=AW
TMWB$Excess=Excess
TMWB$dP=dP
rm(list=c("AW","dP","ET", "Excess"))
detach(TMWB) 
TMWB$Qpred=NA
TMWB$Qpred[1]=0
TMWB$S=NA
TMWB$S[1]=0
attach(TMWB)
fcres=.2
for (t in 2:length(date)){
  S[t]=S[t-1]+Excess[t]
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
}
TMWB$S=S
TMWB$Qpred=Qpred 

detach(TMWB) 
rm(list=c("Qpred","S"))

NSE=function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs,
                                                         na.rm=TRUE))^2, na.rm=TRUE))
}
NSE(TMWB$Qmm,TMWB$Qpred)
