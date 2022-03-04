options(repos ="http://cran.us.r-project.org")  # required to get latest libs
# Installing the packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster,EcoHydRology,rnoaa,curl,httr,ggplot2,data.table)
rm(list=objects())

###Source TMWB model
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
##source CNmodel function
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")

# Download a soils dataset for basin based on the WebSoilSurvey method 
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/tsy2zccfzadznnkp5znbsg2p/wss_aoi_2022-03-03_18-28-08.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
#function to get data from USGS 0205551460 
#LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
mysoil=readOGR("wss_aoi_2022-02-11_09-31-24/spatial/soilmu_a_aoi.shp")    
# Explore the mysoil dataset which is returned
mybbox=c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
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
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2021). 
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
#summary(WXData)  #

# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
# Comparing precipitation to the flow out of basin
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0
modeldata[is.na(modeldata)]=0 # A Quick removal of NAs
TMWB=modeldata

CNmodeldf=modeldata


#..................HW1............................
#optimize TMWB
f <- function (x) {
  fcresopt=x[1]
  SFTmpopt=x[2]
  Tlagopt=x[3]
  AWCopt=x[4]
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=fcresopt,SFTmp = SFTmpopt,Tlag =Tlagopt,AWCval = AWCopt)
  1-NSE(TMWBnew$Qmm,TMWBnew$Qpred)  
}
lower <- c(0.1,-5,0,50)
upper <- c(0.5,20,1,350)

## run DEoptim and set a seed first for replicability
set.seed(1234)
DEoptim(f, lower, upper,control = DEoptim.control(itermax=60))
#0.3215601 9.005322 0.39573402 85.27068
TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.3215,SFTmp =9,Tlag = .4,AWCval = 85)
NSE(TMWBnew$Qmm,TMWBnew$Qpred)
#0.3627943
#optimize CNmodel
f <- function (x) {
  CNopt=x[1]
  IaOpt=x[2]
  DAWCOpt=x[3]
  zOpt=x[4]
  fcres_opt=x[5]
  CNmodelnew=CNmodel(CNmodeldf =CNmodeldf,CNavg = CNopt,IaFrac = IaOpt,func_DAWC=DAWCOpt,func_z=zOpt,fnc_fcres=fcres_opt)
  1-NSE(CNmodelnew$Qmm,CNmodelnew$Qpred)  
}

lower <- c(35,.01,0.1,500,0.1)
upper <- c(99,.25,0.35,1000,0.5)

## run DEoptim and set a seed first for replicability
set.seed(1234)
DEoptim(f, lower, upper,control = DEoptim.control(itermax=60))

#97.00911 0.03616312 0.3247503 954.2761 0.3265114

CNmodelnew=CNmodel(CNmodeldf =CNmodeldf,CNavg =97,IaFrac=0.036,fnc_slope=0, 
                   fnc_aspect=0,func_DAWC=.325,func_z=954.28,fnc_fcres=.3265)

NSE(CNmodelnew$Qmm,CNmodelnew$Qpred)
#0.5276314

#plotting
colors <- c("Q_Actual" = "red","Q_CN_NSE_0.527613"="darkgreen","Q_TMWB_NSE_0.3627943"="blue")

ggplot()+
  
  geom_line(data=TMWBnew, aes(x=date,y=Qpred,color='Q_TMWB_NSE_0.3627943'), size=.5) + 
  geom_line(data=TMWBnew, aes(x=date,y=Qmm,color='Q_Actual'), size=.5)+
  geom_line(data=CNmodelnew, aes(x=date,y=Qpred,color='Q_CN_NSE_0.527613'), size=.5)+
  
  
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
  ggtitle('Comparing Actual,TMWB and CN Modeled Discharge at Lick Run,VA')+
  scale_color_manual('Parameters:',values = colors)


#............................................
#...............HW2........................
#............................................

pacman::p_load(lubridate, data.table)
BasinCN_JO=CNmodelnew[(month(CNmodelnew$date) >= 5 & month(CNmodelnew$date) < 11),]


(1000/85-10)*25.4   # our CN estimate in bold
#44.82

(1000/50-10)*25.4   # our CN estimate in bold
#254

## by visually "guestimate" is that S should be somewhere between 
# 45mm and 260mm… repeat plotting until solution covers the 
# largest Qmm vs dP event (upper right hand corner of plot). 
#

# Assuming that (P-Ia) ~ dP, we can visually compare 
attach(BasinCN_JO)

# Now performing a “Calibration” using our method from Lab3 and the NSE
# as the “Objective Function”.  
# Vary S to maximize NSE using Eq. 4 of Lyon 2004 as our predictor of Q
#   Qpred=dP^2/(dP+S)
#
NSE(Qmm,dP^2/(dP+260))
# [1] 0.02375528
NSE(Qmm,dP^2/(dP+45))

#
# Keep iterating until NSE is as high as possible
# best estimate to S (Sest)
#
f <- function (x) {
  Sest=x
  NSE(Qmm,dP^2/(dP+Sest))
}
optimize(f, c(50,500), tol = 0.0001,maximum = TRUE)$maximum
Sest=121.6215

detach(BasinCN_JO)

# We will split into 5 VSA areas represented by 5 TI Classes
nTIclass=5
VSAsol=data.table(WetClass=seq(from=nTIclass,to=1),
                  As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
#
# Inspecting what the previous command gives us,we may need to 
# shift the index of a vector in the VSAsol data frame 
# using the data.table::shift() function.
#
VSAsol 
#
# Now fill in the missing value
#
VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1

VSAsol
#
# Calculate TI Class localized sigma and Curve Number
#
VSAsol[,sigma:=Sest*sSratio]
VSAsol[,CN:=25400/(sigma+254)]
VSAsol

# Initialize the TI Class objects from top to bottom of slope
TIC05=modeldata
TIC04=modeldata
TIC03=modeldata
TIC02=modeldata
TIC01=modeldata
# For TIC05 CNavg=VSAsol$CN[1] and higher fcres
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAsol$CN[1],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.5)

# Scaling reservoir coefficient between the .2-.5 given in class
# routing Qpred to ExcessIn below

TIC04$P=TIC05$Excess+TIC04$P
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAsol$CN[2],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.4)


TIC03$P=TIC04$Excess+TIC03$P
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAsol$CN[3],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)


TIC02$P=TIC03$Excess+TIC02$P
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAsol$CN[4],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)


TIC01$P=TIC02$Excess+TIC01$P
TIC01 = CNmodel(CNmodeldf = TIC01, CNavg=VSAsol$CN[5],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.2)

mean(TIC05$Qpred)
mean(TIC04$Qpred)
mean(TIC03$Qpred)
mean(TIC02$Qpred)
mean(TIC01$Qpred)

colors <- c("TI1" = "red","TI2"="darkgreen","TI3"="blue","TI4"="black","TI5"="purple")

ggplot()+
  
  geom_line(data=TIC01, aes(x=date,y=Qpred,color='TI1'), size=.5) + 
  geom_line(data=TIC02, aes(x=date,y=Qpred,color='TI2'), size=.5)+
  geom_line(data=TIC03, aes(x=date,y=Qpred,color='TI3'), size=.5)+
  geom_line(data=TIC04, aes(x=date,y=Qpred,color='TI4'), size=.5)+
  geom_line(data=TIC05, aes(x=date,y=Qpred,color='TI5'), size=.5)+
  
  
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
  annotate("text", label = "Average Qpred[mm]\nTI01= 4.704631\nTI02=3.927297\nTI03=3.173904\nTI04=2.430394\nTI05=1.612730", x = as.Date("2017-07-01"), y = 175)+
  xlab('Year')+
  ylab('Q[mm]')+
  ggtitle('Comparing Flow for Different TI at Lick Run,VA')+
  scale_color_manual('Parameters:',values = colors)

#............................................................
#.......................HW3..................................
#............................................................

mean(TIC05$AW)
mean(TIC04$AW)
mean(TIC03$AW)
mean(TIC02$AW)
mean(TIC01$AW)

max(TIC05$AW)
max(TIC04$AW)
max(TIC03$AW)
max(TIC02$AW)
max(TIC01$AW)

min(TIC05$AW)
min(TIC04$AW)
min(TIC03$AW)
min(TIC02$AW)
min(TIC01$AW)

colors <- c("TI1" = "red","TI2"="darkgreen","TI3"="blue","TI4"="black","TI5"="purple")

ggplot()+
  
  geom_line(data=TIC01, aes(x=date,y=AW,color='TI1'), size=.5) + 
  geom_line(data=TIC02, aes(x=date,y=AW,color='TI2'), size=.5)+
  geom_line(data=TIC03, aes(x=date,y=AW,color='TI3'), size=.5)+
  geom_line(data=TIC04, aes(x=date,y=AW,color='TI4'), size=.5)+
  geom_line(data=TIC05, aes(x=date,y=AW,color='TI5'), size=.5)+
  
  
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
  ggtitle('Comparing AW for Different TI at Lick Run,VA')+
  scale_color_manual('Parameters:',values = colors)


#............................................................
#.......................HW4..................................
#............................................................

mean(TIC05$ET)
mean(TIC04$ET)
mean(TIC03$ET)
mean(TIC02$ET)
mean(TIC01$ET)

colors <- c("TI1" = "red","TI2"="green","TI3"="blue","TI4"="black","TI5"="purple")

ggplot()+
  
  geom_line(data=TIC01, aes(x=date,y=ET,color='TI1'), size=.5) + 
  geom_line(data=TIC02, aes(x=date,y=ET,color='TI2'), size=.5)+
  geom_line(data=TIC03, aes(x=date,y=ET,color='TI3'), size=.5)+
  geom_line(data=TIC04, aes(x=date,y=ET,color='TI4'), size=.5)+
  geom_line(data=TIC05, aes(x=date,y=ET,color='TI5'), size=.5)+
  
  
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
  ylab('ET[mm]')+
  ggtitle('Comparing ET for Different TI at Lick Run,VA')+
  scale_color_manual('Parameters:',values = colors)




