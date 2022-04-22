rm(list=objects())
setwd("~")
dir.create("Lab12")
setwd("~/Lab12/")
install.packages("EcoHydRology", repos="http://R-Forge.R-project.org")
library(EcoHydRology)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,curl,data.table,multisensi)

pacman::p_load(rnoaa,ggplot2,fastmatch,forcats,dplyr)
#--------------source CN model function from previous lab----------
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")
#SAUGUS RIVER AT SAUGUS IRONWORKS AT SAUGUS, MA
#---------Getting streamflow data from USGS-------------------
myflowgage_id="01102345"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2010-01-01",end_date = "2022-04-19")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
#--------Getting weather data from FillMissWX (noaa) function
WXData=FillMissWX(declat=myflowgage$declat, declon=myflowgage$declon,60,
                  date_min="2010-01-01",
                  date_max="2022-04-19")
#---building modeldata dataframe by merging streamflow and weather data-------
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
#---Running CNmodel function to check the model-------------------------------
CN_DF=CNmodel(CNmodeldf = modeldata,CNavg = 70,
              IaFrac = 0.05,fnc_slope =  0.05,
              fnc_aspect =0.0,func_DAWC = 0.35,
              func_z=1000,fnc_fcres=0.20,
              declat=myflowgage$declat,
              declon=myflowgage$declon)
NSeff(CN_DF$Qmm,CN_DF$Qpred)
#-----------getting ready for optimization using DEoptim----------------------

f <- function (x) {
  CNopt=x[1]
  IaOpt=x[2]
  fnc_slopeOpt=x[3]
  fnc_aspectOpt=x[4]
  func_DAWCOpt=x[5]
  func_zOpt=x[6]
  fnc_fcresOpt=x[7]
  
  CNmodelnew=CNmodel(CNmodeldf =modeldata,
                     CNavg = CNopt,IaFrac = IaOpt,
                     fnc_slope=fnc_slopeOpt,fnc_aspect=fnc_aspectOpt,
                     func_DAWC=func_DAWCOpt,func_z=func_zOpt,
                     fnc_fcres=fnc_fcresOpt,declat=myflowgage$declat,
                     declon=myflowgage$declon)
  mNSE=1-NSeff(CNmodelnew$Qmm,CNmodelnew$Qpred)
  return(mNSE)
}
#            Educated guesses for ranges… or not.        
#            CN,  Ia, Slp, Asp, AWC,Depth,ResCoef
lower <- c(30,0.01, 0.0, 0.0, 0.2,  500,   0.1)
upper <- c(99,0.20, 0.1, 0.1, 0.4, 2000,   0.5)

# run DEoptim and set a seed first for replicability
set.seed(1234)
#
# The cluster will show you 128 cores, but the management software
# will limit you to running on however many cores you requested.


cl <- parallel::makeCluster(16)
outDEoptim=DEoptim(f, lower, upper,
                   DEoptim.control(cluster=cl,strategy = 6,
                                   NP = 16,itermax=100,parallelType = 1,
                                   packages = c("EcoHydRology"),
                                   parVar=c("CNmodel","SnowMelt","PET_fromTemp",
                                            "SoilStorage","NSeff","modeldata","myflowgage","EstCloudiness","transmissivity","AtmosphericEmissivity","SatVaporDensity","SatVaporPressure","Longwave","Solar","PotentialSolar","declination","slopefactor",
                                            "solarangle","NetRad","PTpet","SatVapPresSlope")))

outDEoptim$optim$bestmem
# NSE = 1- bestval
1-outDEoptim$optim$bestval
x=outDEoptim$optim$bestmem
f(x)
names(outDEoptim$member$lower)=c("CNavg","IaFrac","fnc_slope", 
                                 "fnc_aspect","func_DAWC","func_z",
                                 "fnc_fcres")
plot(outDEoptim)
plot(outDEoptim,plot.type="bestvalit")
#plot(outDEoptim,plot.type="bestmemit")
dev.off()
##-------plot NSE vs iteration------------
NSE=data.frame(1-outDEoptim$member$bestvalit)
colnames(NSE)="NSE"
plot(NSE$NSE,type = "l", ylab = "NSE", xlab = "Iteration number")
#----another way to plot all parameters of differnt iteration----
BestMemit=data.frame(outDEoptim$member$bestmemit)
colnames(BestMemit)=c("CNavg","IaFrac","fnc_slope", 
                      "fnc_aspect","func_DAWC","func_z",
                      "fnc_fcres")
#plotting the progression of improvements similar to above
for(i in 1:length(BestMemit)){
  plot(BestMemit[,i], ylab= "",xlab="",ylim=c(lower[i],upper[i]),cex=0.3)
  mtext(side=1, line=2, "iteration", col="black",cex=1)
  mtext(side=2, line=3,colnames(BestMemit)[i], col="black", font=2, cex=1)
}
#---------merging NSE and BestMemit datafram ---------------------------
Sensitivity_df=cbind(BestMemit,NSE)
#---calculate changes between two succesive rows-----------------------------
deltaparams=(tail(Sensitivity_df,-1)-Sensitivity_df[1:(length(Sensitivity_df[,1])-1),])
# 
#----calculate the relative sensitivity of each parameter------
#       Remember this slide 13 from the sensitivity lecture?
shifted_sensitivity_df=Sensitivity_df[1:(length(Sensitivity_df[,1]))-1,]
for (i in colnames(deltaparams)){
  #i="CNavg"
  nam = paste0("junkpar_",i)
  assign(nam,deltaparams[!deltaparams[,fmatch(i,names(deltaparams))]==0,] )
  nam2=paste0("initialval_",i)
  assign(nam2,shifted_sensitivity_df[!deltaparams[,fmatch(i,names(deltaparams))]==0,])
  nam3=paste0("Sr_",i)
  assign(nam3,abs(((get(paste0("junkpar_",i))$NSE)/
                     (get(paste0("junkpar_",i))[fmatch(i,names(deltaparams))]))*
                    (get(paste0("initialval_",i))[fmatch(i,names(deltaparams))]/
                       get(paste0("initialval_",i))$NSE)))
}

# 
# With sr for each params and you can get the mean of them to compare the sr
# Making a dataframe with bbox and violin plots

RelSensi=base::as.data.frame(stack(list(CNavg=Sr_CNavg$CNavg,
                                        IaFrac=Sr_IaFrac$IaFrac,
                                        fnc_slope=Sr_fnc_slope$fnc_slope, 
                                        fnc_aspect=Sr_fnc_aspect$fnc_aspect,
                                        func_DAWC=Sr_func_DAWC$func_DAWC,
                                        func_z=Sr_func_z$func_z,
                                        fnc_fcres=Sr_fnc_fcres$fnc_fcres
)))
names(RelSensi)[2]="Parameter"
# Basic violin plot
p <- ggplot(RelSensi, aes(x=Parameter, y=values, color=Parameter)) + 
  coord_cartesian(ylim = c(0,25))+
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme(legend.position = "None")

p
p + stat_summary(fun.y=mean, geom="point", shape=8, size=4,col="black") +
  xlab("Parameter")+
  ylab("Relative Sensitivity")+
  theme(axis.text = element_text(size =10))+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_text(size =10))+
  theme(axis.title = element_text(face="bold"))

#.............................................................................


pacman::p_load(neuralnet,ggplot2,Metrics,reshape2)



dat=data.frame(Qm3ps=CNmodelnew$Qmm,
               Jday=strptime(CNmodelnew$date,"%Y-%m-%d")$yday+1,Precip=CNmodelnew$P,TMX=CNmodelnew$MaxTemp,TMN=CNmodelnew$MinTemp,
               jdate=julian(CNmodelnew$date))
# Need the previous days flow
# The NN for this example requires the previous days flow… remember this for 
# the homework questions. 
dat$Qm3ps_lag=dplyr::lag(dat$Qm3ps)
dat[is.na(dat)]=0
View(dat)
head(dat)
datmax <- apply(dat, 2, max)
datmin <- apply(dat, 2, min)
# now split the datat into train, validation and test sets.
scaled <- scale(dat, center = datmin, scale = datmax - datmin)

# sub sample to get the test and train set
index <- seq(1, nrow(dat))
train_index <- sample(index, 0.8*length(index))
train <- scaled[train_index,]
test <- scaled[-train_index,]
#===================================================
#             Now, setup and train the neural network!!!
#===================================================
# first set a random seed to a lucky number, but it does give results!

set.seed(2)
# setup and train the neural network
net <- neuralnet(Qm3ps ~ Qm3ps_lag+Precip+jdate+Jday+TMX+TMN,
                 data = train,
                 hidden = c(6, 3), # can adjust here
                 threshold = 0.02, # if you network does not run, 
                 #                         try turning this up just a little
                 stepmax = 1e+05, # can adjust if you don't converge
                 rep = 1,
                 #startweights = NULL,
                 #learningrate.limit = NULL,
                 #learningrate.factor = list(minus = 0.5, plus = 1.2),
                 learningrate = 0.0027,
                 lifesign = "full",
                 #lifesign.step = 1000,
                 algorithm = "backprop",
                 err.fct = 'sse',
                 act.fct = "logistic",
                 linear.output = FALSE,
                 exclude = NULL,
                 constant.weights = NULL)
#===================================================
#                     Plot the Structure of the Neural Net
#===================================================
# View the structure of the neural net
plot(net)
#===================================================
#             Make Predictions to See how the neural network did
#===================================================
# look at predictions from the train and test sets
#  to see how the network did

# make the predictions on the train data
nn_prediction = predict(net, train[
  ,c("Qm3ps_lag","Precip","jdate","Jday","TMX","TMN")])

# un-scale the predictions
unscaled_train_prediction = (nn_prediction * (max(dat$Qm3ps) - min(dat$Qm3ps))) + min(dat$Qm3ps)

train_jdate = (train[,c("jdate")] * (max(dat$jdate) - min(dat$jdate))) + min(dat$jdate)
plot(train_jdate,unscaled_train_prediction)
# But, notice how predictions are created
lines(train_jdate,unscaled_train_prediction,type="l")

# make predictions on the test data
nn_prediction = predict(net, test[, c("Qm3ps_lag","Precip","jdate","Jday","TMX","TMN")])

# un-scale
unscaled_test_prediction = (nn_prediction * (max(dat$Qm3ps) - min(dat$Qm3ps))) + min(dat$Qm3ps)
test_jdate = (test[,c("jdate")] * (max(dat$jdate) - min(dat$jdate))) + min(dat$jdate)
#===================================================
#           Plot the predictions for both train and test
#===================================================
# plot the train and test and predictions
train_fit <- data.frame(NNjdate = train_jdate, NNFlow = unscaled_train_prediction, color="red",data_set = 'train')
test_fit = data.frame(NNjdate = test_jdate, NNFlow  = unscaled_test_prediction, color="green",data_set = 'test')
# Build a results dataframe for the NN
NNFlowModel=rbind(train_fit,test_fit)
NNFlowModel=NNFlowModel[order(NNFlowModel$NNjdate),]
#
# Rebuild a date from Julian Date
NNFlowModel$date=as.Date(NNFlowModel$NNjdate,origin = as.Date("1970-01-01"))

# Neatly combine results into a single dataframe.
NNFlowModel_test=subset(NNFlowModel,data_set=="test")
NNFlowModel_train=subset(NNFlowModel,data_set=="train")
plot(CNmodelnew$date, CNmodelnew$Qmm,type="l",col="black",xlab="Time", ylab="Discharge[mm/day]")
lines(NNFlowModel_train$date,NNFlowModel_train$NNFlow,col="red",type="l")
lines(NNFlowModel_test$date,NNFlowModel_test$NNFlow,col="green",type="l")
# Compare this NN Model to that of Lab12
lines(CNmodelnew$date, CNmodelnew$Qpred,type="l",col="blue")
legend("topright",legend = c("Observed Flow","NN Train Flow","NN Test Flow", "CN_Model Flow"),col = c("black","red","green","blue"),lty = 1:2, cex = 0.6)

# Calculate NSEs for the CN_Model and NNFlowModel
comparedf=merge(CNmodelnew,NNFlowModel)
NN_NSE<- NSeff(comparedf$Qmm, comparedf$NNFlow) 
print(NN_NSE)
CNmodelnew_NSE<- NSeff(CNmodelnew$Qmm, CNmodelnew$Qpred)
print(CNmodelnew_NSE)

#...............................................................................
time_1 <- c("2022-04-23","2022-04-24","2022-04-25" )
Jday1=as.Date(time_1,"%Y-%m-%d")
test1=data.frame(Qm3ps=c(0.035,0.04,0.045),
               Jday=strptime(time_1,"%Y-%m-%d")$yday+1,Precip=c(0.2,0.1,0),TMX=c(14,8,12),TMN=c(6,6,7),
               jdate=julian(Jday1))
test1$Qm3ps_lag=dplyr::lag(test1$Qm3ps)

nn_prediction1 = predict(net, test1[, c("Qm3ps_lag","Precip","jdate","Jday","TMX","TMN")])

# un-scale
unscaled_test_prediction1 = (nn_prediction1 * (max(dat$Qm3ps) - min(dat$Qm3ps))) + min(dat$Qm3ps)
test_jdate1 = (test1[,c("jdate")])
test_fit1 = data.frame(NNjdate = test_jdate1, NNFlow  = unscaled_test_prediction1, color="green",data_set = 'test')
NNFlowModel1=rbind(test_fit1)
NNFlowModel1=NNFlowModel1[order(NNFlowModel1$NNjdate),]
NNFlowModel1$date=as.Date(NNFlowModel1$NNjdate,origin = as.Date("1970-01-01"))
plot(NNFlowModel1$date,NNFlowModel1$NNFlow,xlab="Time", ylab="Discharge[mm/day]",col="green",type="l",format = "%Y-%m-%d")
