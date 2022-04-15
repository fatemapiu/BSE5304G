if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,data.table,multisensi)
# KEEP OPEN AS YOU WILL BE WALKING THROUGH IT FOR LAB	
#vignette("multisensi-vignette")
#
# Let’s get started as normal. 
setwd("~")
objects()
rm(list=objects())
#
dir.create("~/Week11")
setwd("~/Week11/")
list.files(all.files = T)
objects()
#......................................................................................................................................
#..................................................................PET function........................................................
#......................................................................................................................................

PET_fromTemp <- function (Jday, Tmax_C, Tmin_C, lat_radians, AvgT = (Tmax_C + Tmin_C)/2, albedo = 0.18, TerrestEmiss = 0.97, aspect = 0, slope = 0, forest = 0, PTconstant=1.26, AEparams=list(vp=NULL, opt="linear"))
{
  cloudiness <- EstCloudiness(Tmax_C, Tmin_C)
  DailyRad <- NetRad(lat_radians, Jday, Tmax_C, Tmin_C, albedo, forest, slope, aspect, AvgT, cloudiness, TerrestEmiss, AvgT, AEparams=AEparams)
  potentialET <- PTpet(DailyRad, AvgT, PTconstant)
  potentialET[which(potentialET < 0)] <- 0
  potentialET[which(Tmax_C == -999 | Tmin_C == -999)] <- (-999)
  return(potentialET)
}

J <- seq(from = 1, to = 365, by = 5)

PET_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- PET_fromTemp (lat=X$lat[i],
                              Jday=Jday, Tmax=X$Tmax[i], 
                              Tmin=(X$Tmax[i]-X$Trange[i]),
                              X$slope[i],X$aspect[i],X$albedo[i],X$TerrestEmiss[i])
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}




n <- 10
set.seed(1234)
X <- data.frame(Tmax_C = runif(n, min = 5, max = 30), 
                Trange = runif(n, min = 0,max = 16),
                slope = runif(n, min = 0.0, max = 0.2),
                aspect = runif(n, min = 0.0, max = 0.2),
                albedo = runif(n, min = 0.0, max = 1),
                TerrestEmiss = runif(n, min = 0.0, max = 1),
                lat_radians=runif(n, min = 0.0, max = 1.1))
View(X)

Y <- PET_Looped (X,Jday = J)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", ylab = "PET[m]")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
} 

X <- expand.grid(Tmax = c(5,15,25), 
                 Trange = c(2,9,16), 
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0),
                 albedo = c(0.1,.5,1.0),
                 TerrestEmiss=c(0.1,.5,1.0),
                 lat=c(0.1,.77,1.1))

Y <- PET_Looped(X,Jday=J) 
PET_Looped.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE) 

dev.off() # Clean up previous par()
plot(PET_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
#title(xlab = "Days of the Year.")
plot(PET_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
#title(xlab = "Days of the Year.")

PET_Looped.pca <- multisensi(model=PET_Looped, reduction=basis.ACP, scale=FALSE,
                               design.args = list( Tmax = c(5,15,25), 
                                                   Trange = c(2,9,16), 
                                                   slope = c(0.1,0.2,0.3),
                                                   aspect = c(0.1,.5,1.0),
                                                   lat=c(0.1,.77,1.1),
                                                   albedo = c(0.1,.5,1.0),
                                                   TerrestEmiss=c(0.1,.5,1.0)))

summary(PET_Looped.pca, digits = 2)
dev.off()
plot(PET_Looped.pca, graph = 1)
plot(PET_Looped.pca, graph = 2)
plot(PET_Looped.pca, graph = 3)

#.......................................................................................
##sobol2007
m <- 10000
Xb <- data.frame(Tmax = runif(m, min = 5, max = 30), 
                 Trange = runif(m, min = 2,max = 16), 
                 slope = runif(m, min = 0.0, max = 0.2),
                 aspect = runif(m, min = 0.0, max = 0.2),
                 albedo = runif(n, min = 0.0, max = 1),
                 TerrestEmiss = runif(n, min = 0.0, max = 1),
                 lat=runif(m, min = 0.0, max = 1.1))

PET_Looped.seq.sobol <- multisensi(design = sobol2007, model = PET_Looped,
                                     reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                     design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                     analysis.args = list(keep.outputs = FALSE))

print(PET_Looped.seq.sobol, digits = 2)
dev.off()
#plot(NetRad_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)
plot(PET_Looped.seq.sobol, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

#fast99
PET_Looped.seq.fast <- multisensi(design = fast99, model = PET_Looped,
                                    center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                    design.args=list( factors=c("Tmax","Trange","slope","aspect","lat",'albedo','TerrestEmiss'), 
                                                      n=1000, q = "qunif",
                                                      q.arg = list(list(min=5, max=30), 
                                                                   list(min=2, max=16),
                                                                   list(min=0, max=.2),
                                                                   list(min=0, max=.2),
                                                                   list(min = 0.0, max = 1.1),
                                                                   list(min = 0.0, max = 1),
                                                                   list(min = 0.0, max = 1))),
                                    analysis.args=list(keep.outputs=FALSE))

print(PET_Looped.seq.fast,digits=2)

plot(PET_Looped.seq.fast, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)


##..................................................................................................................................................
#...................................................................NetRad()........................................................................
#...................................................................................................................................................
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,data.table,multisensi)

objects()
rm(list=objects())
list.files(all.files = T)
objects()

#Net Rad


J <- seq(from = 1, to = 365, by = 5)

NetRad_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- NetRad(lat=X$lat[i],
                      Jday=Jday, Tx=X$Tx[i], 
                      Tn=(X$Tx[i]-X$Trange[i]), 
                      X$slope[i],X$aspect[i],units="Wm2")
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}
n <- 10
set.seed(1234)
X <- data.frame(Tx = runif(n, min = 5, max = 30), Trange = runif(n, min = 2,max = 16),
                slope = runif(n, min = 0.0, max = 0.2),
                aspect = runif(n, min = 0.0, max = 0.2),
                lat=runif(n, min = 0.0, max = 1.1))
View(X)
#
Y <- NetRad_Looped(X,Jday = J)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", 
     ylab = "Surface Short Wave Rad(W/m^2)")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
}  

NetRad_Looped.seq <- multisensi(model=NetRad_Looped, reduction=NULL, center=FALSE,
                               design.args = list( Tx = c(5,15,25), 
                                                   Trange = c(2,9,16), 
                                                   slope = c(0.1,0.2,0.3),
                                                   aspect = c(0.1,.5,1.0),
                                                   lat=c(0.1,.77,1.1)))

print(NetRad_Looped.seq, digits = 2)

dev.off() # Clean up previous par()
plot(NetRad_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
#title(xlab = "Days of the Year.")
plot(NetRad_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
#title(xlab = "Days of the Year.")

#creating different set for x
X <- expand.grid(Tx = c(5,15,25), 
                 Trange = c(2,9,16), 
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0),
                 lat=c(0.1,.77,1.1))
# Look at our input
#head(X,10)
Y <- NetRad_Looped(X,Jday=J) # can be performed outside R if necessary
# Look at our model output
#head(Y,10)
NetRad_Looped.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE) 

#PCA
NetRad_Looped.pca <- multisensi(model=NetRad_Looped, reduction=basis.ACP, scale=FALSE,
                               design.args = list( Tx = c(5,15,25), 
                                                   Trange = c(2,9,16), 
                                                   slope = c(0.1,0.2,0.3),
                                                   aspect = c(0.1,.5,1.0),
                                                   lat=c(0.1,.77,1.1)))

summary(NetRad_Looped.pca, digits = 2)
#PCA plot
dev.off()
plot(NetRad_Looped.pca, graph = 1)
plot(NetRad_Looped.pca, graph = 2)
plot(NetRad_Looped.pca, graph = 3)

#sobol2007
m <- 10000
Xb <- data.frame(Tx = runif(m, min = 5, max = 30), 
                 Trange = runif(m, min = 2,max = 16), 
                 slope = runif(m, min = 0.0, max = 0.2),
                 aspect = runif(m, min = 0.0, max = 0.2),
                 lat=runif(m, min = 0.0, max = 1.1))

NetRad_Looped.seq.sobol <- multisensi(design = sobol2007, model = NetRad_Looped,
                                     reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                     design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                     analysis.args = list(keep.outputs = FALSE))

print(NetRad_Looped.seq.sobol, digits = 2)
dev.off()
#plot(NetRad_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)
plot(NetRad_Looped.seq.sobol, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

#fast99
NetRad_Looped.seq.fast <- multisensi(design = fast99, model = NetRad_Looped,
                                    center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                    design.args=list( factors=c("Tx","Trange","slope","aspect","lat"), 
                                                      n=1000, q = "qunif",
                                                      q.arg = list(list(min=5, max=30), 
                                                                   list(min=2, max=16),
                                                                   list(min=0, max=.2),
                                                                   list(min=0, max=.2),
                                                                   list(min = 0.0, max = 1.1))),
                                    analysis.args=list(keep.outputs=FALSE))

print(NetRad_Looped.seq.fast,digits=2)
#plot(NetRad_Looped.seq.fast, normalized = TRUE, color = terrain.colors)
plot(NetRad_Looped.seq.fast, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

#................................................................................................................................
#.....................................................................grad........................................................
#................................................................................................................................

if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,data.table,multisensi)

objects()
rm(list=objects())
list.files(all.files = T)
objects()

#SoilStorage(S_avg, field_capacity, soil_water_content, porosity)
SoilStorage(S_avg=120, field_capacity=0.2, soil_water_content=0.1, porosity=0.3)

J <- seq(from = 10, to = 100, by = 5)

Soil_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- SoilStorage(X$S_avg[i],
                      X$field_capacity[i],
                      X$soil_water_content[i],X$porosity[i])
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}
n <- 10
set.seed(1234)
X <- data.frame(S_avg = runif(n, min = 10, max = 400),
                field_capacity = runif(n, min = 0.0, max = 1),
                soil_water_content= runif(n, min = 0.0, max = 1),
                porosity=runif(n, min = 0.0, max = 1))

Soil_Looped.seq <- multisensi(model=Soil_Looped, reduction=NULL, center=FALSE,
                               design.args = list( S_avg = runif(n, min = 150, max = 350),
                                                   field_capacity = runif(n, min = 0.05, max = 0.3),
                                                   soil_water_content= runif(n, min = 0.05, max =0.49),
                                                   porosity=runif(n, min = 0.05, max =0.49)))

print(Soil_Looped.seq, digits = 2)
#S_avg=150; porosity=.05; log(1 - (porosity/(1 - 2.54/(2.381 * S_avg))) - porosity)


dev.off() # Clean up previous par()
plot(Soil_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
#title(xlab = "Days of the Year.")
plot(Soil_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
#title(xlab = "Days of the Year.")



#PCA
Soil_Looped.pca <- multisensi(model=Soil_Looped, reduction=basis.ACP, scale=FALSE,
                               design.args = list(S_avg = runif(n, min = 150, max = 350),
                                                  field_capacity = runif(n, min = 0.05, max = 0.3),
                                                  soil_water_content= runif(n, min = 0.05, max =0.49),
                                                  porosity=runif(n, min = 0.05, max =0.49)))

summary(Soil_Looped.pca, digits = 2)
#PCA plot
dev.off()
plot(Soil_Looped.pca, graph = 1)
plot(Soil_Looped.pca, graph = 2)
plot(Soil_Looped.pca, graph = 3)

dev.off()

#fast99
Soil_Looped.seq.fast <- multisensi(design = fast99, model = Soil_Looped,
                                    center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                    design.args=list( factors=c("S_avg","field_capacity ","soil_water_content","porosity"), 
                                                      n=1000, q = "qunif",
                                                      q.arg = list(list(min=150, max=350), 
                                                                   list(min=0.05, max=0.3),
                                                                   list(min=0.05, max=.49),
                                                                   list(min=0.05, max=.49))),
                                    analysis.args=list(keep.outputs=FALSE))

print(Soil_Looped.seq.fast,digits=2)
#plot(Solar_Looped.seq.fast, normalized = TRUE, color = terrain.colors)
plot(Soil_Looped.seq.fast, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)


