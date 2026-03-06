########Create the kriging model for various lafe stages of ribbonfish#########
##http://rgeostats.free.fr/demo/RGeostats.acoustic.maur.R
#Set the working directionary
setwd("E:/Bycatchmapping/R")
Grid_monsoon<- read.csv("Grid_monsoon.csv")
Grid_winter<- read.csv("Grid_winter.csv")
Grid_summer<- read.csv("Grid_summer.csv")
head(Grid_monsoon)
hist(Grid_monsoon$Ribbon_medium)

#install.packages("RgeoS")
library(rgeos)
vgm()
#filter seasonal data & data cleaning
library(dplyr)
#remove outlier
#outlier removed from ribbonfish catch in ribbon.tot data frame
library(data.table)
outlierReplace = function(dataframe, cols, rows, newValue = NA) {
  if (any(rows)) {
    set(dataframe, rows, cols, newValue)
  }
}
#outlierReplace(summer.b, "Total_other", which(summer.b$Total_other>70), NA)
#boxplot(summer.b$Total_other)
#summer.b<- na.omit(summer.b)

#############Post-monsoon########################
monsoon<-Grid_monsoon%>%dplyr::select(c(Longitude,Latitude,Total_robbon,U.duvaucelii, S.eliptica,N.japanicus,
                                 Crocker,S.pharaonis,P.hamrur,S.inermis,S.tumbil))

winter<-Grid_winter%>%dplyr::select(c(Longitude,Latitude,Total_robbon,U.duvaucelii, S.eliptica,N.japanicus,
                                 Crocker,S.pharaonis,P.hamrur,S.inermis,S.tumbil))
summer<-Grid_summer%>%dplyr::select(c(Longitude,Latitude,Total_robbon,U.duvaucelii, S.eliptica,N.japanicus,
                                      Crocker,S.pharaonis,P.hamrur,S.inermis,S.tumbil))
############Box-Cox transformation############################
# install.packages(MASS)
library(MASS)

########ribbon
#monsoon
a<-boxcox(lm(Total_robbon ~ 1, data = monsoon))
# Transformed data
new_robbon <- sqrt(monsoon$Total_robbon)
# Histogram
hist(new_robbon)
shapiro.test(new_robbon)

lambda <- a$x[which.max(a$y)] 
lambda
new_robbon_exact <- (new_robbon ^ lambda - 1) / lambda
hist(new_robbon_exact)

#winter
a0<-boxcox(lm(Total_robbon ~ 1, data = winter))
# Transformed data
new_robbon0 <- sqrt(winter$Total_robbon)
# Histogram
hist(new_robbon0)
shapiro.test(new_robbon)

lambda <- a0$x[which.max(a0$y)] 
lambda
new_robbon_exact0 <- (new_robbon0 ^ lambda - 1) / lambda
hist(new_robbon_exact0)

#summer
a000<-boxcox(lm(Total_robbon ~ 1, data = summer))
# Transformed data
new_robbon000 <- log(summer$Total_robbon)
# Histogram
hist(new_robbon000)
shapiro.test(new_robbon000)

lambda <- a000$x[which.max(a000$y)] 
lambda
new_robbon_exact000 <- (new_robbon000 ^ lambda - 1) / lambda
hist(new_robbon_exact000)

########U.duvaucelii
#monsoon
a1<-boxcox(lm(U.duvaucelii ~ 1, data = monsoon))
# Transformed data
new_U.duvaucelii <- sqrt(monsoon$U.duvaucelii)
# Histogram
hist(new_U.duvaucelii)
shapiro.test(new_U.duvaucelii_exact)

lambda <- a1$x[which.max(a1$y)] 
lambda
new_U.duvaucelii_exact <- (new_U.duvaucelii ^ lambda - 1) / lambda
hist(new_U.duvaucelii_exact)

#winter
U.duvaucelii11<-winter$U.duvaucelii+2.12
a11<-boxcox(lm(U.duvaucelii11 ~ 1, data = winter))
# Transformed data
new_U.duvaucelii11 <- sqrt(U.duvaucelii11)
# Histogram
hist(new_U.duvaucelii11)
shapiro.test(new_U.duvaucelii_exact)

lambda <- a11$x[which.max(a11$y)] 
lambda
new_U.duvaucelii_exact11 <- (new_U.duvaucelii11 ^ lambda - 1) / lambda
hist(new_U.duvaucelii_exact11)

#summer
U.duvaucelii111<-summer$U.duvaucelii+0.98
a111<-boxcox(lm(U.duvaucelii111 ~ 1, data = summer))
# Transformed data
new_U.duvaucelii111 <- sqrt(U.duvaucelii111)
# Histogram
hist(new_U.duvaucelii111)
shapiro.test(new_U.duvaucelii111)

lambda <- a111$x[which.max(a111$y)] 
lambda
new_U.duvaucelii_exact111 <- (new_U.duvaucelii111 ^ lambda - 1) / lambda
hist(new_U.duvaucelii_exact111)

######S.eliptica
#monsoon
S.eliptica1<- monsoon$S.eliptica+0.09
a2<-boxcox(lm(S.eliptica1 ~ 1, data = monsoon))
# Transformed data
new_S.eliptica <- sqrt(S.eliptica1)
# Histogram
hist(new_S.eliptica)
shapiro.test(new_S.eliptica)

lambda <- a2$x[which.max(a2$y)] 
lambda
new_S.eliptica_exact <- (new_S.eliptica ^ lambda - 1) / lambda
hist(new_S.eliptica_exact)

#winter
S.eliptica11<- winter$S.eliptica+0.30
a22<-boxcox(lm(S.eliptica11 ~ 1, data = winter))
# Transformed data
new_S.eliptica11 <- sqrt(S.eliptica11)
# Histogram
hist(new_S.eliptica11)
shapiro.test(new_S.eliptica)

lambda <- a22$x[which.max(a22$y)] 
lambda
new_S.eliptica_exact22 <- (new_S.eliptica11 ^ lambda - 1) / lambda
hist(new_S.eliptica_exact22)

#summer
S.eliptica111<- summer$S.eliptica+0.50
a222<-boxcox(lm(S.eliptica111 ~ 1, data = summer))
# Transformed data
new_S.eliptica222 <- sqrt(S.eliptica111)
# Histogram
hist(new_S.eliptica222)
shapiro.test(new_S.eliptica222)

lambda <- a222$x[which.max(a222$y)] 
lambda
new_S.eliptica_exact222 <- (new_S.eliptica222 ^ lambda - 1) / lambda
hist(new_S.eliptica_exact222)

#######N.japanicus
#monsoon
N.japanicus1<- monsoon$N.japanicus+0.22
a3<-boxcox(lm(N.japanicus1 ~ 1, data = monsoon))
# Transformed data
new_N.japanicus <- sqrt(N.japanicus1)
# Histogram
hist(new_N.japanicus)
shapiro.test(new_N.japanicus)

lambda <- a3$x[which.max(a3$y)] 
lambda
new_N.japanicus_exact <- (new_N.japanicus ^ lambda - 1) / lambda
hist(new_N.japanicus_exact)

#winter
N.japanicus11<- winter$N.japanicus+1.11
a33<-boxcox(lm(N.japanicus11 ~ 1, data = winter))
# Transformed data
new_N.japanicus33 <- sqrt(N.japanicus11)
# Histogram
hist(new_N.japanicus33)
shapiro.test(new_N.japanicus33)

lambda <- a33$x[which.max(a33$y)] 
lambda
new_N.japanicus_exact33 <- (new_N.japanicus33 ^ lambda - 1) / lambda
hist(new_N.japanicus_exact33)

#summer
N.japanicus111<- summer$N.japanicus+0.29
a333<-boxcox(lm(N.japanicus111 ~ 1, data = summer))
# Transformed data
new_N.japanicus333 <- sqrt(N.japanicus111)
# Histogram
hist(new_N.japanicus333)
shapiro.test(new_N.japanicus333)

lambda <- a333$x[which.max(a333$y)] 
lambda
new_N.japanicus_exact333 <- (new_N.japanicus333 ^ lambda - 1) / lambda
hist(new_N.japanicus_exact333)

########Crocker
#monsoon
Crocker1<- monsoon$Crocker+0.38
a4<-boxcox(lm(Crocker1 ~ 1, data = monsoon))
# Transformed data
new_Crocker <- sqrt(Crocker1)
# Histogram
hist(new_Crocker)
shapiro.test()

lambda <- a4$x[which.max(a4$y)] 
lambda
new_Crocker_exact <- (new_Crocker ^ lambda - 1) / lambda
hist(new_Crocker_exact)

#winter
Crocker11<- winter$Crocker+0.58
a44<-boxcox(lm(Crocker11 ~ 1, data = winter))
# Transformed data
new_Crocker44 <- sqrt(Crocker11)
# Histogram
hist(new_Crocker44)
shapiro.test(new_Crocker44)

lambda <- a44$x[which.max(a44$y)] 
lambda
new_Crocker_exact44 <- (new_Crocker44 ^ lambda - 1) / lambda
hist(new_Crocker_exact44)

#summer
Crocker111<- summer$Crocker+1.09
a444<-boxcox(lm(Crocker111 ~ 1, data = summer))
# Transformed data
new_Crocker444 <- sqrt(Crocker111)
# Histogram
hist(new_Crocker444)
shapiro.test(new_Crocker444)

lambda <- a444$x[which.max(a444$y)] 
lambda
new_Crocker_exact444 <- (new_Crocker444 ^ lambda - 1) / lambda
hist(new_Crocker_exact444)


#############S.pharaonis
#monsoon
S.pharaonis1<- monsoon$S.pharaonis+0.08
a5<-boxcox(lm(S.pharaonis1 ~ 1, data = monsoon))
# Transformed data
new_S.pharaonis <- sqrt(S.pharaonis1)
# Histogram
hist(new_S.pharaonis)
shapiro.test(new_S.pharaonis)

lambda <- a5$x[which.max(a5$y)] 
lambda
new_S.pharaonis_exact <- (new_S.pharaonis ^ lambda - 1) / lambda
hist(new_S.pharaonis_exact)

#winter
S.pharaonis11<- winter$S.pharaonis+0.25
a55<-boxcox(lm(S.pharaonis11 ~ 1, data = winter))
# Transformed data
new_S.pharaonis55 <- sqrt(S.pharaonis11)
# Histogram
hist(new_S.pharaonis55)
shapiro.test(new_S.pharaonis55)

lambda <- a55$x[which.max(a55$y)] 
lambda
new_S.pharaonis_exact55 <- (new_S.pharaonis55 ^ lambda - 1) / lambda
hist(new_S.pharaonis_exact55)

#summer
S.pharaonis111<- summer$S.pharaonis+0.23
a555<-boxcox(lm(S.pharaonis111 ~ 1, data = summer))
# Transformed data
new_S.pharaonis555 <- sqrt(S.pharaonis111)
# Histogram
hist(new_S.pharaonis555)
shapiro.test(new_S.pharaonis555)

lambda <- a555$x[which.max(a555$y)] 
lambda
new_S.pharaonis_exact555 <- (new_S.pharaonis555 ^ lambda - 1) / lambda
hist(new_S.pharaonis_exact555)

################P.hamrur
#monsoon
P.hamrur1<- monsoon$P.hamrur+0.42
a6<-boxcox(lm(P.hamrur1 ~ 1, data = monsoon))
# Transformed data
new_P.hamrur <- sqrt(P.hamrur1)
# Histogram
hist(new_P.hamrur)
shapiro.test(new_P.hamrur)

lambda <- a6$x[which.max(a6$y)] 
lambda
new_P.hamrur_exact <- (new_P.hamrur ^ lambda - 1) / lambda
hist(new_P.hamrur_exact)

#winter
P.hamrur11<- winter$P.hamrur+0.29
a66<-boxcox(lm(P.hamrur11 ~ 1, data = winter))
# Transformed data
new_P.hamrur66 <- sqrt(P.hamrur11)
# Histogram
hist(new_P.hamrur66)
shapiro.test(new_P.hamrur66)

lambda <- a66$x[which.max(a66$y)] 
lambda
new_P.hamrur_exact66 <- (new_P.hamrur66 ^ lambda - 1) / lambda
hist(new_P.hamrur_exact66)

#summer
P.hamrur111<- summer$P.hamrur+0.27
a666<-boxcox(lm(P.hamrur111 ~ 1, data = summer))
# Transformed data
new_P.hamrur666 <- sqrt(P.hamrur111)
# Histogram
hist(new_P.hamrur666)
shapiro.test(new_P.hamrur666)

lambda <- a666$x[which.max(a666$y)] 
lambda
new_P.hamrur_exact666 <- (new_P.hamrur666 ^ lambda - 1) / lambda
hist(new_P.hamrur_exact666)

################S.inermis
#monsoon
S.inermis1<- monsoon$S.inermis+0.21
a7<-boxcox(lm(S.inermis1 ~ 1, data = monsoon))
# Transformed data
new_S.inermis <- sqrt(S.inermis1)
# Histogram
hist(new_S.inermis)
shapiro.test(new_S.inermis)

lambda <- a7$x[which.max(a7$y)] 
lambda
new_S.inermis_exact <- (new_S.inermis ^ lambda - 1) / lambda
hist(new_S.inermis_exact)

#winter
S.inermis11<- winter$S.inermis+0.12
a77<-boxcox(lm(S.inermis11 ~ 1, data = winter))
# Transformed data
new_S.inermis77 <- sqrt(S.inermis11)
# Histogram
hist(new_S.inermis77)
shapiro.test(new_S.inermis77)

lambda <- a77$x[which.max(a77$y)] 
lambda
new_S.inermis_exact77 <- (new_S.inermis77 ^ lambda - 1) / lambda
hist(new_S.inermis_exact77)

#summer
S.inermis111<- summer$S.inermis+0.15
a777<-boxcox(lm(S.inermis111 ~ 1, data = summer))
# Transformed data
new_S.inermis777 <- sqrt(S.inermis111)
# Histogram
hist(new_S.inermis777)
shapiro.test(new_S.inermis777)

lambda <- a777$x[which.max(a777$y)] 
lambda
new_S.inermis_exact777 <- (new_S.inermis777 ^ lambda - 1) / lambda
hist(new_S.inermis_exact777)

###########S.tumbil
#monsoon
S.tumbil1<- monsoon$S.tumbil+0.32
a8<-boxcox(lm(S.tumbil1 ~ 1, data = monsoon))
# Transformed data
new_S.tumbil <- sqrt(S.tumbil1)
# Histogram
hist(new_S.tumbil)
shapiro.test(new_S.tumbil)

lambda <- a8$x[which.max(a8$y)] 
lambda
new_S.tumbil_exact <- (new_S.tumbil ^ lambda - 1) / lambda
hist(new_S.tumbil_exact)

#winter
S.tumbil11<- winter$S.tumbil+0.42
a88<-boxcox(lm(S.tumbil11 ~ 1, data = winter))
# Transformed data
new_S.tumbil88 <- sqrt(S.tumbil11)
# Histogram
hist(new_S.tumbil88)
shapiro.test(new_S.tumbil88)

lambda <- a88$x[which.max(a88$y)] 
lambda
new_S.tumbil_exact88 <- (new_S.tumbil88 ^ lambda - 1) / lambda
hist(new_S.tumbil_exact88)

#summer
S.tumbil111<- summer$S.tumbil+0.17
a888<-boxcox(lm(S.tumbil111 ~ 1, data = summer))
# Transformed data
new_S.tumbil888 <- sqrt(S.tumbil111)
# Histogram
hist(new_S.tumbil888)
shapiro.test(new_S.tumbil888)

lambda <- a888$x[which.max(a888$y)] 
lambda
new_S.tumbil_exact888 <- (new_S.tumbil888 ^ lambda - 1) / lambda
hist(new_S.tumbil_exact888)


#####save data in csv###########
monsoon_box<- data.frame(monsoon$Longitude,monsoon$Latitude,new_robbon_exact,new_U.duvaucelii_exact,new_S.eliptica_exact,new_N.japanicus_exact,new_Crocker_exact,
                         new_S.pharaonis_exact,new_P.hamrur_exact,new_S.inermis_exact,new_S.tumbil_exact)
View(monsoon_box)
colnames(monsoon_box) <- c("Longitude","Latitude","T.lepturus","U.duvaucelii","S.elliptica","N.japonicus","O.cuveri",
                           "S.pharaonis","P.hamrur","S.inermis","S.tumbil")
#write.csv(monsoon_box,"monsoon_boxcox.csv")

wintre_box<- data.frame(winter$Longitude,winter$Latitude,new_robbon_exact0,new_U.duvaucelii_exact11,
                        new_S.eliptica_exact22,new_N.japanicus_exact33,new_Crocker_exact44,
                        new_S.pharaonis_exact55,new_P.hamrur_exact66,new_S.inermis_exact77,
                        new_S.tumbil_exact88)
colnames(wintre_box) <- c("Longitude","Latitude","T.lepturus","U.duvaucelii","S.elliptica","N.japonicus","O.cuveri",
                           "S.pharaonis","P.hamrur","S.inermis","S.tumbil")
#write.csv(wintre_box,"winter_boxcox.csv")

summer_box<- data.frame(summer$Longitude,summer$Latitude,new_robbon_exact000,new_U.duvaucelii_exact111,
                        new_S.eliptica_exact222,new_N.japanicus_exact333,new_Crocker_exact444,
                        new_S.pharaonis_exact555,new_P.hamrur_exact666,new_S.inermis_exact777,
                        new_S.tumbil_exact888)

colnames(summer_box) <- c("Longitude","Latitude","T.lepturus","U.duvaucelii","S.elliptica","N.japonicus","O.cuveri",
                          "S.pharaonis","P.hamrur","S.inermis","S.tumbil")
#write.csv(summer_box,"summer_boxcox.csv")

######combine seaonal data##########

#install.packages("Jmisc")
library(Jmisc)
addCol <-
  function(x,...,value){
    arguments<-list(...)
    
    if (length(arguments)>0){
      if (!missing(value))
        arguments<-c(arguments,value)
      return(addCol(x,value=arguments))
    }
    
    if (any(sapply(value,length)!=1))
      stop("Element to be added must be singleton.")
    
    cbind(x,repRow(value,nrow(x)))
  }

monsoon_box1<-addCol(monsoon_box,season=1000)
wintre_box1<-addCol(wintre_box, season=2000)
summer_box1<-addCol(summer_box,season=3000)
dim(monsoon_box1)

Boxcox_data<-rbind(monsoon_box1,wintre_box1)
Boxcox_data<-rbind(Boxcox_data,summer_box1)
#save Boxcox_data
#write.csv(Boxcox_data,"Boxcox_data.csv")


#box-cox using geoR
library(geoR)
T.box = boxcoxfit(S.eliptica1, lambda2=TRUE)
T.box
lambda = T.box$lambda[1]
lambda2 = T.box$lambda[2]
if(lambda==0){T.norm = log(S.eliptica1 + lambda2)}
if(lambda!=0){T.norm = (S.eliptica1 ^ lambda - 1) / lambda}
hist(T.norm, col="gray")
View(T.norm)


#Setup boarder file by extracting study area shapefile to csv
#first conver the line to points and points are saved in csv
library(sf)
AOI <- st_read("D:/Bycatch/studyarea/Study_area.shp")
pnts<-st_cast(AOI, "POINT")
#st_write(pnts,"AOI.csv",layer_options="GEOMETRY=AS_XY")
AOI<- read.csv("E:/Bycatchmapping/R/AOI.csv")

pnts<-st_cast(india, "POINT")
#st_write(pnts,"E:/Bycatchmapping/R/INDIA.csv",layer_options="GEOMETRY=AS_XY")

india<-st_read("E:/Bycatchmapping/R/Coast and depth contours/Coast.shp")
m50 <- st_read("E:/Bycatchmapping/R/Coast and depth contours/50m_contour.shp")
m100<-st_read("E:/Bycatchmapping/R/Coast and depth contours/100m_Contour.shp")
m200<-st_read("E:/Bycatchmapping/R/Coast and depth contours/200m_Contour.shp")

india<-st_geometry_type(india)
st_crs(india)
box= c(xmin=67, xmax=74, ymin=15.6, ymax=23)
plot(st_crop(india,box),col = "red",lwd=2)

plot(india, lwd=2,ylim= c(15,24))
plot(m50, add=TRUE,col = "red", lwd=2,ylim= c(15,24))
plot(m100,add=TRUE,col = "yellow", lwd=2,ylim= c(15,24))
plot(m200, add=TRUE,col = "green", lwd=2,ylim= c(15,24))

#####Install geoR for variogram analysis
#install.packages("geoR")
library(geoR)
#install.packages(c("Hmisc", "fields","fBasics","plot3D","viridis"))
library(Hmisc)
library(fields)
library(fBasics)
library(plot3D)
library(viridis)

setwd("E:/Bycatchmapping/monsoon1")
#Create a geodata object
View(monsoon_box)
tumbil.w.geo <- as.geodata(wintre_box,
                            coords.col=1:2, data.col=11,
                            rep.data.action = "first")
summary(elliptica.s.geo)
plot()

# Initial name to save files
nome_analise = 'S.tumbil'
X= monsoon_box$Longitude
Y= monsoon_box$Latitude
Z= monsoon_box$T.lepturus
#  Statistical analysis ----------------------------------------------------
estb=basicStats(Z, ci=0.95); estb
tks=ks.test(Z, 'pnorm', mean=mean(Z), sd=sd(Z)); tks
shapiro.test(Z)
jarqueberaTest(Z)
boxplot(Z)

# Starting geostatistics ------------------------------------------------
# Maximum distance
mx=max(X)
my=max(Y)
mdta=sqrt(mx^2+my^2);mdta; mdta/3

# Generating and plotting the semivariogram -------------------------------------
plot(variog(ribbon.m.geo))
mdt=10; nlg=9 # maximum distance and number of logs (no. points come within said distance mdt)
svgteorico=variog(ribbon.m.geo, max.dist=mdt, uvec=nlg)
svgteorico1=variog(ribbon.m.geo, max.dist=mdt, uvec=nlg, direction=pi/8)
svgteorico2=variog(ribbon.m.geo, max.dist=mdt, uvec=nlg, direction=pi/4)
svgteorico3=variog(ribbon.m.geo, max.dist=mdt, uvec=nlg, direction=pi/2)
par(mfrow=c(2,2))
plot(svgteorico); plot(svgteorico1, main=expression(pi/8)); plot(svgteorico2, main=expression(pi/4)); plot(svgteorico3, main=expression(pi/2)); layout(1)

h=svgteorico$u
v=svgteorico$v
npar=svgteorico$n
ltmx=(max(h)+0.4*max(h))
ltmy=(max(v)+0.4*max(v))
plot(svgteorico, las=1, xaxs='i', yaxs='i',pch=16, col='red', ylim=c(0,ltmy), xlim=c(0,ltmx))
textxy(h,v,npar, cex=0.7)
sm=data.frame(h,v,npar, row.names=NULL);sm

# Adjusting the theoretical semivariogram --------------------------------
esferico=variofit(svgteorico, cov.model='sph', max.dist=max(h), messages=F)
exponencial=variofit(svgteorico, cov.model='exponential', max.dist=max(h), messages=F)
gaussiano=variofit(svgteorico,cov.model='gaussian', max.dist=max(h),messages=F)
X11(); sentimento=eyefit(svgteorico,silent=F)
stp = unlist(sentimento)
sigmasq=as.numeric(stp[2]); phi=as.numeric(stp[3]); tausq=as.numeric(stp[4])



# Plotting all the semivariogram and its adjustments --------------------------gaussiano=variofit(svgteorico, cov.model='gaussian', max.dist=max(h),messages=F)
#jpeg("Semivariograms_monsoon.m.jpg", res = 600,height = 7,width = 8,units = "in")
par(mfrow=c(2,2),mar=c(5,5,2,2))
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Spherical')
lines.variomodel(esferico, col='red' , lwd=2, lty=1)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Exponential')
lines.variomodel(exponencial, col='red' , lwd=2, lty=1)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Gaussian')
lines.variomodel(gaussiano, col='red' , lwd=2, lty=1)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Manual')
lines(sentimento, col='red' , lwd=2, lty=1)
#dev.copy(pdf, paste(nome_analise, '_Semivariogramas.pdf', sep = ''),  width=8, height=6)
dev.off()
layout(1)

#jpeg("Semivariograms_winter.jpg", res = 600,height = 7,width = 8,units = "in")
#par(mfrow=c(1,1),mar=c(5,5,2,2))
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black', ylim=c(0,ltmy), xlim=c(0,9),
     ylab='Variance', xlab='Distance (degree)', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5, 
     main=substitute(paste(italic('Saurida tumbil'))))
lines(sentimento, col='red' , lwd=2, lty=1)
dev.off()

#vario.mon.ribbon<-svgteorico

# Parameters adjusted for models ------------------------------------ 
c0.esf=esferico$nugget; c.esf=esferico$cov.pars[1]; a.esf=esferico$cov.pars[2]
c0.exp=exponencial$nugget; c.exp=exponencial$cov.pars[1]; a.exp=exponencial$cov.pars[2]
c0.gau=gaussiano$nugget; c.gau=gaussiano$cov.pars[1]; a.gau=gaussiano$cov.pars[2]
C0.snt = sigmasq; C.snt = tausq; a.snt = phi
ide.esf = 100*(c0.esf/(c0.esf+c.esf)); ide.exp = 100*(c0.exp/(c0.exp+c.exp)); ide.gau = 100*(c0.gau/(c0.gau+c.gau)); ide.sent = 100*(C.snt/(C0.snt+C.snt))


# Cross-validation of adjustments -------------------------------------------
#change data
cv.esf=xvalid(tumbil.w.geo, model=esferico); zsco.esf=cv.esf$std.error
cv.exp=xvalid(tumbil.w.geo, model=exponencial); zsco.exp=cv.exp$std.error
cv.gau=xvalid(tumbil.w.geo, model=gaussiano); zsco.gau=cv.gau$std.error
cv.sent=xvalid(tumbil.w.geo, model=sentimento); zsco.sent=cv.sent$std.error

# Reduced error mean and variance --------------------------------------
jkmed.esf=round(mean(zsco.esf),5); jkvar.esf=round(var(zsco.esf),5)
jkmed.exp=round(mean(zsco.exp),5); jkvar.exp=round(var(zsco.exp),5)
jkmed.gau=round(mean(zsco.gau),5); jkvar.gau=round(var(zsco.gau),5)
jkmed.sent=round(mean(zsco.sent),5); jkvar.sent=round(var(zsco.sent),5)


# R2 (R square) kriging-------------------------------------------------------------
r2k.esf = cor(cv.esf$predicted,cv.esf$data)^2
r2k.exp = cor(cv.exp$predicted,cv.exp$data)^2
r2k.gau = cor(cv.gau$predicted,cv.gau$data)^2
r2k.sent = cor(cv.sent$predicted,cv.sent$data)^2


# Summary of analysis to choose the best model -----------------------
modelos=c('Esf', 'Exp','Gau','Sent')
m.jk=rbind(jkmed.esf,jkmed.exp, jkmed.gau, jkmed.sent)
v.jk=rbind(jkvar.esf, jkvar.exp, jkvar.gau, jkvar.sent)
c0.smfit=rbind(c0.esf, c0.exp, c0.gau, tausq)
c.smfit=rbind(c.esf, c.exp, c.gau, sigmasq)
a.smfit=rbind(a.esf, a.exp, a.gau, phi)
ide = rbind(ide.esf, ide.exp, ide.gau, ide.sent)
r2k = rbind(r2k.esf, r2k.exp, r2k.gau, r2k.sent)

resumo=data.frame(row.names=modelos, c0.smfit, c.smfit, a.smfit, ide, m.jk, v.jk, r2k)
resumo
#write.table(resumo, 'winter.model selection.csv', sep = 'T')

# Choosing the best model ------------------------------------------------
# esferico; exponencial; gaussiano; sentimento
smfit = sentimento

# Generating the grid for interpolations ---------------------------------------
ndiv = 400 # Interpolation interval size
x.range <- as.integer(c(66,74))
y.range <- as.integer(c(15,24))
grid.map=expand.grid(x=seq(from=x.range[1], to=x.range[2], by= (x.range[2] - x.range[1])/ndiv),
                     y=seq(from=y.range[1], to=y.range[2], by=(y.range[2] - y.range[1])/ndiv))
#plot(grid.map)
grid.map
# Loading border boundaries with contour line in text form ---------------------------------------------
lmt=read.table('E:/Bycatchmapping/R/AOI.txt', h=T)
dlmt=read.geodata('E:/Bycatchmapping/R/AOI.txt', h=T, coords.col=1:2, data.col=NULL,rep.data.action = "none")

# Doing Kriging --------------------------------------------------------
krg=krige.conv(tumbil.w.geo, locations=grid.map, krige=krige.control(obj.model=smfit), borders=dlmt)

# Contour chart ----------------------------------------------------
contour(krg, f=T, col=viridis(14), nlevels=10)

# Settings for saving the semivariogram ------------------------------
#jpeg("Semivario_winter.jpg", res = 600,height = 4,width = 5,units = "in")
par(family='', las=1, xaxs='i', yaxs='i', mar=c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
plot(svgteorico, type='p',pch=19,col='black' , ylim=c(0,(max(sm$v)+max(sm$v)*0.15)), 
     xlim=c(0,(max(sm$h)+max(sm$h)*0.15)),main=substitute(paste(italic('Saurida tumbil'))), 
     ylab='', xlab='Distance (degree)')#, yaxt='n')
#vr=var(Z)
abline(v=NULL, #h=vr, 
       lty=2, lwd=1, untf=3)
lines(smfit, col='black' , lwd=2, lty=1)
mtext(text=expression("Semivariance ("*gamma*")"), side=2, line=3, las=3, cex=1.5)
#dev.copy(pdf, paste(nome_analise, 'semivari_monsoon.s.pdf', sep = ''), width=8, height=6)
dev.off()

####krigging name
#mon.ribbon.krig<-krg
#mon.ribbon.vari<-svgteorico
#mon.ribbon.fit<-smfit
#win.ribbon.krig<-krg
#win.ribbon.vari<-svgteorico
#win.ribbon.fit<-smfit
#sum.ribbon.krig<-krg
#sum.ribbon.vari<-svgteorico
#sum.ribbon.fit<-smfit
#mon.duv.krig<-krg
#mon.duv.vari<-svgteorico
#mon.duv.fit<-smfit
#win.duv.krig<-krg
#win.duv.vari<-svgteorico
#win.duv.fit<-smfit
#sum.duv.krig<-krg
#sum.duv.vari<-svgteorico
#sum.duv.fit<-smfit
#mon.ell.krig<-krg
#mon.ell.vari<-svgteorico
#mon.ell.fit<-smfit
#win.ell.krig<-krg
#win.ell.vari<-svgteorico
#win.ell.fit<-smfit
#sum.ell.krig<-krg
#sum.ell.vari<-svgteorico
#sum.ell.fit<-smfit
#mon.jap.krig<-krg
#mon.jap.vari<-svgteorico
#mon.jap.fit<-smfit
#win.jap.krig<-krg
#win.jap.vari<-svgteorico
#win.jap.fit<-smfit
#sum.jap.krig<-krg
#sum.jap.vari<-svgteorico
#sum.jap.fit<-smfit
#mon.cuv.krig<-krg
#mon.cuv.vari<-svgteorico
#mon.cuv.fit<-smfit
#win.cuv.krig<-krg
#win.cuv.vari<-svgteorico
#win.cuv.fit<-smfit
#sum.cuv.krig<-krg
#sum.cuv.vari<-svgteorico
#sum.cuv.fit<-smfit
#mon.pha.krig<-krg
#mon.pha.vari<-svgteorico
#mon.pha.fit<-smfit
#win.pha.krig<-krg
#win.pha.vari<-svgteorico
#win.pha.fit<-smfit
#sum.pha.krig<-krg
#sum.pha.vari<-svgteorico
#sum.pha.fit<-smfit
#mon.ham.krig<-krg
#mon.ham.vari<-svgteorico
#mon.ham.fit<-smfit
#win.ham.krig<-krg
#win.ham.vari<-svgteorico
#win.ham.fit<-smfit
#sum.ham.krig<-krg
#sum.ham.vari<-svgteorico
#sum.ham.fit<-smfit
#mon.ine.krig<-krg
#mon.ine.vari<-svgteorico
#mon.ine.fit<-smfit
#win.ine.krig<-krg
#win.ine.vari<-svgteorico
#win.ine.fit<-smfit
#sum.ine.krig<-krg
#sum.ine.vari<-svgteorico
#sum.ine.fit<-smfit
#mon.tum.krig<-krg
#mon.tum.vari<-svgteorico
#mon.tum.fit<-smfit
#win.tum.krig<-krg
#win.tum.vari<-svgteorico
#win.tum.fit<-smfit
#sum.tum.krig<-krg
#sum.tum.vari<-svgteorico
#sum.tum.fit<-smfit


##Ploting all the variograms in one layout
jpeg("All_variograms fit.jpg", res = 600,height = 7,width = 7,units = "in")
par(mfrow=c(3,3),mar=c(5,5,2,2))
plot(mon.ribbon.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='Variance', xlab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Trichiurus lepturus'))))
lines(mon.ribbon.fit, col='red' , lwd=2, lty=1)
lines(win.ribbon.fit, col='green' , lwd=2, lty=1)
lines(sum.ribbon.fit, col='blue' , lwd=2, lty=1)

plot(mon.duv.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='', xlab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Uroteuthis duvaucelii'))))
lines(mon.duv.fit, col='red' , lwd=2, lty=1)
lines(win.duv.fit, col='green' , lwd=2, lty=1)
lines(sum.duv.fit, col='blue' , lwd=2, lty=1)

plot(mon.ell.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='', xlab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Sepia elliptica'))))
lines(mon.ell.fit, col='red' , lwd=2, lty=1)
lines(win.ell.fit, col='green' , lwd=2, lty=1)
lines(sum.ell.fit, col='blue' , lwd=2, lty=1)

plot(mon.jap.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='Variance', xlab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Nemipterus japonicus'))))
lines(mon.jap.fit, col='red' , lwd=2, lty=1)
lines(win.jap.fit, col='green' , lwd=2, lty=1)
lines(sum.jap.fit, col='blue' , lwd=2, lty=1)

plot(mon.cuv.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='', xlab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Otolithes cuvieri'))))
lines(mon.cuv.fit, col='red' , lwd=2, lty=1)
lines(win.cuv.fit, col='green' , lwd=2, lty=1)
lines(sum.cuv.fit, col='blue' , lwd=2, lty=1)

plot(mon.pha.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='', xlab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Sepia pharaonis'))))
lines(mon.pha.fit, col='red' , lwd=2, lty=1)
lines(win.pha.fit, col='green' , lwd=2, lty=1)
lines(sum.pha.fit, col='blue' , lwd=2, lty=1)

plot(mon.ham.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='Variance', xlab='Distance (degree)', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Priacanthus hamrur'))))
lines(mon.ham.fit, col='red' , lwd=2, lty=1)
lines(win.ham.fit, col='green' , lwd=2, lty=1)
lines(sum.ham.fit, col='blue' , lwd=2, lty=1)

plot(mon.ine.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='', xlab='Distance (degree)', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Sepiella inermis'))))
lines(mon.ine.fit, col='red' , lwd=2, lty=1)
lines(win.ine.fit, col='green' , lwd=2, lty=1)
lines(sum.ine.fit, col='blue' , lwd=2, lty=1)

plot(mon.tum.vari, las=1, type='p',pch=19, cex=1.4, col='white', ylim=c(0,1), xlim=c(0,9),
     ylab='', xlab='Distance (degree)', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.5,cex.main=1.5, 
     main=substitute(paste(italic('Saurida tumbil'))))
lines(mon.tum.fit, col='red' , lwd=2, lty=1)
lines(win.tum.fit, col='green' , lwd=2, lty=1)
lines(sum.tum.fit, col='blue' , lwd=2, lty=1)
dev.off()

# Map settings for predicted distribution-----------------------------------------------
numlevel=4
# Caption strings
sl=seq(min(krg$predict), max(krg$predict),by=(max(krg$predict)-min(krg$predict))/numlevel)
sll=formatC(sl, digits=2, format='f', decimal.mark = ".")
zlia = range(sl)
# Map
#jpeg("Prediction_winter.jpg", res = 600,height = 3.5,width = 4,units = "in")
par(family='', las=1, xaxs='i', yaxs='i', mar=c(5,5,2,6.5), cex.axis = 1, cex.lab = 1, cex.main=0.9)
plot(X,Y, type='n',xlab='Longitude (E)', ylab='Latitude (N)', #yaxt = 'n',
     main=substitute(paste(italic('B9. Saurida tumbil'))), xlim = c(66.5,74), ylim =c(15.6,23.1))
#image(krg, add=T, col=rev(viridis(numlevel)), zlim=zlia)

#colkey(side = 4, length = 0.7, padj = 0.5, shift = -0.1, dist = -0.01, add = T, cex.axis = 0.9, cex.clab = 1,
       #at = sl, mgp = c(1, 0.3, 0), tcl = -0.2, clim=zlia, labels = sll, col = rev(viridis(numlevel)),
       #side.clab = 1, line.clab = -14, clab = '              sqrt(CPUE)')
image(krg, add=T, col=rev(viridis(25)), zlim=range(krg$predict))

colkey(side = 4, length = 0.6, padj = 0.5, shift = -0.1, dist = -0.43, add = T, cex.axis = 0.8, cex.clab = 0.8,
       at = sl,mgp = c(1, 0.3, 0), tcl = -0.2, clim=range(krg$predict), labels = sll, col = rev(viridis(25)),
       side.clab = 1, line.clab = -7.5, clab = '         log(kg/h)')

plot(st_crop(india,box), add=TRUE,col = "black", lwd=1,xlim = c(66,74), ylim =c(15,24))
plot(st_crop(m50,box), add=TRUE,col = "red", lwd=0.8,xlim = x.range, ylim = y.range)
plot(st_crop(m100, box),add=TRUE,col = "brown", lwd=0.8,xlim = x.range, ylim = y.range)
plot(st_crop(m200, box), add=TRUE,col = "blue", lwd=0.8,xlim = x.range, ylim = y.range)
dev.off()



# Standard deviation map
esc.var=matrix(krg$krige.var, ncol=1)
esc.desv=matrix(sqrt(krg$krige.var), ncol=1)
#numlevel=4

# legend setting
slsd=seq(min(esc.desv), max(esc.desv),by=(max(esc.desv)-min(esc.desv))/numlevel)
sllsd=formatC(slsd, digits=2, format='f', decimal.mark = ".")
zlsd = range(slsd)

# Map
#jpeg("Variance_monsoon.s_man.jpg", res = 600,height = 7,width = 7,units = "in")
par(family='serif', las=1, xaxs='i', yaxs='i', mar=c(5,5,2,8), cex.axis = 1.5, cex.lab = 1.5)
plot(X,Y, type='n', xlab='Longitude (E)', ylab='Latitude (N)', #yaxt = 'n',
     xlim = c(66.5,74), ylim =c(15.6,23.1))
#axis(side = 2, at = seq(min(Y), max(Y), (max(Y)-min(Y))/5), las = 3)
image(krg, val=(krg$krige.var)^0.5, add=T, col=rev(viridis(numlevel)), zlim=c(min(esc.desv),max(esc.desv)))
#lines(lmt, lwd=2)
#box()
#points(dgeo$X,dgeo$Y, pch='+', col='red')

colkey(side = 4, length = 0.7, padj = 0.5, shift = 0, dist = 0.05, add = T, cex.axis = 0.9, cex.clab = 1,
       at = slsd, mgp = c(1, 0.3, 0), tcl = -0.2, clim=zlsd, labels = sllsd, col = rev(viridis(numlevel)),
       side.clab = 2, line.clab = -4.2, clab = 'CPUE')
plot(st_crop(india,box), add=TRUE,col = "black", lwd=2,xlim = c(66,74), ylim =c(15,24))
plot(st_crop(m50,box), add=TRUE,col = "red", lwd=2,xlim = x.range, ylim = y.range)
plot(st_crop(m100, box),add=TRUE,col = "brown", lwd=2,xlim = x.range, ylim = y.range)
plot(st_crop(m200, box), add=TRUE,col = "blue", lwd=2,xlim = x.range, ylim = y.range)
#dev.copy(pdf, paste(nome_analise, '_variance.map._monsoon.s.pdf', sep = ''), width=8, height=8)
dev.off()


##################Centres of gravity and inertia plot#################
##creatre db data
head(monsoon_box)
projec.toggle(0)
db.mon = db.create(monsoon_box)
db.win = db.create(wintre_box)
db.sum = db.create(summer_box)
#We can check the contents of the db.mon by typing any of the following commands:
db.print(db.mon)
print(db.mon)
db.mon
#add locators
db.mon = db.locate(db.mon,2:3,"x")
db.mon
db.mon = db.locate(db.mon,4:12,"z")
db.mon
db.win
db.win = db.locate(db.win,2:3,"x")
db.win = db.locate(db.win,4:12,"z")
db.sum = db.locate(db.sum,2:3,"x")
db.sum = db.locate(db.sum,4:12,"z")

india.map<-read.csv("E:/Bycatchmapping/R/INDIA.csv")

db.map<-db.create(india.map)
db.map
db.map<-db.locate(db.map,c(2,3),"x")

#install.packages("mapdata")
library(maps)
plot(db.mon,pch=19,col="blue")

#projec.define(projection="mean",db=db.mon)
#plot(db.mon,pch=19,col="blue")

map("worldHires",add=T,fill=T,col=8)
SI.cgi(db.win, name = db.getname(db.win,"z",7), flag.plot=TRUE,
       flag.inertia=T, flag.ellipse=F, col="blue",lwd=2)

###save cgi plot
setwd("E:/Bycatchmapping/R")

#jpeg("Saurida tumbil_cgi.jpg", res = 600,height = 5,width = 4,units = "in")
plot(db.mon,pch=19,col="white",title="",xlab='Longitude (E)', ylab='Latitude (N)', 
     main=substitute(paste(italic('Saurida tumbil'))), cex.axis = 1, cex.lab = 1, cex.main=0.9)
map("worldHires",add=T,fill=T,col=8)
SI.cgi(db.mon, name = db.getname(db.mon,"z",9), flag.plot=TRUE,
       flag.inertia=TRUE, col="red",lwd=2)
SI.cgi(db.win, name = db.getname(db.win,"z",9), flag.plot=TRUE,
       flag.inertia=TRUE, col="green",lwd=2)
SI.cgi(db.sum, name = db.getname(db.sum,"z",9), flag.plot=TRUE,
       flag.inertia=TRUE, col="blue",lwd=2)
dev.off()

##############Empirical Orthogonal Functions (EOF)####################
## converting kriging data to SpatialGridDataFrame
GT.pr <- points2grid(SpatialPoints(as.matrix(grid.map)))
GT.pr
ind <- .geoR_inout(grid.map, dlmt)

reorder <- as.vector(matrix(1:nrow(grid.map), nc=slot(GT.pr, "cells.dim")[2])[,slot(GT.pr, "cells.dim")[2]:1])

df.pr <- data.frame(predict=rep(NA, nrow(grid.map)), krige.var=rep(NA, nrow(grid.map)))
df.pr[ind,] <- as.data.frame(sum.ine.krig[1:2])
SGDF.ine.s <- SpatialGridDataFrame(grid=GT.pr, data=df.pr[reorder,])
SGDF.ine.s
na.omit(SGDF.ine.s)
rast.ine.s<-raster(SGDF.ine.s)
rast.ine.s
#saving raster image
#writeRaster(rast.ine.s$predict,"E:/Bycatchmapping/grid/inermis_s.tif")
##image plot of the objects
library(maps)
library(mapdata)
library(raster)
par(mfrow=c(1,1), mar=c(1.5,1.5,1,0.5), mgp=c(1.8,0,0))
raster::plot(rast.tum.w$predict)
maps::map("worldHires",add=T,fill=T,col=8)


##name SGDF and rater files for each species
#SGDF.tum.m #rast.tum.m
#SGDF.tum.w #rast.tum.w
#SGDF.tum.s #rast.tum.s
#SGDF.rib.m #rast.rib.m
#SGDF.rib.w #rast.rib.w
#SGDF.rib.s #rast.rib.s
#SGDF.duv.m #rast.duv.m
#SGDF.duv.w #rast.duv.w
#SGDF.duv.s #rast.duv.s
#SGDF.ell.m #rast.ell.m
#SGDF.ell.w #rast.ell.w
#SGDF.ell.s #rast.ell.s
#SGDF.jap.m #rast.jap.m
#SGDF.jap.w #rast.jap.w
#SGDF.jap.s #rast.jap.s
#SGDF.cuv.m #rast.cuv.m
#SGDF.cuv.w #rast.cuv.w
#SGDF.cuv.s #rast.cuv.s
#SGDF.pha.m #rast.pha.m
#SGDF.pha.w #rast.pha.w
#SGDF.pha.s #rast.pha.s
#SGDF.ham.m #rast.ham.m
#SGDF.ham.w #rast.ham.w
#SGDF.ham.s #rast.ham.s
#SGDF.ine.m #rast.ine.m
#SGDF.ine.w #rast.ine.w
#SGDF.ine.s #rast.ine.s

#create raster stack
rib.layer<-stack(rast.rib.m,rast.rib.w,rast.rib.s)
rib.layer
duv.layer<-stack(rast.duv.m,rast.duv.w,rast.duv.s)
duv.layer
ell.layer<-stack(rast.ell.m,rast.ell.w,rast.ell.s)
ell.layer
jap.layer<-stack(rast.jap.m,rast.jap.w,rast.jap.s)
jap.layer
cuv.layer<-stack(rast.cuv.m,rast.cuv.w,rast.cuv.s)
cuv.layer
pha.layer<-stack(rast.pha.m,rast.pha.w,rast.pha.s)
pha.layer
ham.layer<-stack(rast.ham.m,rast.ham.w,rast.ham.s)
ham.layer
ine.layer<-stack(rast.ine.m,rast.ine.w,rast.ine.s)
ine.layer
tum.layer<-stack(rast.tum.m,rast.tum.w,rast.tum.s)
tum.layer

#creat mean layer from raster stack
rib.layer.mean<-mean(rib.layer)
rib.layer.mean
duv.layer.mean<-mean(duv.layer)
duv.layer.mean
ell.layer.mean<-mean(ell.layer)
ell.layer.mean
jap.layer.mean<-mean(jap.layer)
jap.layer.mean
cuv.layer.mean<-mean(cuv.layer)
cuv.layer.mean
pha.layer.mean<-mean(pha.layer)
pha.layer.mean
ham.layer.mean<-mean(ham.layer)
ham.layer.mean
ine.layer.mean<-mean(ine.layer)
ine.layer.mean
tum.layer.mean<-mean(tum.layer)
tum.layer.mean

raster::plot(rib.layer.mean)
maps::map("worldHires",add=T,fill=T,col=8)

#saving mean layer image
#writeRaster(tum.layer.mean,"E:/Bycatchmapping/grid/tum.layer.mean.tif")

#conver raster data to text for EOF analysis
setwd("E:/Bycatchmapping/EOF")
library("spacetime")
library("tidyr")
ine.txt<-rasterToPoints(ine.layer)
#ine.txt<-write.csv(ine.txt,"ine.txt.csv")
#ine.txt<-read.csv("ine.txt.csv")
colnames(ine.txt)<-c("n","lon","lat","V1","V2","V3")
ine.lonlat<-ine.txt[,2:3]
ine.data<-ine.txt[,4:6]
head(ine.lonlat)

## Put data into space-wide form
Z <- t(ine.data)
dim(Z)
## First find the matrix we need to subtract:
spat_mean <- apply(ine.data, 1, mean)
nT <- ncol(ine.data)
## Then subtract and standardize:
Zspat_detrend <- Z - outer(rep(1, nT), spat_mean)
Zt <- 1/sqrt(nT - 1)*Zspat_detrend
#Finally, to carry out the SVD we run
E <- svd(Zt)
#The matrix V contains the EOFs in space-wide format
V <- E$v
colnames(E$v) <- paste0("EOF", 1:ncol(ine.data)) # label columns
EOFs <- cbind(ine.lonlat, E$v)
head(EOFs)
# convert U to data frame# add a time field# put columns (except time)# into long-table format with# EOF-PC as key-value pair
TS <- data.frame(E$u)%>% mutate(t = 1:nrow(E$u))%>% gather(EOF, PC, -t) 
TS$nPC <- TS$PC * sqrt(nT-1)
TS
##naming eof
#ribbon.eof<-EOFs
#duv.eof<-EOFs
#ell.eof<-EOFs
#jap.eof<-EOFs
#cuv.eof<-EOFs
#pha.eof<-EOFs
#ham.eof<-EOFs
#ine.eof<-EOFs
#tumbil.eof<-EOFs


#saving EOF image
#jpeg("tumbil_eof1.jpg", res = 600,height = 5,width = 4.5,units = "in")
gplot<-ggplot(tumbil.eof) + geom_tile(aes(x = lon, y = lat, fill = EOF1)) +
  scale_fill_viridis_c(direction = -1)+
  theme_bw() +
  xlab("Longitude (E)") + ylab("Latitude (N)")+
  geom_sf(data = india, col="black")+ 
  geom_sf(data = m50, col="red")+geom_sf(data = m100, col="brown")+geom_sf(data = m200, col="blue")+
  labs(title = substitute(paste(italic('Saurida tumbil'))),fill="log(kg/h)")+ theme(plot.title = element_text(hjust = 0.5))+
  xlim(66.5,74)+ylim(15.6,23.1)+
  coord_sf(xlim = c(66.5,74), ylim = c(15.6,23.1),expand = FALSE)
gplot+theme(legend.position = c(0.15,0.2))
dev.off()


#saving EOF image with reduced frame size
#jpeg("tumbil_eof.jpg", res = 600,height = 4,width = 3.5,units = "in")
gplot<-ggplot(tumbil.eof) + geom_tile(aes(x = lon, y = lat, fill = EOF1)) +
  scale_fill_viridis_c(direction = -1)+
  theme_bw() +
  xlab("Longitude (E)") + ylab("Latitude (N)")+
  geom_sf(data = india, col="black")+ 
  geom_sf(data = m50, col="red")+geom_sf(data = m100, col="brown")+geom_sf(data = m200, col="blue")+
  labs(title = substitute(paste(italic('Saurida tumbil'))),fill="log(kg/h)")+ theme(plot.title = element_text(hjust = 0.5))+
  xlim(66.5,74)+ylim(15.6,23.1)+
  coord_sf(xlim = c(66.5,74), ylim = c(15.6,23.1),expand = FALSE)
gplot+theme(legend.position = c(0.18,0.28))
dev.off()
citation("RGeostats")
