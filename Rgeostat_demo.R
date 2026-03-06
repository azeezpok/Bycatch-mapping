#################
## ICES CRR Handbook of geostatistics in R for fisheries and marine ecology
##
## This code computes variograms on hake densities for a series of trawl surveys
## in the gulf of Lion, France. The data were supplied by Ifremer
##
## Author: N.Bez, IRD
#################
# Clean workspace
#rm(list=ls(all=TRUE))
# Load libraries
#install.packages("mapdata")
library(RGeostats)
library(mapdata)
# Inactivate any previous projection
projec.toggle(0)
# Load data
rg.load(filename="Demo.hake.med.db.data", objname="db.data")
# Data presentation
plot(db.sel(db.data,YEAR==1996),zmin=0.001,pch.low=3,cex.low=0.25,las
     =1,pch=21,col=1, inches=5,title="Hake - 1996",asp=1)
map("worldHires",add=T)
# Define a simple projection based on the cosine of the mean latitude
projec.define(projection="mean",db=db.data)

# Evaluate the distance lag (in projected units i.e. n.mi.)
# and the number of lags by clicking on points
plot(db.sel(db.data,YEAR==1996),zmin=0.001,pch.low=3,cex.low=0.25,las
     =1,pch=21,col=1,inches=5,title="Hake - 1996",asp=1)
worldHires <- map("worldHires",plot=F,xlim=c(3,5),ylim=c(42,44))
lines(projec.operate(worldHires $x,worldHires $y))

# click on the graph at two adjacent sample locations
# distance lag
lag <- dist.digit()
lag <- signif(lag,2)
lag

# click on the graph at two sample locations separated by the greatest distance
# number of lag
diagonal <- dist.digit()
nlag <- ceiling((diagonal/2)/lag)
nlag
# Compute and represent omnidirectional variogram
vg.data <- vario.calc(db.sel(db.data,YEAR==1996), lag=lag,nlag=nlag)
# Edit the results
vg.data
# Plot the results
plot(vg.data,las=1,xlab="Distance (n.mi.)")
plot(vg.data,npairdw=T,inches=0.1,las=1,xlab="Distance (n.mi.)")


# Compute annual variograms and superimpose them
for(i in unique(db.data[,"YEAR"])){
  vg.data <- vario.calc(db.sel(db.data,YEAR==i), lag=lag,nlag=nlag)
  plot(vg.data,npairdw=T,inches=0.05,col=rgb(0,0,0,0.25),add=!(i==1996),
       las=1,xlab="Distance (n.mi.)",ylim=c(0,1e+08))
}
# superimpose several omnidirectional standardized variograms
for(i in unique(db.data[,"YEAR"])){
  vg.data <- vario.calc(db.sel(db.data,YEAR==i), lag=lag,nlag=nlag)
  plot(vg.data,npairdw=T,inches=0.1,col=rgb(0,0,0,0.25),add=!(i==1996),
  flag.norm=T,las=1,xlab="Distance (n.mi.)",ylim=c(0,2))
}

# Same computations but for the log-transformation of the hake density
# This transformation can reduce the fluctuations and facilitate the capture of
# a structure, however it does not allow to go from the structure of the log to
# the structure of the raw variable.
for(i in unique(db.data[,"YEAR"])){
  vg.data <- vario.calc(db.sel(db.add(db.data,z1=log(1+MERLMER)),YEAR==i),
                        lag=lag,nlag=nlag)
  plot(vg.data,npairdw=T,inches=0.1,col=rgb(0,0,0,0.25),add=!(i==1996),
  flag.norm=T,las=1,xlab="Distance (n.mi.)",ylim=c(0,2))
}

# Standardized the density by the annual standard deviation and create a new file.
# This should not be done with R function var() or sd() which computes
# the variance in (n-1)
# The YEAR is then attributed the locator "code" for selecting pairs of points of
# similar years from the same year
db.data.std <- db.data
for(i in unique(db.data[,"YEAR"])){
  sel <- db.data.std[,3]==i
  sd.year <- sqrt(mean(db.data.std[,11][sel]^2) - mean(db.data.std[,11][sel])^2)
  db.data.std[,11][sel] <- db.data.std[,11][sel]/sd.year
}
db.data.std <- db.locate(db.data.std,3,"code")

# Compute annual variograms (which are normalized because of the stan dardization of
# the density values)
for(i in unique(db.data[,"YEAR"])){
  vg.data <- vario.calc(db.sel(db.data.std,YEAR==i),lag=lag,nlag=nlag)
  plot(vg.data,npairdw=T,inches=0.1,col=rgb(0,0,0,0.25),add=!(i==1996),
  las=1,xlab="Distance (n.mi.)",ylim=c(0,2))
}
# Compute the mean annual variogram
# Pairs are retained if their codes are the same i.e. if their difference is smaller # or equal than 0
vg.data.std <- vario.calc(db.data.std, lag=5,nlag=15,opt.code=1,tolcode=0)
plot(vg.data.std,npairdw=T,inches=0.1,las=1,add=T,col=2,lwd=2)
vario.average(vg.data.std)
sim<- simtub(model=vg.data.std,dbout=grid.db,nbsim=1,nbtuba=1000)

library(gmGeostats)
TurningBands(nsim = 1, nBands = 1000)
?simtub
###########################Mapping by kriging with a variogram#################
## ICES CRR Handbook of geostatistics in R for fisheries and marine ecology
##
## This code performs variography and mapping by kriging for a fisheries
## acoustic survey on anchovy in the bay of Biscay, France
## The data were supplied by Ifremer
##
## Author: P.Petitgas, Ifremer
#################
# Clean workspace
#rm(list=ls(all=TRUE))
# Load geostatistical package and others
library(RGeostats)
library(mapdata)
# Inactivate any previous projection
projec.toggle(0)
# Load data
rg.load(filename="Demo.anchovy.bob.2d.db.data",objname="db.data")
rg.load(filename="Demo.anchovy.bob.2d.poly.data",objname="poly.data")
# Area limits
y1lim <- 43.3; y2lim <- 47; x1lim <- -4.5; x2lim <- -1
# Plot data
plot(db.data,name="ENGR.ENC",pch=1,asp=1.2,inches=5,col="black",
     xlim=c(x1lim,x2lim),ylim=c(y1lim,y2lim))
plot(poly.data, add=T, lty=1, density=0)
map("worldHires",add=T,fill=T,col=8)
######### Variography#######
# Define projection
projec.define(projection="mean")
# Mask duplicates (points too close)
db.data <- duplicate(db.data)
# Calculate directional variogram
vg2 <- vario.calc(db.data,lag=c(2,15),dirvect=c(35,145), nlag=c(40,7))
plot(vg2,npairpt=0,npairdw=TRUE,title="",inches=.05)
# omni-directional variogram
vg <- vario.calc(db.data,lag=2,dirvect=NA, nlag=40)
plot(vg,npairpt=0,npairdw=TRUE,title="",inches=.05)
# fit isotropic variogram
vg.mod <- model.auto(vario=vg,struct=melem.name(c(1,3,3)))
vg.mod
######### Kriging#######
# Grid for Kriging
x0 <- -4; y0 <- 43.4; dx <- 0.1;dy <- 0.1; nx <- 30; ny <- 37
db.grid <- db.create(x0=c(x0,y0),dx=c(dx,dy),nx=c(nx,ny))
# Select grid points inside polygon
db.grid <- db.polygon(db.grid,poly.data)
# plot data, grid and polygon
plot(db.grid, xlim=c(x1lim,x2lim),ylim=c(y1lim,y2lim),pch=3, col="red",asp=1.2,
     flag.proj=FALSE)
plot(db.data,pch=20,add=T,col="black",inches=3, flag.proj=FALSE)
map("worldHires",add=T,fill=T,col=8); box()
# neighbourhood
neimov <- neigh.create(ndim=2,type=2,nmini=3,nmaxi=10,radius=25)
# Kriging
kres <- kriging(dbin=db.data,dbout=db.grid, model=vg.mod, neigh=neimov)
# plot kriging results: K.estim
plot(kres,name.image="z1",title="K.estim",col=topo.colors(20), asp=1.2,
     xlim=c(x1lim,x2lim),ylim=c(y1lim,y2lim),pos.legend=5,flag.proj=FALSE)
plot(db.data,pch=20,add=T,col="red",inches=3,flag.proj=FALSE)
map("worldHires",add=T,fill=T,col=8)
# plot kriging results: K.std
plot(kres,name.image="z2",title="K.std",col=topo.colors(20), asp=1.2,
     xlim=c(x1lim,x2lim),ylim=c(y1lim,y2lim),pos.legend=5,flag.proj=FALSE)
plot(db.data,pch=20,add=T,col="black",inches=1.5,flag.proj=FALSE)
map("worldHires",add=T,fill=T,col=8)
# ratio of means kriged.map/data
mean(kres[db.grid[,"sel"],"z1"])/mean(db.data[,"z1"])


##################Centre of gravity, inertia, and isotropy of hake###############
#Pre-requisite
projec.toggle(0)
rg.load("Demo.hake.bob.db.data","db.data")
rg.load("Demo.hake.bob.poly.data","poly.data")
projec.define(projection="mean",db=db.data)
# Compute and plot the inertia, the total abundance, the isotropy,
# the center of gravity, and the coordinates of the axes of inertia.
# Note that intermediate results of the PCA decomposition are provided
# (the eigen values and the eigen vectors).
plot(db.data,title="Centre of gravity and inertia of densities and samples",asp=1,
     xlim=c(-300,150),ylim=c(-200,150),inches=5,
     xlab="Nautical mile",ylab="Nautical mile")
plot(poly.data,col=8,add=T)
SI.cgi(db.data,flag.plot=TRUE,flag.inertia=TRUE,col=2)
citation(RGeostats)



#####################Create db adata#################

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
#plot(db.map,add=T,col = "black", lwd=2,)

#projec.define(projection="mean",db=db.mon)
#plot(db.mon,pch=19,col="blue")

map("worldHires",add=T,fill=T,col=8)

#map('world2Hires', xlim=c(67, 77), ylim = c(15.6,23))

SI.cgi(db.mon,flag.plot=TRUE,flag.inertia=TRUE,col=2)
map("worldHires",add=F,fill=T,col=8)

##plot with >2 variable
# Compute the global index of collocation between age 0 and age 1
SI.gic(db1=db.mon,db2=db.mon,name1="T.lepturus",name2="U.duvaucelii",
       flag.plot=F,flag.inertia=T,asp=1,inches=5,
       xlab="Nautical mile",ylab="Nautical mile",
       col1="red",col2="blue",)

SI.cgi(db.win, name = db.getname(db.win,"z",7), flag.plot=TRUE,
       flag.inertia=T, flag.ellipse=F, col="blue",lwd=2)

setwd("E:/Bycatchmapping/R")

#jpeg("Saurida tumbil_cgi.jpg", res = 600,height = 5,width = 4,units = "in")
#par(family='', las=1, xaxs='i', yaxs='i', mar=c(5,5,2,6.5), cex.axis = 1, cex.lab = 1, cex.main=0.9)
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

###################EOF##############



delete_rows <- which(SSTlandmask == 1)
SSTdata <- SSTdata[-delete_rows, 1:396]

## Put data into space-wide form
Z <- t(SSTdata)
dim(Z)
## First find the matrix we need to subtract:
spat_mean <- apply(SSTdata, 1, mean)
nT <- ncol(SSTdata)
## Then subtract and standardize:
Zspat_detrend <- Z - outer(rep(1, nT), spat_mean)
Zt <- 1/sqrt(nT - 1)*Zspat_detrend
#Finally, to carry out the SVD we run
E <- svd(Zt)
#The matrix V contains the EOFs in space-wide format
V <- E$v
colnames(E$v) <- paste0("EOF", 1:ncol(SSTdata)) # label columns
EOFs <- cbind(SSTlonlat[-delete_rows, ], E$v)
head(EOFs[, 1:6])
# convert U to data frame# add a time field# put columns (except time)# into long-table format with# EOF-PC as key-value pair
TS <- data.frame(E$u) %>% 
  mutate(t = 1:nrow(E$u)) %>% 
  gather(EOF, PC, -t) 
TS$nPC <- TS$PC * sqrt(nT-1)
ggplot(EOFs) + geom_tile(aes(x = lon, y = lat, fill = EOF1)) +
   theme_bw() +
  xlab("Longitude (deg)") + ylab("Latitude (deg)")
head(EOFs)
