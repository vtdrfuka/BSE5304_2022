Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")
dir.create("~/Lab10")
setwd("~/Lab10")
options(repos ="http://cran.r-project.org")  # required to get latest libs
if (!require("pacman")) install.packages("pacman")
#install.packages("EcoHydRology", repos="http://R-Forge.R-project.org")
pacman::p_load(devtools,terra)
devtools::install_github("ropensci/FedData")
pacman::p_load(lubridate,aqp,curl,httr,rnoaa,raster,shapefiles,
               rgdal,elevatr,soilDB,circlize,topmodel,DEoptim,
               FedData,raster,EcoHydRology,ggplot2,data.table)

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

proj4_utm = paste0("+proj=utm +zone=", trunc((180+declon)/6+1), " +datum=WGS84 +units=m +no_defs")

# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"

# Now we will build our proj4strings which define our “Coordinate 
# Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
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
#browseURL(url)
pourpoint=SpatialPoints(gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=gagepoint_utm@coords
# We are familiar with this basin, and know it is mostly to the NW 
searchlength=sqrt(area*4)*1000 
bboxpts=rbind(bboxpts,bboxpts)
bboxpts[2,1]=bboxpts[2,1]-searchlength
bboxpts[2,2]=bboxpts[2,2]+searchlength
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                        z = 12, prj = proj4_utm,src ="aws",expand=1)
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
plot(TIC)
# Make a poly with command line gdal (fast)
# 
system("gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp")
mydemw_poly=readOGR("mydemw_poly_gdal.shp")
plot(mydemw_poly,add=T,border="black")
writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)

#source CNmodel function
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")
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
TIC02$P=TIC01$Qpred+TIC02$P
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAsol$CN[2], 
                declat=declat,declon=declon)
TIC03$P=TIC02$Qpred+TIC03$P
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAsol$CN[3], 
                declat=declat,declon=declon)
TIC04$P=TIC03$Qpred+TIC04$P
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAsol$CN[4], 
                declat=declat,declon=declon)
TIC05$P=TIC04$Qpred+TIC05$P
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAsol$CN[5], 
                declat=declat,declon=declon)

PLossFunc=function(DPTIdf=DPTI01
                   ,tau=9.3  # days
                   ,dt=1     # days time step
                   ,kF=.015  # Table 2
                   ,TICn=1){
# should be in m of water
# tau=9.3  # days
# dt=1     # days time step
# kF=.015  # Table 2
# Initialize MF and DF
DPTIdf$MF=0
DPTIdf$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTIdf$MF[(format(DPTIdf$date,"%j") %in% c(121,213,274))]=5.4*10^-4
# Remember what we say about attaching! 
attach(DPTIdf)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTIdf$MF=MF
DPTIdf$DF=DF
detach(DPTIdf)
rm(list=c("MF","DF")) # Clean up the environment 

# Calculate your Export Coef from Easton et al 2007 Figure 2 using TI Class 5
# and note that the bold gives the TI Class. This figure gives a range of 
# 0 - 520 micrograms/litre 

muTS_TI=(((520-0)/5)*TICn+0) # use VSAsol$TIClass table to 
# [microgram P/litre]
# assign TIC to calc (remember TIC05 is in location 1 in the VSAsol$TIClass 
# table
# Setting range of Soil P values (MS) using the range given on page 7 of 
# Easton et al. 2007 (3.7-18.5mg/kg), then 
# Moore 1993… assume soil density is 2000kg/m^3 of soil
MS_TI=(((18.5-3.7)/5)*TICn+3.7)*2000  # range from Easton et 
# MS_TI [mg P/m^3 soil]
# al. 2007, pg 7. Moore suggests linear relationship
# We will take care of all of TIClass 05 now as will so 
# it makes sense when you repeat for TI Classes 1-4
# You will use muTS_TI01 thru muTS_TI04 to finish this lab
#
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTIdf$muS= muTS_TI*QS^((DPTIdf$Tavg-TR)/10)  # Eq. 5
DPTIdf$DS=(DPTIdf$muS*MS_TI*DPTIdf$Rt)/10^9          # Eq. 4
return(DPTIdf)
#### Above is run for each TIClass
}

TIC01$Rt=TIC01$Qpred/1000
TIC01$Tavg=(TIC01$MaxTemp+TIC01$MinTemp)/2
TIC01=PLossFunc(DPTIdf = TIC01,tau=9.3,dt=1,kF=.015,TICn=1)

TIC02$Rt=TIC02$Qpred/1000
TIC02$Tavg=(TIC02$MaxTemp+TIC02$MinTemp)/2
TIC02=PLossFunc(DPTIdf = TIC02,tau=9.3,dt=1,kF=.015,TICn=2)

TIC03$Rt=TIC03$Qpred/1000
TIC03$Tavg=(TIC03$MaxTemp+TIC03$MinTemp)/2
TIC03=PLossFunc(DPTIdf = TIC03,tau=9.3,dt=1,kF=.015,TICn=3)

TIC04$Rt=TIC04$Qpred/1000
TIC04$Tavg=(TIC04$MaxTemp+TIC04$MinTemp)/2
TIC04=PLossFunc(DPTIdf = TIC04,tau=9.3,dt=1,kF=.015,TICn=4)

TIC05$Rt=TIC05$Qpred/1000
TIC05$Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2
TIC05=PLossFunc(DPTIdf = TIC05,tau=9.3,dt=1,kF=.015,TICn=5)

DPLT=data.frame(date=TIC05$date,
                Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
plot(DPLT$date,DPLT$LB)

TIC01$area=myflowgage$area/5*10^6
TIC02$area=myflowgage$area/5*10^6
TIC03$area=myflowgage$area/5*10^6
TIC04$area=myflowgage$area/5*10^6
TIC05$area=myflowgage$area/5*10^6

DPLT$LT=DPLT$LB +
  TIC01$area*(TIC01$DF + TIC01$DS)+
  TIC02$area*(TIC02$DF + TIC02$DS)+
  TIC03$area*(TIC03$DF + TIC03$DS)+
  TIC04$area*(TIC04$DF + TIC04$DS)+
  TIC05$area*(TIC05$DF + TIC05$DS)
plot(DPLT$date,DPLT$LT, type="l")


mean(TIC03$DS) # kg/m^2/d
mean(c(18.5,3.7))*2000/1000 # [mg/kg] * 2000 [kg/m^3] / 1000[mg/kg] ): [kg/m^3]
.1 * .2 * mean(c(18.5,3.7))*2000/1000 # 20% of 10% of 1m , .02m * kg/m^3 = [kg/m^2]
.444/mean(TIC03$DS)/365



rcp01=mean(TIC01$DF+TIC01$DS)*365
rcp02=mean(TIC02$DF+TIC02$DS)*365
rcp03=mean(TIC03$DF+TIC03$DS)*365
rcp04=mean(TIC04$DF+TIC04$DS)*365
rcp05=mean(TIC05$DF+TIC05$DS)*365

plot(TIC)

m <- c(1,rcp01,
       2,rcp02,
       3,rcp03, 
       4,rcp04, 
       5,rcp05)
rclmat <- matrix(m, ncol=2, byrow=TRUE)
rc <- reclassify(TIC, rclmat)
plot(rc)



