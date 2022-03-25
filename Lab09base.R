Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")
options(repos ="http://cran.r-project.org")  # required to get latest libs
if (!require("pacman")) install.packages("pacman")
pacman::p_load(aqp,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,soilDB,circlize,topmodel,DEoptim)
install.packages("EcoHydRology", repos="http://R-Forge.R-project.org")
pacman::p_load(EcoHydRology)
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

trunc((180+declat)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)

# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"

# Now we will build our proj4strings which define our “Coordinate 
# Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)
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
browseURL(url)
pourpoint=SpatialPoints(gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=gagepoint_utm@coords
# We are familiar with this basin, and know it is mostly to the NW 
searchlength=sqrt(area*4)*1000 
bboxpts
bboxpts=rbind(bboxpts,bboxpts)
bboxpts[2,1]=bboxpts[2,1]-searchlength
bboxpts[2,2]=bboxpts[2,2]+searchlength
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                        z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydem)
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
#
# A series of plots to show all of the components
#
par(mfrow = c(2, 2))
plot(TI)
plot(TIC)
plot(mybasinsca)
plot(mybasinslp)
dev.off()
plot(TIC,col=rainbow(5))
# Make a poly with command line gdal (fast)
# 
system("gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp")
mydemw_poly=readOGR("mydemw_poly_gdal.shp")
plot(mydemw_poly,add=T,border="black")
mydemw_poly
writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)
# We will use this ESRI shape file, a zipped group of it, to download 
# our soil extent from the WebSoilSurvey Website
#zip("mydemw.zip",list.files(pattern="mydemw[:.:]"))
# Download to your local machine mydemw.zip from the "Files" tab
# Open the WebSoilSurvey site to: 
#browseURL("https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx")
# "Creat AOI from a zipped shapefile"
# Open "Download Soils Data" Tab
# "Create Download Link" in lower right hand corner
# Right-Click on download link and "Copy Link Address" and 
# paste into a url object
#url=""
#download.file(url,"wss_aoi.zip")
#unzip("wss_aoi.zip")
# Or, if you want to make your life easier, 
mysoil_ll <- mapunit_geom_by_ll_bbox(bbox(spTransform(mydemw_poly,crs_ll)))
crs(mysoil_ll)=crs_ll
mysoil_utm = spTransform(mysoil_ll,crs_utm)
plot(mysoil_utm,add=T)

# Get the NLCD (USA ONLY)
# Returns a raster
# About time for a break!?!?
pacman::p_load(devtools,terra)
Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")
devtools::install_github("ropensci/FedData")
pacman::p_load(FedData)

NLCD <- get_nlcd(template=TIC, label='TIC',force.redo = TRUE)

# Plot with raster::plot
plot(NLCD)
NLCDutm=projectRaster(NLCD,TIC)
NLCDmask=mask(NLCDutm,mybasinmask)
plot(NLCDmask)


rmysoil_utm=rasterize(mysoil_utm,TIC,field=as.numeric(mysoil_utm$mukey))
unique(rmysoil_utm)
pacman::p_load(circlize)
plot(rmysoil_utm,col=rand_color(length(unique(values(rmysoil_utm)))))
rmysoil_utm=mask(rmysoil_utm,mybasinmask)
plot(rmysoil_utm,col=rand_color(length(unique(values(rmysoil_utm)))))
unique(rmysoil_utm)
#
ratify( TIC*10^9)
unique(ratify(round(mybasinslp*10+1)))
unique(ratify((rmysoil_utm*10^3)))


hru=ratify(TIC*10^9 + (rmysoil_utm*10^3) + round(mybasinslp*10+1))
unique(values(hru))
length(unique(values(hru)))
plot(hru,col=rand_color(length(unique(values(hru)))))


hru_table = levels(hru)[[1]]
origID = hru_table$ID # preserve data order for later reordering
# metadata parameters from a string... this will make more sense
# after the next "head()" command
hru_table$TIclass = as.numeric(substr(sprintf("%10.0f", hru_table$ID), 1,1))
hru_table$mukey = as.numeric(substr(sprintf("%10.0f", hru_table$ID), 2,7))
hru_table$slp = (as.numeric(substr(sprintf("%10.0f", 
                                             hru_table$ID), 8,10))-1)/10
#
# Calculate the area for each unique soil (mukey) X TIClass combination
# using res(raster) for x and y resolution in units of m
# Note that table() function returns the count of the occurrences of
# unique values in the hru raster cells.
hru_table$areaSQKM = as.vector(round(res(hru)[1]*res(hru)[2]*
                                         table(values(hru))/10^6, 3))
#
# To better understand what happened, look at the new hru_table
head(hru_table)
summary(hru_table)

rm("mu2co")
for (mymukey in unique(hru_table$mukey)){
  print(mymukey)
  mukey_statement = format_SQL_in_statement(mymukey)
  q_mu2co = paste("SELECT mukey,cokey FROM component 
           WHERE mukey IN ", mukey_statement, sep="")
  if(!exists("mu2co")){
    mu2co=SDA_query(q_mu2co)} 
  else{
    mu2co=rbind(mu2co,SDA_query(q_mu2co))
  } 
}
View(mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
# cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
# q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r,frag3to10_r  
# FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
# co2ch = SDA_query(q_co2ch)
rm("co2ch")
for (mycokey in unique(mu2co$cokey)){
  print(mycokey)
  cokey_statement = format_SQL_in_statement(mycokey)
  q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r,frag3to10_r FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
  print(q_co2ch)
  if(!exists("co2ch")){
    co2ch=SDA_query(q_co2ch)
  } else{
    try((co2ch=rbind(co2ch,SDA_query(q_co2ch))))
  } 
}
View(co2ch)
rm("co2co")
for (mycokey in unique(mu2co$cokey)){
  print(mycokey)
  cokey_statement = format_SQL_in_statement(mycokey)
  q_co2co = paste("SELECT cokey,slopelenusle_r FROM component WHERE cokey IN ", cokey_statement, sep="")
  print(q_co2co)
  if(!exists("co2co")){
    co2co=SDA_query(q_co2co)} 
  else{
    try((co2co=rbind(co2co,SDA_query(q_co2co))))} 
}
View(co2co)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
mu2ch=merge(mu2ch,co2co)
View(mu2ch)
#
#
#source CNmodel function
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")
pacman::p_load(data.table)
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
TIC02$P=TIC01$Excess+TIC02$P
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAsol$CN[2], 
                declat=declat,declon=declon)
TIC03$P=TIC02$Excess+TIC03$P
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAsol$CN[3], 
                declat=declat,declon=declon)
TIC04$P=TIC03$Excess+TIC04$P
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAsol$CN[4], 
                declat=declat,declon=declon)
TIC05$P=TIC04$Excess+TIC05$P
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAsol$CN[5], 
                declat=declat,declon=declon)

