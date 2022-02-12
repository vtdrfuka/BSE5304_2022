source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/soilwetting.R")
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/soil_wetting_above_capacity.R")
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/soildrying.R")

TMWBmodel=function(TMWB=TMWB,fcres=.25,SFTmp=0,bmlt6=2.5,bmlt12=1,Tlag=.5,AWCval=200){
# Now complete the model… what flows from TopSlope to MidSlope, and down to 
# BotSlope. How will these be connected?

  # notice that there is an Energy Balance based Snow Accumulation 
  # and Melt model in the EcoHydRology package.
  attach(TMWB)
SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = 0, aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                      startingSnowDensity_kg_m3=450)
# Note that the -3 in the above 
detach(TMWB)
TMWB$SNO=SNO_Energy$SnowWaterEq_mm
TMWB$SNOmlt=SNO_Energy$SnowMelt_mm
attach(TMWB)
TMWB$Albedo=.23
TMWB$Albedo[TMWB$SNO>0]=.95
PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,lat_radians = myflowgage$declat*pi/180) * 1000
TMWB$PET=PET
detach(TMWB)
# add in rm
rm(list=c("PET"))


#TMWB$AWC=(0.45-0.15)*1000 #Fld Cap = .45, Wilt Pt = .15, z=1000mm
TMWB$AWC=AWCval
# Oh, this we want to vary some of these around our watershed!
TMWB$dP = 0 # Initializing Net Precipitation
TMWB$ET = 0 # Initializing ET
TMWB$AW = 0 # Initializing AW
TMWB$Excess = 0 # Initializing Excess
# Functions for the Thornthwaite-Mather
#

# Loop to calculate AW and Excess
attach(TMWB)
for (t in 2:length(AW)){
  # This is where ET and Net Precipitation is now calculated
    ET[t] = min (AW[t-1],PET[t])
    ET[t] = (AW[t-1]/AWC[t-1])*PET[t] # New Model
     if(AvgTemp[t] >= SFTmp){
        dP[t] = P[t] - ET[t] + SNOmlt[t]
      }  else {
          dP[t] = ET[t]
        }
  # From here onward, everything is the same as Week2’s lab
   if (dP[t]<=0) {
    values<-soildrying(AW[t-1],dP[t],AWC[t])
  } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
    values<-soilwetting(AW[t-1],dP[t],AWC[t])
  } else {
    values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
  }
  AW[t]<-values[1]
  Excess[t]<-values[2]
}
TMWB$AW=AW
TMWB$Excess=Excess
# TMWB$dP=Excess  # This was in error originally
TMWB$dP=dP
TMWB$ET=ET
detach(TMWB) # IMPORTANT TO DETACH
rm(list=c("AW","Excess","dP","ET"))

TMWB$Qpred=NA
TMWB$Qpred[1]=0
TMWB$S=NA
TMWB$S[1]=0
attach(TMWB)
#fcres=.3        # Oh, this we want to vary in different areas
for (t in 2:length(Qpred)){
  S[t]=S[t-1]+Excess[t]     
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
  print(t)
}
TMWB$S=S
TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING
detach(TMWB) # IMPORTANT TO DETACH
rm(list=c("S","Qpred"))
  return(TMWB)
}
