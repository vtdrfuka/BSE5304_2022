NSE=function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs, na.rm=TRUE))^2, na.rm=TRUE))
}
