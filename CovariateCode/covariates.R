#setwd("/Users/xiaoqixie/Desktop/Mountain Pine Beetle")
library(terra)
library(here)
#locations
map=rast(here("study_area.tif"))
xy=xyFromCell(map,1:ncell(map))
v=extract(map,xy)

xy=xy[which(v==1),]#622286 cells overall
#start with pine density

mz2021=rast(here("mz2021"))
mz2021=disagg(mz2021,4)

namR=list.files(path=here("Infest Green Red"),pattern="Red",full.names=T)
namR=namR[c(-1,-17)]
namG=list.files(path=here("Infest Green Red"),pattern="Green",full.names=T)
namG=namG[c(-12,-13)]
#lapply takes a vector or list
rastR=lapply(namR,rast)
rastG=lapply(namG,rast)
for(i in 1:15){
  rastR[[i]]=resample(rastR[[i]],mz2021,"near")
}
for(i in 1:11){
  rastG[[i]]=resample(rastG[[i]],mz2021,"near")
}
GreenFromRed=as.numeric()
for(i in 2:15){
  e=as.matrix(extract(rastR[[i]],xy))
  e[is.na(e)]=0
  GreenFromRed=cbind(GreenFromRed,e)
}
for(i in 2:15){
  colnames(GreenFromRed)[i-1]=paste("GreenFromRed",2005+i)
}
#detected green
#2007-2020
detectedGreen=as.numeric()
for(i in 1:3){
  e=as.matrix(extract(rastG[[i]],xy))
  e[is.na(e)]=0
  detectedGreen=cbind(detectedGreen,e)
}
detectedGreen=cbind(detectedGreen,0)# 2010
for(i in 4:9){
  e=extract(rastG[[i]],xy)
  e[is.na(e)]=0
  detectedGreen=cbind(detectedGreen,e)
}
detectedGreen=cbind(detectedGreen,0,0)#2017,2018, no detected green
e=extract(rastG[[10]],xy)
e[is.na(e)]=0
detectedGreen=cbind(detectedGreen,e)#2019
detectedGreen=cbind(detectedGreen,0)#2020 none
for(i in 2:15){
  colnames(detectedGreen)[i-1]=paste("DetectedGreen",i+2005)
}
#green trees can raise beetles are green from red + detected green
greenStayToNextYear=detectedGreen+GreenFromRed
#adds controlled trees into response variable
#cut green+ green from red +detected green=total green
namCtr=list.files(path=here("ctrol"),pattern="ctr",full.names=T)
rastCtr=lapply(namCtr,rast)
for(i in 1:15){
  rastCtr[[i]]=resample(rastCtr[[i]],mz2021,"near")
}
#extract each year cut trees
cutGreen=as.numeric()
for(i in 2:15){
  e=as.matrix(extract(rastCtr[[i]],xy))
  e[is.na(e)]=0
  cutGreen=cbind(cutGreen,e)
}
for(i in 2:15){
  colnames(cutGreen)[i-1]=paste("cutGreen",i+2005)
}
#add cut trees into 
CurrentYearInfest=greenStayToNextYear+cutGreen
save(CurrentYearInfest,file="Climate Change/CurrentYearInfest.RData")

########################################################################################
#total volume *(percentage abundance of lodgepole pine +jack pine)=volume of lodgepole pine and jack pine
total_volume=rast(here("NFI layers 2011/Structure_Volume_Total.tif"))
total_volume=aggregate(total_volume,fact=2,method="mean")
jack_percent=rast(here("NFI layers 2011/Species_Pinu_Ban.tif"))
jack_percent=aggregate(jack_percent,fact=2,method="mean")
lodge_percent=rast(here("NFI layers 2011/Species_Pinu_Con.tif"))
lodge_percent=aggregate(lodge_percent,fact=2,method="mean")
mz2021=rast("mz2021")
mz2021=disagg(mz2021,4)
total_volume1=project(total_volume,crs(mz2021),res=500)
total_volume2=resample(total_volume1,mz2021,method="bilinear")
total_volume2=crop(total_volume2,mz2021)

jack_percent1=project(jack_percent,crs(mz2021),res=500)
jack_percent2=resample(jack_percent1,mz2021,method="bilinear")
jack_percent2=crop(jack_percent2,mz2021)

lodge_percent1=project(lodge_percent,crs(mz2021),res=500)
lodge_percent2=resample(lodge_percent1,mz2021,method="bilinear")
lodge_percent2=crop(lodge_percent2,mz2021)

#extract xy and values from cell
xy=xyFromCell(total_volume2,1:ncell(total_volume2))
totalVolume=extract(total_volume2,xy)
jackPercent=extract(jack_percent2,xy)/100
lodgePercent=extract(lodge_percent2,xy)/100

pine_percent_v=totalVolume*(jackPercent+lodgePercent)#volume per hectare

E=function(volume_hectare){
  a=17.6
  delta=0.00527
  f=a*volume_hectare*exp(-delta*volume_hectare)
  return(f)
  
}
###################################################################

expected_volume=E(pine_percent_v)
d=cbind(xy,expected_volume)
library(raster)
r=rasterFromXYZ(d)
crs(r)=crs(mz2021)
writeRaster(r, filename = "pine_density_numbers_hectare_beaudoin_11.tif", format = "GTiff",overwrite=TRUE)
#numbers per hectare
#############################################################
density=rast(here("pine_density_numbers_hectare_beaudoin_11.tif"))
load(here("CurrentYearInfest.RData"))
#current pine density before attack
density2011=as.matrix(extract(density,xy)*500^2/100^2)
density2011=round(density2011)
pinedensity=as.numeric()
for(i in 1:3){
  v=apply(CurrentYearInfest[,c(i:4)],1,sum)
  pinedensity=cbind(pinedensity,density2011+v)
}
pinedensity=cbind(pinedensity,density2011+CurrentYearInfest[,4])#2010
pinedensity=cbind(pinedensity,density2011)#add 2011
v=CurrentYearInfest[,5]
v1=density2011-v
v1[which(v1<0)]=0
pinedensity=cbind(pinedensity,v1)#2012
for(i in 6:13){
  v=apply(CurrentYearInfest[,c(5:i)],1,sum)
  v1= density2011-v
  v1[which(v1<0)]=0
  pinedensity=cbind(pinedensity,v1)
}
colnames(pinedensity)=2007:2020

save(pinedensity,file="Climate Change/pinedensity.RData")

########################################################################
#########################################################################

#what happened for Q
#
hybrid=rast(here("hybridaep1000"))
hybrid=focal(hybrid,matrix(1,3,3),fun=mean,na.policy="only")
hybrid=project(hybrid,method='near',crs(mz2021),res=500)
hybrid=resample(hybrid,mz2021,'near')
hybrid=disagg(hybrid,fact=2)
q=as.matrix(extract(hybrid,xy))

save(q,file="Climate Change/Q.RData")


#still keep management zone even though there is no variation in the data

#elevation
elevation=rast(here('elevationaep'))#100
elevation=aggregate(elevation,fact=5,fun=mean)
elevation=resample(elevation,mz2021,'near')
elevation=focal(elevation,matrix(1,5,5),fun=mean,na.policy="only")
e=as.matrix(extract(elevation,xy))#elevation
length(which(is.na(e)))
elevation=as.numeric(replicate(14,as.matrix(e)))
save(elevation,file="Climate Change/elevation.RData")
#
aspect=rast(here("aspectaep"))#100
aspect=aggregate(aspect,fact=5,fum=mean)
aspect=resample(aspect,mz2021,'near')
aspect=focal(aspect,matrix(1,5,5),fun=mean,na.policy="only")
a=as.matrix(extract(aspect,xy))
rad=a/180*pi
northness=as.numeric(replicate(14,cos(rad)))
easterness=as.numeric(replicate(14,sin(rad)))
save(northness,file="Climate Change/northness.RData")
save(easterness,file="Climate Change/easterness.RData")

#slope
slope=rast(here('slopeaep'))#100
slope=aggregate(slope,fact=5,fun=mean)
slope=resample(slope,mz2021,'near')
slope=focal(slope,matrix(1,5,5),fun=mean,na.policy="only")
output=as.matrix(extract(slope,xy))
Slope=as.numeric(replicate(14,output))
save(Slope,file="Climate Change/slope.RData")

#
#extract within cell green red from 2006 to 2019
GreenFromRed=as.numeric()
for(i in 1:14){
  e=as.matrix(extract(rastR[[i]],xy))
  e[is.na(e)]=0
  GreenFromRed=cbind(GreenFromRed,e)
}
GreenFromRed=as.data.frame(GreenFromRed)
for(i in 1:14){
  colnames(GreenFromRed)[i]=paste("GreenFromRed",2005+i)
}
#detected green
detectedGreen=as.numeric()
for(i in 1:3){
  e=as.matrix(extract(rastG[[i]],xy))
  e[is.na(e)]=0
  detectedGreen=cbind(detectedGreen,e)
}
detectedGreen=cbind(0,detectedGreen,0)#add 2006, 2006 has no green trees, similar for 2010
for(i in 4:9){
  e=as.matrix(extract(rastG[[i]],xy))
  e[is.na(e)]=0
  detectedGreen=cbind(detectedGreen,e)
}
detectedGreen=cbind(detectedGreen,0,0)#2017,2018, no detected green
e=as.matrix(extract(rastG[[10]],xy))
e[is.na(e)]=0
detectedGreen=cbind(detectedGreen,e)#2019
detectedGreen=as.data.frame(detectedGreen)
for(i in 1:14){
  colnames(detectedGreen)[i]=paste("DetectedGreen",i+2005)
}
#green trees can raise beetles are green from red + detectd green
LastYearGreen=detectedGreen+GreenFromRed
save(LastYearGreen,file="Climate Change/LastYearGreen.RData")


#last year neighborhood
#
#last year neighborhood
nei3Red=list()
for(i in 1:14){
  nei3Red[[i]]=focal(rastR[[i]],w=matrix(c(1,1,1,1,0,1,1,1,1),nrow=3,byrow=T),fun=sum,na.rm=T)
}
nei3Green=list()
for(i in 1:11){
  nei3Green[[i]]=focal(rastG[[i]],w=matrix(c(1,1,1,1,0,1,1,1,1),nrow=3,byrow=T),fun=sum,na.rm=T)
}
#make table from raster
nei3=as.numeric()
#2006-2019
for(i in 1:14){
  e=as.matrix(extract(nei3Red[[i]],xy))
  e[is.na(e)]=0
  nei3=cbind(nei3,e)
}

Green=as.numeric()#last year green
for(i in 1:3){
  e=as.matrix(extract(nei3Green[[i]],xy))
  e[is.na(e)]=0
  Green=cbind(Green,e)
}
Green=cbind(0,Green,0)#add 2006, 2006 has no green trees, similar for 2010
for(i in 4:9){
  e=as.matrix(extract(nei3Green[[i]],xy))
  e[is.na(e)]=0
  Green=cbind(Green,e)
}
Green=cbind(Green,0,0)#2017,2018, no detected green
e=as.matrix(extract(nei3Green[[10]],xy))
e[is.na(e)]=0
Green=cbind(Green,e)
#
nei3=nei3+Green
save(nei3,file="Climate Change/last year neighborhood infestation 0.75km.RData")
##################################################################

nei5Red=list()
for(i in 1:14){
  nei5Red[[i]]=focal(rastR[[i]],w=matrix(c(1,1,1,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,1),nrow=5,byrow=T),fun=sum,na.rm=T)
}
nei5Green=list()
for(i in 1:11){
  nei5Green[[i]]=focal(rastG[[i]],w=matrix(c(1,1,1,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,1),nrow=5,byrow=T),fun=sum,na.rm=T)
}
#make table from raster
nei5=as.numeric()
#2006-2019
for(i in 1:14){
  e=as.matrix(extract(nei5Red[[i]],xy))
  e[is.na(e)]=0
  nei5=cbind(nei5,e)
}

Green=as.numeric()#last year green
for(i in 1:3){
  e=as.matrix(extract(nei5Green[[i]],xy))
  e[is.na(e)]=0
  Green=cbind(Green,e)
}
Green=cbind(0,Green,0)#add 2006, 2006 has no green trees, similar for 2010
for(i in 4:9){
  e=extract(nei5Green[[i]],xy)
  e[is.na(e)]=0
  Green=cbind(Green,e)
}
Green=cbind(Green,0,0)#2017,2018, no detected green
e=extract(nei5Green[[10]],xy)
e[is.na(e)]=0
Green=cbind(Green,e)
#
nei5=nei5+Green
#save nei3 data
save(nei5,file="Climate Change/last year neighborhood infestation 1.25km.RData")

loc=xy
nam=list.files(path=here(),pattern="mz",full.names=T)
nam=nam[c(1:14)]
rastm=lapply(nam,rast)
for(i in 1:14){
  rastm[[i]]=resample(rastm[[i]],mz2021,"near")
}
i=1
l=extract(rastm[[i]],loc)
l[is.na(l),1]=0
mz=l
for(i in 2:14){
  l=extract(rastm[[i]],loc)
  l[is.na(l),1]=0
  mz=cbind(mz,l)
}

#leading 1, all others 0
l1<-which(mz[,1]==2)
l2<-which(mz[,1]!=2)
mz[l1,1]<-1
mz[l2,1]<-0
#####
######
l1<-which(mz[,2]==2)
l2<-which(mz[,2]!=2)
mz[l1,2]<-1
mz[l2,2]<-0
###09 1 leading
mz[,3][which(mz[,3]!=1)]=0

l1<-which(mz[,4]==2)
l2<-which(mz[,4]!=2)
mz[l1,4]<-1
mz[l2,4]<-0

l1<-which(mz[,5]==2)
l2<-which(mz[,5]!=2)
mz[l1,5]<-1
mz[l2,5]<-0

l1<-which(mz[,6]==2)
l2<-which(mz[,6]!=2)
mz[l1,6]<-1
mz[l2,6]<-0

l1<-which(mz[,7]==2)
l2<-which(mz[,7]!=2)
mz[l1,7]<-1
mz[l2,7]<-0

l1<-which(mz[,8]==2)
l2<-which(mz[,8]!=2)
mz[l1,8]<-1
mz[l2,8]<-0

#15,16,17,18, 1 leading
for(i in 9:12){
  mz[,i][which(mz[,i]!=1)]=0
}
###
l1<-which(mz[,13]==2)
l2<-which(mz[,13]!=2)
mz[l1,13]<-1
mz[l2,13]<-0


l1<-which(mz[,14]==2)
l2<-which(mz[,14]!=2)
mz[l1,14]<-1
mz[l2,14]<-0


save(mz,file="Climate Change/mz.RData")

r1=rast()
crs(r1)=crs(mz2021)
ext(r1)=c(171965-4001,860465+4001,5442952-4001,6334952+4001)
loc1=crop(mz2021,ext(r1))

xy=xyFromCell(loc1,1:ncell(loc1))#larger than the interested locations
df=xyFromCell(map,1:ncell(map))
v=extract(map,df)
df=df[which(v==1),]#interested locations

namR=list.files(path="Infest Green Red",pattern="Red",full.names=T)
namR=namR[c(-1,-17)]
namG=list.files(path="Infest Green Red",pattern="Green",full.names=T)
namG=namG[c(-12,-13)]
#lapply takes a vector or list
rastR=lapply(namR,rast)
rastG=lapply(namG,rast)
for(i in 1:15){
  rastR[[i]]=resample(rastR[[i]],mz2021,"near")
}
for(i in 1:11){
  rastG[[i]]=resample(rastG[[i]],mz2021,"near")
}

library("dplyr")
radius=function(df1,df2){
  return(sqrt((df1[,1]-df2[1])^2+(df1[,2]-df2[2])^2)/1000)
}#in km

locWithin4km=function(df){
  r=radius(xy,df)
  data=data.frame(cbind(xy,r))
  data=data%>%filter(r<=4)
  return(data)
}

all=read.csv(here("Influence of Zone/All.csv"))
colnames(all)=c("distance","cdf")
all$distance=round(all$distance,5)
all$cdf=round(all$cdf,5)
all[89,]=c(4,1)
library('pbapply')
df=data.frame(df)
xy=data.frame(xy)
distance=pbapply(df,1,locWithin4km,cl=5)
save(distance,file="Climate Change/dsitance.RData")
load("Climate Change/distance.RData")
closestR=function(r){
  R=all$distance
  d1=abs(R-r)
  d2=R-r
  #d3=cbind(d1,d2)
  l=which(d1==min(d1))
  if(l>1){
    if(d2[l]>=0){
      prob=all$cdf[l]-all$cdf[l-1]}
    if(d2[l]<0){
      prob=all$cdf[l+1]-all$cdf[l]}}
  if(l==1){
    prob=all$cdf[l]
    
  }
  return(prob)
}
fx=list()
for(i in 1:dim(df)[1]){
  print(i)
  data=distance[[i]]#difference between two row in all is 0.0145
  pdf=lapply(data$r,closestR)
  pdf=as.numeric(pdf)
  pdf[is.na(pdf)]=0
  fx[[i]]=pdf
}
save(fx,file="Climate Change/fx.RData")
load("Climate Change/fx.RData")
func1=function(xy){
  
  GreenFromRed=as.numeric()
  for(i in 1:15){
    e=as.matrix(extract(rastR[[i]],xy))
    e[is.na(e)]=0
    GreenFromRed=cbind(GreenFromRed,e)
  }
  colnames(GreenFromRed)=2006:2020
  
  #cutGreen=as.numeric()
  #for(i in 1:15){
   # e=as.matrix(extract(rastCtr[[i]],xy))
  #  e[is.na(e)]=0
  #  cutGreen=cbind(cutGreen,e)
  #}
  #colnames(cutGreen)=2006:2020
  
  detectedGreen=as.numeric()
  for(i in 1:3){
    e=as.matrix(extract(rastG[[i]],xy))
    e[is.na(e)]=0
    detectedGreen=cbind(detectedGreen,e)
  }
  detectedGreen=cbind(0,detectedGreen,0)#add 2006, 2006 has no green trees, similar for 2010
  for(i in 4:9){
    e=extract(rastG[[i]],xy)
    e[is.na(e)]=0
    detectedGreen=cbind(detectedGreen,e)
  }
  detectedGreen=cbind(detectedGreen,0,0)#2017,2018, no detected green
  e=extract(rastG[[10]],xy)
  e[is.na(e)]=0
  detectedGreen=cbind(detectedGreen,e)
  detectedGreen=cbind(detectedGreen,0)#2020 none
  colnames(detectedGreen)=2006:2020
  
  Green=GreenFromRed+detectedGreen
  return(Green)
}

neib=list()
for(i in 1:dim(df)[1]){
  print(i)
  df1=as.matrix(distance[[i]])
  neib[[i]]=func1(df1[,1:2])
}
save(neib,file="Climate Change/neib.RData")
load("Climate Change/neib.RData")
influencewithin4km=as.numeric()
for(i in 1:dim(df)[1]){
  print(i)
  d1=fx[[i]]
  d2=neib[[i]]
  d3=as.numeric()
  for(j in 1:15){
    d3[j]=sum(d1*d2[,j])
  }
  influencewithin4km=rbind(influencewithin4km,d3)
}
colmax=apply(influencewithin4km,2,max)
save(influencewithin4km,file="Climate Change/influencewithin4km.RData")

locWithin1.25km=function(df){
  location=xy%>%filter(x>df[1]-1250&x<df[1]+1250&y>df[2]-1250&y<df[2]+1250)
  r=radius(location,df)
  data=data.frame(cbind(location,r))
  data=data%>%filter(r<=1.25)
  return(data)
}
distance1=pbapply(df,1,locWithin1.25km,cl=8)#speed slower than my expectation
save(distance1,file="Climate Change/distance1.RData")
load("Climate Change/distance1.RData")
fx1=list()
for(i in 1:dim(df)[1]){
  print(i)
  data=distance1[[i]]#difference between two row in all is 0.0145
  pdf=lapply(data$r,closestR)
  pdf=as.numeric(pdf)
  pdf[is.na(pdf)]=0
  fx1[[i]]=pdf
}#pdf of each point
#func1,2,3 refer to the draft 1
save(fx1,file="Climate Change/fx1.RData")
neib1=list()
for(i in 1:dim(df)[1]){
  print(i)
  df1=as.matrix(distance1[[i]][,1:2])
  neib1[[i]]=func1(df1)
}
save(neib1,file="Climate Change/neib1.RData")
load("Climate Change/neib1.RData")
influencewithin1.25km=as.numeric()
for(i in 1:dim(df)[1]){
  print(i)
  d1=fx1[[i]]
  d2=neib1[[i]]
  d3=as.numeric()
  for(j in 1:15){
    d3[j]=sum(d1*d2[,j])
  }
  influencewithin1.25km=rbind(influencewithin1.25km,d3)
}
save(influencewithin1.25km,file="Climate Change/influencewithin1.25km.RData")

load("Climate Change/influencewithin1.25km.RData")
influenceexclude1.25km=influencewithin4km-influencewithin1.25km
save(influenceexclude1.25km,file="Climate Change/influenceexclude1.25km.RData")

#age
nam=list.files(path=here("NFI layers 2011"),pattern="Age",full.names=T)
age=rast(nam)
age=aggregate(age,2,mean)
age=project(age,crs(mz2021),res=500,method="near")
age=resample(age,mz2021,'near')
Age2011=round(as.matrix(extract(age,xy)),0)
length(which(is.na(Age2011)))
Age1=Age2011

Age=round(Age2011,0)
for(i in 1:4){
  e=Age2011-i
  e[which(e<0)]=0
  Age=cbind(e,Age)
}
for(i in 1:9){
  e=Age2011+i
  Age=cbind(Age,e)
}
colnames(Age)=2007:2020
save(Age,file="Climate Change/Age.RData")













#######################################################################################
#management zone leading 1 and holding 0
loc=xy
mz2021=rast(here("mz2021"))
mz2021=disagg(mz2021,4)
nam=list.files(path=here(),pattern="mz",full.names=T)
nam=nam[c(-15)]
rastm=lapply(nam,rast)
for(i in 1:14){
  rastm[[i]]=resample(rastm[[i]],mz2021,"near")
}
i=1
l=extract(rastm[[i]],loc)
l[is.na(l),1]=0
mz=l
for(i in 2:14){
  l=extract(rastm[[i]],loc)
  l[is.na(l),1]=0
  mz=cbind(mz,l)
}

#leading 1, all others 0
l1<-which(mz[,1]==2)
l2<-which(mz[,1]!=2)
mz[l1,1]<-1
mz[l2,1]<-0
#####
######
l1<-which(mz[,2]==2)
l2<-which(mz[,2]!=2)
mz[l1,2]<-1
mz[l2,2]<-0
###09 1 leading
mz[,3][which(mz[,3]!=1)]=0

l1<-which(mz[,4]==2)
l2<-which(mz[,4]!=2)
mz[l1,4]<-1
mz[l2,4]<-0

l1<-which(mz[,5]==2)
l2<-which(mz[,5]!=2)
mz[l1,5]<-1
mz[l2,5]<-0

l1<-which(mz[,6]==2)
l2<-which(mz[,6]!=2)
mz[l1,6]<-1
mz[l2,6]<-0

l1<-which(mz[,7]==2)
l2<-which(mz[,7]!=2)
mz[l1,7]<-1
mz[l2,7]<-0

l1<-which(mz[,8]==2)
l2<-which(mz[,8]!=2)
mz[l1,8]<-1
mz[l2,8]<-0

#15,16,17,18,19, 1 leading
for(i in 9:12){
  mz[,i][which(mz[,i]!=1)]=0
}
###
l1<-which(mz[,13]==2)
l2<-which(mz[,13]!=2)
mz[l1,13]<-1
mz[l2,13]<-0


l1<-which(mz[,14]==2)
l2<-which(mz[,14]!=2)
mz[l1,14]<-1
mz[l2,14]<-0

save(mz,file="Climate Change/mz.RData")
apply(mz,2,max)

#










