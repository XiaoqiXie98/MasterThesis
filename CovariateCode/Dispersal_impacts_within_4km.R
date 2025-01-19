library(terra)
library(here)
getwd()#check my current directory

#find the interested locations
#We are interested in the locations with an active management status (which means the cell value is 1 in the 'map' map)
map=rast(here("study_area.tif"))#load the raster map. This map contains two values: 1 (active management), 0 (inactive) 
plot(map)#see how the map looks like
xy=xyFromCell(map,1:ncell(map))#extract all xy coordinates from this map. Coordinates contain for cell value 0 and 1.
v=extract(map,xy)#extract the cell value using the extracted xy coordinates

xy=xy[which(v==1),]#only leave cells with a value 1
print(dim(xy)[1])#how many cells (rows) have a cell value of 1
#Overall, we have 622286 cells

mz2021=rast(here("mz2021"))#load the raster map with values representing management status in 2021. We will use this map to align the position of cells in all other maps. We want cell positions in different maps to be same.
mz2021=disagg(mz2021,4)#The original cell size is 2000 m. We disaggregate it into 500 m by separating each cell to four cells.

#Load maps with detection of red-top tree (last year's infested tree) and green-dead tree (current year's infested tree).
#We did not consider cut and burn data because we only care about the influence of last year's infestations that will affect cuurent year.
#We have such data from 2007 to 2021 for red-top trees. Green-dead trees only observed in some years.
#I already separated maps of red-top trees and green-dead tree.
#In the following few lines, we use pattern="Red" and "Green" to call raster maps for red-top tree and green-dead tree.
#namR and namG give the directory of maps for red-top tree and green-dead tree.
namR=list.files(path=here("Infest Green Red"),pattern="Red",full.names=T)#red-top tree map
namR=namR[c(-1,-17)]#remove the unused directories
namG=list.files(path=here("Infest Green Red"),pattern="Green",full.names=T)#green-dead tree map
namG=namG[c(-12,-13)]#remove the unused directories
#
rastR=lapply(namR,rast)#load maps
rastG=lapply(namG,rast)#load maps
#align the cells in different maps by using management zone 2021 map

rastR_new=lapply(rastR,function(x) resample(x,mz2021,"near"))#resample cells by nearest neighborhood method
rastG_new=lapply(rastG,function(x) resample(x,mz2021,"near"))#resample cells by nearest neighborhood method

#

r1=rast()#create a new raster
crs(r1)=crs(mz2021)#coordinate system to be similar as mz2021.
ext(r1)=c(171965-4001,860465+4001,5442952-4001,6334952+4001)#set the extend of this new raster.
loc1=crop(mz2021,ext(r1))#crop mz2021 map by r1. I only keep the area with extend larger than the interested area by 4km. The reason of doing this is because I only consider a 4km radius neighborhood's impacts.

xy=xyFromCell(loc1,1:ncell(loc1))#xy contains coordinates from area r1.
df=xyFromCell(map,1:ncell(map))#coordinates in the area.
v=extract(map,df)#extract cell value by using coordinates.
df=df[which(v==1),]#df contains coordinates from the interested area with a value of 1.


library("dplyr")
#The function is used to count distance from cell x's center to cell y's center.
radius=function(df1,df2){
  return(sqrt((df1[,1]-df2[1])^2+(df1[,2]-df2[2])^2)/1000)
}#in km

#The function is used to only keep cells within a 4 km radius counting from the center of cells.
locWithin4km=function(df){
  r=radius(xy,df)#count distance from the center of interested cell to adjacent cells.
  data=data.frame(cbind(xy,r))#create table with coordinates and the corresponding distance to the terested cell.
  data=data%>%filter(r<=4)#only keep coordinates within a 4 km radius
  return(data)
}

all=read.csv(here("Influence of Zone/All.csv"))#This table contains the probability of finding beetles at different distances.
colnames(all)=c("distance","cdf")#This table shows the cumulative probability up to a distance from the center of target cell.
all$distance=round(all$distance,5)#only keep five decimals
all$cdf=round(all$cdf,5)#only keep five decimals
all[89,]=c(4,1)#set the last row. When the distance reaches 4 km, the cumulative probability reaches 1 as mentioned in Carroll 2017 report.
library('pbapply')#package used to run function in multi cores.
df=data.frame(df)#set data class
xy=data.frame(xy)#set data class
distance=pbapply(df,1,locWithin4km,cl=5)#all coordinates within a 4 km radius for each cell in the interested area.

#This function is used to find the probability of finding a beetle between distance x and y for each coordinate adjacent to each interested cell.
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

#Now, I estimate the above probability for each coordinate(x,y).
fx=list()
for(i in 1:dim(df)[1]){
  print(i)
  data=distance[[i]]#For each interested cell, data contains neighborhoods within a 4 km radius and the distance to the interested cell.
  pdf=lapply(data$r,closestR)#find the probability of finding a beetle between distance x and y.
  pdf=as.numeric(pdf)#set class
  pdf[is.na(pdf)]=0#for those locations with na, i set the probability to be 0.
  fx[[i]]=pdf
}

#This function is used to estimate the number of trees killed in the last year for a cell with coordinate (x,y).
#This function is very similar to the codes used in last year's infestations. For more details about this function, please chekc the notes in another script named "Last_year_infestations.R".
func1=function(xy){
  
  GreenFromRed=as.numeric()
  for(i in 1:15){
    e=as.matrix(extract(rastR_new[[i]],xy))
    e[is.na(e)]=0
    GreenFromRed=cbind(GreenFromRed,e)
  }
  colnames(GreenFromRed)=2006:2020

  detectedGreen=as.numeric()
  for(i in 1:3){
    e=as.matrix(extract(rastG_new[[i]],xy))
    e[is.na(e)]=0
    detectedGreen=cbind(detectedGreen,e)
  }
  detectedGreen=cbind(0,detectedGreen,0)#add 2006, 2006 has no green trees, similar for 2010
  for(i in 4:9){
    e=extract(rastG_new[[i]],xy)
    e[is.na(e)]=0
    detectedGreen=cbind(detectedGreen,e)
  }
  detectedGreen=cbind(detectedGreen,0,0)#2017,2018, no detected green
  e=extract(rastG_new[[10]],xy)
  e[is.na(e)]=0
  detectedGreen=cbind(detectedGreen,e)
  detectedGreen=cbind(detectedGreen,0)#2020 none
  colnames(detectedGreen)=2006:2020
  
  Green=GreenFromRed+detectedGreen
  return(Green)
}

#Calculate the number of trees killed in the last year for each neighborhood cell.
neib=list()
for(i in 1:dim(df)[1]){
  print(i)
  df1=as.matrix(distance[[i]])#let the coordinates of neighborhood to be one table
  neib[[i]]=func1(df1[,1:2])#estimate the number of infested trees in the adjacent cells based on the coordinates given in df1.
}

#calculate the influence by multiplying the probability with last year's infested trees.
influencewithin4km=as.numeric()
for(i in 1:dim(df)[1]){
  print(i)
  d1=fx[[i]]#This contains the probability corresponding to different distances to the center of target cell.
  d2=neib[[i]]#the number of infested trees 
  d3=as.numeric()
  for(j in 1:15){
    d3[j]=sum(d1*d2[,j])#multiply together to count the influence.
  }
  influencewithin4km=rbind(influencewithin4km,d3)
}
colmax=apply(influencewithin4km,2,max)

#Then, I need to remove the counted influence within 1.25 km radius. The following codes are very similar to the above lines with different considered radius (1.25 km).
#This function is used to estimate coordinates within 1.25 radius.
locWithin1.25km=function(df){
  location=xy%>%filter(x>df[1]-1250&x<df[1]+1250&y>df[2]-1250&y<df[2]+1250)
  r=radius(location,df)
  data=data.frame(cbind(location,r))
  data=data%>%filter(r<=1.25)
  return(data)
}
distance1=pbapply(df,1,locWithin1.25km,cl=8)#employ multicores to estimate coordinates of neighborhoods within 1.25 km radius couting the center of interested cell.

#Find the corresponding probabilities for different distances.
fx1=list()
for(i in 1:dim(df)[1]){
  print(i)
  data=distance1[[i]]#difference between two row in all is 0.0145
  pdf=lapply(data$r,closestR)
  pdf=as.numeric(pdf)
  pdf[is.na(pdf)]=0
  fx1[[i]]=pdf
}#pdf of each point

#Find neighborhoods
neib1=list()
for(i in 1:dim(df)[1]){
  print(i)
  df1=as.matrix(distance1[[i]][,1:2])
  neib1[[i]]=func1(df1)
}

#Count influence within 1.25 km.
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

#count the influence within 4km but excluding within 1.25 km radius as it has been counted by other covariates.
influenceexclude1.25km=influencewithin4km-influencewithin1.25km
save(influenceexclude1.25km,file="Climate Change/influenceexclude1.25km.RData")
