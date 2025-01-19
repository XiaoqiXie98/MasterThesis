library(terra)
library(here)
getwd()#check my current directory

#find the interested locations
#We are interested in the locations with an active management status (which means the cell value is 1 in the 'map' map)
map=rast(here("study_area.tif"))#load the raster map 
plot(map)#see how the map looks like
xy=xyFromCell(map,1:ncell(map))#extract all xy coordinates from this map. Coordinates contain for cell value 0 and 1.
v=extract(map,xy)#extract the cell value using the extracted xy coordinates

xy=xy[which(v==1),]#only leave cells with a value 1
print(dim(xy)[1])#how many cells (rows) have a cell value of 1
#Overall, we have 622286 cells

mz2021=rast(here("mz2021"))#load the raster map with values representing management status in 2021. We will use this map to align the position of cells in all other maps. We want cell positions in different maps to be same
mz2021=disagg(mz2021,4)#The orignal cell size is 2000 m. We disaggregate it into 500 m by separate each cell to four cells.

#In this script, we only count last year's infestations that will affect current year. We ignore the cut and burn trees.
#Load maps with detection of red-top tree (last year's infested tree) and green-dead tree (current year's infested tree).
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
GreenFromRed=as.numeric()
for(i in 1:14){#we extract values from 2007-2020(red-top trees), which is equivalent to 2006-2019(green-dead trees).
  e=as.matrix(extract(rastR_new[[i]],xy))#extract the number of red-top trees for interested positions
  e[is.na(e)]=0#For cells with NA, we set their values to be 0.
  GreenFromRed=cbind(GreenFromRed,e)#add column (add data for that year)
}
GreenFromRed=as.data.frame(GreenFromRed)#set class
for(i in 1:14){
  colnames(GreenFromRed)[i]=paste("GreenFromRed",2005+i)#rename columns
}
#detected green
detectedGreen=as.numeric()
for(i in 1:3){#add 2007,2008,2009 green-dead trees
  e=as.matrix(extract(rastG_new[[i]],xy))#extract cell values
  e[is.na(e)]=0#For those cells with NA, we set NA to be 0
  detectedGreen=cbind(detectedGreen,e)#add a column
}
detectedGreen=cbind(0,detectedGreen,0)#add 2006, 2006 has no green trees, similar for 2010
for(i in 4:9){#add 2011, 2012, 2013, 2014, 2015, 2016 green-dead trees
  e=extract(rastG_new[[i]],xy)#extract cell values
  e[is.na(e)]=0#For those cells with NA, we set NA to be 0
  detectedGreen=cbind(detectedGreen,e)#add a column
}
detectedGreen=cbind(detectedGreen,0,0)#2017,2018, no detected green
e=as.matrix(extract(rastG_new[[10]],xy))#extract values for 2019
e[is.na(e)]=0#For cells with NA, we replace NA by 0.
detectedGreen=cbind(detectedGreen,e)#add 2019
detectedGreen=as.data.frame(detectedGreen)#set class
for(i in 1:14){#rename columns
  colnames(detectedGreen)[i]=paste("DetectedGreen",i+2005)
}
#green trees can raise beetles are green from red + detectd green
LastYearGreen=detectedGreen+GreenFromRed
save(LastYearGreen,file="Climate Change/LastYearGreen.RData")

